// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#ifdef HAVE_LIBCURL
#include <curl/curl.h> // needs to be included before siril.h for Windows
#endif
#include "core/siril_log.h"
#include "core/siril.h"
#include "core/siril_networking.h"
#include "io/local_catalogues.h"
#include "io/siril_catalogues.h"
#include <glib/gstdio.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <healpix_base.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <pointing.h>
#include <regex>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

extern "C" {
#include "core/siril_app_dirs.h"
}

//#define HEALPIX_DEBUG

#define ZENODO_GAIA_XPSAMP_RECORD_ID "YOUR_RECORD_ID_HERE"

extern const char** spcc_mirrors;

// RAII handle for a C FILE* opened through g_fopen(). Unlike std::ifstream /
// std::ofstream / std::filesystem::path, g_fopen() correctly handles UTF-8
// paths that contain non-ASCII characters on Windows (it converts to UTF-16
// and uses the wide CRT, instead of mis-decoding the bytes in the legacy ANSI
// code page). All catalogue/cache file I/O in this file must go through it so
// that SPCC works for users whose data dir contains national characters.
using SirilFilePtr = std::unique_ptr<FILE, int(*)(FILE*)>;
static SirilFilePtr siril_fopen_utf8(const std::string& path, const char* mode) {
    return SirilFilePtr(g_fopen(path.c_str(), mode), fclose);
}

// 64-bit fseek: offsets into the catalogue files can exceed 2 GB, but plain
// fseek() takes a 32-bit long on Windows. (std::ifstream::seekg used a 64-bit
// std::streamoff, so we must preserve that here.)
static int siril_fseek64(FILE* f, gint64 offset, int whence) {
#if defined(_WIN32)
    return _fseeki64(f, offset, whence);
#elif defined(__GLIBC__) || defined(__gnu_hurd__)
    return fseeko64(f, offset, whence);
#else
    return fseeko(f, (off_t)offset, whence);
#endif
}

// Enum for Gaia version designator
enum class GaiaVersion {
    DR1 = 0,
    DR2 = 1,
    EDR3 = 2,
    DR3 = 3,
    DR4 = 4,
    DR5 = 5
};

// Enum for Catalogue Type
enum class CatalogueType {
    Cat_None = 0,
    Astrometric = 1,
    Photometric_XP_Sampled = 2,
    Photometric_XP_Continuous = 3
};

// Struct for the Healpix catalogue header
#pragma pack(push, 1)  // Ensure no padding between members
typedef struct HealpixCatHeader {
    std::string title;          // 48 bytes: Catalogue title
    GaiaVersion gaia_version;   // 1 byte: Gaia version designator
    uint8_t healpix_level;      // 1 byte: Healpix indexing level N
    CatalogueType cat_type;     // 1 byte: Catalogue type
    uint8_t chunked;            // 1 byte to indicate whether the catalogue is chunked or not. Zero = monolithic, non-zero = chunked
    uint8_t chunk_level;        // for chunked catalogues, the Healpix level it waas chunked at
    uint32_t chunk_healpix;     // for chunked catalogues, the chunk-level healpixel covered by this file
    uint32_t chunk_first_healpixel; // for chunked catalogues, the index-level healpixel number of the first healpixel in the catalogue
    uint32_t chunk_last_healpixel; // for chunked catalogues, the index-level healpixel number of the last healpixel in the catalogue
    std::array<uint8_t, 63> spare; // 64 bytes reserved for future use
} HealpixCatHeader;
#pragma pack(pop)

// Structure to represent a range of consecutive source IDs
struct HealPixelRange {
    uint64_t start_id;
    uint64_t end_id;
};

// This function creates a range of consecutive HEALpixels: we can scan these
// continuously in the file, rather than doing a binary search for the start of
// each HEALpixel in the range. It is used for both astrometric and photometric
// catalogues.

static std::vector<HealPixelRange> create_healpixel_ranges(const std::vector<int>& pixels) {
    std::vector<HealPixelRange> ranges;
    if (pixels.empty()) {
        return ranges;
    }

    HealPixelRange current_range{
        static_cast<uint64_t>(pixels[0]),
        static_cast<uint64_t>(pixels[0])
    };

    for (size_t i = 1; i < pixels.size(); ++i) {
        uint64_t next_pixel = static_cast<uint64_t>(pixels[i]);
        if (next_pixel == current_range.end_id + 1) {
            // Extend current range
            current_range.end_id = next_pixel;
        } else {
            // Start new range
            ranges.push_back(current_range);
            current_range = HealPixelRange{next_pixel, next_pixel};
        }
    }
    // Add last range
    ranges.push_back(current_range);
    return ranges;
}

// Function to read the Healpix catalogue header
static HealpixCatHeader read_healpix_cat_header(const std::string& filename, int* error_status) {
    HealpixCatHeader header = {
        "",                     // title: empty string
        static_cast<GaiaVersion>(0),  // gaia_version: zero
        0,                      // healpix_level
        static_cast<CatalogueType>(0), // cat_type
        0,                      // chunked
        0,                      // chunk_level
        0,                      // chunk_healpix
        0,                      // chunk_first_healpixel
        0,                      // chunk_last_healpixel
        {}                      // spare: zero-initialized array
    };

    // Initialize error status
    if (error_status) {
        *error_status = 0;
    }

    // Open file in binary mode (g_fopen handles non-ASCII UTF-8 paths on Windows)
    SirilFilePtr file = siril_fopen_utf8(filename, "rb");
    if (!file) {
        if (error_status) {
            *error_status = -1; // File open error
        }
        return header;
    }
    FILE* fp = file.get();

    // Read the fixed-size string (48 bytes)
    char title_buffer[48] = {0};
    if (fread(title_buffer, 1, 48, fp) != 48) {
        if (error_status) {
            *error_status = -2; // Read error
        }
        return header;
    }
    header.title = std::string(title_buffer, strnlen(title_buffer, 48));

    // Read the remaining POD members (note: only 1 byte is read into the
    // enum/uint8 designator fields, matching the on-disk format)
    bool ok = true;
    ok = ok && fread(&header.gaia_version, 1, 1, fp) == 1;
    ok = ok && fread(&header.healpix_level, 1, 1, fp) == 1;
    ok = ok && fread(&header.cat_type, 1, 1, fp) == 1;
    ok = ok && fread(&header.chunked, 1, 1, fp) == 1;
    ok = ok && fread(&header.chunk_level, 1, 1, fp) == 1;
    ok = ok && fread(&header.chunk_healpix, sizeof(uint32_t), 1, fp) == 1;
    ok = ok && fread(&header.chunk_first_healpixel, sizeof(uint32_t), 1, fp) == 1;
    ok = ok && fread(&header.chunk_last_healpixel, sizeof(uint32_t), 1, fp) == 1;
    ok = ok && fread(&header.spare, 1, 63, fp) == 63;

    if (!ok) {
        if (error_status) {
            *error_status = -2; // Read error
        }
        return header;
    }

    return header;
}

// This function queries a catalogue of packed EntryType objects (the template allows
// use of SourceEntryAstro or SourceEntryXPsamp) for different uses. The catalogue MUST
// be sorted by source_id otherwise the method will fail. Each entry that matches the
// vector of source_id ranges is added to the result vector. Additional filtering (e.g.
// on magnitude) can be done by the caller.

template<typename EntryType>
static std::vector<EntryType> query_catalog(const std::string& filename, std::vector<HealPixelRange>& healpixel_ranges, const HealpixCatHeader& header) {
    constexpr size_t HEADER_SIZE = 128; // Fixed header size in bytes
    std::vector<EntryType> results;
    uint32_t nside = 1 << header.healpix_level; // Calculate NSIDE as 2^level
    uint32_t n_healpixels = 12 * nside * nside; // Total number of Healpixels

    if (header.chunked) {
        uint32_t nside_chunks = 1 << header.chunk_level;
        uint32_t n_chunks = 12 * nside_chunks * nside_chunks;
        n_healpixels /= n_chunks;
        for (auto& range : healpixel_ranges) {
            range.start_id -= header.chunk_first_healpixel;
            range.end_id -= header.chunk_first_healpixel;
        }
    }

    size_t INDEX_SIZE = (n_healpixels) * sizeof(uint32_t);

    // Open the catalog file in binary mode (g_fopen handles non-ASCII UTF-8 paths on Windows)
    SirilFilePtr file = siril_fopen_utf8(filename, "rb");
    if (!file) {
        siril_log_color_message(_("Failed to open file: %s\n"), "red", filename.c_str());
        return results;
    }
    FILE* fp = file.get();
    bool read_error = false;

    // Function to read a single index entry at a specific position
    auto read_index_entry = [fp, &results, &read_error](uint32_t healpixel_id) -> uint32_t {
        uint32_t index_value;
        gint64 pos = (gint64)HEADER_SIZE + (gint64)healpixel_id * sizeof(uint32_t);
        if (siril_fseek64(fp, pos, SEEK_SET) != 0 ||
                fread(&index_value, sizeof(uint32_t), 1, fp) != 1) {
            siril_log_color_message(_("Failed to read catalog index entry.\n"), "red");
            results.clear();
            read_error = true;
            return 0;
        }
        return index_value;
    };

    // Process each range
    for (const auto& range : healpixel_ranges) {
        uint32_t start_healpixel = range.start_id;
        uint32_t end_healpixel = range.end_id;

        if (start_healpixel >= n_healpixels || end_healpixel >= n_healpixels) {
            siril_log_color_message(_("ID range exceeds catalog bounds.\n"), "red");
            results.clear();
            return results;
        }

        // Read the index entries, using previous healpixel's value for start
        uint32_t start_offset = (start_healpixel == 0) ? 0 : read_index_entry(start_healpixel - 1);
        if (start_healpixel != 0 && read_error) {
            results.clear();
            return results;
        }

        uint32_t end_offset = read_index_entry(end_healpixel);
        if (read_error) {
            results.clear();
            return results;
        }

        // Calculate position in the data section
        gint64 data_start_pos = (gint64)HEADER_SIZE + (gint64)INDEX_SIZE + (gint64)start_offset * sizeof(EntryType);

        // Read the required data entries
        size_t num_records = end_offset - start_offset;
        std::vector<EntryType> buffer(num_records);

        if (siril_fseek64(fp, data_start_pos, SEEK_SET) != 0 ||
                fread(buffer.data(), sizeof(EntryType), num_records, fp) != num_records) {
            siril_log_color_message(_("Failed to read data entries.\n"), "red");
            results.clear();
            return results;
        }

        // Move entries to results vector
        results.insert(results.end(), buffer.begin(), buffer.end());
    }

    return results;
}

static bool header_compatible(HealpixCatHeader& a, HealpixCatHeader& b) {
    bool same_version = a.gaia_version == b.gaia_version;
    bool same_chunk_level = a.chunk_level == b.chunk_level;
    bool same_index_level = a.healpix_level == b.healpix_level;
    return same_version && same_chunk_level && same_index_level;
}

// Get Siril data directory and ensure the relevant subdir exists. Returns a
// UTF-8 path (built with GLib so non-ASCII data dirs work on Windows).
static std::string get_or_create_cache_dir() {
    const char *uddir = siril_get_user_data_dir();
    gchar *siril_subdir = g_build_filename(uddir, "spcc-cache", NULL);
    if (g_mkdir_with_parents(siril_subdir, 0755) != 0) {
        siril_debug_print(_("Can't create directory %s, falling back to system temp\n"),
			siril_subdir);
        g_free(siril_subdir);
        return std::string(g_get_tmp_dir());
    }
    std::string result(siril_subdir);
    g_free(siril_subdir);
    return result;
}

#ifdef HAVE_LIBCURL

// Function to read header via HTTP range request using Siril's networking
static HealpixCatHeader read_healpix_cat_header_http_with_curl(CURL* curl, const std::string& url, int* error_status) {
    HealpixCatHeader header = {
        "",
        static_cast<GaiaVersion>(0),
        0,
        static_cast<CatalogueType>(0),
        0,
        0,
        0,
        0,
        0,
        {}
    };

    if (error_status) {
        *error_status = 0;
    }

    // Setup Cache Path
    std::string tempdir = get_or_create_cache_dir();
    std::string filename = url.substr(url.find_last_of("/\\") + 1);
    gchar *cache_path_c = g_build_filename(tempdir.c_str(), (filename + ".header").c_str(), NULL);
    std::string cache_path(cache_path_c);
    g_free(cache_path_c);

    const size_t HEADER_SIZE = 128;

    bool cache_exists = false;

    // Validate the cache file
    GStatBuf st;
    if (g_stat(cache_path.c_str(), &st) == 0) {
        if ((size_t)st.st_size == HEADER_SIZE) {
            siril_debug_print(_("Header exists in cache\n"));
            cache_exists = true;
        } else {
            siril_log_color_message(_("Cache file %s is corrupted or incomplete. Deleting...\n"),
                                "salmon", cache_path.c_str());
            g_remove(cache_path.c_str());
        }
    }

    // Download if it doesn't exist
    if (!cache_exists) {
        siril_debug_print(_("Fetching header to cache: %s\n"), filename.c_str());

        gsize response_length;
        int error;
        char* buffer = fetch_url_range_with_curl((void*)curl, url.c_str(), 0, HEADER_SIZE, &response_length, &error, FALSE);

        if (error || !buffer || response_length < HEADER_SIZE) {
            if (error_status) *error_status = -1;
            free(buffer);
            return {}; // Return empty header
        }

        // Save to disk
        SirilFilePtr out = siril_fopen_utf8(cache_path, "wb");
        if (out) {
            fwrite(buffer, 1, HEADER_SIZE, out.get());
        }
        free(buffer);
    }

    return read_healpix_cat_header(cache_path, error_status);

}

template<typename EntryType>

static std::vector<EntryType> query_catalog_http_with_curl(CURL* curl,
                                                  const std::string& base_url,
                                                  const std::string& filename,
                                                  std::vector<HealPixelRange>& healpixel_ranges,
                                                  const HealpixCatHeader& header) {
    constexpr size_t HEADER_SIZE = 128;
    std::vector<EntryType> results;
    uint32_t nside = 1 << header.healpix_level;
    uint32_t n_healpixels = 12 * nside * nside;

    if (header.chunked) {
        uint32_t nside_chunks = 1 << header.chunk_level;
        uint32_t n_chunks = 12 * nside_chunks * nside_chunks;
        n_healpixels /= n_chunks;

        for (auto& range : healpixel_ranges) {
            range.start_id -= header.chunk_first_healpixel;
            range.end_id -= header.chunk_first_healpixel;
        }
    }

    size_t INDEX_SIZE = n_healpixels * sizeof(uint32_t);
    std::string full_url = base_url + "/" + filename;

    std::string tempdir = get_or_create_cache_dir();
    gchar *cache_path_c = g_build_filename(tempdir.c_str(), (filename + ".index").c_str(), NULL);
    std::string cache_path(cache_path_c);
    g_free(cache_path_c);

    // Load or download index
    std::vector<uint32_t> full_index(n_healpixels);

    bool cache_exists = false;
    GStatBuf st;
    if (g_stat(cache_path.c_str(), &st) == 0) {
        // Check if the size on disk matches our expected index size
        if ((size_t)st.st_size == INDEX_SIZE) {
            cache_exists = true;
        } else {
            siril_log_color_message(_("Cache file %s is corrupted or incomplete. Deleting...\n"),
                                "salmon", cache_path.c_str());
            g_remove(cache_path.c_str());
        }
    }

    // Ensure cache exists; if not, fetch and create it
    if (!cache_exists) {
        siril_debug_print(_("Fetching index to cache\n"));

        gsize response_length;
        int error;

        char* buffer = fetch_url_range_with_curl(
            (void*)curl,
            full_url.c_str(),
            HEADER_SIZE,          // offset
            INDEX_SIZE,           // length
            &response_length,
            &error,
            FALSE
        );

        if (error || !buffer || response_length != INDEX_SIZE) {
            siril_log_color_message(_("Failed to download index via HTTP\n"), "red");
            return results;
        }

        // Write buffer to cache file
        SirilFilePtr out = siril_fopen_utf8(cache_path, "wb");
        if (!out) {
            siril_log_color_message(_("Failed to write index cache file\n"), "red");
            free(buffer);
            return results;
        }

        fwrite(buffer, 1, INDEX_SIZE, out.get());
        out.reset();
        free(buffer);
    }

    // At this point, the cache file is guaranteed to exist.
    siril_debug_print(_("Loading index from cache\n"));

    SirilFilePtr in = siril_fopen_utf8(cache_path, "rb");
    if (!in) {
        siril_log_color_message(_("Failed to read index cache file\n"), "red");
        return results;
    }

    if (fread(full_index.data(), 1, INDEX_SIZE, in.get()) != INDEX_SIZE) {
        siril_log_color_message(_("Failed to read index cache file\n"), "red");
        return results;
    }

    // Process each healpixel range
    siril_debug_print(_("Fetching data\n"));
    for (const auto& range : healpixel_ranges) {
        uint32_t start_healpixel = range.start_id;
        uint32_t end_healpixel = range.end_id;

        if (start_healpixel >= n_healpixels || end_healpixel >= n_healpixels) {
            siril_log_color_message(_("ID range exceeds catalog bounds.\n"), "red");
            results.clear();
            return results;
        }

        // Index lookup from cached index
        uint32_t index_start = (start_healpixel == 0) ? 0 : full_index[start_healpixel - 1];
        uint32_t index_end   = full_index[end_healpixel];

        size_t num_records = index_end - index_start;

	// It's possible that we're fetching a single healpixel in this range 
	// and that the number of data records is zero
	if (num_records == 0) continue;

        size_t data_start_pos =
            HEADER_SIZE + INDEX_SIZE + index_start * sizeof(EntryType);

        size_t data_length = num_records * sizeof(EntryType);

        // Fetch data entries via HTTP
        gsize data_response_length;
        int data_error;

        siril_debug_print(_("Fetching range %zu - %zu\n"), data_start_pos, data_length);
        char* data_buffer = fetch_url_range_with_curl(
            (void*)curl,
            full_url.c_str(),
            data_start_pos,
            data_length,
            &data_response_length,
            &data_error,
            FALSE
        );

        if (data_error || !data_buffer) {
            siril_log_color_message(_("Failed to read data entries via HTTP\n"), "red");
            results.clear();
            return results;
        }

        std::vector<EntryType> buffer(num_records);
        memcpy(buffer.data(), data_buffer, data_response_length);
        free(data_buffer);

        results.insert(results.end(), buffer.begin(), buffer.end());
    }

    return results;
}


// Original template function now wraps the new one with a temporary CURL handle
template<typename EntryType>
static std::vector<EntryType> query_catalog_http(const std::string& base_url,
                                                  const std::string& filename,
                                                  std::vector<HealPixelRange>& healpixel_ranges,
                                                  const HealpixCatHeader& header) {
    CURL* curl = curl_easy_init();
    if (!curl) {
        siril_log_color_message(_("Error initialising CURL handle\n"), "red");
        return std::vector<EntryType>();
    }

    std::vector<EntryType> results = query_catalog_http_with_curl<EntryType>(curl, base_url, filename, healpixel_ranges, header);
    curl_easy_cleanup(curl);
    return results;
}

/**
 * Try fetching from the current mirror, and if it fails, try other mirrors
 * Updates com.spcc_remote_catalogue to the working mirror
 * Returns true if a working mirror was found, false otherwise
 */
static bool try_mirrors_and_update(CURL* curl, const char* path_suffix,
                                   HealpixCatHeader* header_out, int* error_status) {

    // First try the current mirror
    std::string current_url = std::string(com.spcc_remote_catalogue) + "/" + path_suffix;
    *header_out = read_healpix_cat_header_http_with_curl(curl, current_url, error_status);

    if (*error_status == 0) {
        return true;  // Current mirror works
    }

    siril_log_color_message(_("Mirror %s failed, trying alternatives...\n"),
                           "salmon", com.spcc_remote_catalogue);

    // Try each mirror in sequence
    for (int i = 0; spcc_mirrors[i] != NULL; i++) {
        // Skip if this is the current mirror (already tried)
        if (g_strcmp0(spcc_mirrors[i], com.spcc_remote_catalogue) == 0) {
            continue;
        }

        std::string test_url = std::string(spcc_mirrors[i]) + "/" + path_suffix;
        siril_log_message(_("Trying mirror: %s\n"), spcc_mirrors[i]);

        *header_out = read_healpix_cat_header_http_with_curl(curl, test_url, error_status);

        if (*error_status == 0) {
            // This mirror works! Update the global setting
            g_free(com.spcc_remote_catalogue);
            com.spcc_remote_catalogue = g_strdup(spcc_mirrors[i]);
            siril_log_color_message(_("Switched to working mirror: %s\n"),
                                   "green", com.spcc_remote_catalogue);
            return true;
        }
    }

    // All mirrors failed
    siril_log_color_message(_("All catalogue mirrors failed to respond.\n"), "red");
    return false;
}

/**
 * Wrapper for query_catalog_http_with_curl that implements mirror fallback
 */
/**
 * Wrapper for query_catalog_http_with_curl that implements mirror fallback
 * Note: Makes a copy of healpixel_ranges since query_catalog_http_with_curl
 * modifies them in place based on chunk offsets
 */
template<typename EntryType>
static std::vector<EntryType> query_catalog_http_with_fallback(
                                    CURL* curl,
                                    const std::string& filename,
                                    const std::vector<HealPixelRange>& healpixel_ranges,
                                    const HealpixCatHeader& header) {
    std::vector<EntryType> results;

    // Try current mirror first - make a copy since the function modifies ranges
    std::vector<HealPixelRange> ranges_copy = healpixel_ranges;
    std::string base_url(com.spcc_remote_catalogue);
    results = query_catalog_http_with_curl<EntryType>(curl, base_url, filename,
                                                      ranges_copy, header);

    if (!results.empty()) {
        return results;  // Success with current mirror
    }

    // Current mirror failed, try alternatives
    siril_log_color_message(_("Query failed on mirror %s, trying alternatives...\n"),
                           "salmon", com.spcc_remote_catalogue);

    for (int i = 0; spcc_mirrors[i] != NULL; i++) {
        // Skip if this is the current mirror (already tried)
        if (g_strcmp0(spcc_mirrors[i], com.spcc_remote_catalogue) == 0) {
            continue;
        }

        siril_log_message(_("Trying mirror: %s\n"), spcc_mirrors[i]);
        std::string test_base_url(spcc_mirrors[i]);

        // Make a fresh copy of the original ranges for each attempt
        ranges_copy = healpixel_ranges;
        results = query_catalog_http_with_curl<EntryType>(curl, test_base_url, filename,
                                                          ranges_copy, header);

        if (!results.empty()) {
            // This mirror works! Update the global setting
            g_free(com.spcc_remote_catalogue);
            com.spcc_remote_catalogue = g_strdup(spcc_mirrors[i]);
            siril_log_color_message(_("Switched to working mirror: %s\n"),
                                   "green", com.spcc_remote_catalogue);
            return results;
        }
    }

    // All mirrors failed
    siril_log_color_message(_("All catalogue mirrors failed for query.\n"), "red");
    return results;  // Return empty vector
}

#endif


#ifdef HEALPIX_DEBUG
// Function to show entries for a specific healpixel
static void show_healpixel_entries(uint32_t healpixel_id) {
    try {
        std::string filename(com.pref.catalogue_paths[4]);
        // Create a single-element range for the requested healpixel
        std::vector<HealPixelRange> range = {{healpixel_id, healpixel_id}};

        // Read the header first to get the healpix level
        int status = 0;
        HealpixCatHeader header = read_healpix_cat_header(std::string(filename), &status);

        // Query the catalog for this healpixel
        auto entries = query_catalog<SourceEntryAstro>(std::string(filename), range, header);

        // Print each entry with the same formatting as Python
        for (const auto& entry : entries) {
            double ra = entry.ra_scaled * 1.67638063509e-07;
            double dec = entry.dec_scaled * 1.67638063509e-07;
            double mag = entry.mag_scaled / 1000.0;

            std::cout << "Record for healpixid " << healpixel_id
            << ": ra=" << std::fixed << std::setprecision(6) << ra
            << ", dec=" << std::fixed << std::setprecision(5) << dec
            << ", pmra=" << entry.dra_scaled
            << ", pmdex=" << entry.ddec_scaled
            << ", teff=" << entry.teff
            << ", mag=" << std::fixed << std::setprecision(3) << mag
            << std::endl;
        }

        std::cout << "Numrecords = " << entries.size() << std::endl;
        // Print catalog information (matching Python output)
        std::cout << "Catalogue Title: " << header.title << std::endl;
        std::cout << "Gaia Version: " << static_cast<int>(header.gaia_version) << std::endl;
        std::cout << "Healpix Level: " << static_cast<int>(header.healpix_level) << std::endl;
        std::cout << "Catalogue Type: " << static_cast<int>(header.cat_type) << std::endl;


    }
    catch (const std::exception& e) {
        std::cerr << "Error showing healpixel entries: " << e.what() << std::endl;
    }
}
#endif

static int convert_healpix_level(int pixel_index, int from_level, int to_level) {
    if (from_level < to_level) {
        // Going to a finer level - shift left
        // this returns the numerically first fine-level healpixel within the coarser healpixel
        return pixel_index << (2 * (to_level - from_level));
    } else {
        // Going to a coarser level - shift right
        return pixel_index >> (2 * (from_level - to_level));
    }
}

// Assuming convert_healpix_level is defined as:

static std::vector<std::vector<int>> group_pixels_by_chunk(
            const std::vector<int>& chunk_indices,
            const std::vector<int>& pixel_indices,
            int chunk_level,
            int healpix_base_level) {
    // Map chunks to their positions for quick lookup
    std::map<int, size_t> chunk_map;
    for (size_t i = 0; i < chunk_indices.size(); ++i) {
        chunk_map[chunk_indices[i]] = i;
    }

    // Initialize the vector of vectors for each chunk
    std::vector<std::vector<int>> pixels_in_chunks(chunk_indices.size());

    // Iterate over pixel_indices
    for (int pixel : pixel_indices) {
        // Convert pixel index from base level to chunk level
        int chunk_id = convert_healpix_level(pixel, healpix_base_level, chunk_level);

        // Check if the chunk ID exists in the chunk map
        auto it = chunk_map.find(chunk_id);
        if (it != chunk_map.end()) {
            // Add the pixel to the appropriate vector
            pixels_in_chunks[it->second].push_back(pixel);
        }
    }

    return pixels_in_chunks;
}

template <typename T>
static std::vector<T> flatten(const std::vector<std::vector<T>>& nested_vectors) {
    std::vector<T> flattened;

    // Reserve space in advance for optimization
    size_t total_size = 0;
    for (const auto& sub_vector : nested_vectors) {
        total_size += sub_vector.size();
    }
    flattened.reserve(total_size);

    // Insert elements from each sub-vector into the flattened vector
    for (const auto& sub_vector : nested_vectors) {
        flattened.insert(flattened.end(), sub_vector.begin(), sub_vector.end());
    }

    return flattened;
}

static std::string find_matching_cat_file(std::string& path) {
    // Convert the pattern to a regex pattern
    // Replace %u with regex pattern for unsigned integers
    std::string pattern = "siril_cat(\\d+)_healpix(\\d+)_xpsamp_(\\d+)\\.dat";
    std::regex file_regex(pattern);

    // g_dir_open handles non-ASCII (UTF-8) directory paths on Windows
    GError *error = NULL;
    GDir *dir = g_dir_open(path.c_str(), 0, &error);
    if (!dir) {
        g_warning("Error accessing directory: %s", error ? error->message : "(unknown)");
        if (error)
            g_error_free(error);
        return "";
    }

    std::string match;
    const gchar *name;
    while ((name = g_dir_read_name(dir)) != NULL) {
        std::string filename(name);
        std::smatch matches;

        // Check if filename matches our pattern
        if (std::regex_match(filename, matches, file_regex)) {
            gchar *full = g_build_filename(path.c_str(), name, NULL);
            gboolean is_regular = g_file_test(full, G_FILE_TEST_IS_REGULAR);
            g_free(full);
            if (is_regular) {
                match = filename;
                break;
            }
        }
    }
    g_dir_close(dir);

    return match;
}

// This function is the main entry point and is declared extern "C" for ease of
// calling from the Siril C code. For integration into the code it will need to
// return an array of structs for use with SPCC, or possibly (maybe as a separate
// function) an array of deep_star structs for use with existing astrometry
// functions.

extern "C" {
    int local_gaia_xpsamp_available() {
        std::string chunkpath(com.pref.catalogue_paths[5]);
        std::string first_chunk = find_matching_cat_file(chunkpath);
        return (!first_chunk.empty());
    }

    int get_raw_stars_from_local_gaia_astro_catalogue(double ra, double dec, double radius, double limitmag, gboolean phot, deepStarData **stars, uint32_t *nb_stars) {

#ifdef HEALPIX_DEBUG
        show_healpixel_entries(0);
#endif

        radius /= 60.0; // the catalogue radius is in arcmin, we want it in degrees to convert to radians
        siril_debug_print("Search radius: %f deg\n", radius);
        const double DEG_TO_RAD = M_PI / 180.0;
        double radius_rad = radius * DEG_TO_RAD;
        double ra_rad = ra * DEG_TO_RAD;
        double dec_rad = dec * DEG_TO_RAD;
        double theta = M_PI / 2.0 - dec_rad;
        double phi = ra_rad;
        double radius_h = pow(sin(0.5 * radius_rad), 2);
        // Check the correct healpixel level and create our healpix_base
        std::string filename(com.pref.catalogue_paths[4]);
        int status = 0;
        HealpixCatHeader header = read_healpix_cat_header(filename, &status);
        if (status) {
            *stars = nullptr;
            *nb_stars = 0;
            return 1;
        }
        T_Healpix_Base<int> healpix_base(header.healpix_level, NEST);

        pointing point(theta, phi);
        std::vector<int> pixel_indices;
        siril_debug_print("query_disc_inclusive params: theta = %f, phi = %f, rad = %f\n", point.theta, point.phi, radius_rad);

        // Perform the cone search
        healpix_base.query_disc_inclusive(point, radius_rad, pixel_indices);

        // Reduce the pixels to a list of continuous ranges
        std::vector<HealPixelRange> healpixel_ranges = create_healpixel_ranges(pixel_indices);

        // Query the catalogue for matches
        std::vector<SourceEntryAstro> matches = query_catalog<SourceEntryAstro>(filename, healpixel_ranges, header);

        // Check if no matches were found, return 1 in that case
        if (matches.empty()) {
            *nb_stars = 0;
            *stars = nullptr;
            return 1;
        }

        // Filter by magnitude
        double scaled_limitmag = limitmag * 1000.0;
        matches.erase(
            std::remove_if(matches.begin(), matches.end(),
                           [scaled_limitmag](const SourceEntryAstro& entry) {
                               return entry.mag_scaled > scaled_limitmag;
                           }
            ),
            matches.end()
        );

        // Filter by distance from the center of the cone
        double scale = 360.0 / (double) INT32_MAX;
        matches.erase(
            std::remove_if(matches.begin(), matches.end(),
                           [scale, radius_h, ra, dec](const SourceEntryAstro& entry) {
                               return compute_coords_distance_h(ra, dec, (double)entry.ra_scaled * scale, (double)entry.dec_scaled * scale) > radius_h;
                           }
            ),
            matches.end()
        );

        // If we are retrieving stars for PCC, remove any without a valid Teff entry
        if (phot) {
            matches.erase(
                std::remove_if(matches.begin(), matches.end(),
                               [](const SourceEntryAstro& entry) {
                                   return entry.teff <= 0;
                               }
                ),
                matches.end()
            );
        }

        // Check if no matches remain after filtering, again return 1 in that case
        if (matches.empty()) {
            *nb_stars = 0;
            *stars = nullptr;
            return 1;
        }

        // Populate return data
        *nb_stars = matches.size();
        *stars = (deepStarData*) malloc(*nb_stars * sizeof(deepStarData));
        if (*stars == nullptr) {
            return -1;  // Memory allocation failed, return -1
        }

        uint32_t i = 0;
        for (const auto& entry : matches) {
            (*stars)[i++] = (deepStarData) {
                .RA = entry.ra_scaled,
                .Dec = entry.dec_scaled,
                .dRA = entry.dra_scaled,
                .dDec = entry.ddec_scaled,
                .B = static_cast<int16_t>(entry.teff), // this knowingly abuses the struct, .B must be cast back to uint16_t on retrieval
                .V = entry.mag_scaled
            };
        }

        return 0;
    }

    int get_raw_stars_from_local_gaia_xpsampled_catalogue(double ra, double dec, double radius, double limitmag, SourceEntryXPsamp **stars, uint32_t *nb_stars) {
        radius /= 60.0; // the catalogue radius is in arcmin, we want it in degrees to convert to radians
        siril_debug_print("Search radius: %f deg\n", radius);
        const double DEG_TO_RAD = M_PI / 180.0;
        double radius_rad = radius * DEG_TO_RAD;
        double ra_rad = ra * DEG_TO_RAD;
        double dec_rad = dec * DEG_TO_RAD;
        double theta = M_PI / 2.0 - dec_rad;
        double phi = ra_rad;
        double radius_h = pow(sin(0.5 * radius_rad), 2);

        // Read the basic catalogue header from the first available chunk that matches
        // the expected pattern

        std::string chunkpath(com.pref.catalogue_paths[5]);
        std::string first_chunk = find_matching_cat_file(chunkpath);

        if (first_chunk.empty()) {
            siril_log_color_message(_("Error: no chunks detected.\n"), "red");
            *stars = nullptr;
            *nb_stars = 0;
            return 1;
        }
        // Check the correct healpixel level and create our healpix_base
        int status = 0;
        gchar *final_path_c = g_build_filename(chunkpath.c_str(), first_chunk.c_str(), NULL);
        std::string final_path(final_path_c);
        g_free(final_path_c);
        HealpixCatHeader header = read_healpix_cat_header(final_path, &status);
        if (status) {
            *stars = nullptr;
            *nb_stars = 0;
            return 1;
        }

        T_Healpix_Base<int> healpix_base(header.healpix_level, NEST); // For getting the pixel indices to look up

        pointing point(theta, phi);

        // Perform the cone search
        std::vector<int> pixel_indices;
        healpix_base.query_disc_inclusive(point, radius_rad, pixel_indices);

        // Get the set of unique chunks
        std::set<int> chunks;  // Use a set for automatic uniqueness
        for (const int& healpix : pixel_indices) {
            chunks.insert(convert_healpix_level(healpix, header.healpix_level, header.chunk_level));
        }
        std::vector<int> chunk_indices(chunks.begin(), chunks.end());

        // Split the pixel_indices into separate vectors for each chunk
        std::vector<std::vector<int>> chunked_pixel_indices = group_pixels_by_chunk(chunk_indices, pixel_indices, header.chunk_level, header.healpix_level);

        // Initialize the vector of results vectors for each chunk
        std::vector<std::vector<SourceEntryXPsamp>> results_in_chunks(chunk_indices.size());

        // Iterate over the chunk indices and their pixel lists
        bool file_error = false;
        for (size_t i = 0; i < chunk_indices.size(); ++i) {
            const int chunk_id = chunk_indices[i];
            const std::vector<int>& chunk_pixels = chunked_pixel_indices[i];
            // Reduce the pixels to a list of continuous ranges
            std::vector<HealPixelRange> healpixel_ranges = create_healpixel_ranges(chunk_pixels);

            gchar *filename = g_strdup_printf("siril_cat%u_healpix%u_xpsamp_%d.dat", header.chunk_level, header.healpix_level, chunk_id);
            gchar *this_chunk_path_c = g_build_filename(chunkpath.c_str(), filename, NULL);
            std::string this_chunk_path(this_chunk_path_c);
            g_free(this_chunk_path_c);
            g_free(filename);

            if (!g_file_test(this_chunk_path.c_str(), G_FILE_TEST_EXISTS)) {
                siril_log_color_message(_("Chunk file not found: %s\n"), "red", this_chunk_path.c_str());
                file_error = true;
                break;
            }
            // Read this specific header
            int status = 0;
            HealpixCatHeader this_header = read_healpix_cat_header(this_chunk_path, &status);
            if (!header_compatible(header, this_header)) {
                siril_log_color_message(_("Error: catalog header values for chunk %lu are incompatible with previous values. All chunk files must have the same chunk level and indexing level and must represent the same Gaia data release\n"), "red");
                return 1;
            }
            // Query the catalogue for matches
            results_in_chunks[i] = query_catalog<SourceEntryXPsamp>(this_chunk_path, healpixel_ranges, this_header);
        }
        if (file_error)
            return 1;

        std::vector<SourceEntryXPsamp> matches = flatten(results_in_chunks);


        // Check if no matches were found, return 1 in that case
        if (matches.empty()) {
            *nb_stars = 0;
            *stars = nullptr;
            return 1;
        }

        // Filter by magnitude
        double scaled_limitmag = limitmag * 1000.0;
        matches.erase(
            std::remove_if(matches.begin(), matches.end(),
                           [scaled_limitmag](const SourceEntryXPsamp& entry) {
                               return entry.mag_scaled > scaled_limitmag;
                           }
            ),
            matches.end()
        );

        // Filter by distance from the center of the cone
        double scale = 360.0 / (double) INT32_MAX;
        matches.erase(
            std::remove_if(matches.begin(), matches.end(),
                           [scale, radius_h, ra, dec](const SourceEntryXPsamp& entry) {
                               return compute_coords_distance_h(ra, dec, (double)entry.ra_scaled * scale, (double)entry.dec_scaled * scale) > radius_h;
                           }
            ),
            matches.end()
        );

        // Check if no matches remain after filtering, again return 1 in that case
        if (matches.empty()) {
            *nb_stars = 0;
            *stars = nullptr;
            return 1;
        }

        // Populate return data
        *nb_stars = matches.size();
        *stars = (SourceEntryXPsamp*) malloc(*nb_stars * sizeof(SourceEntryXPsamp));
        if (*stars == nullptr) {
            return -1;  // Memory allocation failed, return -1
        }
        std::copy(matches.begin(), matches.end(), *stars);

        return 0;
    }

#ifdef HAVE_LIBCURL
    int get_raw_stars_from_remote_gaia_xpsampled_catalogue(double ra, double dec, double radius,
                                                           double limitmag, SourceEntryXPsamp **stars,
                                                           uint32_t *nb_stars) {
        radius /= 60.0;
        siril_debug_print("Search radius: %f deg\n", radius);
        const double DEG_TO_RAD = M_PI / 180.0;
        double radius_rad = radius * DEG_TO_RAD;
        double ra_rad = ra * DEG_TO_RAD;
        double dec_rad = dec * DEG_TO_RAD;
        double theta = M_PI / 2.0 - dec_rad;
        double phi = ra_rad;
        double radius_h = pow(sin(0.5 * radius_rad), 2);

        // Initialize CURL handle for all requests
        CURL* curl = curl_easy_init();
        if (!curl) {
            siril_log_color_message(_("Error initialising CURL handle\n"), "red");
            *stars = nullptr;
            *nb_stars = 0;
            return 1;
        }

        // Read header from first chunk with mirror fallback
        std::string first_chunk_filename = "siril_cat1_healpix8_xpsamp_0.dat";

        int status = 0;
        HealpixCatHeader header;
        if (!try_mirrors_and_update(curl, first_chunk_filename.c_str(), &header, &status)) {
            curl_easy_cleanup(curl);
            *stars = nullptr;
            *nb_stars = 0;
            return 1;
        }

        T_Healpix_Base<int> healpix_base(header.healpix_level, NEST);
        pointing point(theta, phi);

        // Perform cone search
        std::vector<int> pixel_indices;
        healpix_base.query_disc_inclusive(point, radius_rad, pixel_indices);

        // Get unique chunks
        std::set<int> chunks;
        for (const int& healpix : pixel_indices) {
            chunks.insert(convert_healpix_level(healpix, header.healpix_level, header.chunk_level));
        }
        std::vector<int> chunk_indices(chunks.begin(), chunks.end());

        // Split pixel_indices into separate vectors for each chunk
        std::vector<std::vector<int>> chunked_pixel_indices =
            group_pixels_by_chunk(chunk_indices, pixel_indices, header.chunk_level, header.healpix_level);

        std::vector<std::vector<SourceEntryXPsamp>> results_in_chunks(chunk_indices.size());

        bool file_error = false;
        for (size_t i = 0; i < chunk_indices.size(); ++i) {
            const int chunk_id = chunk_indices[i];
            const std::vector<int>& chunk_pixels = chunked_pixel_indices[i];

            std::vector<HealPixelRange> healpixel_ranges = create_healpixel_ranges(chunk_pixels);

            gchar *filename = g_strdup_printf("siril_cat%u_healpix%u_xpsamp_%d.dat",
                                             header.chunk_level, header.healpix_level, chunk_id);
            std::string chunk_filename(filename);
            g_free(filename);

            // Verify chunk header with fallback - this helps catch if a mirror
            // has incomplete data
            std::string chunk_url = std::string(com.spcc_remote_catalogue) + "/" + chunk_filename;
            int chunk_status = 0;
            HealpixCatHeader this_header = read_healpix_cat_header_http_with_curl(curl, chunk_url, &chunk_status);

            // If initial header read fails, try other mirrors
            if (chunk_status != 0) {
                if (!try_mirrors_and_update(curl, chunk_filename.c_str(), &this_header, &chunk_status)) {
                    siril_log_color_message(_("Failed to read header for chunk %d from any mirror\n"),
                                          "red", chunk_id);
                    file_error = true;
                    break;
                }
            }

            if (!header_compatible(header, this_header)) {
                siril_log_color_message(_("Error: catalog header values for chunk %d are incompatible\n"),
                                       "red", chunk_id);
                file_error = true;
                break;
            }

            // Query with automatic fallback
            results_in_chunks[i] = query_catalog_http_with_fallback<SourceEntryXPsamp>(
                curl, chunk_filename, healpixel_ranges, this_header);

            if (results_in_chunks[i].empty() && !chunk_pixels.empty()) {
                // Query returned nothing but we expected data - this might indicate failure
                siril_log_color_message(_("Warning: No data retrieved for chunk %d\n"),
                                       "salmon", chunk_id);
            }
        }

        // Clean up CURL handle
        curl_easy_cleanup(curl);

        if (file_error) {
            *stars = nullptr;
            *nb_stars = 0;
            return 1;
        }

        std::vector<SourceEntryXPsamp> matches = flatten(results_in_chunks);

        if (matches.empty()) {
            *nb_stars = 0;
            *stars = nullptr;
            return 1;
        }

        // Filter by magnitude
        double scaled_limitmag = limitmag * 1000.0;
        matches.erase(
            std::remove_if(matches.begin(), matches.end(),
                          [scaled_limitmag](const SourceEntryXPsamp& entry) {
                              return entry.mag_scaled > scaled_limitmag;
                          }),
            matches.end()
        );

        // Filter by distance
        double scale = 360.0 / (double) INT32_MAX;
        matches.erase(
            std::remove_if(matches.begin(), matches.end(),
                          [scale, radius_h, ra, dec](const SourceEntryXPsamp& entry) {
                              return compute_coords_distance_h(ra, dec,
                                  (double)entry.ra_scaled * scale,
                                  (double)entry.dec_scaled * scale) > radius_h;
                          }),
            matches.end()
        );

        if (matches.empty()) {
            *nb_stars = 0;
            *stars = nullptr;
            return 1;
        }

        // Populate return data
        *nb_stars = matches.size();
        *stars = (SourceEntryXPsamp*) malloc(*nb_stars * sizeof(SourceEntryXPsamp));
        if (*stars == nullptr) {
            return -1;
        }
        std::copy(matches.begin(), matches.end(), *stars);

        return 0;
    }
#else
    int get_raw_stars_from_remote_gaia_xpsampled_catalogue(double ra, double dec, double radius,
                                                           double limitmag, SourceEntryXPsamp **stars,
                                                           uint32_t *nb_stars) {
        siril_log_color_message(_("Error: Siril was compiled without libcurl support. Remote catalogue access is unavailable.\n"), "red");
        *stars = nullptr;
        *nb_stars = 0;
        return 1;
    }
#endif
}
