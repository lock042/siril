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
#include "io/healpix/fluxcache.h"
#include "io/healpix/healpix_cat.h"
#include "io/healpix/xp_continuous.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <healpix_base.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <pointing.h>
#include <regex>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <optional>
#include <string_view>
#ifndef G_PI
#define G_PI 3.14159265358979323846  /* pi */
#endif

extern "C" {
#include "core/siril_app_dirs.h"
}

//#define HEALPIX_DEBUG

#define ZENODO_GAIA_XPSAMP_RECORD_ID "YOUR_RECORD_ID_HERE"

extern "C" gchar **spcc_mirrors;
extern "C" gchar **spcc_mirrors_xpcts;

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
   char title[48];          // 48 bytes: Catalogue title
   uint8_t  gaia_version;   // 1 byte: Gaia version designator
   uint8_t  healpix_level;      // 1 byte: Healpix indexing level N
   uint8_t  cat_type;     // 1 byte: Catalogue type
   uint8_t  chunked;            // 1 byte to indicate whether the catalogue is chunked or not. Zero = monolithic, non-zero = chunked
   uint8_t  chunk_level;        // for chunked catalogues, the Healpix level it waas chunked at
   uint32_t chunk_healpix;     // for chunked catalogues, the chunk-level healpixel covered by this file
   uint32_t chunk_first_healpixel; // for chunked catalogues, the index-level healpixel number of the first healpixel in the catalogue
   uint32_t chunk_last_healpixel; // for chunked catalogues, the index-level healpixel number of the last healpixel in the catalogue
   uint8_t  spare[63]; // 63 bytes reserved for future use
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

    // Open the catalog file in binary mode
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        siril_log_error(_("Failed to open file: %s\n"), filename.c_str());
        return results;
    }

    // Function to read a single index entry at a specific position
    auto read_index_entry = [&file, &results](uint32_t healpixel_id) -> uint32_t {
        uint32_t index_value;
        size_t pos = HEADER_SIZE + healpixel_id * sizeof(uint32_t);
        file.seekg(pos, std::ios::beg);
        file.read(reinterpret_cast<char*>(&index_value), sizeof(uint32_t));
        if (!file) {
            siril_log_error(_("Failed to read catalog index entry.\n"));
            results.clear();
            return 0;
        }
        return index_value;
    };

    // Process each range
    for (const auto& range : healpixel_ranges) {
        uint32_t start_healpixel = range.start_id;
        uint32_t end_healpixel = range.end_id;

        if (start_healpixel >= n_healpixels || end_healpixel >= n_healpixels) {
            siril_log_error(_("ID range exceeds catalog bounds.\n"));
            results.clear();
            return results;
        }

        // Read the index entries, using previous healpixel's value for start
        uint32_t start_offset = (start_healpixel == 0) ? 0 : read_index_entry(start_healpixel - 1);
        if (start_healpixel != 0 && !file) {
            results.clear();
            return results;
        }

        uint32_t end_offset = read_index_entry(end_healpixel);
        if (!file) {
            results.clear();
            return results;
        }

        // Calculate position in the data section
        size_t data_start_pos = HEADER_SIZE + INDEX_SIZE + start_offset * sizeof(EntryType);

        // Read the required data entries
        size_t num_records = end_offset - start_offset;
        std::vector<EntryType> buffer(num_records);

        file.seekg(data_start_pos, std::ios::beg);
        file.read(reinterpret_cast<char*>(buffer.data()), num_records * sizeof(EntryType));

        if (!file) {
            siril_log_error(_("Failed to read data entries.\n"));
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

// Function to read the Healpix catalogue header
static HealpixCatHeader read_healpix_cat_header(const std::string& filename, int* error_status) {
    HealpixCatHeader header{};

    // Initialize error status
    if (error_status) {
        *error_status = 0;
    }

    // Open file in binary mode
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        if (error_status) {
            *error_status = -1; // File open error
        }
        return header;
    }

    // Read the header
    file.read(reinterpret_cast<char*>(&header), sizeof(header));
    if (!file) {
        *error_status = READ_ERROR;
        return {};
    }

    return header;
}

#ifdef HAVE_LIBCURL
static HealpixCatHeader read_healpix_cat_header_http_with_curl(CURL* curl, const std::string& url, int* error_status) {
    if (error_status) *error_status = 0;

    // Setup Cache Path
    auto tempdir = FluxCache::get_or_create_cache_dir();
    std::string filename = url.substr(url.find_last_of("/\\") + 1);
    std::string cache_path = (tempdir / (filename + ".header")).string();

    const size_t HEADER_SIZE = 128;

    bool cache_exists = false;
    std::error_code ec;

    // Validate the cache file
    if (std::filesystem::exists(cache_path, ec)) {
        if (std::filesystem::file_size(cache_path, ec) == HEADER_SIZE) {
            cache_exists = true;
        } else {
            siril_log_warning(_("Cache file %s is corrupted or incomplete. Deleting...\n"), cache_path.c_str());
            std::filesystem::remove(cache_path, ec);
        }
    }

    // Download if it doesn't exist
    if (!cache_exists) {
        siril_log_debug(_("Fetching header to cache: %s\n"), filename.c_str());

        gsize response_length;
        int error;
        char* buffer = fetch_url_range_with_curl((void*)curl, url.c_str(), 0, HEADER_SIZE, &response_length, &error, FALSE);

        if (error || !buffer || response_length < HEADER_SIZE) {
            if (error_status) *error_status = -1;
            free(buffer);
            return {}; // Return empty header
        }

        // Save to disk
        std::ofstream out(cache_path, std::ios::binary);
        if (out.good()) {
            out.write(buffer, HEADER_SIZE);
            out.close();
        }
        free(buffer);
    }

    return read_healpix_cat_header(cache_path, error_status);
}

template<typename EntryType>

static std::optional<std::vector<EntryType>> query_catalog_http_with_curl(CURL* curl,
                                                  const std::string& base_url,
                                                  const std::string& filename,
                                                  std::vector<HealPixelRange>& healpixel_ranges,
                                                  const HealpixCatHeader& header,
                                                  const char *type_tag) {
    constexpr size_t HEADER_SIZE = 128;
    std::vector<EntryType> results;
    uint32_t nside = 1 << header.healpix_level;
    uint32_t n_healpixels = 12 * nside * nside;

    // Do not mutate caller-owned ranges when converting to chunk-local IDs
    std::vector<HealPixelRange> local_ranges = healpixel_ranges;

    // Convert to chunk local range IDs if necessary
    if (header.chunked) {
        uint32_t nside_chunks = 1 << header.chunk_level;
        uint32_t n_chunks = 12 * nside_chunks * nside_chunks;
        n_healpixels /= n_chunks;

        for (auto& range : local_ranges) {
            range.start_id -= header.chunk_first_healpixel;
            range.end_id -= header.chunk_first_healpixel;
        }
    }

    size_t INDEX_SIZE = n_healpixels * sizeof(uint32_t);
    std::string full_url = base_url + "/" + filename;

    auto tempdir = FluxCache::get_or_create_cache_dir();
    std::string cache_path = (tempdir / (filename + ".index")).string();

    // Load or download index
    std::vector<uint32_t> full_index(n_healpixels);

    bool cache_exists = false;
    std::error_code ec;
    if (std::filesystem::exists(cache_path, ec)) {
        // Check if the size on disk matches our expected index size
        if (std::filesystem::file_size(cache_path, ec) == INDEX_SIZE) {
            cache_exists = true;
        } else {
            siril_log_warning(_("Cache file %s is corrupted or incomplete. Deleting...\n"), cache_path.c_str());
            std::filesystem::remove(cache_path, ec);
        }
    }

    // Ensure cache exists; if not, fetch and create it
    if (!cache_exists) {
        siril_log_debug(_("Fetching index to cache\n"));

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
            siril_log_error(_("Failed to download index via HTTP\n"));
            return std::nullopt;
        }

        // Write buffer to cache file
        std::ofstream out(cache_path, std::ios::binary);
        if (!out.good()) {
            siril_log_error(_("Failed to write index cache file\n"));
            free(buffer);
            return std::nullopt;
        }

        out.write(buffer, INDEX_SIZE);
        out.close();
        free(buffer);
    }

    // At this point, the cache file is guaranteed to exist.
    siril_log_debug(_("Loading index from cache\n"));

    std::ifstream in(cache_path, std::ios::binary);
    if (!in.good()) {
        siril_log_error(_("Failed to read index cache file\n"));
        return std::nullopt;
    }

    in.read(reinterpret_cast<char*>(full_index.data()), INDEX_SIZE);
    in.close();

    // Process each healpixel range
    for (const auto& range : local_ranges) {
        uint32_t start_healpixel = range.start_id;
        uint32_t end_healpixel = range.end_id;

        if (start_healpixel >= n_healpixels || end_healpixel >= n_healpixels) {
            siril_log_error(_("ID range exceeds catalog bounds.\n"));
            results.clear();
            return std::nullopt;
        }

        // Index lookup from cached index
        uint32_t index_start = (start_healpixel == 0) ? 0 : full_index[start_healpixel - 1];
        uint32_t index_end   = full_index[end_healpixel];

        size_t num_records = index_end - index_start;
	char* data_buffer = nullptr;

	// Only fetch via HTTP if there's something to fetch
	if (num_records > 0) {

            size_t data_start_pos =
                HEADER_SIZE + INDEX_SIZE + index_start * sizeof(EntryType);

            size_t data_length = num_records * sizeof(EntryType);

            // Fetch data entries via HTTP
            gsize data_response_length;
            int data_error;

            siril_log_debug(_("Fetching range %zu - %zu\n"), data_start_pos, data_length);
            data_buffer = fetch_url_range_with_curl(
                (void*)curl,
                full_url.c_str(),
                data_start_pos,
                data_length,
                &data_response_length,
                &data_error,
                FALSE
            );

            if (data_error || !data_buffer) {
                siril_log_error(_("Failed to read data entries via HTTP\n"));
                results.clear();
                return std::nullopt;
            }

           std::vector<EntryType> buffer(num_records);
           memcpy(buffer.data(), data_buffer, data_response_length);
           results.insert(results.end(), buffer.begin(), buffer.end());
	}

	// Get our data cache regardless. We always want to write a cache entry even if there's no data
	auto cache = FluxCache::getCache(header.chunk_level, header.healpix_level, header.chunk_healpix, type_tag);
        std::vector<FluxCache::CacheSegment> segments;

	// Unfortunately, the binary data does not contain the healpix id so we need process it
	// to split it into healpixel chunks before we pass to the cache
        uint32_t range_base_record = (range.start_id == 0) ? 0 : full_index[range.start_id - 1];
        for (uint32_t local_hp = range.start_id; local_hp <= range.end_id; ++local_hp) {
            // Get the offsets from the local index
            uint32_t p_start = (local_hp == 0) ? 0 : full_index[local_hp - 1];
            uint32_t p_end   = full_index[local_hp];
            uint32_t count   = p_end - p_start;

            // Convert local_hp to global_hp if the file is chunked
            uint32_t global_hp = local_hp;
            if (header.chunked) {
                global_hp += header.chunk_first_healpixel;
            }

            if (count > 0 && data_buffer) {
                // Calculate byte offset into the data_buffer
                size_t byte_offset = (p_start - range_base_record) * sizeof(EntryType);
                segments.push_back({global_hp, data_buffer + byte_offset, count * sizeof(EntryType)});
            } else {
                segments.push_back({global_hp, nullptr, 0});
	    }
        }

	// Commit to cache
	if (!segments.empty()) {
            siril_log_debug(_("Committing %zu pixels (%zu total records) to SQLite cache\n"),
                              segments.size(), num_records);
            cache->setCacheEntry(segments);
        } else {
            // Shouldn't happen
            siril_log_debug(_("No data found in range %u-%u to cache.\n"), start_healpixel, end_healpixel);
        }

        free(data_buffer);

    }

    return results;
}


// Original template function now wraps the new one with a temporary CURL handle
// TODO: Is this required?
/*
template<typename EntryType>
static std::vector<EntryType> query_catalog_http(const std::string& base_url,
                                                  const std::string& filename,
                                                  std::vector<HealPixelRange>& healpixel_ranges,
                                                  const HealpixCatHeader& header) {
    CURL* curl = curl_easy_init();
    if (!curl) {
        siril_log_error(_("Error initialising CURL handle\n"));
        return std::vector<EntryType>();
    }

    std::vector<EntryType> results = query_catalog_http_with_curl<EntryType>(curl, base_url, filename, healpixel_ranges, header);
    curl_easy_cleanup(curl);
    return results;
}
*/

/**
 * Try fetching from the current mirror, and if it fails, walk the supplied
 * mirror list. Updates *current_mirror_ptr in place to the working mirror so
 * subsequent calls in this session start there. Returns true if any mirror
 * succeeded, false otherwise. Caller-owned globals so the same code can
 * service xpsamp and xpcts catalogue kinds independently.
 */
static bool read_healpix_cat_header_http_with_fallback(CURL* curl, const char* path_suffix,
                                   HealpixCatHeader* header_out, int* error_status,
                                   gchar **mirrors, gchar **current_mirror_ptr) {
    if (!mirrors || !mirrors[0]) {
        *error_status = -1;
        return false;
    }
    std::string current_url = std::string(*current_mirror_ptr) + "/" + path_suffix;
    *header_out = read_healpix_cat_header_http_with_curl(curl, current_url, error_status);
    if (*error_status == 0) return true;

    siril_log_warning(_("Mirror %s failed, trying alternatives...\n"), *current_mirror_ptr);

    for (int i = 0; mirrors[i] != NULL; i++) {
        if (g_strcmp0(mirrors[i], *current_mirror_ptr) == 0) continue;

        std::string test_url = std::string(mirrors[i]) + "/" + path_suffix;
        siril_log_message(_("Trying mirror: %s\n"), mirrors[i]);
        *header_out = read_healpix_cat_header_http_with_curl(curl, test_url, error_status);

        if (*error_status == 0) {
            g_free(*current_mirror_ptr);
            *current_mirror_ptr = g_strdup(mirrors[i]);
            siril_log_info(_("Switched to working mirror: %s\n"), *current_mirror_ptr);
            return true;
        }
    }
    siril_log_error(_("All catalogue mirrors failed to respond.\n"));
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
                                    const HealpixCatHeader& header,
                                    gchar **mirrors, gchar **current_mirror_ptr,
                                    const char *type_tag) {
    if (!mirrors || !mirrors[0]) return {};

    std::vector<HealPixelRange> ranges_copy = healpixel_ranges;
    std::string base_url(*current_mirror_ptr);
    auto response = query_catalog_http_with_curl<EntryType>(curl, base_url, filename,
                                                      ranges_copy, header, type_tag);
    if (response.has_value()) return response.value();

    siril_log_warning(_("Query failed on mirror %s, trying alternatives...\n"), *current_mirror_ptr);

    for (int i = 0; mirrors[i] != NULL; i++) {
        if (g_strcmp0(mirrors[i], *current_mirror_ptr) == 0) continue;

        siril_log_message(_("Trying mirror: %s\n"), mirrors[i]);
        std::string test_base_url(mirrors[i]);
        ranges_copy = healpixel_ranges;
        auto resp = query_catalog_http_with_curl<EntryType>(curl, test_base_url, filename,
                                                            ranges_copy, header, type_tag);
        if (resp.has_value()) {
            g_free(*current_mirror_ptr);
            *current_mirror_ptr = g_strdup(mirrors[i]);
            siril_log_info(_("Switched to working mirror: %s\n"), *current_mirror_ptr);
            return resp.value();
        }
    }

    siril_log_error(_("All catalogue mirrors failed for query.\n"));
    return {};
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
        std::cout << "Catalogue Title: " << std::string_view(header.title, strnlen(header.title, sizeof(header.title))) << std::endl;
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

// Find the first chunk file in `path` whose <type> field matches `type_tag`.
// Pass type_tag = nullptr to accept any photometric type ("xpsamp" or "xpcts").
static std::string find_matching_cat_file(std::string& path, const char *type_tag = "xpsamp") {
    // siril_cat<chunk_level>_healpix<level>_<type>_<N>.dat per the format spec.
    std::string type_re = type_tag ? std::string(type_tag) : std::string("(?:xpsamp|xpcts)");
    std::string pattern = "siril_cat(\\d+)_healpix(\\d+)_" + type_re + "_(\\d+)\\.dat";
    std::regex file_regex(pattern);

    try {
        for (const auto& entry : std::filesystem::directory_iterator(path)) {
            if (!entry.is_regular_file()) continue;
            std::string filename = entry.path().filename().string();
            std::smatch matches;
            if (std::regex_match(filename, matches, file_regex)) return filename;
        }
    } catch (const std::filesystem::filesystem_error& e) {
        g_warning("Error accessing directory: %s", e.what());
        return "";
    }
    return "";
}

// This function is the main entry point and is declared extern "C" for ease of
// calling from the Siril C code. For integration into the code it will need to
// return an array of structs for use with SPCC, or possibly (maybe as a separate
// function) an array of deep_star structs for use with existing astrometry
// functions.

/* Templated local-loader body shared by the xpsamp and xpcts wrappers below.
 * Kept at namespace scope (templates can't have C linkage). EntryType must
 * be a packed struct with int32_t ra_scaled/dec_scaled and int16_t mag_scaled
 * at the start (true for both SourceEntryXPsamp and SourceEntryXPcts).
 * type_tag is the filename token ("xpsamp" or "xpcts"); expected_cat_type
 * is the cat_type byte we require in the on-disk header (2 = xp_sampled,
 * 3 = xp_continuous), per the Siril HEALpix Catalog Format 1.0.0 spec. */
template<typename EntryType>
static int local_gaia_xp_query(double ra, double dec, double radius, double limitmag,
                               const char *type_tag, uint8_t expected_cat_type,
                               EntryType **stars, uint32_t *nb_stars) {
    radius /= 60.0; // arcmin -> deg, then radians below
    siril_log_debug("Search radius: %f deg\n", radius);
    const double DEG_TO_RAD = G_PI / 180.0;
    double radius_rad = radius * DEG_TO_RAD;
    double ra_rad = ra * DEG_TO_RAD;
    double dec_rad = dec * DEG_TO_RAD;
    double theta = G_PI / 2.0 - dec_rad;
    double phi = ra_rad;
    double radius_h = pow(sin(0.5 * radius_rad), 2);

    std::string chunkpath(com.pref.catalogue_paths[5]);
    std::string first_chunk = find_matching_cat_file(chunkpath, type_tag);
    if (first_chunk.empty()) {
        siril_log_error(_("Error: no %s chunks detected.\n"), type_tag);
        *stars = nullptr; *nb_stars = 0;
        return 1;
    }

    int status = 0;
    std::string final_path = (std::filesystem::path(chunkpath) / first_chunk).string();
    HealpixCatHeader header = read_healpix_cat_header(final_path, &status);
    if (status) {
        *stars = nullptr; *nb_stars = 0;
        return 1;
    }
    if (header.cat_type != expected_cat_type) {
        siril_log_error(_("Error: catalogue at %s has cat_type=%d, expected %d for %s.\n"),
                        chunkpath.c_str(), header.cat_type, expected_cat_type, type_tag);
        *stars = nullptr; *nb_stars = 0;
        return 1;
    }

    T_Healpix_Base<int> healpix_base(header.healpix_level, NEST);
    pointing point(theta, phi);
    std::vector<int> pixel_indices;
    healpix_base.query_disc_inclusive(point, radius_rad, pixel_indices);

    std::set<int> chunks;
    for (const int& healpix : pixel_indices)
        chunks.insert(convert_healpix_level(healpix, header.healpix_level, header.chunk_level));
    std::vector<int> chunk_indices(chunks.begin(), chunks.end());

    std::vector<std::vector<int>> chunked_pixel_indices =
        group_pixels_by_chunk(chunk_indices, pixel_indices, header.chunk_level, header.healpix_level);

    std::vector<std::vector<EntryType>> results_in_chunks(chunk_indices.size());

    bool file_error = false;
    for (size_t i = 0; i < chunk_indices.size(); ++i) {
        const int chunk_id = chunk_indices[i];
        const std::vector<int>& chunk_pixels = chunked_pixel_indices[i];
        std::vector<HealPixelRange> healpixel_ranges = create_healpixel_ranges(chunk_pixels);

        gchar *filename = g_strdup_printf("siril_cat%u_healpix%u_%s_%d.dat",
                                          header.chunk_level, header.healpix_level, type_tag, chunk_id);
        std::string chunkfile(filename);
        g_free(filename);
        std::string this_chunk_path = (std::filesystem::path(chunkpath) / chunkfile).string();

        if (!std::filesystem::exists(this_chunk_path)) {
            siril_log_error(_("Chunk file not found: %s\n"), this_chunk_path.c_str());
            file_error = true;
            break;
        }
        int chunk_status = 0;
        HealpixCatHeader this_header = read_healpix_cat_header(this_chunk_path, &chunk_status);
        if (!header_compatible(header, this_header)) {
            siril_log_error(_("Error: catalog header values for chunk %lu are incompatible with previous values. All chunk files must have the same chunk level and indexing level and must represent the same Gaia data release\n"), (unsigned long)i);
            return 1;
        }
        results_in_chunks[i] = query_catalog<EntryType>(this_chunk_path, healpixel_ranges, this_header);
    }
    if (file_error) return 1;

    std::vector<EntryType> matches = flatten(results_in_chunks);
    if (matches.empty()) { *nb_stars = 0; *stars = nullptr; return 1; }

    double scaled_limitmag = limitmag * 1000.0;
    matches.erase(
        std::remove_if(matches.begin(), matches.end(),
                       [scaled_limitmag](const EntryType& entry) {
                           return entry.mag_scaled > scaled_limitmag;
                       }),
        matches.end());

    double scale = 360.0 / (double) INT32_MAX;
    matches.erase(
        std::remove_if(matches.begin(), matches.end(),
                       [scale, radius_h, ra, dec](const EntryType& entry) {
                           return compute_coords_distance_h(ra, dec,
                                                            (double)entry.ra_scaled * scale,
                                                            (double)entry.dec_scaled * scale) > radius_h;
                       }),
        matches.end());

    if (matches.empty()) { *nb_stars = 0; *stars = nullptr; return 1; }

    *nb_stars = matches.size();
    *stars = (EntryType*) malloc(*nb_stars * sizeof(EntryType));
    if (*stars == nullptr) return -1;
    std::copy(matches.begin(), matches.end(), *stars);
    return 0;
}

extern "C" {
    int local_gaia_xpsamp_available() {
        // Back-compat: report 'available' if the photometric dir contains
        // either xpsamp OR xpcts chunks (callers route via the loader, which
        // converts xpcts to xp_sampled at load time).
        std::string chunkpath(com.pref.catalogue_paths[5]);
        std::string first_chunk = find_matching_cat_file(chunkpath, nullptr);
        return (!first_chunk.empty());
    }

    int local_gaia_xpcts_available() {
        std::string chunkpath(com.pref.catalogue_paths[5]);
        std::string first_chunk = find_matching_cat_file(chunkpath, "xpcts");
        return (!first_chunk.empty());
    }

    /* Walk the photometric catalogue directory and report which kind of
     * catalogue it holds. Filename pattern narrows candidates; the cat_type
     * byte in the first chunk's header confirms (per the Siril HEALpix
     * Catalog Format 1.0.0 spec). */
    local_gaia_photo_kind detect_local_gaia_photo_kind(const char *dir) {
        if (!dir || !*dir) return LOCAL_GAIA_PHOTO_NONE;
        std::string chunkpath(dir);
        std::string xpsamp_chunk = find_matching_cat_file(chunkpath, "xpsamp");
        std::string xpcts_chunk  = find_matching_cat_file(chunkpath, "xpcts");
        if (xpsamp_chunk.empty() && xpcts_chunk.empty())
            return LOCAL_GAIA_PHOTO_NONE;
        if (!xpsamp_chunk.empty() && !xpcts_chunk.empty())
            return LOCAL_GAIA_PHOTO_MIXED;

        std::string chosen = xpsamp_chunk.empty() ? xpcts_chunk : xpsamp_chunk;
        std::string final_path = (std::filesystem::path(chunkpath) / chosen).string();
        int status = 0;
        HealpixCatHeader header = read_healpix_cat_header(final_path, &status);
        if (status) return LOCAL_GAIA_PHOTO_BAD;

        if (header.cat_type == 2) return LOCAL_GAIA_PHOTO_XPSAMP;
        if (header.cat_type == 3) return LOCAL_GAIA_PHOTO_XPCTS;
        return LOCAL_GAIA_PHOTO_BAD;
    }

    int get_raw_stars_from_local_gaia_astro_catalogue(double ra, double dec, double radius, double limitmag, gboolean phot, deepStarData **stars, uint32_t *nb_stars) {

#ifdef HEALPIX_DEBUG
        show_healpixel_entries(0);
#endif

        radius /= 60.0; // the catalogue radius is in arcmin, we want it in degrees to convert to radians
        siril_log_debug("Search radius: %f deg\n", radius);
        const double DEG_TO_RAD = G_PI / 180.0;
        double radius_rad = radius * DEG_TO_RAD;
        double ra_rad = ra * DEG_TO_RAD;
        double dec_rad = dec * DEG_TO_RAD;
        double theta = G_PI / 2.0 - dec_rad;
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
        siril_log_debug("query_disc_inclusive params: theta = %f, phi = %f, rad = %f\n", point.theta, point.phi, radius_rad);

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
        return local_gaia_xp_query<SourceEntryXPsamp>(ra, dec, radius, limitmag, "xpsamp", 2, stars, nb_stars);
    }

    int get_raw_stars_from_local_gaia_xpcontinuous_catalogue(double ra, double dec, double radius, double limitmag, SourceEntryXPcts **stars, uint32_t *nb_stars) {
        return local_gaia_xp_query<SourceEntryXPcts>(ra, dec, radius, limitmag, "xpcts", 3, stars, nb_stars);
    }

#ifdef HAVE_LIBCURL
}  /* close extern "C" briefly so the template below can have C++ linkage */

/* Templated outer remote loader: fetches XP records by cone from a chunked
 * remote catalogue. Parameterised on:
 *   EntryType         - SourceEntryXPsamp or SourceEntryXPcts.
 *   type_tag          - "xpsamp" / "xpcts" - filename token AND user-facing
 *                       name in error messages.
 *   expected_cat_type - 2 (xpsamp) or 3 (xpcts) - cat_type byte we require.
 *   mirrors           - NULL-terminated mirror URL list for this kind.
 *   current_mirror_ptr- pointer to the gchar* holding the current preferred
 *                       mirror; updated in place when fallback succeeds.
 * Returns 0 on success, 1 on no-data or unrecoverable mirror failure, -1 on OOM. */
template<typename EntryType>
static int remote_gaia_xp_query(double ra, double dec, double radius, double limitmag,
                                const char *type_tag, uint8_t expected_cat_type,
                                gchar **mirrors, gchar **current_mirror_ptr,
                                EntryType **stars, uint32_t *nb_stars) {
    if (!mirrors || !mirrors[0] || !current_mirror_ptr || !*current_mirror_ptr) {
        siril_log_error(_("No remote %s mirrors are configured.\n"), type_tag);
        *stars = nullptr; *nb_stars = 0;
        return 1;
    }

    radius /= 60.0;
    siril_log_debug("Search radius: %f deg\n", radius);
    const double DEG_TO_RAD = G_PI / 180.0;
    double radius_rad = radius * DEG_TO_RAD;
    double ra_rad = ra * DEG_TO_RAD;
    double dec_rad = dec * DEG_TO_RAD;
    double theta = G_PI / 2.0 - dec_rad;
    double phi = ra_rad;
    double radius_h = pow(sin(0.5 * radius_rad), 2);

    CURL* curl = curl_easy_init();
    if (!curl) {
        siril_log_error(_("Error initialising CURL handle\n"));
        *stars = nullptr; *nb_stars = 0;
        return 1;
    }

    /* The first chunk is conventionally chunk_level=1 / healpix_level=8 / id=0.
     * We re-read the header once more per chunk below to pick up per-chunk values. */
    gchar *first_filename_g = g_strdup_printf("siril_cat1_healpix8_%s_0.dat", type_tag);
    std::string first_chunk_filename(first_filename_g);
    g_free(first_filename_g);

    int status = 0;
    HealpixCatHeader header;
    if (!read_healpix_cat_header_http_with_fallback(curl, first_chunk_filename.c_str(),
                                                    &header, &status,
                                                    mirrors, current_mirror_ptr)) {
        curl_easy_cleanup(curl);
        *stars = nullptr; *nb_stars = 0;
        return 1;
    }
    if (header.cat_type != expected_cat_type) {
        siril_log_error(_("Remote catalogue at %s reports cat_type=%d, expected %d for %s.\n"),
                        *current_mirror_ptr, header.cat_type, expected_cat_type, type_tag);
        curl_easy_cleanup(curl);
        *stars = nullptr; *nb_stars = 0;
        return 1;
    }
    curl_easy_cleanup(curl);

    T_Healpix_Base<int> healpix_base(header.healpix_level, NEST);
    pointing point(theta, phi);

    std::vector<int> pixel_indices;
    healpix_base.query_disc_inclusive(point, radius_rad, pixel_indices);

    std::set<int> chunks;
    for (const int& healpix : pixel_indices)
        chunks.insert(convert_healpix_level(healpix, header.healpix_level, header.chunk_level));
    std::vector<int> chunk_indices(chunks.begin(), chunks.end());

    std::vector<std::vector<int>> chunked_pixel_indices =
        group_pixels_by_chunk(chunk_indices, pixel_indices, header.chunk_level, header.healpix_level);

    std::vector<std::vector<EntryType>> results_in_chunks(chunk_indices.size());

    bool file_error = false;
    #pragma omp parallel for shared(file_error) num_threads(4)
    for (size_t i = 0; i < chunk_indices.size(); ++i) {
        if (file_error) continue;

        const int chunk_id = chunk_indices[i];
        const std::vector<int>& chunk_pixels = chunked_pixel_indices[i];

        auto cache = FluxCache::getCache(header.chunk_level, header.healpix_level, chunk_id, type_tag);
        std::vector<int> remaining_pixels;
        std::vector<EntryType> results_for_this_chunk;

        int cache_count = 0;
        for (int hp : chunk_pixels) {
            std::vector<char> cached_blob;
            if (cache->getCacheEntry(static_cast<uint32_t>(hp), cached_blob)) {
                if (!cached_blob.empty()) {
                    const EntryType* entries = reinterpret_cast<const EntryType*>(cached_blob.data());
                    size_t count = cached_blob.size() / sizeof(EntryType);
                    results_for_this_chunk.insert(results_for_this_chunk.end(), entries, entries + count);
                }
                cache_count++;
            } else {
                remaining_pixels.push_back(hp);
            }
        }

#ifdef HAVE_SQLITE
        siril_log_message(_("Chunk %d: %d healpixels from cache, %zu to fetch from remote\n"),
                          chunk_id, cache_count, remaining_pixels.size());
#endif

        if (remaining_pixels.empty()) {
            results_in_chunks[i] = std::move(results_for_this_chunk);
            continue;
        }

        CURL* thread_curl = curl_easy_init();
        if (!thread_curl) {
            #pragma omp atomic write
            file_error = true;
            continue;
        }

        std::vector<HealPixelRange> healpixel_ranges = create_healpixel_ranges(remaining_pixels);

        gchar *filename = g_strdup_printf("siril_cat%u_healpix%u_%s_%d.dat",
                                         header.chunk_level, header.healpix_level, type_tag, chunk_id);
        std::string chunk_filename(filename);
        g_free(filename);

        std::string chunk_url = std::string(*current_mirror_ptr) + "/" + chunk_filename;
        int chunk_status = 0;
        HealpixCatHeader this_header = read_healpix_cat_header_http_with_curl(thread_curl, chunk_url, &chunk_status);

        if (chunk_status != 0) {
            if (!read_healpix_cat_header_http_with_fallback(thread_curl, chunk_filename.c_str(),
                                                            &this_header, &chunk_status,
                                                            mirrors, current_mirror_ptr)) {
                siril_log_error(_("Failed to read header for chunk %d from any mirror\n"), chunk_id);
                #pragma omp atomic write
                file_error = true;
                curl_easy_cleanup(thread_curl);
                continue;
            }
        }

        if (!header_compatible(header, this_header)) {
            siril_log_error(_("Error: catalog header values for chunk %d are incompatible\n"), chunk_id);
            #pragma omp atomic write
            file_error = true;
            curl_easy_cleanup(thread_curl);
            continue;
        }

        std::vector<EntryType> remote_results = query_catalog_http_with_fallback<EntryType>(
            thread_curl, chunk_filename, healpixel_ranges, this_header,
            mirrors, current_mirror_ptr, type_tag);

        results_for_this_chunk.insert(results_for_this_chunk.end(),
                                      remote_results.begin(), remote_results.end());
        results_in_chunks[i] = std::move(results_for_this_chunk);

        if (results_in_chunks[i].empty() && !chunk_pixels.empty()) {
            siril_log_warning(_("Warning: No data retrieved for chunk %d\n"), chunk_id);
        }

        curl_easy_cleanup(thread_curl);
    }

    if (file_error) {
        *stars = nullptr; *nb_stars = 0;
        return 1;
    }

    std::vector<EntryType> matches = flatten(results_in_chunks);
    if (matches.empty()) { *nb_stars = 0; *stars = nullptr; return 1; }

    double scaled_limitmag = limitmag * 1000.0;
    matches.erase(
        std::remove_if(matches.begin(), matches.end(),
                      [scaled_limitmag](const EntryType& entry) {
                          return entry.mag_scaled > scaled_limitmag;
                      }),
        matches.end());

    double scale = 360.0 / (double) INT32_MAX;
    matches.erase(
        std::remove_if(matches.begin(), matches.end(),
                      [scale, radius_h, ra, dec](const EntryType& entry) {
                          return compute_coords_distance_h(ra, dec,
                              (double)entry.ra_scaled * scale,
                              (double)entry.dec_scaled * scale) > radius_h;
                      }),
        matches.end());

    if (matches.empty()) { *nb_stars = 0; *stars = nullptr; return 1; }

    *nb_stars = matches.size();
    *stars = (EntryType*) malloc(*nb_stars * sizeof(EntryType));
    if (*stars == nullptr) return -1;
    std::copy(matches.begin(), matches.end(), *stars);
    return 0;
}

extern "C" {
    int get_raw_stars_from_remote_gaia_xpsampled_catalogue(double ra, double dec, double radius,
                                                           double limitmag, SourceEntryXPsamp **stars,
                                                           uint32_t *nb_stars) {
        return remote_gaia_xp_query<SourceEntryXPsamp>(ra, dec, radius, limitmag,
                                                       "xpsamp", 2,
                                                       spcc_mirrors, &com.spcc_remote_catalogue,
                                                       stars, nb_stars);
    }
#else
    int get_raw_stars_from_remote_gaia_xpsampled_catalogue(double ra, double dec, double radius,
                                                           double limitmag, SourceEntryXPsamp **stars,
                                                           uint32_t *nb_stars) {
        siril_log_error(_("Error: Siril was compiled without libcurl support. Remote catalogue access is unavailable.\n"));
        *stars = nullptr;
        *nb_stars = 0;
        return 1;
    }
#endif

    /* Remote xp_continuous loader. Becomes operational the moment a non-empty
     * spcc_mirrors_xpcts list is provided to initialize_spcc_mirrors() — until
     * then the templated body returns "no mirrors configured". */
#ifdef HAVE_LIBCURL
    int get_raw_stars_from_remote_gaia_xpcontinuous_catalogue(double ra, double dec, double radius,
                                                              double limitmag, SourceEntryXPcts **stars,
                                                              uint32_t *nb_stars) {
        return remote_gaia_xp_query<SourceEntryXPcts>(ra, dec, radius, limitmag,
                                                      "xpcts", 3,
                                                      spcc_mirrors_xpcts, &com.spcc_remote_catalogue_xpcts,
                                                      stars, nb_stars);
    }
#else
    int get_raw_stars_from_remote_gaia_xpcontinuous_catalogue(double ra, double dec, double radius,
                                                              double limitmag, SourceEntryXPcts **stars,
                                                              uint32_t *nb_stars) {
        (void)ra; (void)dec; (void)radius; (void)limitmag;
        siril_log_error(_("Error: Siril was compiled without libcurl support. Remote catalogue access is unavailable.\n"));
        *stars = nullptr; *nb_stars = 0;
        return 1;
    }
#endif
}
