// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#include "core/siril_log.h"
#include "core/siril.h"
#include "core/siril_networking.h"
#include "io/local_catalogues.h"
#include "io/siril_catalogues.h"
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
#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

#define HEALPIX_DEBUG

#define ZENODO_GAIA_XPSAMP_RECORD_ID "YOUR_RECORD_ID_HERE"

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
        siril_log_color_message(_("Failed to open file: %s\n"), "red", filename.c_str());
        return results;
    }

    // Function to read a single index entry at a specific position
    auto read_index_entry = [&file, &results](uint32_t healpixel_id) -> uint32_t {
        uint32_t index_value;
        size_t pos = HEADER_SIZE + healpixel_id * sizeof(uint32_t);
        file.seekg(pos, std::ios::beg);
        file.read(reinterpret_cast<char*>(&index_value), sizeof(uint32_t));
        if (!file) {
            siril_log_color_message(_("Failed to read catalog index entry.\n"), "red");
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
            siril_log_color_message(_("ID range exceeds catalog bounds.\n"), "red");
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

    gsize response_length;
    int error;
    char* buffer = fetch_url_range_with_curl(curl, url.c_str(), 0, 128, &response_length, &error, FALSE);

    if (error || !buffer || response_length < 128) {
        if (error_status) {
            *error_status = -1;
        }
        free(buffer);
        return header;
    }

    try {
        size_t offset = 0;

        // Read title (48 bytes)
        char title_buffer[48] = {0};
        memcpy(title_buffer, buffer + offset, 48);
        header.title = std::string(title_buffer, strnlen(title_buffer, 48));
        offset += 48;

        // Read remaining fields
        memcpy(&header.gaia_version, buffer + offset, 1); offset += 1;
        memcpy(&header.healpix_level, buffer + offset, 1); offset += 1;
        memcpy(&header.cat_type, buffer + offset, 1); offset += 1;
        memcpy(&header.chunked, buffer + offset, 1); offset += 1;
        memcpy(&header.chunk_level, buffer + offset, 1); offset += 1;
        memcpy(&header.chunk_healpix, buffer + offset, sizeof(uint32_t)); offset += sizeof(uint32_t);
        memcpy(&header.chunk_first_healpixel, buffer + offset, sizeof(uint32_t)); offset += sizeof(uint32_t);
        memcpy(&header.chunk_last_healpixel, buffer + offset, sizeof(uint32_t)); offset += sizeof(uint32_t);
        memcpy(&header.spare, buffer + offset, 63);

        free(buffer);
    }
    catch (const std::exception&) {
        if (error_status) {
            *error_status = -3;
        }
        free(buffer);
        return header;
    }

    return header;
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

    // Function to read index entries via HTTP range request using provided CURL handle
    auto read_index_entries = [curl, &full_url](uint32_t start_healpixel, uint32_t end_healpixel)
        -> std::pair<bool, std::vector<uint32_t>> {
        size_t start_pos = HEADER_SIZE + start_healpixel * sizeof(uint32_t);
        size_t length = (end_healpixel - start_healpixel + 1) * sizeof(uint32_t);

        gsize response_length;
        int error;
        char* buffer = fetch_url_range_with_curl(curl, full_url.c_str(), start_pos, length,
                                      &response_length, &error, FALSE);

        if (error || !buffer) {
            return {false, {}};
        }

        std::vector<uint32_t> indices(response_length / sizeof(uint32_t));
        memcpy(indices.data(), buffer, response_length);
        free(buffer);

        return {true, indices};
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

        // Read the index entries
        uint32_t index_start = (start_healpixel == 0) ? 0 : start_healpixel - 1;
        auto [success, indices] = read_index_entries(index_start, end_healpixel);

        if (!success) {
            siril_log_color_message(_("Failed to read index entries via HTTP\n"), "red");
            results.clear();
            return results;
        }

        uint32_t start_offset = (start_healpixel == 0) ? 0 : indices[0];
        uint32_t end_offset = indices[indices.size() - 1];

        // Calculate position in the data section
        size_t data_start_pos = HEADER_SIZE + INDEX_SIZE + start_offset * sizeof(EntryType);
        size_t num_records = end_offset - start_offset;
        size_t data_length = num_records * sizeof(EntryType);

        // Read the data entries via HTTP range request using provided CURL handle
        gsize data_response_length;
        int data_error;
        char* data_buffer = fetch_url_range_with_curl(curl, full_url.c_str(), data_start_pos, data_length,
                                           &data_response_length, &data_error, FALSE);

        if (data_error || !data_buffer) {
            siril_log_color_message(_("Failed to read data entries via HTTP\n"), "red");
            results.clear();
            return results;
        }

        // Convert buffer to EntryType vector
        std::vector<EntryType> buffer(num_records);
        memcpy(buffer.data(), data_buffer, data_response_length);
        free(data_buffer);

        // Move entries to results vector
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

#endif

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

    // Open file in binary mode
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        if (error_status) {
            *error_status = -1; // File open error
        }
        return header;
    }

    try {
        // Read the fixed-size string (48 bytes)
        char title_buffer[48] = {0};
        file.read(title_buffer, 48);
        if (file.fail()) {
            if (error_status) {
                *error_status = -2; // Read error
            }
            return header;
        }
        header.title = std::string(title_buffer, strnlen(title_buffer, 48));

        // Read the remaining POD members
        file.read(reinterpret_cast<char*>(&header.gaia_version), 1);
        file.read(reinterpret_cast<char*>(&header.healpix_level), 1);
        file.read(reinterpret_cast<char*>(&header.cat_type), 1);
        file.read(reinterpret_cast<char*>(&header.chunked), 1);
        file.read(reinterpret_cast<char*>(&header.chunk_level), 1);
        file.read(reinterpret_cast<char*>(&header.chunk_healpix), sizeof(uint32_t));
        file.read(reinterpret_cast<char*>(&header.chunk_first_healpixel), sizeof(uint32_t));
        file.read(reinterpret_cast<char*>(&header.chunk_last_healpixel), sizeof(uint32_t));
        file.read(reinterpret_cast<char*>(&header.spare), 63);

        if (file.fail()) {
            if (error_status) {
                *error_status = -2; // Read error
            }
            return header;
        }
    }
    catch (const std::exception&) {
        if (error_status) {
            *error_status = -3; // Exception during read
        }
        return header;
    }

    return header;
}


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

    try {
        // Iterate through directory entries
        for (const auto& entry : std::filesystem::directory_iterator(path)) {
            if (!entry.is_regular_file()) {
                continue;
            }

            std::string filename = entry.path().filename().string();
            std::smatch matches;

            // Check if filename matches our pattern
            if (std::regex_match(filename, matches, file_regex)) {
                return filename;
            }
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
        std::string final_path = (std::filesystem::path(chunkpath) / first_chunk).string();
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
            std::string chunkfile(filename);
            g_free(filename);
            std::string this_chunk_path = (std::filesystem::path(chunkpath) / chunkfile).string();

            if (!std::filesystem::exists(this_chunk_path)) {
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

        // Construct base URL
        std::string base_url(com.spcc_remote_catalogue);

        // Initialize CURL handle for all requests
        CURL* curl = curl_easy_init();
        if (!curl) {
            siril_log_color_message(_("Error initialising CURL handle\n"), "red");
            *stars = nullptr;
            *nb_stars = 0;
            return 1;
        }

        // Read header from first chunk to get catalog parameters
        std::string first_chunk_filename = "siril_cat1_healpix8_xpsamp_0.dat";
        std::string first_chunk_url = base_url + "/" + first_chunk_filename;

        int status = 0;
        HealpixCatHeader header = read_healpix_cat_header_http_with_curl(curl, first_chunk_url, &status);
        if (status) {
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

            std::string chunk_url = base_url + "/" + chunk_filename;

            int chunk_status = 0;
            HealpixCatHeader this_header = read_healpix_cat_header_http_with_curl(curl, chunk_url, &chunk_status);

            if (chunk_status) {
                siril_log_color_message(_("Failed to read header from: %s\n"), "red", chunk_url.c_str());
                file_error = true;
                break;
            }

            if (!header_compatible(header, this_header)) {
                siril_log_color_message(_("Error: catalog header values for chunk %d are incompatible\n"), "red", chunk_id);
                file_error = true;
                break;
            }

            results_in_chunks[i] = query_catalog_http_with_curl<SourceEntryXPsamp>(curl, base_url, chunk_filename,
                                                                         healpixel_ranges, this_header);
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
