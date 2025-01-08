// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#include <healpix_base.h>
#include <pointing.h>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <utility>
#include "core/siril.h"
#include "core/siril_log.h"
#include "io/local_catalogues.h"
#include "io/siril_catalogues.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

// Enum for Gaia version designator
enum class GaiaVersion {
    DR1 = 1,
    DR2 = 2,
    EDR3 = 3,
    DR3 = 4,
    DR4 = 5,
    DR5 = 6
};

// Enum for Catalogue Type
enum class CatalogueType {
    Astrometric = 1,
    Photometric_XP_Sampled = 2,
    Photometric_XP_Continuous = 3
};

// Struct for the Healpix catalogue header
typedef struct HealpixCatHeader {
    std::string title;          // 48 bytes: Catalogue title
    GaiaVersion gaia_version;   // 1 byte: Gaia version designator
    uint8_t healpix_level;      // 1 byte: Healpix indexing level N
    CatalogueType cat_type;     // 1 byte: Catalogue type
    std::array<uint8_t, 77> spare; // 77 bytes reserved
} HealpixCatHeader;

// Structure to represent a range of consecutive source IDs
struct HealPixelRange {
    uint64_t start_id;
    uint64_t end_id;
};

// This function creates a range of consecutive HEALpixels: we can scan these
// continuously in the file, rather than doing a binary search for the start of
// each HEALpixel in the range. It is used for both astrometric and photometric
// catalogues.

std::vector<HealPixelRange> create_healpixel_ranges(const std::vector<int>& pixels) {
    std::vector<HealPixelRange> ranges;
    if (pixels.empty()) {
        return ranges;
    }

    const uint64_t PIXEL_MULTIPLIER = 34359738368ULL; // 2^35
    HealPixelRange current_range{
        static_cast<uint64_t>(pixels[0]) * PIXEL_MULTIPLIER,
        (static_cast<uint64_t>(pixels[0]) + 1) * PIXEL_MULTIPLIER - 1
    };

    for (size_t i = 1; i < pixels.size(); ++i) {
        uint64_t next_start = static_cast<uint64_t>(pixels[i]) * PIXEL_MULTIPLIER;
        uint64_t next_end = next_start + PIXEL_MULTIPLIER - 1;

        if (next_start == current_range.end_id + 1) {
            // Extend current range
            current_range.end_id = next_end;
        } else {
            // Start new range
            ranges.push_back(current_range);
            current_range = HealPixelRange{next_start, next_end};
        }
    }
    // Add last range
    ranges.push_back(current_range);

    return ranges;
}

// This function queries a catalogue of packed EntryType objects (the template allows
// use of SourceEntryAstro or SourceEntryPhoto) for different uses. The catalogue MUST
// be sorted by source_id otherwise the method will fail. Each entry that matches the
// vector of source_id ranges is added to the result vector. Additional filtering (e.g.
// on magnitude) can be done by the caller.

template<typename EntryType>
std::vector<EntryType> query_catalog(const std::string& filename, const std::vector<HealPixelRange>& healpixel_ranges, uint8_t healpix_level) {
    constexpr size_t HEADER_SIZE = 128; // Fixed header size in bytes
    std::vector<EntryType> results;

    try {
        // Open the catalog file in binary mode
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            siril_log_message(_("Failed to open file: %s\n"), "red", filename.c_str());
            return results;
        }

        uint32_t nside = 1 << healpix_level; // Calculate NSIDE as 2^level
        uint32_t n_healpixels = 12 * nside * nside; // Total number of Healpixels

        // Function to read a single index entry at a specific position
        auto read_index_entry = [&file](uint32_t healpixel_id) -> uint32_t {
            uint32_t index_value;
            file.seekg(HEADER_SIZE + healpixel_id * sizeof(uint32_t), std::ios::beg);
            file.read(reinterpret_cast<char*>(&index_value), sizeof(uint32_t));
            if (!file) {
                throw std::runtime_error("Failed to read catalog index entry.");
            }
            return index_value;
        };

        // Process each range
        for (const auto& range : healpixel_ranges) {
            uint32_t start_healpixel = range.start_id;
            uint32_t end_healpixel = range.end_id;

            if (start_healpixel >= n_healpixels || end_healpixel >= n_healpixels) {
                siril_log_message(_("ID range exceeds catalog bounds."));
                return results;
            }

            // Read only the necessary index entries
            uint32_t start_offset = read_index_entry(start_healpixel);
            uint32_t end_offset = read_index_entry(end_healpixel + 1); // Read the next entry to get exclusive end

            // Calculate position in the data section
            size_t INDEX_SIZE = n_healpixels * sizeof(uint32_t);
            size_t data_start_pos = HEADER_SIZE + INDEX_SIZE + start_offset * sizeof(EntryType);

            // Read the required data entries
            size_t num_records = end_offset - start_offset;
            std::vector<EntryType> buffer(num_records);

            file.seekg(data_start_pos, std::ios::beg);
            file.read(reinterpret_cast<char*>(buffer.data()), num_records * sizeof(EntryType));
            if (!file) {
                throw std::runtime_error(_("Failed to read data entries."));
            }

            // Move entries to results vector
            results.insert(results.end(), buffer.begin(), buffer.end());
        }
    }
    catch (const std::exception& e) {
        siril_log_color_message(e.what(), "red");
        results.clear();  // Ensure vector is empty in case of error
    }

    return results;
}

// Function to read the Healpix catalogue header
HealpixCatHeader read_healpix_cat_header(const std::string& filename) {
    HealpixCatHeader header;

    // Open file in binary mode
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    // Read 48-byte title
    char title_buffer[48] = {0};
    file.read(title_buffer, 48);
    if (!file) {
        throw std::runtime_error("Error reading title from file.");
    }
    header.title = std::string(title_buffer);

    // Read 1-byte Gaia version
    uint8_t gaia_version_raw;
    file.read(reinterpret_cast<char*>(&gaia_version_raw), 1);
    if (!file) {
        throw std::runtime_error("Error reading Gaia version from file.");
    }
    header.gaia_version = static_cast<GaiaVersion>(gaia_version_raw);

    // Read 1-byte Healpix level
    file.read(reinterpret_cast<char*>(&header.healpix_level), 1);
    if (!file) {
        throw std::runtime_error("Error reading Healpix level from file.");
    }

    // Read 1-byte catalogue type
    uint8_t cat_type_raw;
    file.read(reinterpret_cast<char*>(&cat_type_raw), 1);
    if (!file) {
        throw std::runtime_error("Error reading catalogue type from file.");
    }
    header.cat_type = static_cast<CatalogueType>(cat_type_raw);

    // Read 77 spare bytes
    file.read(reinterpret_cast<char*>(header.spare.data()), 77);
    if (!file) {
        throw std::runtime_error("Error reading spare bytes from file.");
    }

    file.close();
    return header;
}

// This function is the main entry point and is declared extern "C" for ease of
// calling from the Siril C code. For integration into the code it will need to
// return an array of structs for use with SPCC, or possibly (maybe as a separate
// function) an array of deep_star structs for use with existing astrometry
// functions.

extern "C" {
    int get_raw_stars_from_local_gaia_astro_catalogue(double ra, double dec, double radius, double limitmag, deepStarData **stars, uint32_t *nb_stars) {
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
        HealpixCatHeader header = read_healpix_cat_header(filename);
        unsigned int nsides = 1 << header.healpix_level;
        T_Healpix_Base<int> healpix_base(nsides, NEST, SET_NSIDE);

        pointing point(theta, phi);
        std::vector<int> pixel_indices;
        siril_debug_print("query_disc_inclusive params: theta = %f, phi = %f, rad = %f\n", point.theta, point.phi, radius_rad);

        // Perform the cone search
        healpix_base.query_disc_inclusive(point, radius_rad, pixel_indices);

        // Reduce the pixels to a list of continuous ranges
        std::vector<HealPixelRange> healpixel_ranges = create_healpixel_ranges(pixel_indices);

        // Query the catalogue for matches
        std::vector<SourceEntryAstro> matches = query_catalog<SourceEntryAstro>(filename, healpixel_ranges, header.healpix_level);

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
        matches.erase(
            std::remove_if(matches.begin(), matches.end(),
                           [radius_h, ra, dec](const SourceEntryAstro& entry) {
                               return compute_coords_distance_h(ra, dec, (double)entry.ra_scaled * 0.000001, (double)entry.dec_scaled * .00001) > radius_h;
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

    int get_raw_stars_from_local_gaia_photo_catalogue(double ra, double dec, double radius, double limitmag, SourceEntryPhoto **stars, uint32_t *nb_stars) {
        radius /= 60.0; // the catalogue radius is in arcmin, we want it in degrees to convert to radians
        const double DEG_TO_RAD = M_PI / 180.0;
        double radius_rad = radius * DEG_TO_RAD;
        double ra_rad = ra * DEG_TO_RAD;
        double dec_rad = dec * DEG_TO_RAD;
        double theta = M_PI / 2.0 - dec_rad;
        double phi = ra_rad;

        // Check the correct healpixel level and create our healpix_base
        std::string filename(com.pref.catalogue_paths[5]);
        HealpixCatHeader header = read_healpix_cat_header(filename);
        unsigned int nsides = 1 << header.healpix_level;
        T_Healpix_Base<int> healpix_base(nsides, NEST, SET_NSIDE);

        pointing point(theta, phi);
        std::vector<int> pixel_indices;
        healpix_base.query_disc_inclusive(point, radius_rad, pixel_indices);

        std::vector<HealPixelRange> healpixel_ranges = create_healpixel_ranges(pixel_indices);

        // Query the catalog for matches
        std::vector<SourceEntryPhoto> matches = query_catalog<SourceEntryPhoto>(filename, healpixel_ranges, header.healpix_level);

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
                           [scaled_limitmag](const SourceEntryPhoto& entry) {
                               return entry.mag_scaled > scaled_limitmag;
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

        *nb_stars = matches.size();
        *stars = (SourceEntryPhoto*)malloc(matches.size() * sizeof(SourceEntryPhoto));
        if (*stars == nullptr) {
            return -1;  // Memory allocation failed, return -1
        }
        std::memcpy(*stars, matches.data(), matches.size() * sizeof(SourceEntryPhoto));
        return 0;
    }
}
