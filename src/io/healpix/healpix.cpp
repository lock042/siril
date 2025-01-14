// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#include <iomanip>
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

//#define HEALPIX_DEBUG

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
    uint8_t chunk_level;        // for chunked catalogues, the Healpix level it waas chunked at
    uint32_t chunk_healpix;     // for chunked catalogues, the chunk-level healpixel covered by this file
    uint32_t chunk_first_healpixel; // for chunked catalogues, the index-level healpixel number of the first healpixel in the catalogue
    uint32_t chunk_last_healpixel; // for chunked catalogues, the index-level healpixel number of the last healpixel in the catalogue
    std::array<uint8_t, 64> spare; // 64 bytes reserved for future use
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

std::vector<HealPixelRange> create_healpixel_ranges(const std::vector<int>& pixels) {
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
// use of SourceEntryAstro or SourceEntryPhoto) for different uses. The catalogue MUST
// be sorted by source_id otherwise the method will fail. Each entry that matches the
// vector of source_id ranges is added to the result vector. Additional filtering (e.g.
// on magnitude) can be done by the caller.

template<typename EntryType>
std::vector<EntryType> query_catalog(const std::string& filename, const std::vector<HealPixelRange>& healpixel_ranges, uint8_t healpix_level) {
    constexpr size_t HEADER_SIZE = 128; // Fixed header size in bytes
    std::vector<EntryType> results;
    uint32_t nside = 1 << healpix_level; // Calculate NSIDE as 2^level
    uint32_t n_healpixels = 12 * nside * nside; // Total number of Healpixels
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
            siril_log_color_message(_("Failed to read catalog index entry."), "red");
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
            siril_log_color_message(_("ID range exceeds catalog bounds."), "red");
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
            siril_log_color_message(_("Failed to read data entries."), "red");
            results.clear();
            return results;
        }

        // Move entries to results vector
        results.insert(results.end(), buffer.begin(), buffer.end());
    }

    return results;
}

// Error codes
enum HealpixHeaderReadError {
    HEALPIX_SUCCESS = 0,
    HEALPIX_FILE_OPEN_ERROR = -1,
    HEALPIX_READ_TITLE_ERROR = -2,
    HEALPIX_READ_VERSION_ERROR = -3,
    HEALPIX_READ_LEVEL_ERROR = -4,
    HEALPIX_READ_TYPE_ERROR = -5,
    HEALPIX_READ_SPARE_ERROR = -6
};

// Function to read the Healpix catalogue header
HealpixCatHeader read_healpix_cat_header(const std::string& filename, int* error_status) {
    HealpixCatHeader header;
    *error_status = HEALPIX_SUCCESS;  // Initialize to success

    // Open file in binary mode
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        *error_status = HEALPIX_FILE_OPEN_ERROR;
        return header;
    }

    // Read 48-byte title
    char title_buffer[48] = {0};
    file.read(title_buffer, 48);
    if (!file) {
        *error_status = HEALPIX_READ_TITLE_ERROR;
        file.close();
        return header;
    }
    header.title = std::string(title_buffer);

    // Read 1-byte Gaia version
    uint8_t gaia_version_raw;
    file.read(reinterpret_cast<char*>(&gaia_version_raw), 1);
    if (!file) {
        *error_status = HEALPIX_READ_VERSION_ERROR;
        file.close();
        return header;
    }
    header.gaia_version = static_cast<GaiaVersion>(gaia_version_raw);

    // Read 1-byte Healpix level
    file.read(reinterpret_cast<char*>(&header.healpix_level), 1);
    if (!file) {
        *error_status = HEALPIX_READ_LEVEL_ERROR;
        file.close();
        return header;
    }

    // Read 1-byte catalogue type
    uint8_t cat_type_raw;
    file.read(reinterpret_cast<char*>(&cat_type_raw), 1);
    if (!file) {
        *error_status = HEALPIX_READ_TYPE_ERROR;
        file.close();
        return header;
    }
    header.cat_type = static_cast<CatalogueType>(cat_type_raw);

    // Read 77 spare bytes
    file.read(reinterpret_cast<char*>(header.spare.data()), 64);
    if (!file) {
        *error_status = HEALPIX_READ_SPARE_ERROR;
        file.close();
        return header;
    }

    file.close();
    return header;
}
#ifdef HEALPIX_DEBUG
// Function to show entries for a specific healpixel
void show_healpixel_entries(uint32_t healpixel_id) {
    try {
        std::string filename(com.pref.catalogue_paths[4]);
        // Create a single-element range for the requested healpixel
        std::vector<HealPixelRange> range = {{healpixel_id, healpixel_id}};

        // Read the header first to get the healpix level
        HealpixCatHeader header = read_healpix_cat_header(filename);

        // Query the catalog for this healpixel
        auto entries = query_catalog<SourceEntryAstro>(filename, range, header.healpix_level);

        // Print each entry with the same formatting as Python
        for (const auto& entry : entries) {
            double ra = entry.ra_scaled / 1000000.0;
            double dec = entry.dec_scaled / 100000.0;
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

// This function is the main entry point and is declared extern "C" for ease of
// calling from the Siril C code. For integration into the code it will need to
// return an array of structs for use with SPCC, or possibly (maybe as a separate
// function) an array of deep_star structs for use with existing astrometry
// functions.

extern "C" {
    int get_raw_stars_from_local_gaia_astro_catalogue(double ra, double dec, double radius, double limitmag, gboolean phot, deepStarData **stars, uint32_t *nb_stars) {
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

}
