// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#include <healpix_base.h>
#include <pointing.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <utility>
#include "io/local_catalogues.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

// Structure to represent a range of consecutive source IDs
struct SourceIdRange {
    uint64_t start_id;
    uint64_t end_id;
};

// This function creates a range of consecutive HEALpixels: we can scan these
// continuously in the file, rather than doing a binary search for the start of
// each HEALpixel in the range. It is used for both astrometric and photometric
// catalogues.

std::vector<SourceIdRange> create_sourceid_ranges(const std::vector<int>& pixels) {
    std::vector<SourceIdRange> ranges;
    if (pixels.empty()) {
        return ranges;
    }

    const uint64_t PIXEL_MULTIPLIER = 34359738368ULL; // 2^35
    SourceIdRange current_range{
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
            current_range = SourceIdRange{next_start, next_end};
        }
    }
    // Add last range
    ranges.push_back(current_range);

    return ranges;
}

// This function finds the first record that could be in the range.

template<typename EntryType>
size_t find_first_record(std::ifstream& infile, uint64_t target_id, size_t file_size) {
    static constexpr size_t RECORD_SIZE = sizeof(EntryType);
    size_t record_count = file_size / RECORD_SIZE;
    size_t left = 0;
    size_t right = record_count;

    while (left < right) {
        size_t mid = left + (right - left) / 2;
        EntryType entry;
        infile.seekg(mid * RECORD_SIZE);
        infile.read(reinterpret_cast<char*>(&entry), RECORD_SIZE);
        if (entry.source_id < target_id) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    if (left < record_count) {
        EntryType entry;
        infile.seekg(left * RECORD_SIZE);
        infile.read(reinterpret_cast<char*>(&entry), RECORD_SIZE);
    }
    return left;
}

// This function queries a catalogue of packed EntryType objects (the template allows
// use of SourceEntryAstro or SourceEntryPhoto) for different uses. The catalogue MUST
// be sorted by source_id otherwise the method will fail. Each entry that matches the
// vector of source_id ranges is added to the result vector. Additional filtering (e.g.
// on magnitude) can be done by the caller.

template<typename EntryType>
std::vector<EntryType> query_catalog(const std::string& filename, const std::vector<SourceIdRange>& id_ranges) {
    std::vector<EntryType> results;
    static constexpr size_t RECORD_SIZE = sizeof(EntryType);
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    infile.seekg(0, std::ios::end);
    size_t file_size = infile.tellg();
    infile.seekg(0, std::ios::beg);

    static constexpr size_t BUFFER_SIZE = 512;
    std::vector<EntryType> buffer(BUFFER_SIZE);

    for (const auto& range : id_ranges) {
        size_t current_pos = find_first_record<EntryType>(infile, range.start_id, file_size);
        if (current_pos * RECORD_SIZE >= file_size) {
            continue;
        }

        while (current_pos * RECORD_SIZE < file_size) {
            size_t records_to_read = std::min(BUFFER_SIZE,
                                              (file_size - current_pos * RECORD_SIZE) / RECORD_SIZE);
            infile.seekg(current_pos * RECORD_SIZE);
            infile.read(reinterpret_cast<char*>(buffer.data()),
                        records_to_read * RECORD_SIZE);

            for (size_t i = 0; i < records_to_read; ++i) {
                if (buffer[i].source_id > range.end_id) {
                    break;
                }
                if (buffer[i].source_id >= range.start_id) {
                    results.push_back(buffer[i]);
                }
            }
            current_pos += records_to_read;
        }
    }
    std::cerr << "\nTotal results found: " << results.size() << std::endl;
    return results;
}

// This function is the main entry point and is declared extern "C" for ease of
// calling from the Siril C code. For integration into the code it will need to
// return an array of structs for use with SPCC, or possibly (maybe as a separate
// function) an array of deep_star structs for use with existing astrometry
// functions.

extern "C" {
    int get_raw_stars_from_local_gaia_astro_catalogue(double ra, double dec, double radius, double limitmag, deepStarData **stars, uint32_t *nb_stars) {
        const char* catalogue_file = "sourceid-data-m15.bin"; // TODO: add this catalogue to com.pref and UI elements to set it
        const double DEG_TO_RAD = M_PI / 180.0;
        double radius_rad = radius * DEG_TO_RAD;
        T_Healpix_Base<int> healpix_base(4096, NEST, SET_NSIDE);

        double ra_rad = ra * DEG_TO_RAD;
        double dec_rad = dec * DEG_TO_RAD;
        double theta = M_PI / 2.0 - dec_rad;
        double phi = ra_rad;

        pointing point(theta, phi);
        std::vector<int> pixel_indices;
        healpix_base.query_disc_inclusive(point, radius_rad, pixel_indices);

        std::vector<SourceIdRange> id_ranges = create_sourceid_ranges(pixel_indices);

        // Convert char* to std::string for query_catalog
        std::string filename(catalogue_file);
        std::vector<SourceEntryAstro> matches = query_catalog<SourceEntryAstro>(filename, id_ranges);

        // Filter by magnitude
        double scaled_limitmag = limitmag * 1000.0;
        matches.erase(
            std::remove_if(matches.begin(), matches.end(),
                [scaled_limitmag](const SourceEntryAstro& entry) {
                    return entry.mag_scaled < scaled_limitmag;
                }
            ),
            matches.end()
        );
        *nb_stars = matches.size();
        uint32_t i = 0;
        *stars = (deepStarData*) malloc(*nb_stars * sizeof(deepStarData));
        if (*stars == nullptr) {
            return -1;  // Memory allocation failed
        }
        for (const auto& entry : matches) {
            (*stars)[i++] = (deepStarData) {
                .RA = entry.ra_scaled,
                .Dec = entry.dec_scaled,
                .dRA = entry.dra_scaled,
                .dDec = entry.ddec_scaled,
                .B = entry.mag_scaled,
                // deepStarData.V is not used with the Gaia catalogue
                .V = 0
            };
        }
        return 0;
    }

    int get_raw_stars_from_local_gaia_photo_catalogue(double ra, double dec, double radius, double limitmag, SourceEntryPhoto **stars, uint32_t *nb_stars) {
        const char* catalogue_file = "sourceid-data-m15.bin"; // TODO: add this catalogue to com.pref and UI elements to set it
        const double DEG_TO_RAD = M_PI / 180.0;
        double radius_rad = radius * DEG_TO_RAD;
        T_Healpix_Base<int> healpix_base(4096, NEST, SET_NSIDE);

        double ra_rad = ra * DEG_TO_RAD;
        double dec_rad = dec * DEG_TO_RAD;
        double theta = M_PI / 2.0 - dec_rad;
        double phi = ra_rad;

        pointing point(theta, phi);
        std::vector<int> pixel_indices;
        healpix_base.query_disc_inclusive(point, radius_rad, pixel_indices);

        std::vector<SourceIdRange> id_ranges = create_sourceid_ranges(pixel_indices);

        // Convert char* to std::string for query_catalog
        std::string filename(catalogue_file);
        std::vector<SourceEntryPhoto> matches = query_catalog<SourceEntryPhoto>(filename, id_ranges);

        // Filter by magnitude
        double scaled_limitmag = limitmag * 1000.0;
        matches.erase(
            std::remove_if(matches.begin(), matches.end(),
                           [scaled_limitmag](const SourceEntryPhoto& entry) {
                               return entry.mag_scaled < scaled_limitmag;
                           }
            ),
            matches.end()
        );
        *nb_stars = matches.size();
        *stars = (SourceEntryPhoto*)malloc(matches.size() * sizeof(SourceEntryPhoto));
        if (*stars == nullptr) {
            return -1;  // Memory allocation failed
        }
        std::memcpy(*stars, matches.data(), matches.size() * sizeof(SourceEntryPhoto));
        return 0;
    }
}
