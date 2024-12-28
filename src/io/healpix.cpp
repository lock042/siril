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
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <utility>

// Source data structure with packed attributes to match file format
#pragma pack(push, 1)  // Ensure no padding between members
struct SourceEntry {
    uint64_t source_id;  // 8 bytes
    int32_t ra_scaled;   // 4 bytes
    int32_t dec_scaled;  // 4 bytes
    int16_t mag_scaled; // 2 bytes
    uint8_t fexpo;      // 1 byte
    int16_t flux[343];  // 686 bytes: xp_sampled flux values
};
#pragma pack(pop)

// Structure to represent a range of consecutive source IDs
struct SourceIdRange {
    uint64_t start_id;
    uint64_t end_id;
};

// This function creates a range of consecutive HEALpixels: we can scan these
// continuously in the file, rather than doing a binary search for the start of
// each HEALpixel in the range

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

// This function finds the first record that could be in the range

size_t find_first_record(std::ifstream& infile, uint64_t target_id, size_t file_size) {
    static constexpr size_t RECORD_SIZE = sizeof(SourceEntry);
    size_t record_count = file_size / RECORD_SIZE;
    size_t left = 0;
    size_t right = record_count;

//    std::cerr << "Binary search for target_id: " << target_id << std::endl;
//    std::cerr << "File has " << record_count << " records" << std::endl;

    while (left < right) {
        size_t mid = left + (right - left) / 2;

        SourceEntry entry;
        infile.seekg(mid * RECORD_SIZE);
        infile.read(reinterpret_cast<char*>(&entry), RECORD_SIZE);

//        std::cerr << "  At position " << mid << " found source_id: " << entry.source_id
//        << " (left=" << left << ", right=" << right << ")" << std::endl;

        if (entry.source_id < target_id) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
//    std::cerr << "Binary search returned position: " << left << std::endl;

    // Verify the position we found
    if (left < record_count) {
        SourceEntry entry;
        infile.seekg(left * RECORD_SIZE);
        infile.read(reinterpret_cast<char*>(&entry), RECORD_SIZE);
//        std::cerr << "First record at returned position has source_id: " << entry.source_id << std::endl;
    }

    return left;
}

// This function queries a catalogue of packed SourceEntry objects. The catalogue MUST
// be sorted by source_id otherwise the method will fail. Each entry that matches the
// vector of source_id ranges is added to the result vector. Additional filtering (e.g.
// on magnitude) would be easy to add, or could be dealt with by the caller, when
// integrated into Siril.

std::vector<SourceEntry> query_catalog(const std::string& filename, const std::vector<SourceIdRange>& id_ranges) {
    std::vector<SourceEntry> results;
    static constexpr size_t RECORD_SIZE = sizeof(SourceEntry);

    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    // Get file size
    infile.seekg(0, std::ios::end);
    size_t file_size = infile.tellg();
    infile.seekg(0, std::ios::beg);

    std::cerr << "\nFile size: " << file_size << " bytes" << std::endl;
    std::cerr << "Record size: " << RECORD_SIZE << " bytes" << std::endl;
    std::cerr << "Total records: " << (file_size / RECORD_SIZE) << std::endl;
    // Read and print first few records for verification
/*    std::cerr << "\nFirst few records in file:" << std::endl;
    for (int i = 0; i < 5; i++) {
        SourceEntry entry;
        infile.seekg(i * RECORD_SIZE);
        infile.read(reinterpret_cast<char*>(&entry), RECORD_SIZE);
        int healpix = static_cast<int>(entry.source_id / 34359738368ULL);
        std::cerr << "Record " << i << ": healpix=" << healpix
        << " source_id=" << entry.source_id
        << " ra=" << (entry.ra_scaled/1000000.0)
        << " dec=" << (entry.dec_scaled/100000.0)
        << std::endl;
    }*/

    // Buffer for reading records
    static constexpr size_t BUFFER_SIZE = 512;
    std::vector<SourceEntry> buffer(BUFFER_SIZE);

    // Process each source_id range
    for (const auto& range : id_ranges) {
//        std::cerr << "\nProcessing range: " << range.start_id << " to " << range.end_id << std::endl;

        // Find starting position
        size_t current_pos = find_first_record(infile, range.start_id, file_size);

        if (current_pos * RECORD_SIZE >= file_size) {
//            std::cerr << "Starting position " << current_pos << " is beyond file size" << std::endl;
            continue;
        }

        // Read and process records until we pass the end of the range
        bool found_any_in_range = false;
        size_t records_checked = 0;

        while (current_pos * RECORD_SIZE < file_size) {
            size_t records_to_read = std::min(BUFFER_SIZE,
                        (file_size - current_pos * RECORD_SIZE) / RECORD_SIZE);

            infile.seekg(current_pos * RECORD_SIZE);
            infile.read(reinterpret_cast<char*>(buffer.data()),
                        records_to_read * RECORD_SIZE);

//            std::cerr << "Read " << records_to_read << " records starting at position "
//            << current_pos << std::endl;


            bool range_complete = false;
            for (size_t i = 0; i < records_to_read; ++i) {
                records_checked++;

/*                if (records_checked <= 10 || buffer[i].source_id >= range.start_id) {
                    std::cerr << "Checking record: source_id=" << buffer[i].source_id;
                    std::cerr << " (range: " << range.start_id << " to " << range.end_id << ")";
                    std::cerr << " at file position " << (current_pos + i) << std::endl;
                }*/

                // If we've passed the end of the range, we're done
                if (buffer[i].source_id > range.end_id) {

//                    std::cerr << "Passed end of range at source_id: " << buffer[i].source_id << std::endl;
                    range_complete = true;

                    break;
                }

                // If we're in the range, add to results
                if (buffer[i].source_id >= range.start_id) {
                    found_any_in_range = true;
//                    std::cerr << "Found matching record: " << buffer[i].source_id << std::endl;

                    results.push_back(buffer[i]);
                }
            }

            if (range_complete) {
//                std::cerr << "Completed range after checking " << records_checked << " records" << std::endl;
                break;
            }

            current_pos += records_to_read;
        }

        if (!found_any_in_range) {
//            std::cerr << "No records found in range after checking " << records_checked << " records" << std::endl;
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
    void conesearch(double ra, double dec, double radius, const char* catalogue_file) {
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
        std::vector<SourceEntry> matches = query_catalog(filename, id_ranges);

        for (const auto& entry : matches) {
            std::cout << entry.source_id << ","
            << (entry.ra_scaled/1000000.0) << ","
            << (entry.dec_scaled/100000.0) << ","
            << (entry.mag_scaled / 1000.0) << ","
            << (static_cast<int>(entry.fexpo)) << "\n";
        }
    }
}

// main() function for testng the code base from the commandline prior to
// integrating it into Siril
/*
int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: conesearch {ra} {dec} {radius}" << std::endl;
        return 1;
    }

    double ra = std::stod(argv[1]);
    double dec = std::stod(argv[2]);
    double radius = std::stod(argv[3]);

    const char* catalog_file = "sourceid-data-m15.bin";

    try {
        conesearch(ra, dec, radius, catalog_file);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
*/
