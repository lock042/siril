// Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "core/siril.h"
#include "core/siril_log.h"

extern "C" {
#include "io/healpix/healpix_image.h"

// Forward-declare the two helpers we need from algos/siril_wcs.h. Pulling in
// that header would also drag in <wcslib.h>, which conflicts with healpix's
// own typedefs (both define `int64`).
gboolean has_wcs(fits *fit);
void pix2wcs(fits *fit, double pixel_x, double pixel_y, double *world_x, double *world_y);
}

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include <healpix_base.h>
#include <pointing.h>
#include <rangeset.h>

#ifndef G_PI
#define G_PI 3.14159265358979323846
#endif

namespace {

constexpr int LEVEL_LOW = 1;   // Nside = 2
constexpr int LEVEL_HIGH = 8;  // Nside = 256

// Convert RA/Dec in degrees to a healpix pointing (theta=colatitude, phi=longitude in radians).
pointing radec_to_pointing(double ra_deg, double dec_deg) {
    const double deg2rad = G_PI / 180.0;
    double theta = (90.0 - dec_deg) * deg2rad;
    double phi = ra_deg * deg2rad;
    while (phi < 0.0) phi += 2.0 * G_PI;
    while (phi >= 2.0 * G_PI) phi -= 2.0 * G_PI;
    if (theta < 0.0) theta = 0.0;
    if (theta > G_PI) theta = G_PI;
    return pointing(theta, phi);
}

} // namespace

extern "C" int log_image_healpixels(fits *fit) {
    if (!fit) {
        siril_log_error(_("No image loaded\n"));
        return 1;
    }
    if (!has_wcs(fit)) {
        siril_log_error(_("This command only works on plate solved images\n"));
        return 1;
    }

    const double w = static_cast<double>(fit->rx);
    const double h = static_cast<double>(fit->ry);
    if (w <= 0.0 || h <= 0.0) {
        siril_log_error(_("Image has no usable dimensions\n"));
        return 1;
    }

    // Four corner pixel centres traversed counter-clockwise in pixel space.
    // pix2wcs applies SIP distortion if present.
    const double px[4] = { 0.0,      w - 1.0, w - 1.0, 0.0     };
    const double py[4] = { 0.0,      0.0,     h - 1.0, h - 1.0 };

    std::vector<pointing> vertices;
    vertices.reserve(4);
    for (int i = 0; i < 4; ++i) {
        double ra = 0.0, dec = 0.0;
        pix2wcs(fit, px[i], py[i], &ra, &dec);
        if (std::isnan(ra) || std::isnan(dec)) {
            siril_log_error(_("Failed to project image corner (%.0f, %.0f) to world coordinates\n"),
                            px[i], py[i]);
            return 1;
        }
        vertices.push_back(radec_to_pointing(ra, dec));
    }

    std::vector<int> low_pixels;
    std::vector<int> high_pixels;
    try {
        T_Healpix_Base<int> base_low(LEVEL_LOW, NEST);
        T_Healpix_Base<int> base_high(LEVEL_HIGH, NEST);

        rangeset<int> set_low, set_high;
        // fact=4 is the typical conservative choice for inclusive polygon queries
        // in the NESTED scheme (must be a power of two).
        base_low.query_polygon_inclusive(vertices, set_low, 4);
        base_high.query_polygon_inclusive(vertices, set_high, 4);
        set_low.toVector(low_pixels);
        set_high.toVector(high_pixels);
    } catch (const std::exception &e) {
        siril_log_error(_("HEALPix query failed: %s\n"), e.what());
        return 1;
    }

    // In NESTED ordering, the level-1 parent of a level-N pixel is obtained by
    // shifting right by 2 * (N - 1) bits.
    const int shift = 2 * (LEVEL_HIGH - LEVEL_LOW);
    std::map<int, std::vector<int>> grouped;
    for (int p : high_pixels)
        grouped[p >> shift].push_back(p);

    // Ensure every level-1 pixel from the level-1 query is shown, even if its
    // group of level-8 children is empty (rare, but possible for very small
    // footprints when the two inclusive queries disagree at the margins).
    for (int p : low_pixels)
        grouped.emplace(p, std::vector<int>());

    siril_log_message(_("HEALPix coverage for image footprint (NESTED ordering)\n"));
    siril_log_message(_("Level 1 (Nside=2): %u pixel(s); level 8 (Nside=256): %u pixel(s)\n"),
                      (unsigned) grouped.size(), (unsigned) high_pixels.size());

    for (const auto &kv : grouped) {
        int parent = kv.first;
        const std::vector<int> &children = kv.second;
        siril_log_message(_("  Level 1 pixel %d (%u level 8 children):\n"),
                          parent, (unsigned) children.size());
        const size_t per_line = 8;
        for (size_t i = 0; i < children.size(); i += per_line) {
            std::string line = "    ";
            for (size_t j = i; j < std::min(i + per_line, children.size()); ++j) {
                char buf[32];
                snprintf(buf, sizeof(buf), "%d", children[j]);
                if (j > i) line += ", ";
                line += buf;
            }
            line += "\n";
            siril_log_message("%s", line.c_str());
        }
    }

    return 0;
}
