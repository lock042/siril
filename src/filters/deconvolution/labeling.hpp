#pragma once

#include <algorithm>
#include <functional>
#include <map>

#include "image.hpp"

namespace labeling {

    // compute the connected components of an image
    // using 8-connected neighbors
template <typename T>
int labels(img_t<int>& labels, const img_t<T>& img) {
    labels.resize(img);
    labels.set_value(0);
    std::vector<int> equiv;
    equiv.push_back(0);
    int nblabels = 0;

    // Parallel region for the main loop
#ifdef _OPENMP
    #pragma omp parallel num_threads(com.fftw_max_thread)
    {
#endif
        // Each thread will have its own local nblabels
        int local_nblabels = 0;
        std::vector<int> local_equiv;
        local_equiv.push_back(0);

        // Parallelize the outermost loop
#ifdef _OPENMP
        #pragma omp for schedule(dynamic) collapse(2)
#endif
        for (int d = 0; d < img.d; d++) {
            for (int y = 0; y < img.h; y++) {
                for (int x = 0; x < img.w; x++) {
                    T val = img(x, y, d);
                    if (!val) {
                        continue;
                    }
                    int tl = y > 0 && x > 0 ? labels(x-1, y-1, d) : 0;
                    int t = y > 0 ? labels(x, y-1, d) : 0;
                    int tr = y > 0 && x < img.w-1 ? labels(x+1, y-1, d) : 0;
                    int l = x > 0 ? labels(x-1, y, d) : 0;

                    if (tl + t + tr + l == 0) {
                        local_nblabels++;
                        labels(x, y, d) = local_nblabels;
                        local_equiv.push_back(local_nblabels);
                        continue;
                    }

                    int max = std::max({tl, t, tr, l});
                    labels(x, y, d) = max;

                    if (tl && tl != max)
                        local_equiv[tl] = max;
                    if (t && t != max)
                        local_equiv[t] = max;
                    if (tr && tr != max)
                        local_equiv[tr] = max;
                    if (l && l != max)
                        local_equiv[l] = max;
                }
            }
        }

        // Combine local results
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            int offset = nblabels;
            nblabels += local_nblabels;

            // Adjust local labels and merge local_equiv into equiv
            for (int i = 1; i <= local_nblabels; i++) {
                equiv.push_back(local_equiv[i] + offset);
            }
        }
#ifdef _OPENMP
    }
#endif
    // Assign one label per connected component
    std::function<int(int)> getroot = [&](int l) {
        while (l != equiv[l]) {
            equiv[l] = equiv[equiv[l]];  // Path compression
            l = equiv[l];
        }
        return l;
    };

#ifdef _OPENMP
    #pragma omp parallel for num_threads(com.fftw_max_thread)
#endif
    for (int i = 0; i < labels.size; i++) {
        labels[i] = getroot(labels[i]);
    }

    return nblabels;
}


    template <typename T>
    std::map<int, T> sum(const img_t<int>& labels, const img_t<T>& img) {
        std::map<int, T> sum;
        for (int i = 0; i < img.size; i++) {
            sum[labels[i]] += img[i];
        }
        return sum;
    }

};


