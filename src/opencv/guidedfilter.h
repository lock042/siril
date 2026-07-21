#ifndef GUIDED_FILTER_H
#define GUIDED_FILTER_H

/* Only cv::Mat is needed here.  The opencv.hpp umbrella pulled in
 * stitching.hpp, whose exposure_compensate.hpp emits -Woverloaded-virtual
 * warnings on every include with current GCC/OpenCV. */
#include <opencv2/core.hpp>

class GuidedFilterImpl;

class GuidedFilter
{
public:
    GuidedFilter(const cv::Mat &I, int r, double eps);
    ~GuidedFilter();

    cv::Mat filter(const cv::Mat &p, int depth = -1) const;

private:
    GuidedFilterImpl *impl_;
};

cv::Mat guidedFilter(const cv::Mat &I, const cv::Mat &p, int r, double eps, int depth = -1);

#endif
