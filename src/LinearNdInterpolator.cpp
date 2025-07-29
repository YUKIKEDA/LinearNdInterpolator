#include "LinearNdInterpolator.h"
#include <stdexcept>
#include <cmath>
#include <limits>

LinearNdInterpolator::LinearNdInterpolator(
    const std::vector<std::vector<double>>& points, 
    const std::vector<double>& values
)
    : points_(points), values_(values) {
    
    // Input validation
    if (points.empty() || values.empty()) {
        throw std::invalid_argument("Points and values cannot be empty");
    }
    
    if (points.size() != values.size()) {
        throw std::invalid_argument("Number of points and values must match");
    }
    
    num_points_ = points.size();
    dimension_ = points[0].size();
    
    // Dimension consistency check
    for (const auto& point : points) {
        if (point.size() != dimension_) {
            throw std::invalid_argument("All points must have the same dimension");
        }
    }
    
    if (dimension_ == 0) {
        throw std::invalid_argument("Points must have at least one dimension");
    }
}

LinearNdInterpolator::~LinearNdInterpolator() {
    // No special cleanup needed currently
}

double LinearNdInterpolator::interpolate(const std::vector<double>& point) const {
    // Input point dimension check
    if (point.size() != dimension_) {
        throw std::invalid_argument("Query point dimension does not match interpolator dimension");
    }
    
    // Temporary implementation: return nearest neighbor value
    // In the future, implement Delaunay triangulation and barycentric coordinate interpolation
    
    if (num_points_ == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    size_t nearest_index = 0;
    double min_distance_sq = std::numeric_limits<double>::max();
    
    for (size_t i = 0; i < num_points_; ++i) {
        double distance_sq = 0.0;
        for (size_t j = 0; j < dimension_; ++j) {
            double diff = point[j] - points_[i][j];
            distance_sq += diff * diff;
        }
        
        if (distance_sq < min_distance_sq) {
            min_distance_sq = distance_sq;
            nearest_index = i;
        }
    }
    
    return values_[nearest_index];
}