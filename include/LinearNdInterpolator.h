#pragma once

#include <vector>

/**
 * @brief N-dimensional linear interpolation class
 * 
 * Provides N-dimensional linear interpolation of scattered data based on Qhull library.
 * Uses Delaunay triangulation to create simplices and performs linear interpolation
 * using barycentric coordinates. Returns nearest neighbor values for points outside convex hull.
 */
class LinearNdInterpolator {
public:
    /**
     * @brief Constructor
     * @param points Interpolation point coordinates (list of N-dimensional points)
     * @param values Values corresponding to each point
     */
    LinearNdInterpolator(const std::vector<std::vector<double>>& points, 
                        const std::vector<double>& values);

    /**
     * @brief Destructor
     */
    ~LinearNdInterpolator();

    /**
     * @brief Calculate interpolated value at specified point
     * @param point Point coordinates where interpolation is performed
     * @return Interpolated value (nearest neighbor value if outside convex hull)
     */
    double interpolate(const std::vector<double>& point) const;

private:
    // Implementation will be added in the future
    std::vector<std::vector<double>> points_;
    std::vector<double> values_;
    size_t dimension_;
    size_t num_points_;
};