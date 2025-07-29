#pragma once

#include <vector>
#include <memory>

// Forward declarations for Qhull
namespace orgQhull {
    class Qhull;
    class QhullFacet;
    class QhullPoint;
}

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
    // Data storage
    std::vector<std::vector<double>> points_;
    std::vector<double> values_;
    size_t dimension_;
    size_t num_points_;
    
    // Qhull triangulation
    std::unique_ptr<orgQhull::Qhull> qhull_;
    std::vector<double> flat_points_;  // Flattened point coordinates for Qhull
    
    // Private methods
    void setupTriangulation();
    double findNearestNeighborValue(const std::vector<double>& point) const;
    bool findContainingSimplex(const std::vector<double>& point, 
                              orgQhull::QhullFacet& facet) const;
    std::vector<double> calculateBarycentricCoordinates(
        const std::vector<double>& point, 
        const orgQhull::QhullFacet& facet) const;
    double interpolateInSimplex(const std::vector<double>& barycentricCoords,
                               const orgQhull::QhullFacet& facet) const;
    std::vector<double> solveLinearSystem(std::vector<std::vector<double>>& matrix, 
                                         std::vector<double>& rhs) const;
    size_t findPointIndex(const orgQhull::QhullPoint& vertex_point) const;
    bool isPointInSimplex(const std::vector<double>& point,
                         const orgQhull::QhullFacet& facet) const;
    std::vector<double> calculateBarycentricCoordinatesForTest(
        const std::vector<double>& point, 
        const orgQhull::QhullFacet& facet) const;
    std::vector<double> solveLinearSystemForTest(std::vector<std::vector<double>>& matrix, 
                                                std::vector<double>& rhs) const;
};