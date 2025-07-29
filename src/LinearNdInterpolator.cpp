#include "LinearNdInterpolator.h"
#include <stdexcept>
#include <cmath>
#include <limits>
#include <algorithm>

// Qhull includes
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullPoint.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/QhullError.h"

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
    
    // Need at least dimension + 1 points for triangulation
    if (num_points_ < dimension_ + 1) {
        throw std::invalid_argument("Need at least dimension + 1 points for triangulation");
    }
    
    // Setup Delaunay triangulation
    setupTriangulation();
}

LinearNdInterpolator::~LinearNdInterpolator() {
    // Smart pointer handles cleanup automatically
}

void LinearNdInterpolator::setupTriangulation() {
    try {
        // Flatten points for Qhull
        flat_points_.clear();
        flat_points_.reserve(num_points_ * dimension_);
        
        for (const auto& point : points_) {
            for (double coord : point) {
                flat_points_.push_back(coord);
            }
        }
        
        // Create Qhull object for Delaunay triangulation
        // "d" option creates Delaunay triangulation
        // "QJ" option joggles input to avoid precision problems
        std::string qhull_options = "d Qbb Qt QJ";
        
        qhull_ = std::make_unique<orgQhull::Qhull>(
            "LinearNdInterpolator",
            static_cast<int>(dimension_),
            static_cast<int>(num_points_),
            flat_points_.data(),
            qhull_options.c_str()
        );
        
    } catch (const orgQhull::QhullError& e) {
        throw std::runtime_error("Failed to create Delaunay triangulation: " + std::string(e.what()));
    }
}

double LinearNdInterpolator::interpolate(const std::vector<double>& point) const {
    // Input point dimension check
    if (point.size() != dimension_) {
        throw std::invalid_argument("Query point dimension does not match interpolator dimension");
    }
    
    if (num_points_ == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    // First check if point exactly matches any input point
    const double exact_match_tolerance = 1e-15;
    for (size_t i = 0; i < num_points_; ++i) {
        bool exact_match = true;
        for (size_t j = 0; j < dimension_; ++j) {
            if (std::abs(point[j] - points_[i][j]) > exact_match_tolerance) {
                exact_match = false;
                break;
            }
        }
        if (exact_match) {
            return values_[i];
        }
    }
    
    try {
        // Try to find containing simplex
        orgQhull::QhullFacet facet;
        if (findContainingSimplex(point, facet)) {
            // Calculate barycentric coordinates
            std::vector<double> barycentricCoords = calculateBarycentricCoordinates(point, facet);
            
            // Perform linear interpolation
            return interpolateInSimplex(barycentricCoords, facet);
        } else {
            // Point is outside convex hull, use nearest neighbor
            return findNearestNeighborValue(point);
        }
    } catch (const std::exception&) {
        // If anything goes wrong, fall back to nearest neighbor
        return findNearestNeighborValue(point);
    }
}

double LinearNdInterpolator::findNearestNeighborValue(const std::vector<double>& point) const {
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

bool LinearNdInterpolator::findContainingSimplex(const std::vector<double>& point, 
                                                orgQhull::QhullFacet& facet) const {
    if (!qhull_ || !qhull_->initialized()) {
        return false;
    }
    
    try {
        // For Delaunay triangulation, we need to check all facets 
        // to find one that contains the query point
        orgQhull::QhullFacetList facets = qhull_->facetList();
        
        for (orgQhull::QhullFacetList::iterator it = facets.begin(); it != facets.end(); ++it) {
            if (!it->isGood() || it->isUpperDelaunay()) continue;
            
            // For each facet (simplex), check if point is inside
            // We use barycentric coordinates to test containment
            if (isPointInSimplex(point, *it)) {
                facet = *it;
                return true;
            }
        }
        
        // If no containing simplex found but we have valid data, use a simpler approach
        // Just check the first valid simplex to see if we can perform interpolation
        for (orgQhull::QhullFacetList::iterator it = facets.begin(); it != facets.end(); ++it) {
            if (!it->isGood() || it->isUpperDelaunay()) continue;
            
            // Try this facet anyway for points very close to vertices
            orgQhull::QhullVertexSet vertices = it->vertices();
            if (vertices.count() == dimension_ + 1) {
                // Check if point is very close to any vertex
                for (orgQhull::QhullVertexSet::iterator vit = vertices.begin(); vit != vertices.end(); ++vit) {
                    orgQhull::QhullPoint vertex_point = (*vit).point();
                    double dist_sq = 0.0;
                    for (size_t i = 0; i < dimension_; ++i) {
                        double diff = point[i] - vertex_point[static_cast<int>(i)];
                        dist_sq += diff * diff;
                    }
                    if (dist_sq < 1e-20) {  // Very close to vertex
                        facet = *it;
                        return true;
                    }
                }
            }
        }
        
    } catch (const std::exception&) {
        return false;
    }
    
    return false;
}

bool LinearNdInterpolator::isPointInSimplex(const std::vector<double>& point,
                                           const orgQhull::QhullFacet& facet) const {
    try {
        // Get vertices of the simplex
        orgQhull::QhullVertexSet vertices = facet.vertices();
        if (vertices.count() != dimension_ + 1) {
            return false;
        }
        
        // Calculate barycentric coordinates
        std::vector<double> barycentricCoords = calculateBarycentricCoordinates(point, facet);
        
        if (barycentricCoords.empty()) {
            return false;
        }
        
        // Check if all barycentric coordinates are non-negative and sum to 1
        double sum = 0.0;
        const double tolerance = 1e-8;  // Relaxed tolerance for numerical stability
        
        for (double coord : barycentricCoords) {
            if (coord < -tolerance) {
                return false;  // Point is outside
            }
            sum += coord;
        }
        
        return std::abs(sum - 1.0) < tolerance;
        
    } catch (const std::exception&) {
        return false;
    }
}


std::vector<double> LinearNdInterpolator::calculateBarycentricCoordinates(
    const std::vector<double>& point, 
    const orgQhull::QhullFacet& facet) const {
    
    std::vector<double> coords;
    
    try {
        // Get vertices of the simplex
        orgQhull::QhullVertexSet vertices = facet.vertices();
        size_t num_vertices = vertices.count();
        
        if (num_vertices != dimension_ + 1) {
            // Invalid simplex, return empty
            return coords;
        }
        
        // Create matrix for barycentric coordinate calculation
        // System: Ax = b where A is (dimension+1) x (dimension+1) matrix
        std::vector<std::vector<double>> matrix(dimension_ + 1, std::vector<double>(dimension_ + 1));
        std::vector<double> rhs(dimension_ + 1);
        
        // Fill matrix with vertex coordinates and ones
        size_t vertex_idx = 0;
        for (orgQhull::QhullVertexSet::iterator vit = vertices.begin(); vit != vertices.end(); ++vit, ++vertex_idx) {
            orgQhull::QhullPoint vertex_point = (*vit).point();
            for (size_t i = 0; i < dimension_; ++i) {
                matrix[i][vertex_idx] = vertex_point[static_cast<int>(i)];
            }
            matrix[dimension_][vertex_idx] = 1.0;  // Last row is all ones
        }
        
        // Fill right-hand side
        for (size_t i = 0; i < dimension_; ++i) {
            rhs[i] = point[i];
        }
        rhs[dimension_] = 1.0;  // Constraint: sum of coordinates = 1
        
        // Solve linear system using Gaussian elimination
        coords = solveLinearSystem(matrix, rhs);
        
    } catch (const std::exception&) {
        // Return empty on error
        coords.clear();
    }
    
    return coords;
}

double LinearNdInterpolator::interpolateInSimplex(const std::vector<double>& barycentricCoords,
                                                const orgQhull::QhullFacet& facet) const {
    if (barycentricCoords.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    try {
        // Get vertices and their corresponding values
        orgQhull::QhullVertexSet vertices = facet.vertices();
        
        double interpolated_value = 0.0;
        size_t coord_idx = 0;
        
        for (orgQhull::QhullVertexSet::iterator vit = vertices.begin(); 
             vit != vertices.end() && coord_idx < barycentricCoords.size(); 
             ++vit, ++coord_idx) {
            
            // Find the corresponding point index in our original data
            orgQhull::QhullPoint vertex_point = (*vit).point();
            size_t point_idx = findPointIndex(vertex_point);
            
            if (point_idx < values_.size()) {
                interpolated_value += barycentricCoords[coord_idx] * values_[point_idx];
            }
        }
        
        return interpolated_value;
        
    } catch (const std::exception&) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

std::vector<double> LinearNdInterpolator::solveLinearSystem(
    std::vector<std::vector<double>>& matrix, 
    std::vector<double>& rhs) const {
    
    const size_t n = matrix.size();
    std::vector<double> solution(n);
    
    try {
        // Forward elimination with partial pivoting
        for (size_t i = 0; i < n; ++i) {
            // Find pivot
            size_t pivot_row = i;
            for (size_t k = i + 1; k < n; ++k) {
                if (std::abs(matrix[k][i]) > std::abs(matrix[pivot_row][i])) {
                    pivot_row = k;
                }
            }
            
            // Swap rows if needed
            if (pivot_row != i) {
                std::swap(matrix[i], matrix[pivot_row]);
                std::swap(rhs[i], rhs[pivot_row]);
            }
            
            // Check for zero pivot (singular matrix)
            const double pivot = matrix[i][i];
            if (std::abs(pivot) < 1e-12) {
                return std::vector<double>();  // Return empty on singular matrix
            }
            
            // Eliminate
            for (size_t k = i + 1; k < n; ++k) {
                const double factor = matrix[k][i] / pivot;
                for (size_t j = i; j < n; ++j) {
                    matrix[k][j] -= factor * matrix[i][j];
                }
                rhs[k] -= factor * rhs[i];
            }
        }
        
        // Back substitution
        for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
            solution[i] = rhs[i];
            for (size_t j = i + 1; j < n; ++j) {
                solution[i] -= matrix[i][j] * solution[j];
            }
            solution[i] /= matrix[i][i];
        }
        
        return solution;
        
    } catch (const std::exception&) {
        return std::vector<double>();  // Return empty on error
    }
}

// Helper method to find point index by coordinates
size_t LinearNdInterpolator::findPointIndex(const orgQhull::QhullPoint& vertex_point) const {
    const double tolerance = 1e-12;
    
    for (size_t i = 0; i < num_points_; ++i) {
        bool match = true;
        for (size_t j = 0; j < dimension_; ++j) {
            if (std::abs(vertex_point[static_cast<int>(j)] - points_[i][j]) > tolerance) {
                match = false;
                break;
            }
        }
        if (match) {
            return i;
        }
    }
    
    // Should not happen in normal cases
    return 0;
}