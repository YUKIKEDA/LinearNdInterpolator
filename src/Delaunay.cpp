#include "Delaunay.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>

Delaunay::Delaunay(const std::vector<std::vector<double>>& points) {
    // SciPyに基づく入力検証
    if (points.empty()) {
        throw std::invalid_argument("Input points cannot be empty.");
    }
    
    if (points[0].empty()) {
        throw std::invalid_argument("Input points must have at least one dimension.");
    }
    
    const size_t n_points = points.size();
    const size_t n_dims = points[0].size();
    
    // 最小次元数の確認（SciPyでは2次元以上が必要）
    if (n_dims < 2) {
        throw std::invalid_argument("Points must have at least 2 dimensions for triangulation.");
    }
    
    // 最小点数の確認（n次元では最低n+1個の点が必要）
    if (n_points < n_dims + 1) {
        throw std::invalid_argument(
            "Need at least " + std::to_string(n_dims + 1) + 
            " points for " + std::to_string(n_dims) + "-dimensional triangulation.");
    }
    
    // 次元の一貫性とNaN/無限大値のチェック
    for (size_t i = 0; i < n_points; ++i) {
        if (points[i].size() != n_dims) {
            throw std::invalid_argument("All points must have the same dimension.");
        }
        
        for (size_t j = 0; j < n_dims; ++j) {
            if (!std::isfinite(points[i][j])) {
                throw std::invalid_argument("Input points contain NaN or infinite values.");
            }
        }
    }

    // 点データを平坦化し、内部で保持
    points_ = points;
    std::vector<double> flat_points(n_points * n_dims);
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < n_dims; ++j) {
            flat_points[i * n_dims + j] = points[i][j];
        }
    }

    // SciPyと同等のQhullオプション設定
    // d: Delaunay三角分割, Qbb: バウンディングボックス計算, Qc: 共面点保持
    // Qz: 無限遠点を追加（共円・共球問題対応）, Q12: 広角を許可, Qt: 三角分割出力
    std::string options = "d Qbb Qc Qz Q12 Qt";
    if (n_dims >= 5) {
        options += " Qx";
    }

    try {
        qhull_ = std::make_unique<orgQhull::Qhull>();
        qhull_->runQhull("", n_dims, n_points, flat_points.data(), options.c_str());
        
        if (qhull_->qhullStatus() != 0) {
            throw std::runtime_error(
                "Qhull triangulation failed with status: " + 
                std::to_string(qhull_->qhullStatus()));
        }
    } catch (const orgQhull::QhullError &e) {
        throw std::runtime_error("Qhull error: " + std::string(e.what()));
    }
}

Delaunay::~Delaunay() = default;

orgQhull::Qhull* Delaunay::getQhull() const {
    return qhull_.get();
}

int Delaunay::findSimplex(const std::vector<double>& point) const {
    if (!qhull_ || qhull_->qhullStatus() != 0) {
        return -1;
    }
    
    const size_t n_dims = point.size();
    
    // QhullのfacetListを使用して単体を検索
    auto facetList = qhull_->facetList();
    int simplex_id = 0;
    
    for (auto facet : facetList) {
        if (facet.isUpperDelaunay()) {
            continue; // 上半分のDelaunay面は無視
        }
        
        // 点が単体内部にあるかチェック
        bool inside = true;
        auto vertices = facet.vertices();
        
        if (vertices.size() != n_dims + 1) {
            continue; // n次元では n+1 個の頂点が必要
        }
        
        // 重心座標で内部判定（簡略版）
        std::vector<double> barycentric = calculateBarycentricCoordinates(point, simplex_id);
        
        for (double coord : barycentric) {
            if (coord < -1e-10) { // 数値誤差を考慮した許容範囲
                inside = false;
                break;
            }
        }
        
        if (inside) {
            return simplex_id;
        }
        
        simplex_id++;
    }
    
    return -1; // 凸包外
}

std::vector<double> Delaunay::calculateBarycentricCoordinates(
    const std::vector<double>& point, int simplex_id) const {
    
    if (!qhull_ || qhull_->qhullStatus() != 0) {
        return {};
    }
    
    const size_t n_dims = point.size();
    std::vector<double> barycentric(n_dims + 1, 0.0);
    
    // Get the specific simplex vertices
    auto simplices = getSimplices();
    if (simplex_id < 0 || simplex_id >= static_cast<int>(simplices.size())) {
        return {};
    }
    
    const auto& simplex = simplices[simplex_id];
    if (simplex.size() != n_dims + 1) {
        return {};
    }
    
    // Build the transformation matrix A and vector b
    // A * lambda = point - vertex[0], where lambda are the last n_dims barycentric coordinates
    std::vector<std::vector<double>> A(n_dims, std::vector<double>(n_dims));
    std::vector<double> b(n_dims);
    
    // Use first vertex as origin
    int first_vertex_idx = simplex[0];
    if (first_vertex_idx < 0 || first_vertex_idx >= static_cast<int>(points_.size())) {
        return {};
    }
    
    const auto& first_vertex = points_[first_vertex_idx];
    
    // Build matrix A: columns are (vertex[i] - vertex[0]) for i = 1, 2, ..., n_dims
    for (size_t i = 0; i < n_dims; ++i) {
        int vertex_idx = simplex[i + 1];
        if (vertex_idx < 0 || vertex_idx >= static_cast<int>(points_.size())) {
            return {};
        }
        
        const auto& vertex = points_[vertex_idx];
        for (size_t j = 0; j < n_dims; ++j) {
            A[j][i] = vertex[j] - first_vertex[j];
        }
    }
    
    // Build vector b: point - vertex[0]
    for (size_t j = 0; j < n_dims; ++j) {
        b[j] = point[j] - first_vertex[j];
    }
    
    // Solve linear system A * lambda = b
    std::vector<double> lambda = solveLinearSystem(A, b);
    
    // Convert to barycentric coordinates
    // barycentric[0] = 1 - sum(lambda)
    // barycentric[i] = lambda[i-1] for i = 1, ..., n_dims
    double sum = 0.0;
    for (size_t i = 0; i < n_dims; ++i) {
        barycentric[i + 1] = lambda[i];
        sum += lambda[i];
    }
    barycentric[0] = 1.0 - sum;
    
    return barycentric;
}

std::vector<std::vector<int>> Delaunay::getSimplices() const {
    std::vector<std::vector<int>> simplices;
    
    if (!qhull_ || qhull_->qhullStatus() != 0) {
        return simplices;
    }
    
    auto facetList = qhull_->facetList();
    
    for (auto facet : facetList) {
        if (facet.isUpperDelaunay()) {
            continue;
        }
        
        auto vertices = facet.vertices();
        std::vector<int> simplex;
        
        for (auto vertex : vertices) {
            // Get vertex coordinates
            auto vertex_point = vertex.point();
            std::vector<double> vertex_coords(vertex_point.coordinates(), 
                                            vertex_point.coordinates() + vertex_point.dimension());
            
            // Find matching point in original points_ array
            int point_index = -1;
            for (size_t i = 0; i < points_.size(); ++i) {
                bool match = true;
                for (size_t j = 0; j < points_[i].size(); ++j) {
                    if (std::abs(points_[i][j] - vertex_coords[j]) > 1e-12) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    point_index = static_cast<int>(i);
                    break;
                }
            }
            
            // Add valid indices only
            if (point_index >= 0) {
                simplex.push_back(point_index);
            }
        }
        
        // Only add simplices with correct number of vertices
        if (simplex.size() == points_[0].size() + 1) {
            simplices.push_back(simplex);
        }
    }
    
    return simplices;
}

// ヘルパーメソッド: 線形システムを解く
std::vector<double> Delaunay::solveLinearSystem(
    const std::vector<std::vector<double>>& A, 
    const std::vector<double>& b) const {
    
    const size_t n = A.size();
    std::vector<double> x(n, 0.0);
    
    if (n == 0 || A[0].size() != n || b.size() != n) {
        return x;
    }
    
    // コピーを作成（ガウス消去法で変更するため）
    std::vector<std::vector<double>> mat(n, std::vector<double>(n + 1));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            mat[i][j] = A[i][j];
        }
        mat[i][n] = b[i];
    }
    
    // 前進消去
    for (size_t i = 0; i < n; ++i) {
        // ピボット選択
        size_t pivot_row = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::abs(mat[k][i]) > std::abs(mat[pivot_row][i])) {
                pivot_row = k;
            }
        }
        
        if (pivot_row != i) {
            std::swap(mat[i], mat[pivot_row]);
        }
        
        // 特異行列のチェック
        if (std::abs(mat[i][i]) < 1e-12) {
            return x; // 解けない場合は零ベクトルを返す
        }
        
        // 消去
        for (size_t k = i + 1; k < n; ++k) {
            double factor = mat[k][i] / mat[i][i];
            for (size_t j = i; j <= n; ++j) {
                mat[k][j] -= factor * mat[i][j];
            }
        }
    }
    
    // 後退代入
    for (int i = n - 1; i >= 0; --i) {
        x[i] = mat[i][n];
        for (size_t j = i + 1; j < n; ++j) {
            x[i] -= mat[i][j] * x[j];
        }
        x[i] /= mat[i][i];
    }
    
    return x;
}
