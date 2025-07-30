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

    // 点データを平坦化
    std::vector<double> flat_points(n_points * n_dims);
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < n_dims; ++j) {
            flat_points[i * n_dims + j] = points[i][j];
        }
    }

    // SciPyと同等のQhullオプション設定
    std::string options = "Qbb Qc Qz Q12";
    if (n_dims >= 5) {
        options += " Qx";
    }
    // Qt（三角分割出力）オプションを追加
    options += " Qt";

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
    
    auto facetList = qhull_->facetList();
    int current_id = 0;
    
    for (auto facet : facetList) {
        if (facet.isUpperDelaunay()) {
            continue;
        }
        
        if (current_id == simplex_id) {
            auto vertices = facet.vertices();
            
            if (vertices.size() == n_dims + 1) {
                // 線形システムを解いて重心座標を計算
                // Ax = b の形で、A は (n_dims x n_dims) 行列
                std::vector<std::vector<double>> A(n_dims, std::vector<double>(n_dims));
                std::vector<double> b(n_dims);
                
                auto vertex_iter = vertices.begin();
                auto first_vertex = *vertex_iter;
                ++vertex_iter;
                
                // 最初の頂点を基準にして相対座標を計算
                for (size_t i = 0; i < n_dims; ++i) {
                    auto vertex = *vertex_iter;
                    for (size_t j = 0; j < n_dims; ++j) {
                        A[j][i] = vertex.point()[j] - first_vertex.point()[j];
                    }
                    ++vertex_iter;
                }
                
                // 右辺ベクトル b を設定
                for (size_t j = 0; j < n_dims; ++j) {
                    b[j] = point[j] - first_vertex.point()[j];
                }
                
                // ガウス消去法で解く（簡単な実装）
                std::vector<double> solution = solveLinearSystem(A, b);
                
                // 重心座標を設定
                double sum = 0.0;
                for (size_t i = 0; i < n_dims; ++i) {
                    barycentric[i + 1] = solution[i];
                    sum += solution[i];
                }
                barycentric[0] = 1.0 - sum;
            }
            break;
        }
        current_id++;
    }
    
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
            simplex.push_back(vertex.id());
        }
        
        if (!simplex.empty()) {
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
