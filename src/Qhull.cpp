#include "Qhull.h"
#include <map>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace qhull {

// ---------------------------------------------------------------------------
// Qhull クラスの実装
// ---------------------------------------------------------------------------

Qhull::Qhull(
    const std::string& command, 
    const std::vector<std::vector<double>>& points, 
    const std::string& options
) 
    : points_(points), command_(command), options_(options), computed_(false) 
{
    
    if (points_.empty()) {
        throw std::invalid_argument("Points array cannot be empty");
    }
    
    ndim_ = points_[0].size();
    if (ndim_ < 2) {
        throw std::invalid_argument("Need at least 2-D data");
    }
    
    // Check for NaN values
    for (const auto& point : points_) {
        for (double coord : point) {
            if (std::isnan(coord)) {
                throw std::invalid_argument("Points cannot contain NaN");
            }
        }
    }
    
    // Check if this is Delaunay mode
    is_delaunay_ = (command == "d" || command == "v");
    
    // Qhullコンテキストの初期化
    qh_zero(&qh_qh, nullptr);
}

Qhull::~Qhull() {
    std::cout << "DEBUG: Qhull destructor called, computed_=" << computed_ << std::endl;
    try {
        // 一時的にクリーンアップを無効にして問題を特定
        std::cout << "DEBUG: Qhull cleanup skipped for debugging (memory leak potential)" << std::endl;
        /*
        if (computed_ && qh_qh.hull_dim > 0) {
            std::cout << "DEBUG: About to call qh_freeqhull..." << std::endl;
            // Qhullの状態をチェックしてからクリーンアップ
            if (qh_qh.facet_list) {
                qh_freeqhull(&qh_qh, qh_ALL);
                qh_memfreeshort(&qh_qh, nullptr, nullptr);
                std::cout << "DEBUG: Qhull cleanup completed" << std::endl;
            } else {
                std::cout << "DEBUG: Qhull cleanup skipped - invalid state" << std::endl;
            }
        } else {
            std::cout << "DEBUG: Qhull cleanup skipped - not computed or invalid" << std::endl;
        }
        */
    } catch (...) {
        std::cout << "DEBUG: Exception in Qhull destructor" << std::endl;
    }
    std::cout << "DEBUG: Qhull destructor completed" << std::endl;
}

void Qhull::triangulate() {
    if (computed_) return;
    
    try {
        // 点群データをQhull形式に変換
        int numpoints = static_cast<int>(points_.size());
        int dim = static_cast<int>(ndim_);
        
        // フラットな配列に変換
        std::vector<coordT> qpoints(numpoints * dim);
        for (int i = 0; i < numpoints; ++i) {
            for (int j = 0; j < dim; ++j) {
                qpoints[i * dim + j] = static_cast<coordT>(points_[i][j]);
            }
        }
        
        // Qhullオプションの準備
        std::string qhull_command;
        if (is_delaunay_) {
            qhull_command = "qhull d " + options_ + " Qt";  // delaunay + triangulated output
        } else {
            qhull_command = "qhull " + options_;
        }
        
        // Qhullの実行
        int exitcode = qh_new_qhull(
            &qh_qh, dim, numpoints, qpoints.data(), 
            False, const_cast<char*>(qhull_command.c_str()), 
            nullptr, stderr
        );
        
        if (exitcode) {
            throw std::runtime_error("Qhull failed with exit code: " + std::to_string(exitcode));
        }
        
        // 三角分割の実行
        qh_triangulate(&qh_qh);
        
        computed_ = true;
        
    } catch (const std::exception& e) {
        throw;
    }
}

std::pair<double, double> Qhull::getParaboloidShiftScale() const {
    if (!computed_) {
        throw std::runtime_error("Triangulation must be computed before accessing paraboloid parameters");
    }
    
    double paraboloid_scale;
    double paraboloid_shift;
    
    // SciPyの実装に基づく
    if (qh_qh.SCALElast) {
        paraboloid_scale = qh_qh.last_newhigh / (qh_qh.last_high - qh_qh.last_low);
        paraboloid_shift = -qh_qh.last_low * paraboloid_scale;
    } else {
        paraboloid_scale = 1.0;
        paraboloid_shift = 0.0;
    }
    
    return {paraboloid_shift, paraboloid_scale};
}

SimplexFacetResult Qhull::getSimplexFacetArray() const {
    if (!computed_) {
        throw std::runtime_error("Triangulation must be computed before accessing results");
    }
    
    std::vector<std::vector<int>> simplices;
    std::vector<std::vector<int>> neighbors;
    std::vector<std::vector<double>> equations;
    std::vector<std::vector<int>> coplanar;
    std::vector<bool> good;
    
    try {
        facetT *facet;
        vertexT *vertex;
        int facet_ndim = static_cast<int>(ndim_);
        
        if (is_delaunay_) {
            facet_ndim += 1;  // Delaunayの場合は次元を1つ増やす
        }
        
        // ファセットIDマップを作成
        std::vector<int> id_map(qh_qh.facet_id, -1);
        
        // 有効なファセットをカウントしてIDマップを作成
        int j = 0;
        for (facet = qh_qh.facet_list; facet && facet->next; facet = facet->next) {
            if (!is_delaunay_ || facet->upperdelaunay == qh_qh.UPPERdelaunay) {
                id_map[facet->id] = j;
                j++;
            }
        }
        
        int num_facets = j;
        
        // 結果の配列を初期化
        simplices.resize(num_facets);
        neighbors.resize(num_facets);
        equations.resize(num_facets);
        good.resize(num_facets);
        
        // ファセット情報を取得
        j = 0;
        for (facet = qh_qh.facet_list; facet && facet->next; facet = facet->next) {
            if (is_delaunay_ && facet->upperdelaunay != qh_qh.UPPERdelaunay) {
                continue;
            }
            
            // 頂点情報を保存
            simplices[j].resize(facet_ndim);
            int i = 0;
            setT *vertices = facet->vertices;
            for (int k = 0; k < qh_setsize(const_cast<qhT*>(&qh_qh), vertices) && i < facet_ndim; k++) {
                vertex = static_cast<vertexT*>(vertices->e[k].p);
                int point_id = qh_pointid(const_cast<qhT*>(&qh_qh), vertex->point);
                simplices[j][i] = point_id;
                i++;
            }
            
            // 隣接情報を保存
            neighbors[j].resize(facet_ndim);
            i = 0;
            setT *neighs = facet->neighbors;
            for (int k = 0; k < qh_setsize(const_cast<qhT*>(&qh_qh), neighs) && i < facet_ndim; k++) {
                facetT *neighbor = static_cast<facetT*>(neighs->e[k].p);
                neighbors[j][i] = (neighbor && neighbor->id < id_map.size() && id_map[neighbor->id] >= 0) ? id_map[neighbor->id] : -1;
                i++;
            }
            
            // 方程式情報を保存
            equations[j].resize(facet_ndim + 1);
            for (i = 0; i < facet_ndim; ++i) {
                equations[j][i] = facet->normal[i];
            }
            equations[j][facet_ndim] = facet->offset;
            
            good[j] = true;
            j++;
        }
        
        // coplanar情報（簡略化）
        coplanar.resize(num_facets);
        
    } catch (const std::exception& e) {
        throw;
    }
    
    return std::make_tuple(simplices, neighbors, equations, coplanar, good);
}

} // namespace qhull
