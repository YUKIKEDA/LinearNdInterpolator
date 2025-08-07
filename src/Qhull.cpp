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
) : points_(points), command_(command), options_(options), initialized_(false), computed_(false) 
{
    // 入力データの検証
    // 点の配列が空でないことを確認
    if (points_.empty()) {
        throw std::invalid_argument("Points array cannot be empty");
    }

    // 次元数を取得し、最低2次元であることを確認
    ndim_ = points_[0].size();
    if (ndim_ < 2) {
        throw std::invalid_argument("Need at least 2-D data");
    }

    // 全ての点がNaNでないことを確認
    for (const auto& point : points_) {
        for (double coord : point) {
            if (std::isnan(coord)) {
                throw std::invalid_argument("Points cannot contain NaN");
            }
        }
    }
    
    // Check if this is Delaunay mode
    is_delaunay_ = (command == "d" || command == "v");
    
    // SciPy準拠: コンストラクタでqh_new_qhullを実行（_Qhull.__init__相当）
    qh_zero(&qh_qh, nullptr);
    
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
        qhull_command = "qhull d " + options_;
    } else {
        qhull_command = "qhull " + options_;
    }
    
    // SciPy準拠: required_options ("Qt") を追加
    qhull_command += " Qt";
    
    // SciPy準拠: qh_new_qhullの実行（_Qhull.__init__相当）
    int exitcode = qh_new_qhull(
        &qh_qh, dim, numpoints, qpoints.data(), 
        False, const_cast<char*>(qhull_command.c_str()), 
        nullptr, stderr
    );
    
    if (exitcode) {
        throw std::runtime_error("Qhull failed with exit code: " + std::to_string(exitcode));
    }
    
    initialized_ = true;
}

Qhull::~Qhull() {
    try {
        // SciPy準拠: _Qhull.__dealloc__() に対応するクリーンアップ処理
        if (initialized_) {
            int curlong, totlong;
            
            // SciPy準拠: qh_freeqhull, qh_memfreeshortの順で実行
            qh_freeqhull(&qh_qh, qh_ALL);
            qh_memfreeshort(&qh_qh, &curlong, &totlong);
            
            // SciPy準拠: メモリリークのチェック（デストラクタでは警告のみ）
            if (curlong != 0 || totlong != 0) {
                std::cerr << "qhull: did not free " << totlong << " bytes (" << curlong << " pieces)" << std::endl;
            }
            
            initialized_ = false;
            computed_ = false;
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception in Qhull destructor: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception in Qhull destructor" << std::endl;
    }
}

void Qhull::triangulate() {
    if (computed_) return;
    
    try {
        // SciPy準拠: qh_triangulateのみを実行（qh_new_qhullはコンストラクタで実行済み）
        qh_triangulate(&qh_qh);
        
        computed_ = true;
        
    } catch (const std::exception& e) {
        throw;
    }
}

std::pair<double, double> Qhull::getParaboloidShiftScale() const {
    if (!initialized_) {
        throw std::runtime_error("Qhull must be initialized before accessing paraboloid parameters");
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
        facetT *neighbor;
        vertexT *vertex;
        pointT *point;
        int i, j, ipoint, ncoplanar;
        int facet_ndim = static_cast<int>(ndim_);
        
        // SciPy準拠: _is_halfspaces処理（今回は非対応）
        // if (_is_halfspaces) { facet_ndim = ndim_ - 1; }
        
        if (is_delaunay_) {
            facet_ndim += 1;  // Delaunayの場合は次元を1つ増やす
        }
        
        // SciPy準拠: id_mapの作成と初期化 (lines 610-615)
        std::vector<int> id_map(qh_qh.facet_id, -1);
        
        // SciPy準拠: 第1段階 - ファセット数のカウントとid_mapの作成 (lines 617-632)
        facet = qh_qh.facet_list;
        j = 0;
        while (facet && facet->next) {
            if (!is_delaunay_ || facet->upperdelaunay == qh_qh.UPPERdelaunay) {
                // SciPy準拠: simplicialチェック (lines 621-627)
                if (!facet->simplicial && (
                    qh_setsize(const_cast<qhT*>(&qh_qh), facet->vertices) != facet_ndim ||
                    qh_setsize(const_cast<qhT*>(&qh_qh), facet->neighbors) != facet_ndim)) {
                    throw std::runtime_error(
                        "non-simplicial facet encountered: " + 
                        std::to_string(qh_setsize(const_cast<qhT*>(&qh_qh), facet->vertices)) + " vertices");
                }
                
                id_map[facet->id] = j;
                j++;
            }
            facet = facet->next;
        }
        
        int num_facets = j;
        
        // SciPy準拠: 出力配列の割り当て (lines 635-638)
        simplices.resize(num_facets);
        for (int k = 0; k < num_facets; ++k) {
            simplices[k].resize(facet_ndim);
        }
        good.resize(num_facets, false);
        neighbors.resize(num_facets);
        for (int k = 0; k < num_facets; ++k) {
            neighbors[k].resize(facet_ndim);
        }
        equations.resize(num_facets);
        for (int k = 0; k < num_facets; ++k) {
            equations[k].resize(facet_ndim + 1);
        }
        
        // SciPy準拠: coplanar配列の初期化 (lines 640-642)
        ncoplanar = 0;
        coplanar.reserve(10);  // 初期サイズ10で予約
        
        // SciPy準拠: 第2段階 - ファセット情報の取得 (lines 645-712)
        facet = qh_qh.facet_list;
        j = 0;
        while (facet && facet->next) {
            if (is_delaunay_ && facet->upperdelaunay != qh_qh.UPPERdelaunay) {
                facet = facet->next;
                continue;
            }
            
            // SciPy準拠: 3次元でのorientation処理 (lines 656-671)
            unsigned int lower_bound = 0;
            if (is_delaunay_ && facet->toporient == qh_ORIENTclock && facet_ndim == 3) {
                // 反時計回りを維持するために最初と2番目のインデックスを交換
                for (i = 0; i < 2; ++i) {
                    unsigned int swapped_index = 1 ^ i;
                    
                    // 頂点情報を保存
                    vertex = static_cast<vertexT*>(facet->vertices->e[i].p);
                    ipoint = qh_pointid(const_cast<qhT*>(&qh_qh), vertex->point);
                    simplices[j][swapped_index] = ipoint;
                    
                    // 隣接情報を保存
                    neighbor = static_cast<facetT*>(facet->neighbors->e[i].p);
                    neighbors[j][swapped_index] = id_map[neighbor->id];
                }
                lower_bound = 2;
            }
            
            // SciPy準拠: 残りの頂点と隣接情報を保存 (lines 673-681)
            for (i = lower_bound; i < static_cast<unsigned int>(facet_ndim); ++i) {
                // 頂点情報を保存
                vertex = static_cast<vertexT*>(facet->vertices->e[i].p);
                ipoint = qh_pointid(const_cast<qhT*>(&qh_qh), vertex->point);
                simplices[j][i] = ipoint;
                
                // 隣接情報を保存
                neighbor = static_cast<facetT*>(facet->neighbors->e[i].p);
                neighbors[j][i] = id_map[neighbor->id];
            }
            
            // SciPy準拠: 単体方程式情報を保存 (lines 684-686)
            for (i = 0; i < facet_ndim; ++i) {
                equations[j][i] = facet->normal[i];
            }
            equations[j][facet_ndim] = facet->offset;
            
            // SciPy準拠: coplanar情報を保存 (lines 688-705)
            if (facet->coplanarset) {
                for (i = 0; i < qh_setsize(const_cast<qhT*>(&qh_qh), facet->coplanarset); ++i) {
                    point = static_cast<pointT*>(facet->coplanarset->e[i].p);
                    
                    // SciPy準拠: qh_nearvertex相当の処理
                    double dist;
                    vertex = qh_nearvertex(const_cast<qhT*>(&qh_qh), facet, point, &dist);
                    
                    // 動的配列拡張
                    if (ncoplanar >= static_cast<int>(coplanar.size())) {
                        coplanar.resize(2 * ncoplanar + 1, std::vector<int>(3));
                    } else if (coplanar.size() <= static_cast<size_t>(ncoplanar)) {
                        coplanar.resize(ncoplanar + 1, std::vector<int>(3));
                    }
                    
                    coplanar[ncoplanar] = {
                        qh_pointid(const_cast<qhT*>(&qh_qh), point),
                        id_map[facet->id],
                        qh_pointid(const_cast<qhT*>(&qh_qh), vertex->point)
                    };
                    ncoplanar++;
                }
            }
            
            // SciPy準拠: good情報を保存 (line 708)
            good[j] = facet->good;
            
            j++;
            facet = facet->next;
        }
        
        // coplanar配列のサイズを実際の要素数に調整
        coplanar.resize(ncoplanar);
        
    } catch (const std::exception& e) {
        throw;
    }
    
    return std::make_tuple(simplices, neighbors, equations, coplanar, good);
}

} // namespace qhull
