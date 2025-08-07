#include "Delaunay.h"
#include <stdexcept>
#include <iostream>
#include <limits>
#include <algorithm>
#include <cmath>

namespace qhull {

// ---------------------------------------------------------------------------
// Delaunay クラスの実装
// ---------------------------------------------------------------------------

Delaunay::Delaunay(const std::vector<std::vector<double>>& input_points) 
{
    // --------------------------------------------------
    // Delaunayクラスの __init__ 相当の処理
    // --------------------------------------------------
    points_ = input_points;
    npoints_ = points_.size();
    ndim_ = npoints_ > 0 ? points_[0].size() : 0;

    std::string qhull_options = "Qbb Qc Qz Q12";
    if (ndim_ >= 5) {
        qhull_options += " Qx";
    }
    
    //HACK: 要修正

    // Scipy: `required_options=b"Qt"` を追加
    std::string full_qhull_options = qhull_options + " Qt";

    // Scipy: `qhull = _Qhull(b"d", points, ...)`
    // `b"d"`はドロネー三角形分割(delaunay)を意味するコマンド。
    Qhull qhull_runner("d", points_, full_qhull_options);

    // --------------------------------------------------
    // _QhullUserクラスの __init__ 相当の処理
    // --------------------------------------------------

    // Qhullの実行結果をメンバ変数に格納する。
    update(qhull_runner);
}

const std::vector<std::vector<int>>& Delaunay::get_simplices() const {
    return simplices_;
}

int Delaunay::findSimplex(
    std::vector<double>& barycentric_coords,
    const std::vector<double>& point,
    int& start_simplex_hint,
    double eps,
    double eps_broad
) const noexcept 
{
    // 点が凸包の外側にある場合は-1を返す
    if (isPointFullyOutside(point, eps)) {
        return -1;
    }

    // 単体が存在しない場合は-1を返す
    if (nsimplex_ <= 0) {
        return -1;
    }

    // 検索開始位置のヒントを設定
    int isimplex = start_simplex_hint;
    // 検索開始位置のヒントが無効な場合は0から開始
    if (isimplex < 0 || isimplex >= static_cast<int>(nsimplex_)) {
        isimplex = 0;
    }

    // 点をパラボロイド変換
    std::vector<double> lifted_point;
    liftPoint(point, lifted_point);

    // Walk the tessellation searching for a facet with a positive planar distance
    double best_dist = distplane(isimplex, lifted_point);

    bool changed = true;
    while (changed) {
        // 正の距離を持つ単体が見つかったらウォークを終了
        if (best_dist > 0) {
            break; 
        }

        changed = false;

        for (size_t k = 0; k < ndim_ + 1; ++k) {
            // SciPy互換: 単体isimplexのk番目の隣接単体を取得
            // SciPy: ineigh = d.neighbors[(ndim+1)*isimplex + k] (1次元配列での線形アクセス)
            // C++:   neighbors_[isimplex][k] (2次元vectorでの直接アクセス)
            // 両方とも機能的に同等で、同じ隣接単体情報にアクセスする
            int ineigh = neighbors_[isimplex][k];
            if (ineigh == -1) {
                continue;
            }
            
            double dist = distplane(ineigh, lifted_point);

            if (dist > best_dist + eps * (1.0 + std::abs(best_dist))) {
                isimplex = ineigh;
                best_dist = dist;
                changed = true;
            }
        }
    }

    // We should now be somewhere near the simplex containing the point,
    // locate it with a directed search
    start_simplex_hint = isimplex;
    return findSimplexDirected(barycentric_coords, point, start_simplex_hint, eps, eps_broad);
}

void Delaunay::update(Qhull& qhull) 
{
    // ---------------------------------
    // --- `Delaunay._update` の処理 ---
    // ---------------------------------
    // 1. 三角形分割の実行
    qhull.triangulate(); //HACK: 要確認

    // 2. パラボロイド変換パラメータの取得
    auto ps = qhull.getParaboloidShiftScale(); //HACK: 要確認
    paraboloid_shift_ = ps.first;
    paraboloid_scale_ = ps.second;

    // 3. 幾何データの抽出
    std::tie(simplices_, neighbors_, equations_, coplanar_, good_) = qhull.getSimplexFacetArray(); //HACK: 要確認

    // 4. 単体数を保存
    nsimplex_ = simplices_.size();

    // 5. 遅延計算フラグの初期化
    // C++では、遅延評価されるポインタをnullptrに、フラグをfalseに設定する
    transform_computed_ = false;
    vertex_to_simplex_computed_ = false;
    vertex_neighbor_vertices_computed_ = false;

    // -----------------------------------
    // --- `_QhullUser._update` の処理 ---
    // -----------------------------------
    // 6. 最後に共通属性を計算する
    calculateBounds();
}

void Delaunay::calculateBounds() {
    if (npoints_ == 0) {
        return;
    }

    // 初期値を最初の点で設定
    min_bound_ = points_[0];
    max_bound_ = points_[0];

    // 2番目以降の点で比較して更新
    for (size_t i = 1; i < npoints_; ++i) {
        for (size_t j = 0; j < ndim_; ++j) {
            min_bound_[j] = std::min(min_bound_[j], points_[i][j]);
            max_bound_[j] = std::max(max_bound_[j], points_[i][j]);
        }
    }
}

bool Delaunay::isPointFullyOutside(const std::vector<double>& point, double eps) const noexcept 
{
    for (size_t j = 0; j < ndim_; ++j) {
        // 点の座標が境界の最小値より小さいか、最大値より大きい場合
        if (point[j] < min_bound_[j] - eps || point[j] > max_bound_[j] + eps) {
            return true;
        }
    }
    return false;
}

void Delaunay::liftPoint(const std::vector<double>& point, std::vector<double>& lifted_point) const noexcept 
{
    lifted_point.resize(ndim_ + 1);
    double sum_sq = 0;
    for (size_t j = 0; j < ndim_; ++j) {
        lifted_point[j] = point[j];
        sum_sq += point[j] * point[j];
    }
    // パラボロイド変換の計算
    lifted_point[ndim_] = sum_sq * paraboloid_scale_ + paraboloid_shift_;
}

double Delaunay::distplane(int simplex_index, const std::vector<double>& lifted_point) const noexcept 
{
    // Scipyの `d.equations` は (nsimplex, ndim+2) のフラットな配列。
    // C++では `std::vector<std::vector<double>>` なのでアクセスが簡単。
    const auto& eq = equations_[simplex_index]; //HACK: equations_はどこで初期化されている？
    double dist = eq[ndim_ + 1]; // オフセット項
    for (size_t k = 0; k < ndim_ + 1; ++k) {
        dist += eq[k] * lifted_point[k];
    }
    
    return dist;
}

int Delaunay::findSimplexDirected(
    std::vector<double>& barycentric_coords,
    const std::vector<double>& point,
    int& start_simplex_hint,
    double eps, 
    double eps_broad) const noexcept 
{
    int isimplex = start_simplex_hint;

    if (isimplex < 0 || isimplex >= static_cast<int>(nsimplex_)) {
        isimplex = 0;
    }

    // Scipy: for cycle_k in range(1 + d.nsimplex//4):
    // 無限ループを防ぐための最大試行回数
    const int max_cycles = 1 + nsimplex_ / 4;
    int cycle_k;
    
    for (cycle_k = 0; cycle_k < max_cycles; ++cycle_k) {
        
        if (isimplex == -1) {
            break;
        }

        // TODO: 座標変換行列の計算を行う
        // NOTE: `transform`行列の計算が必要です。
        // Scipyでは、この行列はDelaunayオブジェクトの属性として
        // 遅延評価で計算されます。ここではダミーの行列を使います。
        std::vector<double> transform_matrix_for_simplex; // ダミー

        int inside_status = 1; // 1: inside, -1: hopped, 0: outside/degenerate

        for (size_t k = 0; k < ndim_ + 1; ++k) {
            barycentricCoordinateSingle(transform_matrix_for_simplex, point, barycentric_coords, k);

            if (barycentric_coords[k] < -eps) {
                // SciPy互換: directed search中の隣接単体への移動
                // SciPy: m = d.neighbors[(ndim+1)*isimplex + k]
                int m = neighbors_[isimplex][k];
                if (m == -1) {
                    start_simplex_hint = isimplex;
                    return -1; // 凸包の外
                }

                isimplex = m;
                inside_status = -1; // 隣へホップした
                break;

            } else if (barycentric_coords[k] <= 1.0 + eps) {
                // OK, まだ内部にいる可能性がある
            } else {
                // 範囲外 or NaN (縮退した単体)
                inside_status = 0;
            }
        }

        if (inside_status == -1) {
            continue; // 次のサイクルへ

        } else if (inside_status == 1) {
            break; // 正しい単体を発見！ループを抜ける

        } else {
            // 縮退などで失敗 -> 総当たりへ
            isimplex = findSimplexBruteforce(barycentric_coords, point, eps, eps_broad);
            break;
        }
    }

    // for-else: ループが正常に終了しなかった（収束しなかった）場合
    // C++では `cycle_k == max_cycles` で判定
    if (cycle_k == max_cycles && isimplex != -1) {
        isimplex = findSimplexBruteforce(barycentric_coords, point, eps, eps_broad);
    }

    start_simplex_hint = isimplex;
    return isimplex;
}

void Delaunay::barycentricCoordinateSingle(
    const std::vector<double>& transform_matrix_for_simplex,
    const std::vector<double>& point,
    std::vector<double>& barycentric_coords,
    int coord_index
) const noexcept 
{
    const int i = coord_index;

    // Scipy: if i == ndim:
    if (i == static_cast<int>(ndim_)) {
        // c[ndim] = 1.0 - sum(c[0...ndim-1])
        barycentric_coords[ndim_] = 1.0;

        for (size_t j = 0; j < ndim_; ++j) {
            barycentric_coords[ndim_] -= barycentric_coords[j];
        }

    } else {
        // c[i] = sum_j( T_inv[i,j] * (x[j] - r[j]) )
        barycentric_coords[i] = 0.0;
        for (size_t j = 0; j < ndim_; ++j) {
            // transform[ndim*i + j] は T_inv[i,j]
            // transform[ndim*ndim + j] は r[j]
            barycentric_coords[i] += transform_matrix_for_simplex[ndim_ * i + j] * (point[j] - transform_matrix_for_simplex[ndim_ * ndim_ + j]);
        }
    }
}

int Delaunay::findSimplexBruteforce(
    std::vector<double>& barycentric_coords,
    const std::vector<double>& point,
    double eps, double eps_broad
) const noexcept 
{
    // if _is_point_fully_outside(d, x, eps): return -1
    if (isPointFullyOutside(point, eps)) { // eps_broadではなくepsを使うのがScipyのコード通り
        return -1;
    }

    // for isimplex in range(d.nsimplex):
    for (int isimplex = 0; isimplex < static_cast<int>(nsimplex_); ++isimplex) {

        // TODO: 座標変換行列の計算を行う
        // NOTE: `this->transform` は (nsimplex, ndim+1, ndim) の3D配列と見なせる
        // フラットな `std::vector<double>` としてメンバに保持する必要がある。
        // `transform_matrix_for_simplex` はその一部を指すビューまたはコピー。
        const auto& transform_matrix_for_simplex = get_transform_for_simplex(isimplex); // ヘルパー関数を仮定

        // if transform[0] == transform[0]:
        if (!std::isnan(transform_matrix_for_simplex[0])) {
            // transform is valid (non-nan)
            if (barycentricInside(transform_matrix_for_simplex, point, barycentric_coords, eps)) {
                return isimplex;
            }

        } else {
            // transform is invalid (nan, implying degenerate simplex)
            // check neighbors
            for (size_t k = 0; k < ndim_ + 1; ++k) {
                // SciPy互換: 縮退単体の隣接単体を確認
                int ineighbor = neighbors_[isimplex][k];
                
                if (ineighbor == -1) {
                    continue;
                }

                // TODO: 座標変換行列の計算を行う
                const auto& neighbor_transform = get_transform_for_simplex(ineighbor);

                if (std::isnan(neighbor_transform[0])) {
                    continue; // another bad simplex
                }

                barycentricCoordinates(neighbor_transform, point, barycentric_coords);

                bool inside_neighbor = true;
                for (size_t m = 0; m < ndim_ + 1; ++m) {
                    // SciPy互換: 隣接単体が元の単体を参照しているかチェック
                    if (neighbors_[ineighbor][m] == isimplex) {
                        
                        // allow extra leeway towards isimplex
                        if (!(barycentric_coords[m] >= -eps_broad && barycentric_coords[m] <= 1.0 + eps)) {
                            inside_neighbor = false;
                            break;
                        }

                    } else {
                        // normal check
                        if (!(barycentric_coords[m] >= -eps && barycentric_coords[m] <= 1.0 + eps)) {
                            inside_neighbor = false;
                            break;
                        }
                    }
                }
                if (inside_neighbor) {
                    return ineighbor;
                }
            }
        }
    }
    return -1;
}

bool Delaunay::barycentricInside(
    const std::vector<double>& transform_matrix,
    const std::vector<double>& point,
    std::vector<double>& barycentric_coords,
    double eps
) const noexcept 
{
    barycentric_coords.resize(ndim_ + 1);

    // c[ndim] = 1.0
    barycentric_coords[ndim_] = 1.0;

    // for i in range(ndim):
    for (size_t i = 0; i < ndim_; ++i) {
        // c[i] = 0
        barycentric_coords[i] = 0.0;
        // for j in range(ndim):
        for (size_t j = 0; j < ndim_; ++j) {
            barycentric_coords[i] += transform_matrix[ndim_ * i + j] * (point[j] - transform_matrix[ndim_ * ndim_ + j]);
        }
        // c[ndim] -= c[i]
        barycentric_coords[ndim_] -= barycentric_coords[i];

        // if not (-eps <= c[i] <= 1 + eps): return 0
        if (!(barycentric_coords[i] >= -eps && barycentric_coords[i] <= 1.0 + eps)) {
            return false;
        }
    }

    // if not (-eps <= c[ndim] <= 1 + eps): return 0
    if (!(barycentric_coords[ndim_] >= -eps && barycentric_coords[ndim_] <= 1.0 + eps)) {
        return false;
    }

    // return 1
    return true;
}

void Delaunay::barycentricCoordinates(const std::vector<double>& transform_matrix,
    const std::vector<double>& point,
    std::vector<double>& barycentric_coords
) const noexcept {

    barycentric_coords.resize(ndim_ + 1);

    // c[ndim] = 1.0
    barycentric_coords[ndim_] = 1.0;

    // for i in range(ndim):
    for (size_t i = 0; i < ndim_; ++i) {
        // c[i] = 0
        barycentric_coords[i] = 0.0;
        // for j in range(ndim):
        for (size_t j = 0; j < ndim_; ++j) {
            // c[i] += transform[ndim*i + j] * (x[j] - transform[ndim*ndim + j])
            barycentric_coords[i] += transform_matrix[ndim_ * i + j] * (point[j] - transform_matrix[ndim_ * ndim_ + j]);
        }
        // c[ndim] -= c[i]
        barycentric_coords[ndim_] -= barycentric_coords[i];
    }
}

std::vector<double> Delaunay::get_transform_for_simplex(int simplex_index) const {
    if (!transform_computed_) {
        calculateTransformMatrix();
    }
    
    const size_t simplex_size = (ndim_ + 1) * ndim_;
    const size_t start_offset = simplex_index * simplex_size;
    
    return std::vector<double>(
        transform_.begin() + start_offset, 
        transform_.begin() + start_offset + simplex_size);
}

void Delaunay::calculateTransformMatrix() const {
    std::cout << "DEBUG: calculateTransformMatrix() start" << std::endl;
    if (transform_computed_) return;
    
    const size_t nsimplex = simplices_.size();
    std::cout << "DEBUG: nsimplex=" << nsimplex << ", ndim_=" << ndim_ << std::endl;
    
    if (nsimplex == 0) {
        std::cout << "DEBUG: No simplices, skipping transform calculation" << std::endl;
        const_cast<bool&>(transform_computed_) = true;
        return;
    }
    
    const size_t transform_size = nsimplex * (ndim_ + 1) * ndim_;
    const_cast<std::vector<double>&>(transform_).resize(transform_size, 0.0);
    
    // 各単体について変換行列を計算
    for (size_t isimplex = 0; isimplex < nsimplex; ++isimplex) {
        std::cout << "DEBUG: Processing simplex " << isimplex << std::endl;
        const auto& simplex = simplices_[isimplex];
        
        if (simplex.size() != ndim_ + 1) {
            std::cout << "DEBUG: Invalid simplex size: " << simplex.size() << " (expected " << (ndim_+1) << ")" << std::endl;
            continue;
        }
        
        // 頂点インデックスの境界チェック
        bool valid_simplex = true;
        for (size_t k = 0; k < simplex.size(); ++k) {
            if (simplex[k] >= points_.size()) {
                std::cout << "DEBUG: Invalid vertex index: " << simplex[k] << " (max: " << points_.size()-1 << ")" << std::endl;
                valid_simplex = false;
                break;
            }
        }
        if (!valid_simplex) continue;
        
        // SciPy式: T_ij = (r_j - r_n)_i の計算
        // T行列を構築（最後の頂点r_nを基準とする）
        std::vector<std::vector<double>> T(ndim_, std::vector<double>(ndim_, 0.0));
        const auto& r_n = points_[simplex[ndim_]]; // 最後の頂点
        std::cout << "DEBUG: r_n = (" << r_n[0] << ", " << r_n[1] << ")" << std::endl;
        
        for (size_t j = 0; j < ndim_; ++j) {
            const auto& r_j = points_[simplex[j]]; // j番目の頂点
            std::cout << "DEBUG: r_" << j << " = (" << r_j[0] << ", " << r_j[1] << ")" << std::endl;
            for (size_t i = 0; i < ndim_; ++i) {
                T[i][j] = r_j[i] - r_n[i]; // T_ij = (r_j - r_n)_i
            }
        }
        
        std::cout << "DEBUG: About to call invertMatrix..." << std::endl;
        // T行列の逆行列を計算（LU分解使用）
        std::vector<std::vector<double>> T_inv;
        try {
            T_inv = invertMatrix(T);
            std::cout << "DEBUG: invertMatrix completed successfully" << std::endl;
        } catch (const std::exception& e) {
            std::cout << "DEBUG: invertMatrix failed: " << e.what() << std::endl;
            continue;
        }
        
        // transform_配列への書き込み（SciPy形式）
        // Tinvs[i,:ndim,:ndim] = T^-1, Tinvs[i,ndim,:] = r_n
        const size_t base_offset = isimplex * (ndim_ + 1) * ndim_;
        
        // T^-1部分の書き込み
        for (size_t row = 0; row < ndim_; ++row) {
            for (size_t col = 0; col < ndim_; ++col) {
                const_cast<std::vector<double>&>(transform_)[base_offset + row * ndim_ + col] = T_inv[row][col];
            }
        }
        
        // r_n部分の書き込み
        for (size_t col = 0; col < ndim_; ++col) {
            const_cast<std::vector<double>&>(transform_)[base_offset + ndim_ * ndim_ + col] = r_n[col];
        }
        
        std::cout << "DEBUG: Simplex " << isimplex << " processing completed" << std::endl;
    }
    
    std::cout << "DEBUG: calculateTransformMatrix() end" << std::endl;
    const_cast<bool&>(transform_computed_) = true;
}

std::vector<std::vector<double>> Delaunay::invertMatrix(const std::vector<std::vector<double>>& matrix) const {
    const size_t n = matrix.size();
    
    // 単位行列を作成
    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        result[i][i] = 1.0;
    }
    
    // 作業用の拡大行列を作成 [A|I]
    std::vector<std::vector<double>> augmented(n, std::vector<double>(2 * n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmented[i][j] = matrix[i][j];
            augmented[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // Gauss-Jordan elimination
    for (size_t i = 0; i < n; ++i) {
        // ピボット要素を1にスケーリング
        double pivot = augmented[i][i];
        if (std::abs(pivot) < 1e-14) {
            throw std::runtime_error("Matrix is singular, cannot compute inverse");
        }
        
        for (size_t j = 0; j < 2 * n; ++j) {
            augmented[i][j] /= pivot;
        }
        
        // 他の行を0にする
        for (size_t k = 0; k < n; ++k) {
            if (k != i) {
                double factor = augmented[k][i];
                for (size_t j = 0; j < 2 * n; ++j) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
    }
    
    // 逆行列部分を抽出
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result[i][j] = augmented[i][j + n];
        }
    }
    
    return result;
}

} // namespace qhull
