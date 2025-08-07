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

const std::vector<std::vector<int>>& Delaunay::getSimplices() const {
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
    // 重心座標配列のサイズを(ndim+1)に設定（シンプレックスの頂点数）
    barycentric_coords.resize(ndim_ + 1);
    
    // 点が凸包の境界ボックス外側にある場合は早期リターン
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

    // 点をパラボロイド変換（n次元→n+1次元への射影）
    std::vector<double> lifted_point;
    liftPoint(point, lifted_point);

    // 探索開始シンプレックスから点までの距離を計算
    double best_dist = distplane(isimplex, lifted_point);

    // 隣接シンプレックスを辿りながら最適なシンプレックスを探索
    bool changed = true;
    while (changed) {
        // 正の距離を持つシンプレックスが見つかったら探索終了
        // （点がシンプレックスの内部または境界上にある）
        if (best_dist > 0) {
            break; 
        }

        // この反復で変更があったかどうかのフラグをリセット
        changed = false;

        // 現在のシンプレックスの全ての隣接シンプレックスを調査
        for (size_t k = 0; k < ndim_ + 1; ++k) {
            
            // 配列境界の安全性チェック
            if (isimplex >= static_cast<int>(neighbors_.size())) {
                continue;
            }
            
            // 隣接関係配列の境界チェック
            if (k >= neighbors_[isimplex].size()) {
                continue;
            }
            
            // SciPy互換: 単体isimplexのk番目の隣接単体を取得
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

    qhull.triangulate();

    // 2. パラボロイド変換パラメータの取得
    auto ps = qhull.getParaboloidShiftScale();
    paraboloid_shift_ = ps.first;
    paraboloid_scale_ = ps.second;

    // 3. 幾何データの抽出
    std::tie(simplices_, neighbors_, equations_, coplanar_, good_) = qhull.getSimplexFacetArray();

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
    
    if (points_.empty() || points_[0].size() != ndim_) {
        return;
    }

    // 初期値を最初の点で設定
    min_bound_ = points_[0];
    max_bound_ = points_[0];
    
    // 2番目以降の点で比較して更新
    for (size_t i = 1; i < npoints_; ++i) {
        if (i >= points_.size() || points_[i].size() != ndim_) {
            continue;
        }
        
        for (size_t j = 0; j < ndim_; ++j) {
            min_bound_[j] = std::min(min_bound_[j], points_[i][j]);
            max_bound_[j] = std::max(max_bound_[j], points_[i][j]);
        }
    }    
}

bool Delaunay::isPointFullyOutside(const std::vector<double>& point, double eps) const noexcept 
{    
    // 入力パラメータの次元数整合性チェック
    if (point.size() != ndim_ || min_bound_.size() != ndim_ || max_bound_.size() != ndim_) {
        return true; // 安全のためtrueを返す（次元不整合は外側扱い）
    }
    
    try {
        // 各次元について境界内外判定を実行
        for (size_t j = 0; j < ndim_; ++j) {            
            // 配列アクセスの安全性チェック（二重チェック）
            if (j >= point.size() || j >= min_bound_.size() || j >= max_bound_.size()) {
                return true;
            }
            
            // 各次元の座標値と境界値を取得
            double point_val = point[j];
            double min_val = min_bound_[j];
            double max_val = max_bound_[j];
            
            // 点の座標が境界の最小値より小さいか、最大値より大きい場合
            // eps分のマージンを考慮して判定（数値誤差対応）
            if (point_val < min_val - eps || point_val > max_val + eps) {
                return true; // この次元で境界外なので全体として境界外
            }
        }
    } catch (...) {
        // 例外発生時は安全のため境界外として扱う
        return true;
    }
    
    // 全ての次元で境界内であれば境界内として判定
    return false;
}

void Delaunay::liftPoint(const std::vector<double>& point, std::vector<double>& lifted_point) const noexcept 
{
    // lifted_pointを(ndim+1)次元に拡張（最後の次元はパラボロイド座標用）
    lifted_point.resize(ndim_ + 1);
    
    // 点の各座標の二乗和を計算（パラボロイド変換のため）
    double sum_sq = 0;
    for (size_t j = 0; j < ndim_; ++j) {
        // 元の座標をそのままコピー
        lifted_point[j] = point[j];
        // 各座標の二乗を累積
        sum_sq += point[j] * point[j];
    }
    // パラボロイド変換の計算: z = scale * ||x||^2 + shift
    // この変換により、n次元の点を(n+1)次元のパラボロイド面上の点に射影
    lifted_point[ndim_] = sum_sq * paraboloid_scale_ + paraboloid_shift_;
}

double Delaunay::distplane(int simplex_index, const std::vector<double>& lifted_point) const noexcept 
{
    // 境界チェック
    if (simplex_index < 0 || simplex_index >= static_cast<int>(equations_.size())) {
        return 0.0; // 安全な値を返す
    }
    
    const auto& eq = equations_[simplex_index];
    
    if (eq.size() < ndim_ + 2) {
        return 0.0; // 不正なサイズの場合は安全な値を返す
    }
    
    // 超平面方程式の定数項（オフセット項）で距離計算を初期化
    double dist = eq[ndim_ + 1]; // オフセット項
    
    // 超平面方程式: a₀x₀ + a₁x₁ + ... + aₙxₙ + aₙ₊₁ = 0 の計算
    for (size_t k = 0; k < ndim_ + 1; ++k) {
        if (k < lifted_point.size() && k < eq.size()) {
            // 各次元の係数と座標の積を累積
            dist += eq[k] * lifted_point[k];
        }
    }
    
    // 点から超平面までの符号付き距離を返す
    // 正の値: 超平面の法線方向側, 負の値: 反対側
    return dist;
}

int Delaunay::findSimplexDirected(
    std::vector<double>& barycentric_coords,
    const std::vector<double>& point,
    int& start_simplex_hint,
    double eps, 
    double eps_broad) const noexcept 
{
    // 探索開始点として与えられたヒントを使用
    int isimplex = start_simplex_hint;

    if (isimplex < 0 || isimplex >= static_cast<int>(nsimplex_)) {
        isimplex = 0;  // 無効なヒントの場合は0番目の単体から開始
    }

    // Scipy: for cycle_k in range(1 + d.nsimplex//4):
    // 無限ループを防ぐための最大試行回数（単体数の1/4 + 1回）
    // これはSciPyの実装と同じ制限値
    const int max_cycles = 1 + nsimplex_ / 4;
    int cycle_k;
    
    // Directed search（方向探索）アルゴリズムのメインループ
    for (cycle_k = 0; cycle_k < max_cycles; ++cycle_k) {
        
        // 単体が見つからない場合（-1）は探索終了
        if (isimplex == -1) {
            break;
        }

        // 現在の単体に対する重心座標変換行列を取得
        std::vector<double> transform_matrix_for_simplex;
        
        try {
            transform_matrix_for_simplex = getBarycentricTransforms(isimplex);
        } catch (const std::exception& e) {
            // 特異行列の場合は総当たり検索にフォールバック
            // 単体が縮退している可能性がある
            return findSimplexBruteforce(barycentric_coords, point, eps, eps_broad);
        }

        // 単体内部判定のステータス
        int inside_status = 1; // 1: inside, -1: hopped, 0: outside/degenerate

        // 各頂点に対する重心座標を計算し、点が単体内部にあるかチェック
        for (size_t k = 0; k < ndim_ + 1; ++k) {
            // k番目の重心座標を計算
            barycentricCoordinateSingle(transform_matrix_for_simplex, point, barycentric_coords, k);

            // 重心座標が負の場合、点はk番目の面の外側にある
            if (barycentric_coords[k] < -eps) {
                // 境界チェック：隣接単体のデータが有効か確認
                if (isimplex >= static_cast<int>(neighbors_.size()) || 
                    k >= neighbors_[isimplex].size()) {
                    start_simplex_hint = isimplex;
                    return -1;  // データ構造エラー
                }
                
                // SciPy互換: directed search中の隣接単体への移動
                // k番目の面に隣接する単体に移動
                int m = neighbors_[isimplex][k];
                if (m == -1) {
                    start_simplex_hint = isimplex;
                    return -1; // 凸包の外（隣接単体が存在しない）
                }

                // 隣接単体に移動
                isimplex = m;
                inside_status = -1; // 隣へホップした
                break;  // 他の座標は計算せずに次のサイクルへ

            } else if (barycentric_coords[k] <= 1.0 + eps) {
                // OK, 重心座標が有効範囲内（まだ内部にいる可能性がある）
                // 何もしない、次の座標をチェック
            } else {
                // 範囲外 or NaN (縮退した単体)
                // 重心座標が1を大きく超える場合は数値的に不安定
                inside_status = 0;
            }
        }

        // ステータスに応じて次の処理を決定
        if (inside_status == -1) {
            continue; // 隣接単体に移動したので次のサイクルへ

        } else if (inside_status == 1) {
            // 全ての重心座標が有効範囲内
            break; // 正しい単体を発見！ループを抜ける

        } else {
            // 縮退などで失敗 -> 総当たりへ
            // 数値的に不安定な場合は安全な総当たり検索を使用
            isimplex = findSimplexBruteforce(barycentric_coords, point, eps, eps_broad);
            break;
        }
    }

    // for-else: ループが正常に終了しなかった（収束しなかった）場合
    // C++では `cycle_k == max_cycles` で判定
    // SciPyのPythonコードのfor-else文に相当する処理
    if (cycle_k == max_cycles && isimplex != -1) {
        // 最大試行回数に達しても収束しなかった場合は総当たり検索
        isimplex = findSimplexBruteforce(barycentric_coords, point, eps, eps_broad);
    }

    // 次回の探索のヒントとして現在の単体を保存
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
    // 最後の重心座標は他の座標から計算される（重心座標の和は1という制約を利用）
    if (i == static_cast<int>(ndim_)) {
        // c[ndim] = 1.0 - sum(c[0...ndim-1])
        // 最後の座標 = 1 - (他の全ての座標の合計)
        // これにより、全ての重心座標の和が1になることが保証される
        barycentric_coords[ndim_] = 1.0;

        // 既に計算された他の重心座標を減算
        for (size_t j = 0; j < ndim_; ++j) {
            barycentric_coords[ndim_] -= barycentric_coords[j];
        }

    } else {
        // c[i] = sum_j( T_inv[i,j] * (x[j] - r[j]) )
        // i番目の重心座標を変換行列を使って直接計算
        // T_inv: 単体の逆変換行列、r: 参照点（通常は単体の最初の頂点）
        barycentric_coords[i] = 0.0;
        
        // 各次元について変換を適用
        for (size_t j = 0; j < ndim_; ++j) {
            // transform_matrix_for_simplex の構造:
            // - [0 ... ndim*ndim-1]: 逆変換行列 T_inv (ndim x ndim)
            // - [ndim*ndim ... ndim*ndim+ndim-1]: 参照点 r (ndim要素)
            
            // transform[ndim*i + j] は T_inv[i,j] - 逆変換行列の(i,j)成分
            // transform[ndim*ndim + j] は r[j] - 参照点のj成分
            barycentric_coords[i] 
                += transform_matrix_for_simplex[ndim_ * i + j] * (point[j] - transform_matrix_for_simplex[ndim_ * ndim_ + j]);
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
    // 点が凸包の外部にある場合は早期リターン
    if (isPointFullyOutside(point, eps)) { // eps_broadではなくepsを使うのがScipyのコード通り
        return -1;
    }

    // for isimplex in range(d.nsimplex):
    // 全ての単体を順番にチェックする総当たり検索
    for (int isimplex = 0; isimplex < static_cast<int>(nsimplex_); ++isimplex) {

        // NOTE: `this->transform` は (nsimplex, ndim+1, ndim) の3D配列と見なせる
        // フラットな `std::vector<double>` としてメンバに保持する必要がある。
        // `transform_matrix_for_simplex` はその一部を指すビューまたはコピー。
        std::vector<double> transform_matrix_for_simplex;
        try {
            // 現在の単体の変換行列を取得
            transform_matrix_for_simplex = getBarycentricTransforms(isimplex); // ヘルパー関数を仮定
        } catch (const std::exception& e) {
            continue; // この単体をスキップして次へ
        }

        // if transform[0] == transform[0]:
        // 変換行列が有効（NaNでない）かチェック
        if (!std::isnan(transform_matrix_for_simplex[0])) {
            // transform is valid (non-nan)
            // 通常の単体：重心座標を計算して内部判定
            if (barycentricInside(transform_matrix_for_simplex, point, barycentric_coords, eps)) {
                return isimplex;
            }

        } else {
            // transform is invalid (nan, implying degenerate simplex)
            // 縮退した単体：隣接単体を調べる特別処理
            // check neighbors
            for (size_t k = 0; k < ndim_ + 1; ++k) {
                // 境界チェック
                if (isimplex >= static_cast<int>(neighbors_.size()) || 
                    k >= neighbors_[isimplex].size()) {
                    continue;
                }
                
                // SciPy互換: 縮退単体の隣接単体を確認
                int ineighbor = neighbors_[isimplex][k];
                
                if (ineighbor == -1) {
                    continue; // 隣接単体が存在しない
                }

                // 隣接単体の変換行列を取得
                std::vector<double> neighbor_transform;
                try {
                    neighbor_transform = getBarycentricTransforms(ineighbor);
                } catch (const std::exception& e) {
                    continue;
                }

                if (std::isnan(neighbor_transform[0])) {
                    continue; // another bad simplex
                }

                // 隣接単体での重心座標を計算
                barycentricCoordinates(neighbor_transform, point, barycentric_coords);

                // 隣接単体内部判定（特別な許容値を使用）
                bool inside_neighbor = true;
                for (size_t m = 0; m < ndim_ + 1; ++m) {
                    // 境界チェック
                    if (ineighbor >= static_cast<int>(neighbors_.size()) || 
                        m >= neighbors_[ineighbor].size()) {
                        inside_neighbor = false;
                        break;
                    }
                    
                    // SciPy互換: 隣接単体が元の単体を参照しているかチェック
                    if (neighbors_[ineighbor][m] == isimplex) {
                        // 元の縮退単体に向かう方向では緩い許容値を使用
                        // allow extra leeway towards isimplex
                        if (!(barycentric_coords[m] >= -eps_broad && barycentric_coords[m] <= 1.0 + eps)) {
                            inside_neighbor = false;
                            break;
                        }

                    } else {
                        // 通常の方向では標準の許容値を使用
                        // normal check
                        if (!(barycentric_coords[m] >= -eps && barycentric_coords[m] <= 1.0 + eps)) {
                            inside_neighbor = false;
                            break;
                        }
                    }
                }
                if (inside_neighbor) {
                    return ineighbor; // 隣接単体が答え
                }
            }
        }
    }
    return -1; // どの単体にも含まれない
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
    // 最後の重心座標を1で初期化（他の座標を引いて最終値を求める）
    barycentric_coords[ndim_] = 1.0;

    // for i in range(ndim):
    // 各次元の重心座標を順次計算し、同時に有効性をチェック
    for (size_t i = 0; i < ndim_; ++i) {
        // c[i] = 0
        barycentric_coords[i] = 0.0;
        
        // for j in range(ndim):
        // 変換行列を使用してi番目の重心座標を計算
        for (size_t j = 0; j < ndim_; ++j) {
            // transform_matrix[ndim*i + j]: 変換行列T_inv[i,j]
            // transform_matrix[ndim*ndim + j]: 参照点r[j]
            barycentric_coords[i] += transform_matrix[ndim_ * i + j] * (point[j] - transform_matrix[ndim_ * ndim_ + j]);
        }
        
        // c[ndim] -= c[i]
        // 最後の座標から現在の座標を減算（重心座標の和=1制約を維持）
        barycentric_coords[ndim_] -= barycentric_coords[i];

        // if not (-eps <= c[i] <= 1 + eps): return 0
        // 早期リターン：i番目の座標が有効範囲外なら即座にfalseを返す
        if (!(barycentric_coords[i] >= -eps && barycentric_coords[i] <= 1.0 + eps)) {
            return false;
        }
    }

    // if not (-eps <= c[ndim] <= 1 + eps): return 0
    // 最後の重心座標も有効範囲内かチェック
    if (!(barycentric_coords[ndim_] >= -eps && barycentric_coords[ndim_] <= 1.0 + eps)) {
        return false;
    }

    // return 1
    // 全ての重心座標が有効範囲内：点は単体内部にある
    return true;
}

void Delaunay::barycentricCoordinates(const std::vector<double>& transform_matrix,
    const std::vector<double>& point,
    std::vector<double>& barycentric_coords
) const noexcept {

    barycentric_coords.resize(ndim_ + 1);

    // c[ndim] = 1.0
    // 最後の重心座標を1で初期化（他の座標を減算して最終値を求める）
    barycentric_coords[ndim_] = 1.0;

    // for i in range(ndim):
    // 各次元の重心座標を順次計算
    for (size_t i = 0; i < ndim_; ++i) {
        // c[i] = 0
        barycentric_coords[i] = 0.0;
        
        // for j in range(ndim):
        // 変換行列を使用してi番目の重心座標を計算
        for (size_t j = 0; j < ndim_; ++j) {
            // c[i] += transform[ndim*i + j] * (x[j] - transform[ndim*ndim + j])
            // transform_matrix[ndim*i + j]: 変換行列T_inv[i,j]
            // transform_matrix[ndim*ndim + j]: 参照点r[j]
            barycentric_coords[i] += transform_matrix[ndim_ * i + j] * (point[j] - transform_matrix[ndim_ * ndim_ + j]);
        }
        
        // c[ndim] -= c[i]
        // 重心座標の和=1制約を満たすため、最後の座標から現在の座標を減算
        barycentric_coords[ndim_] -= barycentric_coords[i];
    }
}

std::vector<double> Delaunay::getBarycentricTransforms(int simplex_index) const {
    if (simplex_index < 0 || simplex_index >= static_cast<int>(nsimplex_)) {
        // 無効なインデックスの場合はNaNで満たされた配列を返す
        const size_t simplex_size = (ndim_ + 1) * ndim_;
        return std::vector<double>(simplex_size, std::numeric_limits<double>::quiet_NaN());
    }
    
    if (!transform_computed_) {
        calculateTransformMatrix();
    }
    
    const size_t simplex_size = (ndim_ + 1) * ndim_;
    const size_t start_offset = simplex_index * simplex_size;
    
    // 範囲チェック
    if (start_offset + simplex_size > transform_.size()) {
        // 範囲外の場合はNaNで満たされた配列を返す
        return std::vector<double>(simplex_size, std::numeric_limits<double>::quiet_NaN());
    }
    
    return std::vector<double>(
        transform_.begin() + start_offset, 
        transform_.begin() + start_offset + simplex_size);
}

void Delaunay::calculateTransformMatrix() const {
    if (transform_computed_) return;
    
    const size_t nsimplex = simplices_.size();
    
    if (nsimplex == 0 || ndim_ == 0) {
        const_cast<bool&>(transform_computed_) = true;
        return;
    }
    
    try {
        const size_t transform_size = nsimplex * (ndim_ + 1) * ndim_;
        const_cast<std::vector<double>&>(transform_).resize(transform_size, std::numeric_limits<double>::quiet_NaN());
        
        // 各単体について変換行列を計算
        for (size_t isimplex = 0; isimplex < nsimplex; ++isimplex) {
            
            if (isimplex >= simplices_.size()) {
                continue;
            }
            
            const auto& simplex = simplices_[isimplex];
            
            if (simplex.size() != ndim_ + 1) {
                continue;
            }
            
            // 頂点インデックスの境界チェック
            bool valid_simplex = true;
            for (size_t k = 0; k < simplex.size(); ++k) {
                if (simplex[k] < 0 || simplex[k] >= static_cast<int>(points_.size())) {
                    valid_simplex = false;
                    break;
                }
            }
            if (!valid_simplex) continue;
            
            // SciPy式: T_ij = (r_j - r_n)_i の計算
            // T行列を構築（最後の頂点r_nを基準とする）
            std::vector<std::vector<double>> T(ndim_, std::vector<double>(ndim_, 0.0));
            const auto& r_n = points_[simplex[ndim_]]; // 最後の頂点
            
            if (r_n.size() != ndim_) {
                continue;
            }
            
            for (size_t j = 0; j < ndim_; ++j) {
                const auto& r_j = points_[simplex[j]]; // j番目の頂点
                
                if (r_j.size() != ndim_) {
                    valid_simplex = false;
                    break;
                }
                
                for (size_t i = 0; i < ndim_; ++i) {
                    T[i][j] = r_j[i] - r_n[i]; // T_ij = (r_j - r_n)_i
                }
            }
            
            if (!valid_simplex) continue;
            
            // T行列の逆行列を計算（LU分解使用）
            std::vector<std::vector<double>> T_inv;
            try {
                T_inv = invertMatrix(T);
            } catch (const std::exception& e) {
                // 特異行列の場合は、Pseudoinverse（擬似逆行列）を試行する
                // 退化したsimplexでも、可能な限り妥当な変換行列を作成
                try {
                    T_inv = pseudoInverseMatrix(T);
                } catch (const std::exception& pe) {
                    // 最後の手段として、単位行列ベースのフォールバック
                    T_inv = std::vector<std::vector<double>>(ndim_, std::vector<double>(ndim_, 0.0));
                    for (size_t i = 0; i < ndim_; ++i) {
                        T_inv[i][i] = 1.0; // 単位行列
                    }
                }
            }
            
            // transform_配列への書き込み（SciPy形式）
            // Tinvs[i,:ndim,:ndim] = T^-1, Tinvs[i,ndim,:] = r_n
            const size_t base_offset = isimplex * (ndim_ + 1) * ndim_;
            
            if (base_offset + (ndim_ + 1) * ndim_ > transform_.size()) {
                continue;
            }
            
            // T^-1部分の書き込み
            for (size_t row = 0; row < ndim_; ++row) {
                for (size_t col = 0; col < ndim_; ++col) {
                    if (row < T_inv.size() && col < T_inv[row].size()) {
                        const_cast<std::vector<double>&>(transform_)[base_offset + row * ndim_ + col] = T_inv[row][col];
                    }
                }
            }
            
            // r_n部分の書き込み
            for (size_t col = 0; col < ndim_; ++col) {
                if (col < r_n.size()) {
                    const_cast<std::vector<double>&>(transform_)[base_offset + ndim_ * ndim_ + col] = r_n[col];
                }
            }            
        }
        
    } catch (const std::exception& e) {
    }
    
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

std::vector<std::vector<double>> Delaunay::pseudoInverseMatrix(const std::vector<std::vector<double>>& matrix) const {
    // 簡易的な擬似逆行列の実装
    // A^+ = (A^T * A)^-1 * A^T (左擬似逆行列)
    // ただし、(A^T * A)が特異の場合は Moore-Penrose 擬似逆行列を近似
    
    const size_t m = matrix.size();        // 行数
    const size_t n = matrix[0].size();     // 列数
        
    // A^T (転置行列) を計算
    std::vector<std::vector<double>> AT(n, std::vector<double>(m, 0.0));
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            AT[j][i] = matrix[i][j];
        }
    }
    
    // A^T * A を計算
    std::vector<std::vector<double>> ATA(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < m; ++k) {
                ATA[i][j] += AT[i][k] * matrix[k][j];
            }
        }
    }
    
    // (A^T * A)の逆行列を計算（できれば）
    std::vector<std::vector<double>> ATA_inv;
    try {
        ATA_inv = invertMatrix(ATA);
    } catch (const std::exception& e) {
        // 正則化: A^T * A + λI を計算
        const double lambda = 1e-6; // 正則化パラメータ
        for (size_t i = 0; i < n; ++i) {
            ATA[i][i] += lambda;
        }
        try {
            ATA_inv = invertMatrix(ATA);
        } catch (const std::exception& e2) {
            // 最小限の解として、対角要素のみ設定
            ATA_inv = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
            for (size_t i = 0; i < n; ++i) {
                ATA_inv[i][i] = 1.0 / std::max(1e-12, std::abs(ATA[i][i]));
            }
        }
    }
    
    // A^+ = (A^T * A)^-1 * A^T を計算
    std::vector<std::vector<double>> pseudo_inv(n, std::vector<double>(m, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            for (size_t k = 0; k < n; ++k) {
                pseudo_inv[i][j] += ATA_inv[i][k] * AT[k][j];
            }
        }
    }
    
    return pseudo_inv;
}

} // namespace qhull
