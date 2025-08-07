#include <vector>
#include <iostream>
#include "Qhull.h"

namespace qhull {

/**
 * @class Delaunay
 * @brief N次元点群に対するドロネー三角形分割を実行するクラス
 * 
 * このクラスはQhullライブラリを使用して、任意の次元の点群に対してドロネー三角形分割を実行します。
 * SciPyの `scipy.spatial.Delaunay` クラスの機能とAPIを忠実に再現しており、
 * 科学計算やデータ解析における空間データの処理に適用できます。
 * 
 * @note このクラスはSciPyの `scipy.spatial.Delaunay` v1.x系の動作を再現します。
 * @note 入力点群は全て同じ次元数を持つ必要があります。
 * @note Qhullライブラリに依存するため、適切にリンクされている必要があります。
 * 
 */
class Delaunay {
public:
    /**
     * @brief ドロネー三角形分割を実行するコンストラクタ
     * 
     * 指定された点群に対してQhullライブラリを使用してドロネー三角形分割を実行します。
     * SciPyの `scipy.spatial.Delaunay` クラスの実装に基づいており、
     * `qhull_options=None` および `incremental=False` の場合の動作を再現します。
     * 
     * @param input_points 三角形分割を行う点群。各点は座標のベクトルとして表現され、
     *                     全ての点は同じ次元数を持つ必要があります。
     * 
     * @note Qhullオプション:
     *       - 基本オプション: "Qbb Qc Qz Q12"
     *       - 5次元以上の場合: "Qx" を追加
     *       - 必須オプション: "Qt" (三角形分割用)
     * 
     * @throws std::runtime_error Qhullの実行に失敗した場合
     * @throws std::invalid_argument 入力点群が無効な場合（空、次元不一致など）
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1870
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1597
     * 
     */
    Delaunay(const std::vector<std::vector<double>>& input_points) 
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

    /**
     * @brief ドロネー三角形分割で得られた単体（simplex）の頂点インデックス配列を取得する
     * 
     * 各単体は入力点群のインデックスで表現された頂点の集合です。
     * N次元空間において、各単体はN+1個の頂点を持ちます（例：2次元では三角形、3次元では四面体）。
     * 
     * @return 単体の配列への const 参照。各要素は単体を構成する頂点のインデックス配列
     *         - simplices_[i] は i 番目の単体
     *         - simplices_[i][j] は i 番目の単体の j 番目の頂点のインデックス（入力点群での位置）
     * 
     * @note この配列は三角形分割の実行時に計算され、以降は変更されません
     * @note 戻り値は const 参照のため、呼び出し元での変更はできません
     * @note SciPyの `scipy.spatial.Delaunay.simplices` プロパティに相当します
     * 
     */
    const std::vector<std::vector<int>>& get_simplices() const {
        return simplices_;
    }

    /**
     * @brief 指定された点を含む単体（simplex）を検索し、重心座標を計算する
     * 
     * 指定された点がどの単体に含まれるかを効率的に検索し、
     * 同時にその点の重心座標（barycentric coordinates）を計算します。
     * 
     * 検索アルゴリズムは2段階で動作します：
     * 1. ウォーキング段階：隣接する単体を順次探索して目標点に近い単体を見つける
     * 2. 指向性検索段階：より精密な検索で実際に点を含む単体を特定する
     * 
     * @param[out] barycentric_coords 計算された重心座標が格納されるベクトル
     *                               点が単体内にある場合、全要素の合計は1.0になります
     *                               サイズは (ndim+1) に設定されます
     * @param[in] point 検索対象の点の座標（ndim次元）
     * @param[in,out] start_simplex_hint 検索開始位置のヒント
     *                                  入力：検索開始する単体のインデックス（-1の場合は0から開始）
     *                                  出力：実際に検索を開始した単体のインデックス
     * @param[in] eps 数値計算の許容誤差（通常は1e-12程度）
     *                重心座標の計算や境界判定で使用されます
     * @param[in] eps_broad 粗い境界判定の許容誤差（通常はepsより大きい値）
     *                      点が明らかに凸包の外側にある場合の早期判定に使用されます
     * 
     * @return 点を含む単体のインデックス（0以上）。点が三角形分割の外側にある場合は-1
     * 
     * @note このメソッドはSciPyの内部実装 `_find_simplex` のロジックを忠実に再現しています
     * @note 点が単体の境界上にある場合、その境界を共有する単体のいずれかが返される可能性があります
     * @note メソッドはnoexcept指定されており、例外を投げません
     * @note 三角形分割が空の場合（nsimplex <= 0）は-1を返します
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1474
     * 
     */
    int findSimplex(
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


private:
    /**
     * @brief Qhullの計算結果からDelaunayオブジェクトの内部状態を更新する
     * 
     * このメソッドはSciPyの `Delaunay._update` メソッドと `_QhullUser._update` メソッドの
     * 機能を統合したものです。Qhullライブラリによる三角形分割の実行結果を取得し、
     * 必要なデータ構造を構築してメンバ変数に格納します。
     * 
     * ## 実行される処理の詳細
     * 
     * ### Delaunay._update 相当の処理:
     * 1. **三角形分割の実行**: `qhull.triangulate()` を呼び出してドロネー三角形分割を実行
     * 2. **パラボロイド変換パラメータの取得**: 高次元での数値安定性のための変換係数を取得
     * 3. **幾何データの抽出**: 単体、隣接関係、方程式、共面点、有効フラグを一括取得
     * 4. **単体数の記録**: 生成された単体の総数を記録
     * 5. **遅延計算フラグの初期化**: 変換行列や頂点関係などの重い計算を必要時まで遅延
     * 
     * ### _QhullUser._update 相当の処理:
     * 6. **境界計算**: 点群の各次元における最小値・最大値を計算
     * 
     * @param qhull 三角形分割が設定済みのQhullオブジェクトの参照。
     *              このオブジェクトは適切に初期化され、点群データが設定されている必要があります。
     * 
     * @pre qhullオブジェクトは有効な点群データで初期化されている
     * @pre qhullオブジェクトのコマンドは "d" (delaunay) に設定されている
     * 
     * @post すべての基本的な幾何データ（単体、隣接関係等）がメンバ変数に格納される
     * @post 遅延計算フラグがfalseに初期化される
     * @post 点群の境界情報が計算される
     * 
     * @throws std::runtime_error Qhullの三角形分割実行に失敗した場合
     * @throws std::bad_alloc メモリ不足の場合
     * 
     * @note このメソッドはコンストラクタから呼び出されるため、直接呼び出す必要はありません
     * @note SciPyとの互換性のため、処理順序と内容を忠実に再現しています
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1895
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1625
     */
    void update(Qhull& qhull) 
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
    
    /**
     * @brief 全ての点の座標における各次元の最小値と最大値を計算する
     * 
     * このメソッドは、入力された点群から各次元における境界（bounding box）を計算し、
     * min_bound_とmax_bound_メンバ変数に結果を格納します。
     * 
     * @details 
     * - 点が存在しない場合（npoints_ == 0）は何も実行しません
     * - 最初の点を初期値として設定し、残りの点と比較して更新します
     * - 各次元について独立して最小値・最大値を求めます
     * 
     * @note このメソッドはDelaunay三角分割の前処理として使用され、
     *       数値的安定性の向上に寄与します
     */
    void calculateBounds() {
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

    /**
     * @brief 指定された点が凸包の境界ボックスの外側に完全に位置するかを判定する
     * 
     * このメソッドは、指定された点が三角形分割によって形成される凸包の
     * 軸並行境界ボックス（axis-aligned bounding box）の外側に完全に位置するかを
     * 効率的に判定します。これは、より重い幾何学的計算を実行する前の
     * 早期除外（early rejection）に使用されます。
     * 
     * 各次元において、点の座標が境界の最小値より小さいか、
     * 最大値より大きい場合（許容誤差を含む）、その点は境界ボックスの外側にあると判定されます。
     * いずれかの次元で外側にある場合、点は完全に外側にあると見なされます。
     * 
     * @param[in] point 判定対象の点の座標（ndim次元）
     * @param[in] eps 境界判定の許容誤差
     *                境界から eps だけ外側までは「内側」として扱われます
     * 
     * @return true: 点が境界ボックスの外側に完全に位置する場合
     * @return false: 点が境界ボックス内または境界上（許容誤差内）にある場合
     * 
     * @note このメソッドはSciPyの `_is_point_fully_outside` 関数と同等の機能を提供します
     * @note 境界ボックスは min_bound_ と max_bound_ によって定義されます
     * @note 計算複雑度は O(ndim) で、非常に高速です
     * @note メソッドはnoexcept指定されており、例外を投げません
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1302
     * 
     */
    bool isPointFullyOutside(const std::vector<double>& point, double eps) const noexcept 
    {
        for (size_t j = 0; j < ndim_; ++j) {
            // 点の座標が境界の最小値より小さいか、最大値より大きい場合
            if (point[j] < min_bound_[j] - eps || point[j] > max_bound_[j] + eps) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief 点をパラボロイド（放物面）に持ち上げる変換を行う
     * 
     * このメソッドは、n次元空間の点を(n+1)次元のパラボロイド上の点に変換します。
     * この変換は、ドロネー三角形分割において凸包の計算を効率化するための
     * 標準的な手法で、「リフティング変換」と呼ばれます。
     * 
     * ## 変換の数学的定義
     * 
     * n次元の点 `x = (x₀, x₁, ..., xₙ₋₁)` を(n+1)次元の点 `z` に変換：
     * 
     * ```
     * z[i] = x[i]                                    (i = 0, 1, ..., n-1)
     * z[n] = (Σᵢ x[i]²) × paraboloid_scale + paraboloid_shift
     * ```
     * 
     * この変換により、元の空間での点の関係が高次元空間での凸包の性質として
     * 表現され、複雑な幾何学的判定を線形代数の問題に帰着できます。
     * 
     * ## パラメータの役割
     * 
     * - **paraboloid_scale**: 数値安定性のためのスケーリング係数
     * - **paraboloid_shift**: 数値精度向上のためのオフセット値
     * 
     * これらの値は三角形分割の初期化時にQhullライブラリによって計算されます。
     * 
     * @param[in] point 変換対象のn次元点の座標
     * @param[out] lifted_point 変換結果の(n+1)次元点が格納されるベクトル
     *                         サイズは自動的に(ndim+1)に調整されます
     * 
     * @note このメソッドはSciPyの `_lift_point` 関数と完全に同等の機能を提供します
     * @note 変換は可逆ではありません（高次元から元の次元への逆変換は一意ではない）
     * @note 計算複雑度は O(ndim) です
     * @note メソッドはnoexcept指定されており、例外を投げません
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1277C11-L1277C22
     * 
     */
    void liftPoint(const std::vector<double>& point, std::vector<double>& lifted_point) const noexcept 
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

    /**
     * @brief リフトされた点と単体に対応する超平面との符号付き距離を計算する
     * 
     * このメソッドは、パラボロイド変換によってリフトされた点と、
     * 指定された単体に対応する超平面との符号付き距離を計算します。
     * この距離は、点が単体の「上側」にあるか「下側」にあるかを判定するために使用され、
     * ドロネー三角形分割における点位置判定の核心となる計算です。
     * 
     * ## 数学的定義
     * 
     * 単体の超平面方程式を `ax₀ + bx₁ + ... + cx_{n} + d = 0` とすると、
     * 点 `p = (p₀, p₁, ..., p_{n})` に対する符号付き距離は：
     * 
     * ```
     * distance = a×p₀ + b×p₁ + ... + c×p_{n} + d
     * ```
     * 
     * ## 距離の解釈
     * 
     * - **正の値**: 点が超平面の「上側」（法線ベクトル方向）にある
     * - **負の値**: 点が超平面の「下側」にある  
     * - **ゼロ**: 点が超平面上にある
     * 
     * ドロネー三角形分割では、正の距離を持つ単体が目標点を含む可能性が高く、
     * 単体探索アルゴリズムの方向決定に使用されます。
     * 
     * ## データ構造との対応
     * 
     * - `equations_[simplex_index]` は長さ (ndim+2) の配列
     * - インデックス 0 〜 ndim: 超平面の法線ベクトル成分
     * - インデックス (ndim+1): 定数項（オフセット）
     * 
     * @param[in] simplex_index 距離を計算する単体のインデックス（0 ≤ index < nsimplex）
     * @param[in] lifted_point パラボロイド変換済みの(ndim+1)次元点の座標
     * 
     * @return 符号付き距離値。正の値は点が超平面の上側にあることを示す
     * 
     * @note このメソッドはSciPyの `_distplane` 関数と完全に同等の機能を提供します
     * @note 計算複雑度は O(ndim) です
     * @note メソッドはnoexcept指定されており、例外を投げません
     * @note 入力パラメータの妥当性チェックは行われません（パフォーマンス重視）
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1286
     * 
     */
    double distplane(int simplex_index, const std::vector<double>& lifted_point) const noexcept 
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

    /**
     * @brief 有向探索（directed walk）で点を含む単体を見つける
     * 
     * このメソッドはSciPyの `_find_simplex_directed` 関数と完全に対応しており、
     * 三角分割内を効率的に歩き回って指定された点を含む単体を特定します。
     * 
     * ## アルゴリズムの概要
     * 
     * 1. **有向歩行**: 開始単体から重心座標を計算し、負の座標があればその方向の隣接単体へ移動
     * 2. **収束判定**: 全ての重心座標が有効範囲内にあれば、その単体が解
     * 3. **フォールバック**: 収束しない場合は総当たり探索に切り替え
     * 
     * ## 重心座標による移動判定
     * 
     * - `barycentric_coords[k] < -eps`: k番目の隣接単体方向へ移動
     * - `barycentric_coords[k] <= 1 + eps`: 現在の単体内に留まる可能性
     * - それ以外: 縮退した単体として総当たり探索へ
     * 
     * @param barycentric_coords [out] 計算された重心座標を格納する配列
     *                                 サイズは (ndim+1) である必要がある
     * @param point [in] 探索対象の点の座標（ndim次元）
     * @param start_simplex_hint [in/out] 探索開始単体のインデックス。
     *                                   無効な値の場合は0に設定される。
     *                                   探索終了時には最終的に確認した単体が設定される
     * @param eps [in] 重心座標の数値許容誤差（通常の境界判定用）
     * @param eps_broad [in] より緩い数値許容誤差（縮退単体処理用）
     * 
     * @return 点を含む単体のインデックス。点が三角分割の外部にある場合は -1
     * 
     * @pre point.size() == ndim_
     * @pre barycentric_coords.size() >= ndim_ + 1
     * @pre eps >= 0 && eps_broad >= eps
     * 
     * @post start_simplex_hint には最終的に確認した単体インデックスが設定される
     * @post 戻り値が >= 0 の場合、barycentric_coords には有効な重心座標が格納される
     * 
     * @note 最大反復回数は `1 + nsimplex_/4` に制限され、収束しない場合は
     *       自動的に総当たり探索（_find_simplex_bruteforce）にフォールバックする
     * @note このメソッドはSciPyの実装と同じアルゴリズムを使用しているため、
     *       同じ入力に対して同じ結果を返すことが保証される
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1373
     * 
     */
    int findSimplexDirected(
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
    
    /**
     * @brief 単一の重心座標を計算する。
     *        Scipyの `_barycentric_coordinate_single` に完全対応。
     * @param transform_matrix_for_simplex [in] (ndim, ndim+1) の変換行列（フラットな配列）。
     *        transform[0..ndim*ndim-1] は逆行列T_inv、
     *        transform[ndim*ndim..] はベクトルr。
     * @param point [in] ターゲットの点x。
     * @param barycentric_coords [in/out] 重心座標c。i < ndim の成分が計算済みである必要がある。
     * @param coord_index [in] 計算する重心座標のインデックスi。
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1238
     * 
     */
    void barycentricCoordinateSingle(
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

    /**
     * @brief 総当たりで単体を検索する。
     * 
     * 与えられた点がどの単体（simplex）に含まれるかを総当たりで検索します。
     * SciPyの `_find_simplex_bruteforce` アルゴリズムに完全対応した実装です。
     * 
     * @details 
     * このアルゴリズムは以下の手順で動作します：
     * 1. まず、点が三角分割の境界外にあるかチェック
     * 2. 各単体について順次チェック
     * 3. 有効な単体（非退化）の場合、重心座標を計算して内部判定
     * 4. 退化した単体の場合、隣接する単体も考慮して判定
     * 
     * 退化した単体の処理では、隣接単体との境界で緩い許容範囲（eps_broad）を適用し、
     * 数値誤差による誤判定を防ぎます。
     * 
     * @param[out] barycentric_coords 見つかった単体における点の重心座標
     *                               ndim+1個の要素を持つベクトル
     * @param[in] point 検索対象の点座標（ndim次元）
     * @param[in] eps 通常の許容誤差。重心座標の境界判定に使用
     * @param[in] eps_broad 退化単体の隣接境界での緩い許容誤差
     * 
     * @return 点を含む単体のインデックス。見つからない場合は-1
     * 
     * @note 
     * - この関数はtransformメンバ変数が事前に計算されている必要があります
     * - O(n)の計算量を持つため、大量の点に対してはfindSimplexを使用することを推奨
     * - SciPyとの完全互換性を保つため、数値計算の詳細まで同一の実装です
     * 
     * @see findSimplex() より高速な探索アルゴリズム
     * @see isPointFullyOutside() 境界外判定
     * @see _barycentric_inside() 重心座標による内部判定
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1315
     * 
     */
    int findSimplexBruteforce(
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

    /**
     * @brief 点が単体の内部にあるか（重心座標がすべて範囲内か）をチェックする。
     *        Scipyの `_barycentric_inside` の最適化されたループに完全対応。
     * @param transform_matrix [in] (ndim, ndim+1) の変換行列（フラットな配列）。
     * @param point [in] ターゲットの点x。
     * @param barycentric_coords [out] 計算された重心座標cが格納される。
     * @param eps [in] 許容誤差。
     * @return 点が内部にあればtrue, そうでなければfalse。
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1216
     * 
     */
    bool barycentricInside(
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

    /**
     * @brief 全ての重心座標を計算する。
     *        Scipyの `_barycentric_coordinates` の最適化されたループに完全対応。
     * @param transform_matrix [in] (ndim, ndim+1) の変換行列（フラットな配列）。
     * @param point [in] ターゲットの点x。
     * @param barycentric_coords [out] 計算された重心座標cが格納される。
     * 
     * @sa https://github.com/scipy/scipy/blob/8bd39fe64fde804faf28dc29d7e33000bfe45cd1/scipy/spatial/_qhull.pyx#L1258
     * 
     */
    void barycentricCoordinates(const std::vector<double>& transform_matrix,
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

    /**
     * @brief `transform`行列から特定の単体の部分行列へのビューを返すヘルパー。
     *        パフォーマンスのため、実際にはポインタやspanを返すのが望ましい。
     *        ここでは簡単のためコピーを返す。
     */
    std::vector<double> get_transform_for_simplex(int simplex_index) const {
        const size_t simplex_size = (ndim_ + 1) * ndim_;
        const size_t start_offset = simplex_index * simplex_size;
        
        // ダミー実装: transformが計算されていない場合はダミーを返す
        if (transform_.empty()) {
            return std::vector<double>(simplex_size, 0.0);
        }
        
        return std::vector<double>(
            transform_.begin() + start_offset, 
            transform_.begin() + start_offset + simplex_size);
    }

    // --- Private Member Variables ---
    
    // ==========================================
    // `_QhullUser` 由来の属性
    // ==========================================
    
    /**
     * @brief 入力点群の座標データ
     * 
     * 三角形分割を実行する元となる点群の座標を格納します。
     * 各要素 `points_[i]` は i 番目の点の座標ベクトルを表し、
     * 全ての点は同じ次元数 (ndim_) を持ちます。
     * 
     * - `points_[i][j]` = i番目の点のj次元目の座標
     * - サイズ: `npoints_ × ndim_`
     * 
     * SciPyの `scipy.spatial.Delaunay.points` 属性に対応します。
     */
    std::vector<std::vector<double>> points_;
    
    /**
     * @brief 入力点群の次元数
     * 
     * 各点の座標ベクトルの要素数を表します。
     * 2次元の場合は2、3次元の場合は3となります。
     * 
     * 値は最初の点 `points_[0].size()` から決定され、
     * 全ての点が同じ次元数を持つ必要があります。
     * 
     * SciPyの `scipy.spatial.Delaunay.ndim` 属性に対応します。
     */
    size_t ndim_;
    
    /**
     * @brief 入力点群の点の総数
     * 
     * `points_` ベクトルのサイズと等価で、
     * 三角形分割に使用される点の総数を表します。
     * 
     * 値は `points_.size()` から決定されます。
     * 
     * SciPyの `scipy.spatial.Delaunay.npoints` 属性に対応します。
     */
    size_t npoints_;
    
    /**
     * @brief 各次元における点群の最小値
     * 
     * 全ての入力点の各次元における最小座標値を格納します。
     * 境界ボックス（bounding box）の下限を定義し、
     * 点位置判定の早期除外に使用されます。
     * 
     * - `min_bound_[j]` = j次元目の最小座標値
     * - サイズ: `ndim_`
     * 
     * `calculateBounds()` メソッドにより計算されます。
     * SciPyの `scipy.spatial.Delaunay.min_bound` 属性に対応します。
     */
    std::vector<double> min_bound_;
    
    /**
     * @brief 各次元における点群の最大値
     * 
     * 全ての入力点の各次元における最大座標値を格納します。
     * 境界ボックス（bounding box）の上限を定義し、
     * 点位置判定の早期除外に使用されます。
     * 
     * - `max_bound_[j]` = j次元目の最大座標値  
     * - サイズ: `ndim_`
     * 
     * `calculateBounds()` メソッドにより計算されます。
     * SciPyの `scipy.spatial.Delaunay.max_bound` 属性に対応します。
     */
    std::vector<double> max_bound_;

    // ==========================================
    // `Delaunay` 由来の属性
    // ==========================================
    
    /**
     * @brief パラボロイド変換のスケーリング係数
     * 
     * ドロネー三角形分割において、n次元点を(n+1)次元パラボロイドに
     * 変換する際のスケーリング係数です。数値安定性の向上に使用されます。
     * 
     * 変換式: `z[n] = (Σᵢ x[i]²) × paraboloid_scale + paraboloid_shift`
     * 
     * Qhullライブラリによって自動計算され、
     * `liftPoint()` メソッドで使用されます。
     * 
     * SciPyの `scipy.spatial.Delaunay.paraboloid_scale` 属性に対応します。
     */
    double paraboloid_scale_;
    
    /**
     * @brief パラボロイド変換のオフセット値
     * 
     * ドロネー三角形分割において、n次元点を(n+1)次元パラボロイドに
     * 変換する際のオフセット値です。数値精度の向上に使用されます。
     * 
     * 変換式: `z[n] = (Σᵢ x[i]²) × paraboloid_scale + paraboloid_shift`
     * 
     * Qhullライブラリによって自動計算され、
     * `liftPoint()` メソッドで使用されます。
     * 
     * SciPyの `scipy.spatial.Delaunay.paraboloid_shift` 属性に対応します。
     */
    double paraboloid_shift_;
    
    /**
     * @brief 三角形分割によって生成された単体（simplex）の頂点インデックス配列
     * 
     * 各単体は入力点群のインデックスで表現された(ndim+1)個の頂点を持ちます。
     * - 2次元: 三角形（3頂点）
     * - 3次元: 四面体（4頂点）
     * - N次元: N+1頂点
     * 
     * - `simplices_[i]` = i番目の単体
     * - `simplices_[i][j]` = i番目の単体のj番目の頂点のインデックス
     * - サイズ: `nsimplex_ × (ndim_+1)`
     * 
     * Qhullライブラリによって計算され、
     * `get_simplices()` メソッドで取得できます。
     * 
     * SciPyの `scipy.spatial.Delaunay.simplices` 属性に対応します。
     */
    std::vector<std::vector<int>> simplices_;
    
    /**
     * @brief 各単体の隣接単体インデックス配列
     * 
     * 各単体に対して、その面を共有する隣接単体のインデックスを格納します。
     * 境界にある面（凸包の境界）の隣接は -1 で表現されます。
     * 
     * - `neighbors_[i]` = i番目の単体の隣接単体配列
     * - `neighbors_[i][j]` = i番目の単体のj番目の面の隣接単体インデックス
     * - サイズ: `nsimplex_ × (ndim_+1)`
     * 
     * 値 -1 は境界面（隣接する単体が存在しない）を示します。
     * 
     * 単体探索アルゴリズム（`findSimplex()`）で使用され、
     * 効率的な歩行探索を可能にします。
     * 
     * SciPyの `scipy.spatial.Delaunay.neighbors` 属性に対応します。
     */
    std::vector<std::vector<int>> neighbors_;
    
    /**
     * @brief 各単体に対応する超平面の方程式係数
     * 
     * パラボロイド変換後の(ndim+1)次元空間において、
     * 各単体が定義する超平面の方程式 `ax₀ + bx₁ + ... + cx_n + d = 0` の
     * 係数 (a, b, ..., c, d) を格納します。
     * 
     * - `equations_[i]` = i番目の単体の超平面方程式係数
     * - `equations_[i][j]` = j番目の係数（j=0..ndim: 法線成分、j=ndim+1: 定数項）
     * - サイズ: `nsimplex_ × (ndim_+2)`
     * 
     * 点と単体の位置関係判定（`distplane()`）で使用され、
     * 単体探索の核心となる計算に必要です。
     * 
     * SciPyの `scipy.spatial.Delaunay.equations` 属性に対応します。
     */
    std::vector<std::vector<double>> equations_;
    
    /**
     * @brief 各単体に共面な点のインデックス配列
     * 
     * 退化した単体（体積が0の単体）に共面な入力点のインデックスを格納します。
     * 通常のドロネー三角形分割では使用されませんが、
     * 特殊なケースでの完全性のため保持されます。
     * 
     * - `coplanar_[i]` = i番目の単体に共面な点のインデックス配列
     * - サイズ: `nsimplex_ × 可変長`
     * 
     * ほとんどの場合は空配列となります。
     * 
     * SciPyの `scipy.spatial.Delaunay.coplanar` 属性に対応します。
     */
    std::vector<std::vector<int>> coplanar_;
    
    /**
     * @brief 各単体の有効性フラグ
     * 
     * 各単体が有効（非退化）かどうかを示すブール値を格納します。
     * 退化した単体は数値計算で使用できないため、
     * フラグによって区別されます。
     * 
     * - `good_[i]` = i番目の単体が有効かどうか
     * - サイズ: `nsimplex_`
     * 
     * `true`: 有効な単体（正常な体積を持つ）
     * `false`: 退化した単体（体積が0または数値的に不安定）
     * 
     * 単体探索や座標変換で、有効な単体のみを使用するよう制御されます。
     * 
     * SciPyの `scipy.spatial.Delaunay.good` 属性に対応します。
     */
    std::vector<bool> good_;
    
    /**
     * @brief 生成された単体の総数
     * 
     * 三角形分割によって生成された単体（simplex）の総数を表します。
     * `simplices_.size()` と等価です。
     * 
     * この値はアルゴリズムのループ制御や
     * メモリ確保のサイズ計算に使用されます。
     * 
     * SciPyの `scipy.spatial.Delaunay.nsimplex` 属性に対応します。
     */
    size_t nsimplex_;

    // ==========================================
    // 遅延評価される属性の制御フラグ
    // ==========================================
    
    /**
     * @brief 座標変換行列の計算済みフラグ
     * 
     * 重心座標計算に必要な座標変換行列（`transform_`）が
     * 計算済みかどうかを示すフラグです。
     * 
     * `false`: 未計算（必要時に計算される）
     * `true`: 計算済み（`transform_` が有効）
     * 
     * 変換行列は計算コストが高いため、実際に必要になるまで
     * 計算を遅延させる遅延評価パターンを採用しています。
     * 
     * SciPyの lazy evaluation システムに対応します。
     */
    bool transform_computed_;
    
    /**
     * @brief 頂点-単体対応の計算済みフラグ
     * 
     * 各頂点がどの単体に属するかの対応関係が
     * 計算済みかどうかを示すフラグです。
     * 
     * `false`: 未計算（必要時に計算される）
     * `true`: 計算済み（対応データ構造が有効）
     * 
     * この情報は特定の高度な操作で必要になりますが、
     * 基本的な単体探索では使用されないため遅延評価されます。
     * 
     * SciPyの lazy evaluation システムに対応します。
     */
    bool vertex_to_simplex_computed_;
    
    /**
     * @brief 頂点-隣接頂点の計算済みフラグ
     * 
     * 各頂点に隣接する頂点の一覧が
     * 計算済みかどうかを示すフラグです。
     * 
     * `false`: 未計算（必要時に計算される）
     * `true`: 計算済み（隣接関係データが有効）
     * 
     * この情報はグラフ的な解析で必要になりますが、
     * 基本的な三角形分割操作では使用されないため遅延評価されます。
     * 
     * SciPyの lazy evaluation システムに対応します。
     */
    bool vertex_neighbor_vertices_computed_;

    // ==========================================
    // 遅延評価される重い計算データ
    // ==========================================
    
    /**
     * @brief 重心座標計算用の座標変換行列
     * 
     * 各単体において点の重心座標を計算するための変換行列を格納します。
     * 3次元配列 `[nsimplex_][ndim_+1][ndim_]` をフラット配列として保持します。
     * 
     * ## データ構造
     * - サイズ: `nsimplex_ × (ndim_+1) × ndim_`
     * - `transform_[i*(ndim_+1)*ndim_ + j*ndim_ + k]` = 
     *   i番目の単体のj行k列要素
     * 
     * ## 変換行列の構成
     * 各単体の変換行列は以下の形式を持ちます：
     * ```
     * [ T_inv_0  ]     [ r_0 ]
     * [ T_inv_1  ]  +  [ r_1 ]
     * [   ...    ]     [ ... ]
     * [ T_inv_n  ]     [ r_n ]
     * ```
     * 
     * - `T_inv`: ndim×ndim の逆変換行列
     * - `r`: ndim要素のオフセットベクトル
     * 
     * ## 計算コスト
     * この行列の計算は非常に重いため（O(nsimplex × ndim³)）、
     * 実際に重心座標計算が必要になるまで遅延されます。
     * 
     * `transform_computed_` フラグで計算状態を管理し、
     * `barycentricCoordinates()` 系メソッドで使用されます。
     * 
     * SciPyの `scipy.spatial.Delaunay.transform` 属性に対応します。
     */
    std::vector<double> transform_;
};

} // namespace qhull