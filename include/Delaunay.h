#pragma once

#include <vector>
#include <memory>

// Qhullの前方宣言
namespace orgQhull {
    class Qhull;
}

/**
 * @brief Delaunay三角分割を実行するクラス
 * 
 * このクラスはQhullライブラリを使用してDelaunay三角分割を計算し、
 * LinearNdInterpolatorクラスで使用される補間機能を提供します。
 * N次元空間での点群に対してDelaunay三角分割を実行し、
 * 任意の点がどのsimplexに属するか、および重心座標の計算などを行います。
 */
class Delaunay {
public:
    /**
     * @brief Delaunayクラスのコンストラクタ
     * 
     * 与えられた点群に対してDelaunay三角分割を実行します。
     * 点群は N次元空間の点の集合として与えられ、各点は座標値のベクターとして表現されます。
     * 
     * @param points N次元空間の点群。各要素は1つの点の座標を表すdoubleのベクター
     * @throws std::runtime_error Qhullライブラリでエラーが発生した場合
     * @throws std::invalid_argument 点群が空、または次元が統一されていない場合
     */
    explicit Delaunay(const std::vector<std::vector<double>>& points);
    
    /**
     * @brief デストラクタ
     * 
     * Qhullリソースを適切に解放します。
     */
    ~Delaunay();

    /**
     * @brief コピーコンストラクタ（削除済み）
     * 
     * Qhullリソースの複雑な管理を避けるため、コピーは禁止されています。
     */
    Delaunay(const Delaunay&) = delete;
    
    /**
     * @brief コピー代入演算子（削除済み）
     * 
     * Qhullリソースの複雑な管理を避けるため、コピー代入は禁止されています。
     */
    Delaunay& operator=(const Delaunay&) = delete;
    
    /**
     * @brief ムーブコンストラクタ（削除済み）
     * 
     * Qhullリソースの複雑な管理を避けるため、ムーブは禁止されています。
     */
    Delaunay(Delaunay&&) = delete;
    
    /**
     * @brief ムーブ代入演算子（削除済み）
     * 
     * Qhullリソースの複雑な管理を避けるため、ムーブ代入は禁止されています。
     */
    Delaunay& operator=(Delaunay&&) = delete;
    
    /**
     * @brief 指定された点を含むsimplexのIDを検索（SciPy _find_simplex準拠）
     * 
     * 与えられた点がどのsimplex（N次元三角形）に含まれるかを判定します。
     * SciPyの内部C実装 _find_simplex と同等のインターフェースを提供します。
     * 
     * @param barycentric_coords 出力用重心座標（計算結果が格納される）
     * @param point 検索対象の点の座標（N次元ベクター）
     * @param start_simplex 検索開始simplex（参照渡しで更新される）
     * @param eps 基本許容誤差（SciPy準拠: 100*DBL_EPSILON）
     * @param eps_broad broad検索用許容誤差（SciPy準拠: sqrt(DBL_EPSILON)）
     * @return 点を含むsimplexのID（見つからない場合は-1）
     * @throws std::invalid_argument pointの次元が構築時の点群と一致しない場合
     */
    int findSimplex(std::vector<double>& barycentric_coords,
                    const std::vector<double>& point, 
                    int& start_simplex, 
                    double eps, 
                    double eps_broad) const;
    
    // 古いcalculateBarycentricCoordinatesメソッドは削除済み
    // SciPy準拠のcalculateBarycentricCoordinatesWithTransformを使用
    
    /**
     * @brief 全simplexの頂点インデックス情報を取得
     * 
     * Delaunay三角分割で生成された全simplexについて、各simplexを構成する
     * 点のインデックス情報を返します。SciPyのDelaunay.simplicesプロパティと同等です。
     * 
     * @return simplexの配列。各要素は1つのsimplexを構成する点のインデックスのベクター
     */
    std::vector<std::vector<int>> getSimplices() const;
    
    /**
     * @brief バリセントリック変換行列を取得
     * 
     * 各simplexに対するバリセントリック座標変換行列を返します。
     * SciPyのDelaunay.transformプロパティと同等の機能を提供します。
     * 変換行列Tは以下の式を満たします: T * c = x - r_n
     * ここで、cはバリセントリック座標、xは点座標、r_nは参照頂点です。
     * 
     * @return 変換行列の3次元配列 [simplex_id][row][col]
     *         形状: (nsimplex, ndim+1, ndim)
     */
    const std::vector<std::vector<std::vector<double>>>& getTransform() const;
    
    /**
     * @brief 変換行列を使用した高速バリセントリック座標計算
     * 
     * 事前計算された変換行列を使用して効率的にバリセントリック座標を計算します。
     * SciPyのアルゴリズムに基づく実装で、O(d²)の計算複雑度を実現します。
     * 
     * @param point バリセントリック座標を計算する点の座標
     * @param simplex_id 対象となるsimplexのID
     * @return バリセントリック座標のベクター
     */
    std::vector<double> calculateBarycentricCoordinatesWithTransform(
        const std::vector<double>& point, int simplex_id) const;
    
    /**
     * @brief 隣接simplexの情報を取得
     * 
     * SciPyのDelaunay.neighborsプロパティと同等の機能を提供します。
     * 各simplexに隣接するsimplexのIDを返します。
     * 
     * @return 隣接関係の2次元配列 [simplex_id][neighbor_index]
     */
    std::vector<std::vector<int>> getNeighbors() const;
    
    
    /**
     * @brief Brute Force単体検索（SciPy準拠フォールバック）
     * 
     * Walking Algorithmが失敗した場合のフォールバック検索手法です。
     * 全simplexを順次チェックして点を含むsimplexを見つけます。
     * 
     * @param point 検索対象の点の座標
     * @return 点を含むsimplexのID（見つからない場合は-1）
     */
    int findSimplexBruteForceWithEps(const std::vector<double>& point,
                                     double eps_broad,
                                     std::vector<double>& barycentric_coords) const;
    
    /**
     * @brief SciPy準拠のDirected Search実装
     */
    int findSimplexDirected(const std::vector<double>& point,
                           int start_simplex,
                           double eps,
                           double eps_broad,
                           std::vector<double>& barycentric_coords) const;
    
    /**
     * @brief Paraboloid上での平面距離計算（SciPy _distplane完全準拠）
     */
    double calculatePlaneDistance(int simplex_id, const double* lifted_point, size_t dim) const;
    
    /**
     * @brief 点をparaboloid上に投影（SciPy _lift_point相当）
     */
    void liftPointToParaboloid(const std::vector<double>& point, std::vector<double>& lifted_point) const;
    
    /**
     * @brief 単一barycentric座標計算（SciPy _barycentric_coordinate_single相当）
     */
    void calculateBarycentricCoordinateSingle(int simplex_id, const std::vector<double>& point, 
                                             std::vector<double>& barycentric, int coordinate_index) const;

private:
    /**
     * @brief Qhullライブラリのインスタンス
     * 
     * Delaunay三角分割の実際の計算を行うQhullオブジェクトです。
     * unique_ptrによりRAIIによる自動リソース管理が行われます。
     */
    std::unique_ptr<orgQhull::Qhull> qhull_;
    
    /**
     * @brief 構築時に渡された元の点データ
     * 
     * 三角分割の基となったN次元空間の点群データを保持します。
     * 各要素は1つの点の座標を表すdoubleのベクターです。
     */
    std::vector<std::vector<double>> points_;
    
    /**
     * @brief バリセントリック変換行列
     * 
     * 各simplexに対する変換行列を保持します。SciPyのtransformプロパティに相当します。
     * 形状: [simplex_id][row][col] = (nsimplex, ndim+1, ndim)
     * 各行列は T * c = x - r_n の関係を満たします。
     */
    mutable std::vector<std::vector<std::vector<double>>> transform_;
    
    /**
     * @brief 変換行列が計算済みかどうかのフラグ
     */
    mutable bool transform_computed_;
    
    /**
     * @brief 隣接simplex情報
     * 
     * 各simplexに隣接するsimplexのIDを保持します。SciPyのneighborsプロパティに相当します。
     * 形状: [simplex_id][neighbor_index] = neighbor_simplex_id
     */
    mutable std::vector<std::vector<int>> neighbors_;
    
    /**
     * @brief 隣接情報が計算済みかどうかのフラグ
     */
    mutable bool neighbors_computed_;
    
    /**
     * @brief 点群の境界ボックス（最小値）
     */
    std::vector<double> min_bound_;
    
    /**
     * @brief 点群の境界ボックス（最大値）
     */
    std::vector<double> max_bound_;
    
    /**
     * @brief SciPy準拠のparaboloid scaling factor
     */
    double paraboloid_scale_;
    
    /**
     * @brief SciPy準拠のparaboloid offset
     */
    double paraboloid_shift_;
    
    
    /**
     * @brief バリセントリック変換行列を計算
     * 
     * 全simplexに対するバリセントリック変換行列を事前計算します。
     * SciPyの_get_barycentric_transformsに相当する機能です。
     */
    void computeBarycentricTransforms() const;
    
    /**
     * @brief 隣接simplex情報を計算
     * 
     * 全simplexに対する隣接関係を事前計算します。
     * SciPyのneighborsプロパティの計算に相当する機能です。
     */
    void computeNeighbors() const;
    
    /**
     * @brief 単一simplexの変換行列を計算
     * 
     * 指定されたsimplexに対する変換行列を計算します。
     * 
     * @param simplex_vertices simplexの頂点座標
     * @return 変換行列 (ndim+1 x ndim)
     */
    std::vector<std::vector<double>> computeSingleTransform(
        const std::vector<std::vector<double>>& simplex_vertices) const;
    
    /**
     * @brief 行列の逆行列を計算（基本版）
     * 
     * Gauss-Jordan消去法を使用して逆行列を計算します。
     * 
     * @param matrix 正方行列
     * @return 逆行列
     * @throws std::runtime_error 行列が特異な場合
     */
    std::vector<std::vector<double>> invertMatrix(
        const std::vector<std::vector<double>>& matrix) const;
    
    /**
     * @brief 行列の特異性チェック（SciPy準拠）
     * 
     * SciPyの条件数チェックに基づく特異行列判定を行います。
     * 機械イプシロンの1000倍を閾値として使用します。
     * 
     * @param matrix チェック対象の正方行列
     * @return true: 特異または数値的に不安定, false: 安定
     */
    bool isMatrixSingular(const std::vector<std::vector<double>>& matrix) const;
    
    /**
     * @brief 高精度逆行列計算（SciPy準拠）
     * 
     * SciPyのLAPACKベースの高精度逆行列計算を模倣します。
     * 数値安定性を重視した実装です。
     * 
     * @param matrix 正方行列
     * @return 逆行列
     * @throws std::runtime_error 計算失敗時
     */
    std::vector<std::vector<double>> invertMatrixRobust(
        const std::vector<std::vector<double>>& matrix) const;
    
    /**
     * @brief 行列式を計算
     * 
     * LU分解を使用して行列式を計算します。
     * 特異性判定で使用されます。
     * 
     * @param matrix 正方行列
     * @return 行列式の値
     */
    double calculateDeterminant(const std::vector<std::vector<double>>& matrix) const;
};
