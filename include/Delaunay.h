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
     * @brief 指定された点を含むsimplexのIDを検索
     * 
     * 与えられた点がどのsimplex（N次元三角形）に含まれるかを判定します。
     * SciPyのDelaunay.find_simplexメソッドと同等の機能を提供します。
     * 
     * @param point 検索対象の点の座標（N次元ベクター）
     * @return 点を含むsimplexのID（見つからない場合は-1）
     * @throws std::invalid_argument pointの次元が構築時の点群と一致しない場合
     */
    int findSimplex(const std::vector<double>& point) const;
    
    /**
     * @brief 指定された点の重心座標を計算
     * 
     * 与えられた点について、指定されたsimplex内での重心座標（barycentric coordinates）を計算します。
     * 重心座標は、simplexの各頂点に対する重みを表し、その合計は1になります。
     * 
     * @param point 重心座標を計算する点の座標
     * @param simplex_id 対象となるsimplexのID
     * @return 重心座標のベクター（simplex頂点数+1の要素数）
     * @throws std::invalid_argument simplex_idが無効、またはpointの次元が不正な場合
     * @throws std::runtime_error 重心座標の計算に失敗した場合
     */
    std::vector<double> calculateBarycentricCoordinates(
        const std::vector<double>& point, int simplex_id) const;
    
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
     * @brief Walking algorithmを使用した効率的なsimplex検索
     * 
     * SciPyのfind_simplexで使用される効率的な検索アルゴリズムです。
     * 開始simplexから隣接simplexをたどって目標点を含むsimplexを見つけます。
     * 
     * @param point 検索対象の点の座標
     * @param start_simplex 検索開始点（省略時は0）
     * @return 点を含むsimplexのID（見つからない場合は-1）
     */
    int findSimplexWalking(const std::vector<double>& point, int start_simplex = 0) const;

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
     * @brief 連立一次方程式を解くヘルパーメソッド
     * 
     * 重心座標の計算などで使用される連立一次方程式 Ax = b を解きます。
     * 
     * @param A 係数行列（N×N行列）
     * @param b 定数ベクター（N要素）
     * @return 解ベクター x（N要素）
     * @throws std::runtime_error 行列が特異、または解が存在しない場合
     */
    std::vector<double> solveLinearSystem(
        const std::vector<std::vector<double>>& A, 
        const std::vector<double>& b) const;
    
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
     * @brief 行列の逆行列を計算
     * 
     * Gauss-Jordan消去法を使用して逆行列を計算します。
     * 
     * @param matrix 正方行列
     * @return 逆行列
     * @throws std::runtime_error 行列が特異な場合
     */
    std::vector<std::vector<double>> invertMatrix(
        const std::vector<std::vector<double>>& matrix) const;
};
