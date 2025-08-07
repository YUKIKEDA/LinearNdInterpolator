/**
 * @file Qhull.h
 * @brief Qhull C APIのC++ラッパークラス
 * 
 * このファイルはQhull計算幾何学ライブラリのC APIをC++でラップし、
 * SciPyの_Qhullクラスと互換性のあるインターフェースを提供します。
 * 主にDelaunay三角分割とVoronoi図の計算に使用されます。
 */

#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <tuple>

// Qhull C APIのインクルード - third_partyからの相対パス
extern "C" {
#include "../third_party/qhull-2020.2/src/libqhull_r/qhull_ra.h"
}

namespace qhull {

/**
 * @brief Qhullの計算結果を表現するための型エイリアス
 * 
 * SimplexFacetArrayの戻り値として使用され、以下の要素を含みます：
 * - simplices: 各単体（三角形、四面体など）の頂点インデックス
 * - neighbors: 各単体の隣接単体のインデックス
 * - equations: 各単体の平面方程式
 * - coplanar: 同一平面上の点のインデックス
 * - good: 各単体が有効かどうかのフラグ
 */
using SimplexFacetResult = std::tuple<
    std::vector<std::vector<int>>,    // simplices
    std::vector<std::vector<int>>,    // neighbors
    std::vector<std::vector<double>>, // equations
    std::vector<std::vector<int>>,    // coplanar
    std::vector<bool>                 // good
>;

/**
 * @class Qhull
 * @brief Qhull計算幾何学ライブラリのC++ラッパークラス
 * 
 * このクラスはQhull C APIをラップし、SciPyの_Qhullクラスと
 * 互換性のあるインターフェースを提供します。主にDelaunay三角分割、
 * Voronoi図、凸包の計算に使用されます。
 * 
 * @note このクラスはQhull 2020.2に基づいて実装されています。
 * @note SciPyとの互換性を保つため、メソッド名はPythonの慣例に従っています。
 */
class Qhull {
private:
    /**
     * Qhullコンテキスト（リエントラント版）
     */
    qhT qh_qh;

    /**
     * 入力点群データ
     */
    std::vector<std::vector<double>> points_;

    /**
     * Qhullコマンド（例: "d", "v"）
     */
    std::string command_;

    /**
     * Qhullオプション文字列
     */
    std::string options_;

    /**
     * 計算が実行済みかどうか
     */
    bool computed_;

    /**
     * Delaunay三角分割モードかどうか
     */
    bool is_delaunay_;

    /**
     * 次元数
     */
    size_t ndim_;
    
public:
    /**
     * @brief コンストラクタ
     * 
     * SciPyの_Qhull.__init__を参考に実装されています。
     * 指定されたコマンドとオプションでQhullを初期化します。
     * 
     * @param command Qhullコマンド（"d": Delaunay, "v": Voronoi, "": convex hull）
     * @param points 入力点群（各点はN次元座標のベクトル）
     * @param options Qhullオプション文字列（例: "Qbb Qc Qz"）
     * 
     * @throws std::invalid_argument 無効なパラメータが渡された場合
     * @throws std::runtime_error Qhullの初期化に失敗した場合
     */
    Qhull(
        const std::string& command, 
        const std::vector<std::vector<double>>& points, 
        const std::string& options
    );
    
    /**
     * @brief デストラクタ
     * 
     * Qhullリソースを適切に解放します。
     */
    ~Qhull();

    /**
     * @brief 三角分割を実行
     * 
     * SciPyの_Qhull.triangulate()に対応するメソッドです。
     * Qhullを実行して三角分割または指定された計算を行います。
     * 
     * @throws std::runtime_error 計算に失敗した場合
     * @throws std::logic_error 既に計算が実行済みの場合
     */
    void triangulate();

    /**
     * @brief 放物面の平行移動とスケールを取得
     * 
     * SciPyの_Qhull.get_paraboloid_shift_scale()に対応するメソッドです。
     * Delaunay三角分割で使用される放物面変換のパラメータを返します。
     * 
     * @return std::pair<double, double> (shift, scale)のペア
     * @throws std::logic_error triangulate()が呼ばれていない場合
     * @throws std::runtime_error Delaunayモードでない場合
     */
    std::pair<double, double> getParaboloidShiftScale() const;

    /**
     * @brief 単体とファセットの配列を取得
     * 
     * SciPyの_Qhull.get_simplex_facet_array()に対応するメソッドです。
     * 計算結果として単体（三角形、四面体など）、隣接関係、
     * 平面方程式などの情報を返します。
     * 
     * @return SimplexFacetResult 計算結果の構造体
     *         - simplices: 各単体の頂点インデックス
     *         - neighbors: 各単体の隣接単体インデックス  
     *         - equations: 各単体の平面方程式係数
     *         - coplanar: 同一平面上の点のインデックス
     *         - good: 各単体が有効（非退化）かどうか
     * 
     * @throws std::logic_error triangulate()が呼ばれていない場合
     * @throws std::runtime_error データの取得に失敗した場合
     */
    SimplexFacetResult getSimplexFacetArray() const;
};

}