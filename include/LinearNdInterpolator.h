#pragma once

#include "Delaunay.h"
#include <vector>
#include <memory>
#include <stdexcept>

/**
 * @brief N次元線形補間クラス
 * 
 * Qhullライブラリに基づく散乱データのN次元線形補間を提供します。
 * Delaunay三角分割を使用して単体を作成し、重心座標を使用して線形補間を実行します。
 * 凸包外の点に対しては最近傍値を返します。
 */
class LinearNdInterpolator {
public:
    /**
     * @brief ベクトル値補間用コンストラクタ
     * 
     * N次元空間の散乱データ点に対してベクトル値の線形補間器を構築します。
     * 各補間点に対して複数次元の値が関連付けられている場合に使用します。
     * 
     * @param points 補間点座標（N次元点のリスト）。各要素は1つの点の座標ベクター
     * @param values 各点に対応する値（M次元ベクトルのリスト）。各要素は1つの点での出力値ベクター
     * @throws std::invalid_argument points が空、または points と values のサイズが一致しない場合
     * @throws std::invalid_argument points の各要素の次元が統一されていない場合
     * @throws std::runtime_error Delaunay三角分割の構築に失敗した場合
     */
    LinearNdInterpolator(
        const std::vector<std::vector<double>>& points, 
        const std::vector<std::vector<double>>& values
    );

    /**
     * @brief スカラー値補間用コンストラクタ
     * 
     * N次元空間の散乱データ点に対してスカラー値の線形補間器を構築します。
     * 各補間点に対して単一のスカラー値が関連付けられている場合に使用します。
     * 
     * @param points 補間点座標（N次元点のリスト）。各要素は1つの点の座標ベクター
     * @param values 各点に対応するスカラー値のリスト
     * @throws std::invalid_argument points が空、または points と values のサイズが一致しない場合
     * @throws std::invalid_argument points の各要素の次元が統一されていない場合
     * @throws std::runtime_error Delaunay三角分割の構築に失敗した場合
     */
    LinearNdInterpolator(
        const std::vector<std::vector<double>>& points, 
        const std::vector<double>& values
    );

    /**
     * @brief デストラクタ
     * 
     * Delaunayオブジェクトとそれに関連するQhullリソースを適切に解放します。
     * unique_ptrによる自動リソース管理により、メモリリークを防止します。
     */
    ~LinearNdInterpolator();

    /**
     * @brief コピーコンストラクタ（削除済み）
     * 
     * Delaunayオブジェクト内のQhullリソースの複雑な管理を避けるため、
     * コピーコンストラクタは禁止されています。
     */
    LinearNdInterpolator(const LinearNdInterpolator&) = delete;
    
    /**
     * @brief コピー代入演算子（削除済み）
     * 
     * Delaunayオブジェクト内のQhullリソースの複雑な管理を避けるため、
     * コピー代入演算子は禁止されています。
     */
    LinearNdInterpolator& operator=(const LinearNdInterpolator&) = delete;
    
    /**
     * @brief ムーブコンストラクタ（削除済み）
     * 
     * Delaunayオブジェクト内のQhullリソースの複雑な管理を避けるため、
     * ムーブコンストラクタは禁止されています。
     */
    LinearNdInterpolator(LinearNdInterpolator&&) = delete;
    
    /**
     * @brief ムーブ代入演算子（削除済み）
     * 
     * Delaunayオブジェクト内のQhullリソースの複雑な管理を避けるため、
     * ムーブ代入演算子は禁止されています。
     */
    LinearNdInterpolator& operator=(LinearNdInterpolator&&) = delete;

    /**
     * @brief 複数点での補間値を一括計算
     * 
     * 指定された複数のクエリ点に対して線形補間を実行し、各点での補間値を返します。
     * 各クエリ点について、それを含むsimplexを特定し、重心座標を用いて線形補間を行います。
     * クエリ点が凸包外にある場合は、NaNの値を返します。
     * 
     * @param query 補間を実行するクエリ点の座標リスト。各要素は1つのN次元点
     * @return 各クエリ点での補間値。ベクトル値補間の場合は各要素がM次元ベクトル
     * @throws std::invalid_argument クエリ点の次元が構築時の点群と一致しない場合
     */
    std::vector<std::vector<double>> interpolate(const std::vector<std::vector<double>>& query) const;
    
    /**
     * @brief 単一点での補間値を計算
     * 
     * 指定された単一のクエリ点に対して線形補間を実行します。
     * SciPy準拠：スカラー値補間の場合は1要素のベクトル、ベクトル値補間の場合はM次元ベクトルを返します。
     * 
     * @param query_point 補間を実行する単一のN次元クエリ点
     * @return 補間値のベクトル。スカラー値補間では[value]、ベクトル値補間では[v0,v1,...,vM-1]
     *         凸包外の場合は適切なサイズのNaNベクトルを返す
     * @throws std::invalid_argument query_pointの次元が構築時の点群と一致しない場合
     */
    std::vector<double> interpolate(const std::vector<double>& query_point) const;

private:
    /**
     * @brief 補間の基準となる点群の座標データ
     * 
     * N次元空間における散乱データ点の座標を格納します。
     * 各要素は1つの点の座標を表すdoubleのベクターで、全ての点は同じ次元数を持ちます。
     * この点群データを基にDelaunay三角分割が構築されます。
     */
    std::vector<std::vector<double>> points_;

    /**
     * @brief 各補間点に対応する値データ
     * 
     * points_の各点に対応する出力値を格納します。
     * スカラー値補間の場合は各要素が1次元ベクター、
     * ベクトル値補間の場合は各要素がM次元ベクターとなります。
     * points_と同じサイズである必要があります。
     */
    std::vector<std::vector<double>> values_;

    /**
     * @brief Delaunay三角分割オブジェクトへのポインタ
     * 
     * 補間計算に必要なDelaunay三角分割機能を提供します。
     * unique_ptrによりRAIIに基づく自動リソース管理が行われ、
     * 内部でQhullライブラリのリソースが管理されます。
     */
    std::unique_ptr<Delaunay> delaunay_;

    /**
     * @brief 補間器の共通初期化処理
     * 
     * @param points 補間点座標（N次元点のリスト）。各要素は同じ次元数を持つ必要があります
     * @param values 各点に対応する値（M次元ベクトルのリスト）。pointsと同じサイズである必要があります
     * @throws std::invalid_argument pointsが空、pointsとvaluesのサイズが不一致、
     *                               または次元の統一性に問題がある場合
     * @throws std::runtime_error Delaunay三角分割の構築に失敗した場合
     */
    void initialize(
        const std::vector<std::vector<double>> &points, 
        const std::vector<std::vector<double>> &values
    );

    /**
     * @brief Delaunay三角分割の構築とセットアップ
     * 
     * @param points 三角分割を構築する点群座標
     * @throws std::runtime_error Qhullライブラリでエラーが発生した場合
     * @throws std::invalid_argument 点群データが三角分割に不適切な場合
     */
    void setupTriangulation(const std::vector<std::vector<double>>& points);
};
