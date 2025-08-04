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
     * @brief コンストラクタ
     * @param points 補間点座標（N次元点のリスト）
     * @param values 各点に対応する値（N次元点のリスト）
     */
    LinearNdInterpolator(
        const std::vector<std::vector<double>>& points, 
        const std::vector<std::vector<double>>& values
    );

    /**
     * @brief コンストラクタ
     * @param points 補間点座標（N次元点のリスト）
     * @param values 各点に対応する値
     */
    LinearNdInterpolator(
        const std::vector<std::vector<double>>& points, 
        const std::vector<double>& values
    );

    /**
     * @brief デストラクタ（Qhullのリソースを解放）
     */
    ~LinearNdInterpolator();

    // コピーとムーブは禁止 (Qhullポインタの管理が複雑になるため)
    LinearNdInterpolator(const LinearNdInterpolator&) = delete;
    LinearNdInterpolator& operator=(const LinearNdInterpolator&) = delete;
    LinearNdInterpolator(LinearNdInterpolator&&) = delete;
    LinearNdInterpolator& operator=(LinearNdInterpolator&&) = delete;

    /**
     * @brief 指定された点での補間値を計算
     * @param query 補間が実行される点座標
     * @return 補間値（凸包外の場合は最近傍値）
     */
    std::vector<std::vector<double>> interpolate(const std::vector<std::vector<double>>& query) const;
    
    /**
     * @brief 単一点での補間値を計算（テスト用オーバーロード）
     * @param query_point 補間が実行される単一点座標
     * @return 補間値（スカラー値、凸包外の場合はNaN）
     */
    double interpolate(const std::vector<double>& query_point) const;

private:
    // **************************************************
    // プライベートフィールド
    // **************************************************

    // 補間点座標（N次元点のリスト）
    std::vector<std::vector<double>> points_;

    // 各点に対応する値（N次元点のリスト）
    std::vector<std::vector<double>> values_;

    // Delaunayオブジェクト
    std::unique_ptr<Delaunay> delaunay_;

    // **************************************************
    // プライベートメソッド
    // **************************************************

    /**
     * @brief 補間器の初期化を行う
     * @param points 補間点座標（N次元点のリスト）
     * @param values 各点に対応する値（N次元点のリスト）
     * 
     * 入力データの検証を行い、メンバ変数に格納します。
     * また、ドロネー三角形分割の設定も行います。
     */
    void initialize(
        const std::vector<std::vector<double>> &points, 
        const std::vector<std::vector<double>> &values
    );

    /**
     * @brief 入力データの形状が有効かどうかを検証する
     * @param points 補間点座標（N次元点のリスト）
     * @param values 各点に対応する値（N次元点のリスト）
     * 
     * 補間点座標と値の数が一致するかどうかをチェックします。
     * 一致しない場合はstd::invalid_argument例外を投げます。
     */
    void throwifInvalidInputDataShape(
        const std::vector<std::vector<double>> &points, 
        const std::vector<std::vector<double>> &values
    ) const;

    /**
     * @brief ドロネー三角形分割をセットアップ
     * @param points 補間点座標
     */
    void setupTriangulation(const std::vector<std::vector<double>>& points);
};
