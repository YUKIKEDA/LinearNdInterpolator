#include "LinearNdInterpolator.h"
#include <stdexcept>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <string>

/**
 * @brief N次元線形補間器のコンストラクタ（多次元値用）
 * 
 * 与えられた補間点座標と対応する多次元値を使用してN次元線形補間器を初期化します。
 * 内部でドロネー三角分割を構築し、線形重心座標補間による高速な補間を可能にします。
 * 
 * @param points 補間点座標の配列。各要素は点の座標を表すベクトル（最低2次元）
 * @param values 各補間点に対応する値の配列。各要素は多次元値を表すベクトル
 * 
 * @throws std::invalid_argument 以下の場合に例外が発生します：
 *   - points と values のサイズが一致しない
 *   - 点の次元数が2未満
 *   - 全ての点の次元が一致しない
 *   - 全ての値の次元が一致しない
 *   - points または values が空
 */
LinearNdInterpolator::LinearNdInterpolator(
    const std::vector<std::vector<double>>& points, 
    const std::vector<std::vector<double>>& values) {
    
    // 補間点座標と値を初期化
    initialize(points, values);
}

/**
 * @brief N次元線形補間器のコンストラクタ（1次元値用）
 * 
 * 与えられた補間点座標と対応する1次元スカラー値を使用してN次元線形補間器を初期化します。
 * 1次元値は内部で2次元形式に変換されてから処理されます。
 * 
 * @param points 補間点座標の配列。各要素は点の座標を表すベクトル（最低2次元）
 * @param values 各補間点に対応するスカラー値の配列
 * 
 * @throws std::invalid_argument 以下の場合に例外が発生します：
 *   - points と values のサイズが一致しない
 *   - 点の次元数が2未満
 *   - 全ての点の次元が一致しない
 *   - points または values が空
 */
LinearNdInterpolator::LinearNdInterpolator(
    const std::vector<std::vector<double>>& points, 
    const std::vector<double>& values) {
    
    // 1次元の場合は、2次元の場合に変換（各値を要素数1のベクトルに変換）
    std::vector<std::vector<double>> values2d(values.size(), std::vector<double>(1));
    for (size_t i = 0; i < values.size(); ++i) {
        values2d[i][0] = values[i];
    }

    // 補間点座標と値を初期化
    initialize(points, values2d);
}

/**
 * @brief N次元線形補間器のデストラクタ
 * 
 * オブジェクトの破棄時にリソースのクリーンアップを行います。
 * 現在はQhullのリソース解放処理は自動的に行われます。
 */
LinearNdInterpolator::~LinearNdInterpolator() {
    // Qhullのリソースを解放
}

/**
 * @brief 補間器の初期化処理
 * 
 * 与えられた補間点座標と値を使用して補間器を初期化します。
 * 入力データの妥当性を検証し、メンバ変数に格納後、ドロネー三角分割を構築します。
 * 
 * @param points 補間点座標の配列。各要素は点の座標を表すベクトル
 * @param values 各補間点に対応する値の配列。各要素は値を表すベクトル
 * 
 * @throws std::invalid_argument 以下の場合に例外が発生します：
 *   - points と values のサイズが不一致
 *   - 点の次元数が2未満
 *   - 点の次元に一貫性がない
 *   - 値の次元に一貫性がない
 *   - points または values が空
 */
void LinearNdInterpolator::initialize(
    const std::vector<std::vector<double>> &points, 
    const std::vector<std::vector<double>> &values)
{
    // **************************************************
    // 入力値のバリデーション
    // **************************************************

    // 補間点座標と値の数が一致しない場合はエラー
    if (points.size() != values.size()) {
        throw std::invalid_argument(
            "Points size (" + 
            std::to_string(points.size()) + 
            ") does not match values size (" + 
            std::to_string(values.size()) + 
            ")"
        );
    }

    // 補間点座標と値が空の場合はエラー
    if (points.empty() || values.empty()) {
        throw std::invalid_argument("Points and values cannot be empty");
    }

    // 補間点座標の次元数が2未満の場合はエラー
    const size_t point_dimension = points[0].size();
    if (point_dimension < 2) {
        throw std::invalid_argument("Points must have at least 2 dimensions");
    }

    // 補間点座標の次元の一貫性チェック(全ての補間点の次元が一致しているか)
    for (const auto& point : points) {
        if (point.size() != point_dimension) {
            throw std::invalid_argument("All points must have the same dimension");
        }
    }

    // 各点に対応する値の次元の一貫性チェック(全ての値の次元が一致しているか)
    const size_t value_dimension = values[0].size();
    for (const auto& value : values) {
        if (value.size() != value_dimension) {
            throw std::invalid_argument("All values must have the same dimension");
        }
    }

    // **************************************************

    // 補間点座標と値をメンバ変数に格納
    points_ = points;
    values_ = values;

    // ドロネー三角形分割を作成
    setupTriangulation(points);
}

/**
 * @brief ドロネー三角分割の設定
 * 
 * 与えられた補間点座標を使用してドロネー三角分割を構築します。
 * 構築されたDelaunayオブジェクトは線形補間の基盤として使用されます。
 * 
 * @param points 補間点座標の配列。各要素は点の座標を表すベクトル
 */
void LinearNdInterpolator::setupTriangulation(
    const std::vector<std::vector<double>> &points)
{
    delaunay_ = std::make_unique<Delaunay>(points);
}

/**
 * @brief 複数点に対する線形補間
 * 
 * 与えられた複数のクエリ点に対して線形補間を実行し、補間値を返します。
 * ドロネー三角分割と重心座標を使用した線形補間により、SciPyと互換性のある結果を提供します。
 * 
 * @param query 補間対象の点座標の配列。各要素は点の座標を表すベクトル
 * 
 * @return 各クエリ点に対応する補間値の配列。各要素は補間値を表すベクトル
 *         凸包外の点についてはNaN値が設定されます
 * 
 * @throws std::invalid_argument クエリ点の次元が補間点の次元と一致しない場合
 * 
 * @note 
 *   - 空のクエリが渡された場合は空の結果配列を返します
 *   - 凸包外の点や計算に失敗した点についてはNaN値が設定されます
 *   - SciPyのLinearNDInterpolatorと同等の動作を提供します
 */
std::vector<std::vector<double>> LinearNdInterpolator::interpolate(
    const std::vector<std::vector<double>> &query) const
{    
    // 入力点の検証
    if (query.empty()) {
        return std::vector<std::vector<double>>();
    }

    // 入力点の次元数チェック (すべてのクエリ点の次元が補間点の次元と一致しているか)
    const size_t expected_dimension = points_[0].size();
    for (const auto& q : query) {
        if (q.size() != expected_dimension) {
            throw std::invalid_argument(
                "Input points dimension (" 
                + std::to_string(q.size()) 
                + ") does not match expected dimension (" 
                + std::to_string(expected_dimension) 
                + ")"
            );
        }
    }

    // 結果を格納するベクトルを初期化
    std::vector<std::vector<double>> results;
    results.reserve(query.size());

    // 値の次元数を取得
    const size_t n_values = values_[0].size();
    
    // 各入力点に対して補間を実行
    for (const auto& q : query) {
        // 単体（シンプレックス）を見つける
        int isimplex = delaunay_->findSimplex(q);
        
        // 結果ベクトルを初期化
        std::vector<double> interpolated_values(n_values, 0.0);
        
        // 凸包外の処理
        if (isimplex == -1) {
            // SciPyと同様に、凸包外の場合はNaNを設定
            for (size_t k = 0; k < n_values; ++k) {
                interpolated_values[k] = std::numeric_limits<double>::quiet_NaN();
            }
            results.push_back(interpolated_values);
            continue;
        }
        
        // SciPy互換の最適化された変換行列を使用
        std::vector<double> barycentric = delaunay_->calculateBarycentricCoordinatesWithTransform(q, isimplex);
        
        // フォールバック: 変換行列が失敗した場合は元の方法を使用
        if (barycentric.empty()) {
            barycentric = delaunay_->calculateBarycentricCoordinates(q, isimplex);
        }
        
        if (barycentric.empty()) {
            // 重心座標の計算に失敗した場合もNaNを設定
            for (size_t k = 0; k < n_values; ++k) {
                interpolated_values[k] = std::numeric_limits<double>::quiet_NaN();
            }
            results.push_back(interpolated_values);
            continue;
        }
        
        // 線形重心座標補間を実行
        auto simplices = delaunay_->getSimplices();
        if (isimplex >= 0 && isimplex < static_cast<int>(simplices.size())) {
            const auto& simplex = simplices[isimplex];
            
            for (size_t j = 0; j < barycentric.size() && j < simplex.size(); ++j) {
                int vertex_index = simplex[j];
                
                // 頂点インデックスの範囲チェック
                if (vertex_index >= 0 && vertex_index < static_cast<int>(values_.size())) {
                    for (size_t k = 0; k < n_values; ++k) {
                        interpolated_values[k] += barycentric[j] * values_[vertex_index][k];
                    }
                }
            }
        } else {
            // 単体インデックスが無効な場合はNaNを設定
            for (size_t k = 0; k < n_values; ++k) {
                interpolated_values[k] = std::numeric_limits<double>::quiet_NaN();
            }
        }
        
        // 結果を格納
        results.push_back(interpolated_values);
    }

    return results;
}

/**
 * @brief 単一点に対する線形補間（スカラー値用）
 * 
 * 与えられた単一のクエリ点に対して線形補間を実行し、スカラー補間値を返します。
 * 内部で複数点版のinterpolateメソッドを呼び出し、結果の最初の値を返します。
 * 
 * @param query_point 補間対象の点座標
 * 
 * @return 補間されたスカラー値。計算に失敗した場合やクエリ点が凸包外の場合はNaNを返します
 */
double LinearNdInterpolator::interpolate(const std::vector<double>& query_point) const {
    // 単一点を複数点形式に変換
    std::vector<std::vector<double>> query = {query_point};
    
    // 既存の複数点interpolateメソッドを呼び出し
    auto results = interpolate(query);
    
    // 結果から最初の値を取得（1次元値の場合）
    if (!results.empty() && !results[0].empty()) {
        return results[0][0];
    }
    
    // エラーの場合はNaNを返す
    return std::numeric_limits<double>::quiet_NaN();
}