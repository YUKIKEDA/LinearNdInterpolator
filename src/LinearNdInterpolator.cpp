#include "LinearNdInterpolator.h"
#include <stdexcept>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>

// Qhull includes
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullPoint.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/QhullError.h"

LinearNdInterpolator::LinearNdInterpolator(
    const std::vector<std::vector<double>>& points, 
    const std::vector<std::vector<double>>& values) {
    
    // 補間点座標と値を初期化
    initialize(points, values);
}

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

LinearNdInterpolator::~LinearNdInterpolator() {
    // Smart pointer handles cleanup automatically
}


// **************************************************
// パブリックメソッド
// **************************************************

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
                "Input points dimension (" + 
                std::to_string(q.size()) + 
                ") does not match expected dimension (" + 
                std::to_string(expected_dimension) + 
                ")"
            );
        }
    }

    // 結果を格納するベクトルを初期化
    std::vector<std::vector<double>> results;
    results.reserve(query.size());

    // **************************************************
    // パフォーマンス最適化のための変数
    // **************************************************
    
    // TODO: 以下の変数を追加
    // - 単体検索の開始位置（前回の結果を利用）
    // - 数値精度のためのepsilon値
    // - 値の次元数
    
    // 値の次元数を取得
    const size_t n_values = values_[0].size();
    
    // **************************************************
    // 各入力点に対して補間を実行
    // **************************************************
    for (const auto& q : query) {
        // **************************************************
        // 1) 単体（シンプレックス）を見つける
        // **************************************************
        int isimplex = delaunay_->findSimplex(q);
        
        // 結果ベクトルを初期化
        std::vector<double> interpolated_values(n_values, 0.0);
        
        // **************************************************
        // 2) 凸包外の処理
        // **************************************************
        if (isimplex == -1) {
            // SciPyと同様に、凸包外の場合はNaNを設定
            for (size_t k = 0; k < n_values; ++k) {
                interpolated_values[k] = std::numeric_limits<double>::quiet_NaN();
            }
            results.push_back(interpolated_values);
            continue;
        }
        
        // **************************************************
        // 3) 重心座標を計算
        // **************************************************
        std::vector<double> barycentric = delaunay_->calculateBarycentricCoordinates(q, isimplex);
        
        if (barycentric.empty()) {
            // 重心座標の計算に失敗した場合もNaNを設定
            for (size_t k = 0; k < n_values; ++k) {
                interpolated_values[k] = std::numeric_limits<double>::quiet_NaN();
            }
            results.push_back(interpolated_values);
            continue;
        }
        
        // **************************************************
        // 4) 線形重心座標補間を実行
        // **************************************************
        auto simplices = delaunay_->getSimplices();
        if (isimplex >= 0 && isimplex < static_cast<int>(simplices.size())) {
            const auto& simplex = simplices[isimplex];
            
            // SciPyの実装に基づく補間計算:
            // for (j = 0; j < ndim+1; j++) {
            //     for (k = 0; k < nvalues; k++) {
            //         m = simplices[isimplex,j];
            //         out[i,k] = out[i,k] + c[j] * values[m,k];
            //     }
            // }
            
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
        
        // **************************************************
        // 5) 結果を格納
        // **************************************************
        results.push_back(interpolated_values);
    }

    return results;
}


// **************************************************
// プライベートメソッド
// **************************************************

void LinearNdInterpolator::initialize(
    const std::vector<std::vector<double>> &points, 
    const std::vector<std::vector<double>> &values)
{
    // **************************************************
    // 入力値のバリデーション
    // **************************************************
    throwifInvalidInputDataShape(points, values);

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

    // 補間点座標と値をメンバ変数に格納
    points_ = points;
    values_ = values;

    // ドロネー三角形分割を作成
    setupTriangulation(points);
}

void LinearNdInterpolator::throwifInvalidInputDataShape(
    const std::vector<std::vector<double>> &points,
    const std::vector<std::vector<double>> &values) const
{
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
}

void LinearNdInterpolator::setupTriangulation(
    const std::vector<std::vector<double>> &points)
{
    delaunay_ = std::make_unique<Delaunay>(points);
}
