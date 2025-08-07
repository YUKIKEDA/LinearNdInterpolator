#include "LinearNdInterpolator.h"
#include <stdexcept>
#include <iostream>
#include <limits>
#include <algorithm>
#include <cmath>

// ---------------------------------------------------------------------------
// LinearNdInterpolator クラスの実装
// ---------------------------------------------------------------------------

LinearNdInterpolator::LinearNdInterpolator(
    const std::vector<std::vector<double>>& input_points,
    const std::vector<double>& input_values)
    : LinearNdInterpolator(input_points, _convert_to_2d(input_values)) 
{
}

LinearNdInterpolator::LinearNdInterpolator(
    const std::vector<std::vector<double>>& input_points,
    const std::vector<std::vector<double>>& input_values) 
{
    points_ = input_points;
    calculateTriangulation();
    setValues(input_values);
}

std::vector<std::vector<double>> LinearNdInterpolator::interpolate(const std::vector<std::vector<double>>& xi) 
{
    // -----------------------------------
    // --- 1. _preprocess_xi 相当の処理 ---
    // -----------------------------------
    
    // C++では引数の型で整形済みなので、
    // Scipyの xi = _ndim_coords_from_arrays(args, ndim=self.points.shape[1])
    // の処理は不要。

    checkInterpolateShape(xi);

    // Scipyの xi = xi.reshape(-1, xi.shape[-1]) 相当の処理
    // C++ではxiは既にこの形状なので、処理は不要。

    // Scipy: return self._scale_x(xi)
    // rescale=Falseなので、スケーリング処理はスキップ。

    // -------------------------------------
    // --- 2. _evaluate_double 相当の処理 ---
    // -------------------------------------
    // Scipy: r = self._evaluate_double(xi)

    return evaluate(xi);
}

std::vector<double> LinearNdInterpolator::interpolate(const std::vector<double>& xi) 
{
    // 単一点を点群形式に変換
    std::vector<std::vector<double>> xi_batch = {xi};
    
    // 複数点版のinterpolateを呼び出し
    auto result_batch = interpolate(xi_batch);
    
    // 結果の最初の要素（単一点の結果）を返す
    return result_batch[0];
}

void LinearNdInterpolator::calculateTriangulation() 
{
    tri_ = std::make_unique<qhull::Delaunay>(points_);
}

void LinearNdInterpolator::setValues(const std::vector<std::vector<double>>& input_values) 
{
    checkInitShape(points_, input_values);
    values_ = input_values;
    fill_value_ = std::numeric_limits<double>::quiet_NaN();
}

void LinearNdInterpolator::checkInitShape(
    const std::vector<std::vector<double>>& p, 
    const std::vector<std::vector<double>>& v) 
{
    if (p.empty()) {
        throw std::invalid_argument("points array cannot be empty");
    }
    
    const size_t ndim = p[0].size();
    if (ndim < 2) {
        throw std::invalid_argument("input data must be at least 2-D");
    }
    
    for(const auto& point_vec : p) {
        if(point_vec.size() != ndim) {
            throw std::invalid_argument("points must have consistent dimensions");
        }
    }

    if (v.size() != p.size()) {
        throw std::invalid_argument("different number of values and points");
    }
}

void LinearNdInterpolator::checkInterpolateShape(const std::vector<std::vector<double>>& xi) const 
{
    if (xi.empty()) {
        return; // 空の入力は許容
    }

    // Scipyの self.points.shape[1] は this->get_points()[0].size() に相当
    const size_t ndim = points_[0].size();
    if (xi[0].size() != ndim) {
        throw std::invalid_argument("number of dimensions in xi does not match x");
    }
}

std::vector<std::vector<double>> LinearNdInterpolator::evaluate(const std::vector<std::vector<double>>& xi) const {
    
    // --------------------------------------------------
    // 変数の初期化
    // --------------------------------------------------
    const auto& local_values = values_;
    const double local_fill_value = fill_value_;
    
    const size_t n_interp_points = xi.size();
    
    if (n_interp_points == 0) {
        return {};
    }
    
    if (!tri_) {
        return {};
    }
    
    const auto& simplices = tri_->get_simplices();
    
    if (simplices.empty()) {
        const size_t ndim = xi[0].size();
        std::vector<std::vector<double>> output(n_interp_points, std::vector<double>(local_values[0].size(), local_fill_value));
        return output;
    }
    
    const size_t ndim = xi[0].size();  // Scipy: ndim = xi.shape[1]

    // Scipy: out = np.empty((xi.shape[0], self.values.shape[1]), dtype=self.values.dtype)
    std::vector<std::vector<double>> output(n_interp_points, std::vector<double>(local_values[0].size()));
    
    // Scipy: nvalues = out.shape[1]
    const size_t n_values = output[0].size();
    
    // Scipy: cdef int start = 0
    int start_simplex_hint = 0;
    
    // Scipy: eps = 100 * DBL_EPSILON
    const double eps = 100.0 * std::numeric_limits<double>::epsilon();
    
    // Scipy: eps_broad = sqrt(DBL_EPSILON)
    const double eps_broad = std::sqrt(std::numeric_limits<double>::epsilon());

    // 重心座標を格納するためのベクトル
    std::vector<double> barycentric_coords;
    barycentric_coords.resize(ndim + 1);  // 事前にサイズを確保

    // --- Scipyのメインループを忠実に再現 ---
    for (size_t i = 0; i < n_interp_points; ++i) {
        const auto& point = xi[i];

        // 1) Find the simplex
        try {
            int isimplex = tri_->findSimplex(
                barycentric_coords,
                point,
                start_simplex_hint,
                eps,
                eps_broad
            );

            // 2) Linear barycentric interpolation
            if (isimplex == -1) {
                for (size_t k = 0; k < n_values; ++k) {
                    output[i][k] = local_fill_value;
                }
                continue;
            }
            
            // HACK: 不要かも
            // 境界チェック
            if (isimplex >= static_cast<int>(simplices.size())) {
                for (size_t k = 0; k < n_values; ++k) {
                    output[i][k] = local_fill_value;
                }
                continue;
            }
            
            std::fill(output[i].begin(), output[i].end(), 0.0);
            
            const auto& simplex = simplices[isimplex];
            if (simplex.size() != ndim + 1) {
                for (size_t k = 0; k < n_values; ++k) {
                    output[i][k] = local_fill_value;
                }
                continue;
            }
            
            for (size_t j = 0; j < ndim + 1; ++j) {
                int m = simplex[j];
                
                // 境界チェック
                if (m >= static_cast<int>(local_values.size()) || m < 0) {
                    for (size_t k = 0; k < n_values; ++k) {
                        output[i][k] = local_fill_value;
                    }
                    break;
                }
                
                if (j >= barycentric_coords.size()) {
                    break;
                }
                
                for (size_t k = 0; k < n_values; ++k) {
                    output[i][k] += barycentric_coords[j] * local_values[m][k];
                }
            }
            
        } catch (const std::exception& e) {
            for (size_t k = 0; k < n_values; ++k) {
                output[i][k] = local_fill_value;
            }
        }
    }
    
    return output;
}

std::vector<std::vector<double>> LinearNdInterpolator::_convert_to_2d(const std::vector<double>& v) {
    std::vector<std::vector<double>> v_2d;
    v_2d.reserve(v.size());
    for (const auto& val : v) {
        v_2d.push_back({val});
    }
    return v_2d;
}

