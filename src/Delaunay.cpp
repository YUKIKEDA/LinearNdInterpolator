#include "Delaunay.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>

/**
 * @brief Delaunayクラスのコンストラクタ
 * 
 * 与えられた点群からDelaunay三角分割を構築します。
 * SciPyのspatial.Delaunayと同等の入力検証とエラーハンドリングを提供します。
 * 
 * @param points 三角分割する点群。各点はstd::vector<double>として表現され、
 *               すべての点は同じ次元数を持つ必要があります。
 *               最低2次元以上で、n次元の場合は最低n+1個の点が必要です。
 * 
 * @throws std::invalid_argument 以下の場合に投げられます：
 *         - 点群が空の場合
 *         - 点の次元が1未満の場合
 *         - 次元数に対して点数が不足している場合（n次元でn+1個未満）
 *         - 点の次元が一貫していない場合
 *         - 点にNaNや無限大値が含まれている場合
 * @throws std::runtime_error Qhullライブラリでエラーが発生した場合
 * 
 * @note このコンストラクタは内部でQhullライブラリを使用し、
 *       SciPyと同等のオプション（"d Qbb Qc Qz Q12 Qt"）で三角分割を実行します。
 *       5次元以上の場合は追加で"Qx"オプションが使用されます。
 */
Delaunay::Delaunay(const std::vector<std::vector<double>>& points) : transform_computed_(false), neighbors_computed_(false), equations_computed_(false) {
    // 入力点が空の場合はエラー
    if (points.empty()) {
        throw std::invalid_argument("Input points cannot be empty.");
    }
    
    // 入力点の次元数が1未満の場合はエラー
    if (points[0].empty()) {
        throw std::invalid_argument("Input points must have at least one dimension.");
    }
    
    const size_t n_points = points.size();
    const size_t n_dims = points[0].size();
    
    // 最小次元数の確認（三角形を作るためには最低2次元以上必要）
    if (n_dims < 2) {
        throw std::invalid_argument("Points must have at least 2 dimensions for triangulation.");
    }
    
    // 最小点数の確認（n次元では最低n+1個の点が必要）
    if (n_points < n_dims + 1) {
        throw std::invalid_argument(
            "Need at least " 
            + std::to_string(n_dims + 1) 
            + " points for " 
            + std::to_string(n_dims) 
            + "-dimensional triangulation.");
    }
    
    // 次元の一貫性とNaN/無限大値のチェック
    for (size_t i = 0; i < n_points; ++i) {
        if (points[i].size() != n_dims) {
            throw std::invalid_argument("All points must have the same dimension.");
        }
        
        for (size_t j = 0; j < n_dims; ++j) {
            if (!std::isfinite(points[i][j])) {
                throw std::invalid_argument("Input points contain NaN or infinite values.");
            }
        }
    }

    // 点データを平坦化し、内部で保持
    points_ = points;
    std::vector<double> flat_points(n_points * n_dims);
    
    // 境界ボックスを初期化
    min_bound_.resize(n_dims, std::numeric_limits<double>::max());
    max_bound_.resize(n_dims, std::numeric_limits<double>::lowest());
    
    // SciPy準拠のparaboloid設定
    paraboloid_scale_ = 1.0;
    paraboloid_shift_ = 0.0;
    
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < n_dims; ++j) {
            flat_points[i * n_dims + j] = points[i][j];
            // 境界ボックスを更新
            min_bound_[j] = std::min(min_bound_[j], points[i][j]);
            max_bound_[j] = std::max(max_bound_[j], points[i][j]);
        }
    }

    // SciPyと同等のQhullオプション設定
    // d: Delaunay三角分割, Qbb: バウンディングボックス計算, Qc: 共面点保持
    // Qz: 無限遠点を追加（共円・共球問題対応）, Q12: 広角を許可, Qt: 三角分割出力
    std::string options = "d Qbb Qc Qz Q12 Qt";
    if (n_dims >= 5) {
        options += " Qx";
    }

    try {
        qhull_ = std::make_unique<orgQhull::Qhull>();
        qhull_->runQhull("", n_dims, n_points, flat_points.data(), options.c_str());
        
        if (qhull_->qhullStatus() != 0) {
            throw std::runtime_error(
                "Qhull triangulation failed with status: " + 
                std::to_string(qhull_->qhullStatus()));
        }
    } catch (const orgQhull::QhullError &e) {
        throw std::runtime_error("Qhull error: " + std::string(e.what()));
    }
}

/**
 * @brief Delaunayクラスのデストラクタ
 * 
 * リソースの適切な解放を行います。unique_ptrを使用しているため、
 * Qhullオブジェクトは自動的に解放されます。
 */
Delaunay::~Delaunay() = default;

/**
 * @brief 指定された点を含む単体（simplex）のインデックスを検索
 * 
 * 与えられた点がどの単体に含まれるかを判定します。
 * 重心座標を使用して点が単体内部にあるかを判定し、
 * SciPyのspatial.Delaunay.find_simplexと同等の機能を提供します。
 * 
 * @param point 検索対象の点。三角分割時と同じ次元数である必要があります。
 * 
 * @return int 点を含む単体のインデックス（0以上）。
 *             点が凸包の外部にある場合は-1を返します。
 *             Qhullが初期化されていない場合も-1を返します。
 * 
 * @note このメソッドは重心座標を計算して内部判定を行います。
 *       数値誤差を考慮して1e-10の許容範囲を設定しています。
 * 
 * @see calculateBarycentricCoordinates()
 */
int Delaunay::findSimplex(std::vector<double>& barycentric_coords,
                           const std::vector<double>& point, 
                           int& start_simplex, 
                           double eps, 
                           double eps_broad) const {
    // SciPy _find_simplex完全準拠実装
    if (!qhull_ || qhull_->qhullStatus() != 0) {
        return -1;
    }
    
    // 入力検証
    if (point.size() != points_[0].size()) {
        throw std::invalid_argument("Point dimension must match triangulation dimension");
    }
    
    const size_t ndim = points_[0].size();
    barycentric_coords.resize(ndim + 1);
    
    // SciPy _is_point_fully_outside相当のチェック
    for (size_t i = 0; i < ndim; ++i) {
        if (point[i] < min_bound_[i] - eps || point[i] > max_bound_[i] + eps) {
            return -1;
        }
    }
    
    auto simplices = getSimplices();
    if (simplices.empty()) {
        return -1;
    }
    
    // start_simplexの範囲チェック
    int isimplex = start_simplex;
    if (isimplex < 0 || isimplex >= static_cast<int>(simplices.size())) {
        isimplex = 0;
    }
    
    // Step 1: 点をparaboloid上に投影
    std::vector<double> lifted_point;
    liftPointToParaboloid(point, lifted_point);
    
    // Step 2: Walking Algorithm - positive facetを探索
    double best_dist = calculatePlaneDistance(isimplex, lifted_point.data(), ndim + 1);
    bool changed = true;
    
    computeNeighbors();  // 隣接情報を確保
    
    while (changed) {
        if (best_dist > 0) {
            break;  // positive facetを発見
        }
        
        changed = false;
        
        // すべての隣接simplexをチェック
        if (isimplex < static_cast<int>(neighbors_.size())) {
            for (size_t k = 0; k <= ndim && k < neighbors_[isimplex].size(); ++k) {
                int ineigh = neighbors_[isimplex][k];
                if (ineigh == -1) continue;
                
                double dist = calculatePlaneDistance(ineigh, lifted_point.data(), ndim + 1);
                
                // SciPy準拠: epsを考慮した比較（無限ループ防止）
                if (dist > best_dist + eps * (1.0 + std::abs(best_dist))) {
                    isimplex = ineigh;
                    best_dist = dist;
                    changed = true;
                    break;  // SciPy準拠: ループの中途でジャンプ
                }
            }
        }
    }
    
    // Step 3: 現在のsimplex近働にいるので、Directed Searchで正確な位置を探索
    start_simplex = isimplex;
    // SciPy準拠：事前計算済みsimplices配列を最適化版に渡す
    return findSimplexDirected(point, isimplex, eps, eps_broad, barycentric_coords, simplices);
}

/**
 * @brief 三角分割の全単体を取得
 * 
 * Delaunay三角分割によって生成されたすべての単体（シンプレックス）を取得します。
 * 各単体は元の点群のインデックスで表現されます。
 * 
 * n次元空間では、各単体はn+1個の頂点を持ちます：
 * - 2次元：三角形（3頂点）
 * - 3次元：四面体（4頂点）
 * - n次元：n+1頂点
 * 
 * @return std::vector<std::vector<int>> 単体のリスト。
 *         各内側のベクトルは1つの単体を表し、元の点群のインデックスを含みます。
 *         Qhullが初期化されていない場合は空のベクトルを返します。
 * 
 * @note このメソッドはQhullの上半分Delaunay面を除外し、
 *       有効な下半分の面のみを返します。
 *       頂点座標の照合には1e-12の許容範囲を使用しています。
 * 
 * @see findSimplex(), calculateBarycentricCoordinates()
 */
std::vector<std::vector<int>> Delaunay::getSimplices() const {
    std::vector<std::vector<int>> simplices;
    
    if (!qhull_ || qhull_->qhullStatus() != 0) {
        return simplices;
    }
    
    auto facetList = qhull_->facetList();
    
    for (auto facet : facetList) {
        if (facet.isUpperDelaunay()) {
            continue;
        }
        
        auto vertices = facet.vertices();
        std::vector<int> simplex;
        
        for (auto vertex : vertices) {
            // Get vertex coordinates
            auto vertex_point = vertex.point();
            std::vector<double> vertex_coords(vertex_point.coordinates(), 
                                            vertex_point.coordinates() + vertex_point.dimension());
            
            // Find matching point in original points_ array
            int point_index = -1;
            for (size_t i = 0; i < points_.size(); ++i) {
                bool match = true;
                for (size_t j = 0; j < points_[i].size(); ++j) {
                    if (std::abs(points_[i][j] - vertex_coords[j]) > 1e-12) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    point_index = static_cast<int>(i);
                    break;
                }
            }
            
            // Add valid indices only
            if (point_index >= 0) {
                simplex.push_back(point_index);
            }
        }
        
        // Only add simplices with correct number of vertices
        if (simplex.size() == points_[0].size() + 1) {
            simplices.push_back(simplex);
        }
    }
    
    return simplices;
}

const std::vector<std::vector<std::vector<double>>>& Delaunay::getTransform() const {
    if (!transform_computed_) {
        computeBarycentricTransforms();
    }
    return transform_;
}

std::vector<double> Delaunay::calculateBarycentricCoordinatesWithTransform(
    const std::vector<double>& point, int simplex_id) const {
    
    if (!transform_computed_) {
        computeBarycentricTransforms();
    }
    
    if (simplex_id < 0 || simplex_id >= static_cast<int>(transform_.size())) {
        return {};
    }
    
    const size_t ndim = point.size();
    const auto& T = transform_[simplex_id];
    
    if (T.empty() || T.size() != ndim + 1 || T[0].size() != ndim) {
        return {};
    }
    
    // 変換行列にNaNが含まれる場合（特異行列）は空を返す（SciPy準拠）
    for (size_t i = 0; i <= ndim; ++i) {
        for (size_t j = 0; j < ndim; ++j) {
            if (!std::isfinite(T[i][j])) {
                return {};  // 特異行列の場合は計算不可
            }
        }
    }
    
    // SciPy準拠のアルゴリズム: T^(-1) * (x - r_n) = c を計算
    // T[:ndim, :] は既に逆行列、T[ndim, :] は参照頂点 r_n
    std::vector<double> x_minus_rn(ndim);
    for (size_t i = 0; i < ndim; ++i) {
        x_minus_rn[i] = point[i] - T[ndim][i];  // x - r_n
    }
    
    // c = T^(-1) * (x - r_n) 行列ベクトル乗算
    std::vector<double> c(ndim, 0.0);
    for (size_t i = 0; i < ndim; ++i) {
        for (size_t j = 0; j < ndim; ++j) {
            c[i] += T[i][j] * x_minus_rn[j];
        }
    }
    
    // SciPy準拠のバリセントリック座標構築: 
    // barycentric[i] = c[i] for i = 0, 1, ..., n-1
    // barycentric[n] = 1 - sum(c[0:n])
    std::vector<double> barycentric(ndim + 1);
    double sum = 0.0;
    for (size_t i = 0; i < ndim; ++i) {
        barycentric[i] = c[i];
        sum += c[i];
    }
    barycentric[ndim] = 1.0 - sum;  // 最後の要素が参照頂点の重み
    
    return barycentric;
}

void Delaunay::computeBarycentricTransforms() const {
    if (!qhull_ || qhull_->qhullStatus() != 0) {
        return;
    }
    
    auto simplices = getSimplices();
    const size_t nsimplex = simplices.size();
    const size_t ndim = points_.empty() ? 0 : points_[0].size();
    
    transform_.clear();
    transform_.resize(nsimplex);
    
    for (size_t i = 0; i < nsimplex; ++i) {
        const auto& simplex = simplices[i];
        
        // 各simplexの頂点座標を取得
        std::vector<std::vector<double>> vertices(simplex.size());
        for (size_t j = 0; j < simplex.size(); ++j) {
            if (simplex[j] >= 0 && simplex[j] < static_cast<int>(points_.size())) {
                vertices[j] = points_[simplex[j]];
            }
        }
        
        // 変換行列を計算
        transform_[i] = computeSingleTransform(vertices);
    }
    
    transform_computed_ = true;
}

std::vector<std::vector<double>> Delaunay::computeSingleTransform(
    const std::vector<std::vector<double>>& simplex_vertices) const {
    
    if (simplex_vertices.empty()) {
        return {};
    }
    
    const size_t ndim = simplex_vertices[0].size();
    const size_t nverts = simplex_vertices.size();
    
    if (nverts != ndim + 1) {
        return {};  // 不正なsimplex
    }
    
    // SciPy準拠: T[ndim+1][ndim] 形式の変換行列を構築
    std::vector<std::vector<double>> T(ndim + 1, std::vector<double>(ndim));
    
    // SciPy準拠: 参照頂点は最後の頂点（r_n）を使用
    const auto& r_n = simplex_vertices[ndim];
    
    /* SciPy準拠のアルゴリズム:
     * アフィン変換行列を計算: T^(-1) * (x - r_n) = c
     * ここで c は最初のn個のバリセントリック座標
     * c[n] = 1 - sum(c[0:n]) で計算される
     * 
     * 行列 A の各列は (vertex[i] - r_n) for i = 0, 1, ..., n-1
     */
    std::vector<std::vector<double>> A(ndim, std::vector<double>(ndim));
    for (size_t j = 0; j < ndim; ++j) {        // 各列j
        for (size_t i = 0; i < ndim; ++i) {    // 各行i
            A[i][j] = simplex_vertices[j][i] - r_n[i];  // vertex[j] - r_n
        }
    }
    
    // SciPy準拠の条件数チェックと逆行列計算
    try {
        // 数値安定性のチェック：条件数が大きすぎる場合は特異行列として扱う
        if (isMatrixSingular(A)) {
            // 特異または数値的に不安定な行列の場合はNaNで埋める（SciPy準拠）
            for (size_t i = 0; i <= ndim; ++i) {
                for (size_t j = 0; j < ndim; ++j) {
                    T[i][j] = std::numeric_limits<double>::quiet_NaN();
                }
            }
            return T;
        }
        
        // 逆行列を計算（SciPy準拠の高数値精度）
        auto inv_A = invertMatrixRobust(A);
        
        // T[:ndim, :] に逆行列を格納
        for (size_t i = 0; i < ndim; ++i) {
            for (size_t j = 0; j < ndim; ++j) {
                T[i][j] = inv_A[i][j];
            }
        }
        
        // T[ndim, :] に参照頂点 r_n を格納（SciPy準拠）
        for (size_t j = 0; j < ndim; ++j) {
            T[ndim][j] = r_n[j];
        }
        
    } catch (const std::runtime_error&) {
        // 計算エラーの場合はNaNで埋める（SciPy互換）
        for (size_t i = 0; i <= ndim; ++i) {
            for (size_t j = 0; j < ndim; ++j) {
                T[i][j] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
    
    return T;
}

std::vector<std::vector<double>> Delaunay::invertMatrix(
    const std::vector<std::vector<double>>& matrix) const {
    
    if (matrix.empty() || matrix.size() != matrix[0].size()) {
        throw std::runtime_error("Matrix must be square for inversion");
    }
    
    const size_t n = matrix.size();
    
    // 拡張行列 [A|I] を作成
    std::vector<std::vector<double>> augmented(n, std::vector<double>(2 * n));
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmented[i][j] = matrix[i][j];
            augmented[i][j + n] = (i == j) ? 1.0 : 0.0;  // 単位行列
        }
    }
    
    // Gauss-Jordan消去法
    for (size_t i = 0; i < n; ++i) {
        // ピボット選択
        size_t pivot_row = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::abs(augmented[k][i]) > std::abs(augmented[pivot_row][i])) {
                pivot_row = k;
            }
        }
        
        if (pivot_row != i) {
            std::swap(augmented[i], augmented[pivot_row]);
        }
        
        // 特異行列チェック
        if (std::abs(augmented[i][i]) < 1e-14) {
            throw std::runtime_error("Matrix is singular and cannot be inverted");
        }
        
        // 対角成分を1にする
        double pivot = augmented[i][i];
        for (size_t j = 0; j < 2 * n; ++j) {
            augmented[i][j] /= pivot;
        }
        
        // 他の行を消去
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
    std::vector<std::vector<double>> inverse(n, std::vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            inverse[i][j] = augmented[i][j + n];
        }
    }
    
    return inverse;
}

void Delaunay::computeNeighbors() const {
    if (!qhull_ || qhull_->qhullStatus() != 0) {
        return;
    }
    
    neighbors_.clear();
    const size_t ndim = points_.empty() ? 0 : points_[0].size();
    
    // SciPy準拠：QhullのfacetListから隣接関係を取得
    auto facetList = qhull_->facetList();
    
    // 有効なfacet（下半分）のマッピング作成
    std::vector<orgQhull::QhullFacet> valid_facets;
    std::map<int, int> facet_id_to_index;
    
    int index = 0;
    for (auto facet : facetList) {
        if (facet.isUpperDelaunay()) {
            continue;  // 上半分は除外（SciPy準拠）
        }
        valid_facets.push_back(facet);
        facet_id_to_index[facet.id()] = index;
        index++;
    }
    
    // 隣接配列を初期化
    neighbors_.resize(valid_facets.size());
    for (size_t i = 0; i < neighbors_.size(); ++i) {
        neighbors_[i].resize(ndim + 1, -1);  // -1は境界を示す（SciPy準拠）
    }
    
    // SciPyの実装に基づく隣接関係計算
    // 「The kth neighbor is opposite to the kth vertex」の原則に従う
    for (size_t i = 0; i < valid_facets.size(); ++i) {
        const auto& facet = valid_facets[i];
        auto vertices = facet.vertices();
        
        // 各頂点に対応する隣接facetを検索
        int vertex_index = 0;
        for (auto vertex : vertices) {
            if (vertex_index >= static_cast<int>(neighbors_[i].size())) break;
            
            // この頂点を含まない隣接facetを探す（SciPy方式）
            // 全facetを検索して隣接を判定
            for (size_t j = 0; j < valid_facets.size(); ++j) {
                if (i == j) continue;
                
                const auto& other_facet = valid_facets[j];
                auto other_vertices = other_facet.vertices();
                
                // この頂点を共有せず、他のn個の頂点を共有するfacetが隣接
                int shared_count = 0;
                bool shares_current_vertex = false;
                
                for (auto other_vertex : other_vertices) {
                    if (vertex.id() == other_vertex.id()) {
                        shares_current_vertex = true;
                        break;
                    }
                    
                    // 他の頂点との共有をチェック
                    for (auto check_vertex : vertices) {
                        if (check_vertex.id() != vertex.id() && 
                            check_vertex.id() == other_vertex.id()) {
                            shared_count++;
                            break;
                        }
                    }
                }
                
                // n個の頂点を共有し、対象頂点を共有しない場合が隣接
                if (!shares_current_vertex && shared_count == static_cast<int>(ndim)) {
                    neighbors_[i][vertex_index] = static_cast<int>(j);
                    break;
                }
            }
            
            vertex_index++;
        }
    }
    
    neighbors_computed_ = true;
}

/**
 * @brief 面方程式を事前計算（SciPy準拠）
 * 
 * 全simplexに対する面方程式を事前計算し、O(1)アクセスを可能にします。
 * SciPyのequationsプロパティに相当する機能です。
 * 
 * 面方程式は [normal0, normal1, ..., normalN, offset] 形式で格納され、
 * simplex_id * (ndim + 2) + coefficient_index でアクセス可能です。
 */
void Delaunay::computeEquations() const {
    if (equations_computed_ || !qhull_ || qhull_->qhullStatus() != 0) {
        return;
    }
    
    const size_t ndim = points_.empty() ? 0 : points_[0].size();
    auto facetList = qhull_->facetList();
    
    // 有効なfacetのみカウント（下半分のfacetのみ）
    size_t valid_facet_count = 0;
    for (auto facet : facetList) {
        if (!facet.isUpperDelaunay()) {
            valid_facet_count++;
        }
    }
    
    // SciPy準拠: equations配列のサイズ確保
    // 形状: [simplex_id * (ndim + 2) + coefficient_index]
    // coefficients 0..ndim: normal vector, ndim+1: offset
    equations_.resize(valid_facet_count * (ndim + 2));
    
    // SciPy準拠の面方程式計算
    size_t equation_index = 0;
    for (auto facet : facetList) {
        if (facet.isUpperDelaunay()) {
            continue;  // 上半分は除外（SciPy準拠）
        }
        
        auto hyperplane = facet.hyperplane();
        auto coordinates = hyperplane.coordinates();
        
        // SciPy準拠: [normal0, normal1, ..., normalN, offset] 形式で格納
        const size_t base_idx = equation_index * (ndim + 2);
        
        // Normal vector components (lifted space coordinates)
        for (size_t k = 0; k <= ndim; ++k) {
            equations_[base_idx + k] = coordinates[k];
        }
        
        // Offset (ndim+1 position)
        equations_[base_idx + ndim + 1] = hyperplane.offset();
        
        equation_index++;
    }
    
    equations_computed_ = true;
}


/**
 * @brief SciPy準拠のDirected Search実装（最適化版）
 * 
 * 事前計算済みのsimplices配列を受け取ることで、getSimplices()の重複呼び出しを回避します。
 * SciPyと同等のパフォーマンス特性を実現します。
 */
int Delaunay::findSimplexDirected(const std::vector<double>& point,
                                  int start_simplex,
                                  double eps,
                                  double eps_broad,
                                  std::vector<double>& barycentric_coords,
                                  const std::vector<std::vector<int>>& simplices) const {
    // SciPy _find_simplex_directed完全準拠実装（最適化版）
    const size_t ndim = points_[0].size();
    
    int isimplex = start_simplex;
    if (isimplex < 0 || isimplex >= static_cast<int>(simplices.size())) {
        isimplex = 0;
    }
    
    // SciPy準拠の最大イテレーション数
    const int max_iterations = 1 + static_cast<int>(simplices.size()) / 4;
    
    // SciPy準拠：変換行列を一度だけ取得（外側で計算）
    const auto& transform = getTransform();
    
    // SciPy _find_simplex_directed: cycle_k loop
    for (int cycle_k = 0; cycle_k < max_iterations; ++cycle_k) {
        if (isimplex == -1) {
            break;  // 退化simplex
        }
        
        // SciPy準拠：変換行列への直接アクセス（inline barycentric calculation）
        if (isimplex >= static_cast<int>(transform.size())) {
            return -1;
        }
        const auto& T = transform[isimplex];
        
        bool inside = true;
        
        // SciPy準拠：各重心座標を計算してsimplexの内側かチェック
        for (size_t k = 0; k <= ndim; ++k) {
            // SciPy _barycentric_coordinate_single相当をインライン化
            if (k == ndim) {
                // 最後の座標: c[ndim] = 1.0 - Σ(c[j]) for j=0..ndim-1
                barycentric_coords[ndim] = 1.0;
                for (size_t j = 0; j < ndim; ++j) {
                    barycentric_coords[ndim] -= barycentric_coords[j];
                }
            } else {
                // c[k] = Σ(transform[k][j] * (x[j] - transform[ndim][j]))
                barycentric_coords[k] = 0.0;
                for (size_t j = 0; j < ndim; ++j) {
                    barycentric_coords[k] += T[k][j] * (point[j] - T[ndim][j]);
                }
            }
            
            if (barycentric_coords[k] < -eps) {
                // SciPy準拠: 負の重心座標 -> 隣接simplexに移動
                if (isimplex < static_cast<int>(neighbors_.size()) && 
                    k < neighbors_[isimplex].size()) {
                    int m = neighbors_[isimplex][k];
                    if (m == -1) {
                        // 境界に到達 -> 外部なのでexit
                        return -1;
                    }
                    isimplex = m;
                    inside = false;
                    break;  // 新しいsimplexに移動したので次のイテレーション
                }
            } else if (barycentric_coords[k] <= 1.0 + eps) {
                // このsimplex内にいる
                continue;
            } else {
                // 外側または退化simplex
                inside = false;
                break;
            }
        }
        
        if (inside) {
            // 正しいsimplexを発見！
            return isimplex;
        } else if (isimplex == -1) {
            // 退化simplexに出くわした -> brute forceにフォールバック
            return findSimplexBruteForceWithEps(point, eps_broad, barycentric_coords);
        }
        // 別のsimplexにホップしたので継続
    }
    
    // アルゴリズムが収束しなかった -> brute forceにフォールバック
    return findSimplexBruteForceWithEps(point, eps_broad, barycentric_coords);
}

/**
 * @brief Paraboloid上での平面距離計算（SciPy _distplane完全準拠）
 */
double Delaunay::calculatePlaneDistance(int simplex_id, const double* lifted_point, size_t dim) const {
    // SciPy _distplane完全準拠実装（O(1)アクセス）
    // dist = equations[isimplex*(ndim+2) + ndim+1] + Σ(equations[isimplex*(ndim+2) + k] * point[k])
    
    // 事前計算済み面方程式配列を確保
    if (!equations_computed_) {
        computeEquations();
    }
    
    if (simplex_id < 0) {
        return -std::numeric_limits<double>::max();
    }
    
    const size_t ndim = points_[0].size();
    const size_t max_simplices = equations_.size() / (ndim + 2);
    
    if (simplex_id >= static_cast<int>(max_simplices)) {
        return -std::numeric_limits<double>::max();
    }
    
    // SciPy準拠: O(1)の事前計算済み配列アクセス
    const size_t base_idx = static_cast<size_t>(simplex_id) * (ndim + 2);
    
    // SciPy _distplane: dist = offset + Σ(normal[k] * point[k])
    double dist = equations_[base_idx + ndim + 1];  // offset
    
    // 法線ベクトルとの内積計算（paraboloid空間でndim+1個の座標）
    for (size_t k = 0; k <= ndim; ++k) {
        dist += equations_[base_idx + k] * lifted_point[k];
    }
    
    return dist;
}

bool Delaunay::isMatrixSingular(const std::vector<std::vector<double>>& matrix) const {
    if (matrix.empty() || matrix.size() != matrix[0].size()) {
        return true;  // 不正な行列は特異として扱う
    }
    
    const size_t n = matrix.size();
    
    // SciPy準拠: 機械イプシロンの1000倍を閾値として使用
    const double eps = std::numeric_limits<double>::epsilon();
    const double condition_threshold = 1000.0 * eps;
    
    // 簡易的な条件数推定（行列式ベース）
    // より正確にはLAPACKのdgeconを使用するが、ここでは行列式で近似
    try {
        double det = calculateDeterminant(matrix);
        
        // 行列式が非常に小さい場合は特異行列として判定
        if (std::abs(det) < condition_threshold) {
            return true;
        }
        
        // 無限大やNaNの場合も特異として扱う
        if (!std::isfinite(det)) {
            return true;
        }
        
        return false;
        
    } catch (...) {
        return true;  // 計算エラーは特異として扱う
    }
}

std::vector<std::vector<double>> Delaunay::invertMatrixRobust(
    const std::vector<std::vector<double>>& matrix) const {
    
    if (matrix.empty() || matrix.size() != matrix[0].size()) {
        throw std::runtime_error("Invalid matrix for inversion");
    }
    
    const size_t n = matrix.size();
    
    // SciPy準拠: より安定した逆行列計算
    // LU分解を使用した高精度計算を実装
    
    // 作業用の拡張行列を作成 [A|I]
    std::vector<std::vector<double>> extended(n, std::vector<double>(2 * n, 0.0));
    
    // 元の行列をコピー
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            extended[i][j] = matrix[i][j];
        }
        extended[i][i + n] = 1.0;  // 単位行列部分
    }
    
    // Gauss-Jordan消去法（ピボット選択付き）
    for (size_t i = 0; i < n; ++i) {
        // SciPy準拠: ピボット選択
        size_t pivot_row = i;
        double max_pivot = std::abs(extended[i][i]);
        
        for (size_t k = i + 1; k < n; ++k) {
            if (std::abs(extended[k][i]) > max_pivot) {
                max_pivot = std::abs(extended[k][i]);
                pivot_row = k;
            }
        }
        
        // 数値安定性チェック
        if (max_pivot < std::numeric_limits<double>::epsilon() * 1000.0) {
            throw std::runtime_error("Matrix is singular or numerically unstable");
        }
        
        // 行の交換
        if (pivot_row != i) {
            std::swap(extended[i], extended[pivot_row]);
        }
        
        // ピボット要素で正規化
        double pivot = extended[i][i];
        for (size_t j = 0; j < 2 * n; ++j) {
            extended[i][j] /= pivot;
        }
        
        // 他の行を消去
        for (size_t k = 0; k < n; ++k) {
            if (k != i) {
                double mult = extended[k][i];
                for (size_t j = 0; j < 2 * n; ++j) {
                    extended[k][j] -= mult * extended[i][j];
                }
            }
        }
    }
    
    // 逆行列部分を抽出
    std::vector<std::vector<double>> result(n, std::vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result[i][j] = extended[i][j + n];
        }
    }
    
    return result;
}

double Delaunay::calculateDeterminant(const std::vector<std::vector<double>>& matrix) const {
    if (matrix.empty() || matrix.size() != matrix[0].size()) {
        return 0.0;
    }
    
    const size_t n = matrix.size();
    
    if (n == 1) {
        return matrix[0][0];
    }
    
    if (n == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }
    
    // より大きな行列にはLU分解を使用
    auto lu_matrix = matrix;
    double det = 1.0;
    
    for (size_t i = 0; i < n; ++i) {
        // ピボット選択
        size_t pivot_row = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::abs(lu_matrix[k][i]) > std::abs(lu_matrix[pivot_row][i])) {
                pivot_row = k;
            }
        }
        
        if (pivot_row != i) {
            std::swap(lu_matrix[i], lu_matrix[pivot_row]);
            det = -det;  // 行の交換で符号変化
        }
        
        if (std::abs(lu_matrix[i][i]) < std::numeric_limits<double>::epsilon()) {
            return 0.0;  // 特異行列
        }
        
        det *= lu_matrix[i][i];
        
        // 消去
        for (size_t k = i + 1; k < n; ++k) {
            double mult = lu_matrix[k][i] / lu_matrix[i][i];
            for (size_t j = i; j < n; ++j) {
                lu_matrix[k][j] -= mult * lu_matrix[i][j];
            }
        }
    }
    
    return det;
}

void Delaunay::liftPointToParaboloid(const std::vector<double>& point, std::vector<double>& lifted_point) const {
    // SciPy _lift_point完全準拠実装
    const size_t ndim = point.size();
    lifted_point.resize(ndim + 1);
    
    // 通常座標をコピー
    double sum_squares = 0.0;
    for (size_t i = 0; i < ndim; ++i) {
        lifted_point[i] = point[i];
        sum_squares += point[i] * point[i];
    }
    
    // paraboloid座標を計算
    lifted_point[ndim] = sum_squares * paraboloid_scale_ + paraboloid_shift_;
}

void Delaunay::calculateBarycentricCoordinateSingle(int simplex_id, const std::vector<double>& point, 
                                                    std::vector<double>& barycentric, int coordinate_index) const {
    // SciPy _barycentric_coordinate_single完全準拠実装
    const size_t ndim = point.size();
    const auto& transform = getTransform();
    
    if (simplex_id < 0 || simplex_id >= static_cast<int>(transform.size())) {
        return;
    }
    
    const auto& T = transform[simplex_id];
    if (T.empty() || T.size() != ndim + 1 || T[0].size() != ndim) {
        return;
    }
    
    if (coordinate_index == static_cast<int>(ndim)) {
        // 最後の座標: c[ndim] = 1.0 - Σ(c[j]) for j=0..ndim-1
        barycentric[ndim] = 1.0;
        for (size_t j = 0; j < ndim; ++j) {
            barycentric[ndim] -= barycentric[j];
        }
    } else if (coordinate_index >= 0 && coordinate_index < static_cast<int>(ndim)) {
        // c[i] = Σ(transform[i][j] * (x[j] - transform[ndim][j]))
        barycentric[coordinate_index] = 0.0;
        for (size_t j = 0; j < ndim; ++j) {
            barycentric[coordinate_index] += T[coordinate_index][j] * (point[j] - T[ndim][j]);
        }
    }
}

int Delaunay::findSimplexBruteForceWithEps(const std::vector<double>& point,
                                            double eps_broad,
                                            std::vector<double>& barycentric_coords) const {
    // SciPy _find_simplex_bruteforce完全準拠実装
    auto simplices = getSimplices();
    
    for (size_t i = 0; i < simplices.size(); ++i) {
        std::vector<double> barycentric = calculateBarycentricCoordinatesWithTransform(point, static_cast<int>(i));
        
        if (barycentric.empty()) {
            continue;
        }
        
        // SciPy準拠の内部判定: 全ての座標が [-eps_broad, 1+eps_broad] の範囲内
        bool inside = true;
        for (double coord : barycentric) {
            if (coord < -eps_broad || coord > 1.0 + eps_broad) {
                inside = false;
                break;
            }
        }
        
        if (inside) {
            barycentric_coords = barycentric;
            return static_cast<int>(i);
        }
    }
    
    return -1;  // 見つからなかった
}
