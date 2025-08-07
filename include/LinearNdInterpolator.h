#include <vector>
#include <iostream>
#include "Delaunay.h"

/**
 * @brief N次元線形補間器クラス
 * 
 * このクラスは、N次元空間における点群データに対して線形補間を実行します。
 * SciPyのLinearNDInterpolatorの機能をC++で実装したものです。
 * 
 * Delaunay三角分割を使用して補間を行い、任意の次元数の点群データを
 * サポートします。補間は重心座標を使用した線形補間で実行されます。
 * 
 * @note このクラスは最低2次元以上のデータを必要とします。
 */
class LinearNDInterpolator {
public:
    
    /**
     * @brief スカラー値を補間する際に使用するコンストラクタ
     * 
     * N次元の点群データと1次元の値配列から線形補間器を構築します。
     * 1次元の値配列は内部的に2次元配列に変換されます。
     * 
     * @param input_points N次元の点群データ。各要素は同じ次元数である必要があります。
     * @param input_values 各点に対応する値の配列。点の数と同じサイズである必要があります。
     * 
     * @throws std::invalid_argument 点群が空、2次元未満、または点と値の数が一致しない場合
     * 
     * @example
     * ```cpp
     * std::vector<std::vector<double>> points = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
     * std::vector<double> values = {1.0, 2.0, 3.0};
     * LinearNDInterpolator interp(points, values);
     * ```
     * 
     * @sa https://github.com/scipy/scipy/blob/a0520a1b6b4b3349b8e1fc4d879aae831ecfb193/scipy/interpolate/_interpnd.pyx#L305
     * @sa https://github.com/scipy/scipy/blob/a0520a1b6b4b3349b8e1fc4d879aae831ecfb193/scipy/interpolate/_interpnd.pyx#L62
     * 
     */
    LinearNDInterpolator(
        const std::vector<std::vector<double>>& input_points,
        const std::vector<double>& input_values)
        : LinearNDInterpolator(input_points, _convert_to_2d(input_values)) 
    {
    }

    
    /**
     * @brief 多次元値配列を使用するコンストラクタ
     * 
     * N次元の点群データと多次元の値配列から線形補間器を構築します。
     * 各点に対して複数の値を持つことができます（例：ベクトル場の補間）。
     * 
     * @param input_points N次元の点群データ。各要素は同じ次元数である必要があります。
     * @param input_values 各点に対応する値の配列。外側の配列サイズは点の数と一致し、
     *                    内側の配列は各点での値のベクトルを表します。
     * 
     * @throws std::invalid_argument 点群が空、2次元未満、または点と値の数が一致しない場合
     * 
     * @example
     * ```cpp
     * std::vector<std::vector<double>> points = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
     * std::vector<std::vector<double>> values = {{1.0, 4.0}, {2.0, 5.0}, {3.0, 6.0}};
     * LinearNDInterpolator interp(points, values);
     * ```
     * 
     * @sa https://github.com/scipy/scipy/blob/a0520a1b6b4b3349b8e1fc4d879aae831ecfb193/scipy/interpolate/_interpnd.pyx#L305
     * @sa https://github.com/scipy/scipy/blob/a0520a1b6b4b3349b8e1fc4d879aae831ecfb193/scipy/interpolate/_interpnd.pyx#L62
     * 
     */
    LinearNDInterpolator(
        const std::vector<std::vector<double>>& input_points,
        const std::vector<std::vector<double>>& input_values) 
    {        
        points_ = input_points;
        calculateTriangulation();
        setValues(input_values);        
    }


    /**
     * @brief 指定された点群に対して線形補間を実行するメインメソッド
     * 
     * N次元空間の任意の点において、Delaunay三角分割と重心座標を使用して
     * 線形補間を実行します。SciPyのLinearNDInterpolator.__call__に対応します。
     * 
     * 補間プロセス：
     * 1. 入力点群の形状検証
     * 2. 各点に対してDelaunay三角分割内の単体（simplex）を検索
     * 3. 重心座標を計算して線形補間を実行
     * 4. 補間領域外の点にはfill_value_を返す
     * 
     * @param xi 補間を実行する点群。各要素は元の点群と同じ次元数である必要があります。
     *           形状: [n_points, n_dimensions]
     * 
     * @return 補間結果の値配列。形状: [n_points, n_values]
     *         各要素は対応する入力点での補間値ベクトルです。
     * 
     * @throws std::invalid_argument 以下の場合に例外を投げます：
     *   - 入力点の次元数が元の点群の次元数と一致しない場合
     * 
     * @note 補間領域外（Delaunay三角分割の凸包外）の点には、
     *       fill_value_（デフォルトはNaN）が返されます。
     * @note この関数は内部的にDelaunay三角分割を遅延初期化します。
     * 
     * @example
     * ```cpp
     * std::vector<std::vector<double>> query_points = {{0.5, 0.5}, {1.5, 1.5}};
     * auto result = interp.interpolate(query_points);
     * 
     * result[0]は点(0.5, 0.5)での補間値
     * result[1]は点(1.5, 1.5)での補間値（領域外の場合はNaN）
     * ```
     * 
     * @sa https://github.com/scipy/scipy/blob/a0520a1b6b4b3349b8e1fc4d879aae831ecfb193/scipy/interpolate/_interpnd.pyx#L154
     * 
     */
    std::vector<std::vector<double>> interpolate(const std::vector<std::vector<double>>& xi) {
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

private:
    /**
     * @brief 補間に使用するN次元の点群データ
     * 
     * 各要素は同じ次元数のベクトルで、補間の基準となる点の座標を表します。
     * コンストラクタで渡されたinput_pointsがそのまま格納されます。
     */
    std::vector<std::vector<double>> points_;
    
    /**
     * @brief 各点に対応する値データ
     * 
     * points_の各点に対応する値を格納します。
     * 1次元値の場合は内部的に2次元配列に変換されて格納されます。
     * values_[i]はpoints_[i]に対応する値ベクトルです。
     */
    std::vector<std::vector<double>> values_;
    
    /**
     * @brief 補間領域外の点に対して返すデフォルト値
     * 
     * Delaunay三角分割の凸包外にある点を補間する際に返される値です。
     * デフォルトではNaN（std::numeric_limits<double>::quiet_NaN()）が設定されます。
     */
    double fill_value_;
    
    /**
     * @brief Delaunay三角分割を管理するオブジェクト
     * 
     * points_から計算されたDelaunay三角分割を格納します。
     * 補間計算で使用される三角形の検索と重心座標の計算に使用されます。
     * unique_ptrで管理され、遅延初期化されます。
     */
    std::unique_ptr<qhull::Delaunay> tri_;

    /**
     * @brief ドロネー三角分割を計算し、内部のDelaunayオブジェクトを初期化します。
     * 
     * このメソッドは与えられた点群（points_）に対してQhullライブラリを使用して
     * ドロネー三角分割を実行し、結果をtri_メンバに格納します。
     * 補間処理の前に呼び出される必要があります。
     * 
     * @sa https://github.com/scipy/scipy/blob/a0520a1b6b4b3349b8e1fc4d879aae831ecfb193/scipy/interpolate/_interpnd.pyx#L309
     * 
     */
    void calculateTriangulation() 
    {
        tri_ = std::make_unique<qhull::Delaunay>(points_);
    }
    
    /**
     * @brief 補間器の値データを更新する関数
     * 
     * SciPyの `_set_values` に対応する関数で、既存の補間器オブジェクトの
     * 値データを新しいデータで置き換えます。点群データは変更されず、
     * 値データのみが更新されます。
     * 
     * この関数を使用することで、同じ点群に対して異なる値データでの
     * 補間を効率的に実行できます。Delaunay三角分割は再計算されないため、
     * 計算コストを削減できます。
     * 
     * @param input_values 新しい値データ。各要素は点群の対応する点の値を表します。
     *                     点の数と同じサイズである必要があります。
     * 
     * @throws std::invalid_argument 以下の場合に例外を投げます：
     *   - 新しい値データのサイズが点群のサイズと一致しない場合
     *   - その他checkInitShapeで検証される条件に違反する場合
     * 
     * @note この関数を呼び出すと、fill_value_がNaNにリセットされます。
     * @note 点群データを変更したい場合は、新しいオブジェクトを作成してください。
     * 
     * @example
     * ```cpp
     * LinearNDInterpolator interp(points, values1);
     * 
     * interp.setValues(values2);
     * ```
     * 
     * @sa https://github.com/scipy/scipy/blob/a0520a1b6b4b3349b8e1fc4d879aae831ecfb193/scipy/interpolate/_interpnd.pyx#L107
     * 
     */
    void setValues(const std::vector<std::vector<double>>& input_values) 
    {
        checkInitShape(points_, input_values);
        values_ = input_values;
        fill_value_ = std::numeric_limits<double>::quiet_NaN();
    }
    
    /**
     * @brief 初期化時の点群データと値データの形状を検証する静的関数
     * 
     * SciPyの `_check_init_shape` に対応する関数で、コンストラクタで渡される
     * 点群データと値データが適切な形状を持っているかを検証します。
     * 
     * 検証項目：
     * - 点群データが空でないこと
     * - 点群データが最低2次元であること
     * - 全ての点が同じ次元数を持つこと
     * - 点の数と値の数が一致すること
     * 
     * @param p 検証対象の点群データ。各要素は同じ次元数のベクトルである必要があります。
     * @param v 検証対象の値データ。点の数と同じサイズである必要があります。
     * 
     * @throws std::invalid_argument 以下の場合に例外を投げます：
     *   - 点群データが空の場合
     *   - 点群データが2次元未満の場合
     *   - 点の次元数が一致しない場合
     *   - 点の数と値の数が一致しない場合
     * 
     * @note この関数はコンストラクタから呼び出され、不正なデータでの
     *       オブジェクト構築を防ぎます。
     * 
     * @sa https://github.com/scipy/scipy/blob/a0520a1b6b4b3349b8e1fc4d879aae831ecfb193/scipy/interpolate/_interpnd.pyx#L206
     * 
     */
    static void checkInitShape(
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

    /**
     * @brief 補間対象の点群データの形状を検証します。
     * 
     * この関数は補間を実行する前に、入力された点群データが正しい次元数を
     * 持っているかを確認します。補間器の構築時に使用された点群データと
     * 同じ次元数である必要があります。
     * 
     * @param xi 補間対象の点群データ。各要素は補間器の点群と同じ次元数の
     *           ベクトルである必要があります。空のベクトルは許容されます。
     * 
     * @throws std::invalid_argument 補間対象の点の次元数が補間器の点群の
     *                              次元数と一致しない場合に例外を投げます。
     * 
     * @note 空の入力（xi.empty() == true）は有効な入力として扱われ、
     *       例外は投げられません。
     * 
     * @sa https://github.com/scipy/scipy/blob/a0520a1b6b4b3349b8e1fc4d879aae831ecfb193/scipy/interpolate/_interpnd.pyx#L134
     * 
     */
    void checkInterpolateShape(const std::vector<std::vector<double>>& xi) const 
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

    /**
     * @brief 線形N次元補間の実際の計算を実行します。
     * 
     * このメソッドはSciPyの `_do_evaluate` メソッドに相当し、与えられた点群に対して
     * 線形重心座標補間を実行します。各補間点について、対応するシンプレックスを
     * 見つけ、重心座標を計算して線形補間を行います。
     * 
     * アルゴリズムの流れ：
     * 1. 各補間点について、それを含むシンプレックスを探索
     * 2. シンプレックス内での重心座標を計算
     * 3. 重心座標を使用して線形補間を実行
     * 4. シンプレックスが見つからない点にはfill_valueを設定
     * 
     * @param xi 補間対象の点群。各要素は補間器と同じ次元数のベクトルです。
     *           空のベクトルの場合は空の結果を返します。
     * 
     * @return 補間結果の2次元ベクトル。各要素は対応する補間点での値のベクトルです。
     *         出力のサイズは [xi.size() × values_[0].size()] となります。
     * 
     * @note このメソッドは以下のSciPyパラメータを忠実に再現します：
     *       - eps = 100 * DBL_EPSILON（数値計算の許容誤差）
     *       - eps_broad = sqrt(DBL_EPSILON)（広範囲探索の許容誤差）
     *       - シンプレックスヒント機能による効率的な探索
     * 
     * @warning このメソッドを呼び出す前に checkInterpolateShape() で
     *          入力データの妥当性を検証することを推奨します。
     * 
     * @sa https://github.com/scipy/scipy/blob/a0520a1b6b4b3349b8e1fc4d879aae831ecfb193/scipy/interpolate/_interpnd.pyx#L320
     * 
     */
    std::vector<std::vector<double>> evaluate(const std::vector<std::vector<double>>& xi) const {
        // --------------------------------------------------
        // 変数の初期化
        // --------------------------------------------------
        const auto& local_values = values_;
        const double local_fill_value = fill_value_;
        const auto& simplices = tri_->get_simplices();
        
        const size_t n_interp_points = xi.size();
        
        if (n_interp_points == 0) {
            return {};
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

        // [NOTE]
        // Scipyの qhull._get_delaunay_info(&info, self.tri, ...) について、
        // C++では、`tri_`オブジェクト自体が必要な情報をすべて保持しており、
        // `find_simplex`はそのメソッドであるため、このデータ準備ステップは不要。
        // `tri_`オブジェクトが `info` 構造体の役割を果たす。

        // --- Scipyのメインループを忠実に再現 ---
        for (size_t i = 0; i < n_interp_points; ++i) {
            const auto& point = xi[i];

            // 1) Find the simplex
            int isimplex = tri_->findSimplex( //HACK: 要確認
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
            
            std::fill(output[i].begin(), output[i].end(), 0.0);
            
            for (size_t j = 0; j < ndim + 1; ++j) {
                int m = simplices[isimplex][j];
                for (size_t k = 0; k < n_values; ++k) {
                    output[i][k] += barycentric_coords[j] * local_values[m][k];
                }
            }
        }
        return output;
    }

    /**
     * @brief 1次元のvaluesを2次元に変換するヘルパー関数。
     */
    static std::vector<std::vector<double>> _convert_to_2d(const std::vector<double>& v) {
        std::vector<std::vector<double>> v_2d;
        v_2d.reserve(v.size());
        for (const auto& val : v) {
            v_2d.push_back({val});
        }
        return v_2d;
    }
};
