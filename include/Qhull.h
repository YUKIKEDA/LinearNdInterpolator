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

// 複数の戻り値を表現するための型エイリアス
using SimplexFacetResult = std::tuple<
    std::vector<std::vector<int>>,    // simplices
    std::vector<std::vector<int>>,    // neighbors
    std::vector<std::vector<double>>, // equations
    std::vector<std::vector<int>>,    // coplanar
    std::vector<bool>                 // good
>;

class Qhull {
private:
    qhT qh_qh; // Qhullコンテキスト
    std::vector<std::vector<double>> points_;
    std::string command_;
    std::string options_;
    bool computed_;
    bool is_delaunay_;
    size_t ndim_;
    
public:
    // コンストラクタ - SciPy _Qhull.__init__を参考に実装
    Qhull(
        const std::string& command, 
        const std::vector<std::vector<double>>& points, 
        const std::string& options);
    
    ~Qhull();

    // SciPy _Qhull.triangulate()に対応
    void triangulate();

    // SciPy _Qhull.get_paraboloid_shift_scale()に対応
    std::pair<double, double> getParaboloidShiftScale() const;

    // SciPy _Qhull.get_simplex_facet_array()に対応
    SimplexFacetResult getSimplexFacetArray() const;
};

}