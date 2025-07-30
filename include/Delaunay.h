#pragma once

#include <vector>
#include <memory>

// Qhullの前方宣言
namespace orgQhull {
    class Qhull;
}

class Delaunay {
public:
    explicit Delaunay(const std::vector<std::vector<double>>& points);
    ~Delaunay();

    // コピーとムーブは禁止
    Delaunay(const Delaunay&) = delete;
    Delaunay& operator=(const Delaunay&) = delete;
    Delaunay(Delaunay&&) = delete;
    Delaunay& operator=(Delaunay&&) = delete;

    // LinearNdInterpolatorがQhullの機能にアクセスするためのゲッター
    orgQhull::Qhull* getQhull() const;
    
    // SciPyと同等の補間支援メソッド
    int findSimplex(const std::vector<double>& point) const;
    std::vector<double> calculateBarycentricCoordinates(
        const std::vector<double>& point, int simplex_id) const;
    std::vector<std::vector<int>> getSimplices() const;

private:
    std::unique_ptr<orgQhull::Qhull> qhull_;
    
    // ヘルパーメソッド
    std::vector<double> solveLinearSystem(
        const std::vector<std::vector<double>>& A, 
        const std::vector<double>& b) const;
};
