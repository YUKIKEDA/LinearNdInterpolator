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

private:
    std::unique_ptr<orgQhull::Qhull> qhull_;
};
