#include "Delaunay.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullError.h"
#include <vector>
#include <string>
#include <iostream>

Delaunay::Delaunay(const std::vector<std::vector<double>>& points) {
    if (points.empty() || points[0].empty()) {
        throw std::invalid_argument("Input points cannot be empty.");
    }
    const size_t n_points = points.size();
    const size_t n_dims = points[0].size();

    std::vector<double> flat_points(n_points * n_dims);
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < n_dims; ++j) {
            flat_points[i * n_dims + j] = points[i][j];
        }
    }

    std::string options = "qhull d Qt Qbb Qc Qz Q12";
    if (n_dims >= 5) {
        options += " Qx";
    }

    try {
        qhull_ = std::make_unique<orgQhull::Qhull>();
        qhull_->runQhull("", n_dims, n_points, flat_points.data(), options.c_str());
        if (qhull_->qhullStatus() != 0) {
            throw std::runtime_error("Qhull failed to run.");
        }
    } catch (const orgQhull::QhullError &e) {
        std::cerr << "Qhull error: " << e.what() << std::endl;
        throw;
    }
}

Delaunay::~Delaunay() = default;

orgQhull::Qhull* Delaunay::getQhull() const {
    return qhull_.get();
}
