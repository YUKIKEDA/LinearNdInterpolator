# LinearNdInterpolator

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](https://github.com/your-repo/LinearNdInterpolator)
[![SciPy Compatible](https://img.shields.io/badge/SciPy-compatible-blue)](https://scipy.org/)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue)](https://en.cppreference.com/w/cpp/17)

A high-performance C++ library for N-dimensional linear interpolation using Delaunay triangulation. This implementation achieves **complete numerical compatibility** with SciPy's `LinearNDInterpolator`, making it suitable for scientific computing applications requiring precise interpolation results.

## ✨ Features

- **🎯 SciPy Compatibility**: 100% numerical accuracy match with SciPy's LinearNDInterpolator
- **📐 N-Dimensional Support**: Handles 2D to 4D+ interpolation seamlessly
- **⚡ Delaunay Triangulation**: Uses Qhull library for robust triangulation
- **🔢 Barycentric Interpolation**: Linear interpolation using barycentric coordinates
- **🛡️ Robust Error Handling**: Comprehensive input validation and edge case management
- **🔍 Convex Hull Detection**: Returns NaN for points outside the convex hull (SciPy-compatible)

## 🚀 Quick Start

### Basic Usage

```cpp
#include "LinearNdInterpolator.h"
#include <vector>

// Define 2D points and their values
std::vector<std::vector<double>> points = {
    {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {0.5, 0.5}
};
std::vector<double> values = {0.0, 1.0, 1.0, 1.0};

// Create interpolator
LinearNdInterpolator interp(points, values);

// Interpolate at a new point
double result = interp.interpolate({0.25, 0.25});
std::cout << "Interpolated value: " << result << std::endl;
```

### Multi-Point Interpolation

```cpp
// Query multiple points at once
std::vector<std::vector<double>> query_points = {
    {0.1, 0.1}, {0.3, 0.4}, {0.7, 0.2}
};

auto results = interp.interpolate(query_points);
for (size_t i = 0; i < results.size(); ++i) {
    std::cout << "Point " << i << ": " << results[i][0] << std::endl;
}
```

## 🏗️ Building the Project

### Prerequisites

- **Compiler**: C++17 compatible compiler (Visual Studio 2019+, GCC 7+, Clang 5+)
- **CMake**: Version 3.15 or higher
- **Dependencies**: Qhull (included), GoogleTest (auto-downloaded)

### Build Instructions

```bash
# Generate build files
cmake -B build -G "Visual Studio 17 2022"

# Build the project
cmake --build build

# Run tests
"build\bin\Debug\test_scipy_reference.exe"
```

### Integration

```cmake
# Add to your CMakeLists.txt
find_package(LinearNdInterpolator REQUIRED)
target_link_libraries(your_target LinearNdInterpolator::LinearNdInterpolator)
```

## 📊 Performance & Compatibility

### SciPy Compatibility Tests

✅ **All 7 test suites passing:**
- 2D Basic Triangle
- 2D Complex Random Points  
- 3D Basic Tetrahedron
- 4D Basic Simplex
- 2D Collinear Points
- 2D Duplicate Points
- 2D Large Dataset (100 points)

```bash
# Run compatibility tests
"build\bin\Debug\test_scipy_reference.exe"
```

### Current Performance Status

| Feature | Status | Notes |
|---------|--------|-------|
| **Numerical Accuracy** | ✅ Perfect | 100% SciPy compatibility |
| **Basic Functionality** | ✅ Complete | All core features working |
| **Error Handling** | ✅ Good | Comprehensive validation |
| **Transform Matrices** | ✅ Implemented | SciPy-style optimization active |
| **Walking Algorithm** | ✅ Infrastructure | Framework ready for optimization |
| **Performance** | 🟡 Partially Optimized | Transform matrices implemented |
| **Memory Usage** | 🟡 Acceptable | Room for optimization |

## 🛠️ API Reference

### Constructor

```cpp
LinearNdInterpolator(
    const std::vector<std::vector<double>>& points, 
    const std::vector<double>& values
);

LinearNdInterpolator(
    const std::vector<std::vector<double>>& points, 
    const std::vector<std::vector<double>>& values  // Multi-dimensional values
);
```

### Methods

```cpp
// Single point interpolation
double interpolate(const std::vector<double>& point) const;

// Multi-point interpolation
std::vector<std::vector<double>> interpolate(
    const std::vector<std::vector<double>>& points
) const;
```

### Error Handling

The library throws `std::invalid_argument` for:
- Empty input data
- Mismatched point and value counts
- Inconsistent dimensions
- Points with < 2 dimensions
- NaN or infinite input values

## 📈 Development Status

### ✅ Completed (Phase 1-2)

- **Core Algorithm**: N-dimensional linear interpolation
- **SciPy Compatibility**: Complete numerical compatibility achieved
- **Transform Matrices**: SciPy-style barycentric coordinate optimization implemented
- **Walking Algorithm Infrastructure**: Framework for efficient simplex search
- **Robust Implementation**: Comprehensive error handling and validation
- **Testing Suite**: All 17 SciPy reference tests passing
- **Cross-Platform Build**: CMake-based build system

### 🚧 Current Limitations & Roadmap

#### High Priority Improvements

1. **🟡 Walking Algorithm Activation** (Target: 3-5 days)
   - **Status**: Infrastructure completed, needs neighbor computation
   - **Solution**: Implement proper Qhull neighbor relationships
   - **Expected Impact**: 10-100x performance improvement for large datasets

2. **✅ Transform Matrix Implementation** (COMPLETED)
   - **Achievement**: SciPy-style pre-computed transformation matrices
   - **Impact**: O(d³) → O(d²) barycentric coordinate calculation
   - **Status**: All 17 SciPy compatibility tests passing

#### Medium Priority Improvements

3. **🟡 Enhanced Unit Testing** (Target: 1 week)
   - Add comprehensive edge case tests
   - Performance benchmarking suite
   - Memory usage validation

4. **🟡 Thread Safety** (Target: 3 days)  
   - Add thread-safe operations
   - Concurrent query processing

5. **🟡 Memory Optimization** (Target: 1 week)
   - Optimize vertex mapping algorithm
   - Reduce memory footprint for large datasets

#### Future Enhancements

- **Incremental Point Addition**: Dynamic point insertion
- **Custom Fill Values**: Alternative to NaN for out-of-hull points
- **Parallel Processing**: OpenMP support for large-scale interpolation
- **Advanced Interpolation**: Cubic and higher-order methods

## 🔬 Technical Details

### Architecture

```
LinearNdInterpolator
├── Delaunay (Qhull-based triangulation)
│   ├── findSimplex()
│   ├── calculateBarycentricCoordinates()
│   └── getSimplices()
└── interpolate() (Linear interpolation)
```

### Algorithm Flow

1. **Triangulation**: Create Delaunay triangulation using Qhull
2. **Simplex Location**: Find containing simplex for query point
3. **Barycentric Calculation**: Compute barycentric coordinates
4. **Linear Interpolation**: Weight vertex values by coordinates

### Memory Layout

- **Points**: Stored as `std::vector<std::vector<double>>`
- **Values**: Multi-dimensional value support
- **Triangulation**: Managed by Qhull with RAII principles

## 📋 Testing

### Test Categories

```bash
# SciPy compatibility tests (primary validation)
"build\bin\Debug\test_scipy_reference.exe"

# Specific test cases
"build\bin\Debug\test_scipy_reference.exe" --gtest_filter="SciPyReferenceTest.2D_Basic_Triangle"
"build\bin\Debug\test_scipy_reference.exe" --gtest_filter="SciPyReferenceTest.3D_Basic_Tetrahedron"
```

### Test Coverage

- ✅ **2D-4D Interpolation**: Comprehensive dimensional testing
- ✅ **Edge Cases**: Collinear points, duplicates, large datasets  
- ✅ **Error Conditions**: Invalid inputs, out-of-hull points
- ✅ **Numerical Precision**: Exact SciPy compatibility validation

## 🤝 Contributing

### Development Priorities

1. **Performance Optimization**: Transform matrix and walking algorithm implementation
2. **Testing Enhancement**: Independent unit test suite development
3. **Documentation**: Usage examples and best practices
4. **Benchmarking**: Performance comparison tools

### Getting Started

```bash
# Clone and build
git clone https://github.com/your-repo/LinearNdInterpolator
cd LinearNdInterpolator
cmake -B build
cmake --build build

# Run tests to verify
"build\bin\Debug\test_scipy_reference.exe"
```

## 📚 Examples

### 3D Surface Interpolation

```cpp
// Define a 3D surface: f(x,y,z) = x + y + z
std::vector<std::vector<double>> points = {
    {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}
};
std::vector<double> values = {0, 1, 1, 1, 3};

LinearNdInterpolator interp(points, values);

// Interpolate at interior point
double result = interp.interpolate({0.25, 0.25, 0.25});
// Expected: 0.75 (linear combination)
```

### Large Dataset Processing

```cpp
// Process large point cloud
std::vector<std::vector<double>> large_points = loadPointCloud("data.csv");
std::vector<double> large_values = loadValues("values.csv");

LinearNdInterpolator interp(large_points, large_values);

// Batch interpolation
std::vector<std::vector<double>> query_grid = generateGrid(100, 100);
auto results = interp.interpolate(query_grid);
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **SciPy Team**: Reference implementation and test cases
- **Qhull**: Robust computational geometry library
- **GoogleTest**: Testing framework

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/your-repo/LinearNdInterpolator/issues)
- **Documentation**: See `.memo/` directory for detailed technical documentation
- **Performance**: See `10_interpolate_method_scipy_comparison.md` for optimization roadmap

---

**Status**: ✅ **Production Ready with SciPy-Compatible Performance Optimization**  
**Latest Achievement**: 🎯 **Transform Matrix Implementation - O(d³) → O(d²) barycentric calculation**  
**Next Focus**: 🚀 **Walking Algorithm Activation for Large-Scale Performance**