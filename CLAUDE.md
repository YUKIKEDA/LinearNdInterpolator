# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

LinearNdInterpolator is a C++ library for N-dimensional linear interpolation using Qhull for Delaunay triangulation. The library performs linear interpolation within convex hulls using barycentric coordinates, and returns nearest neighbor values for points outside the convex hull.

## Architecture

### Core Implementation
- **Main library**: `include/LinearNdInterpolator.h` and `src/LinearNdInterpolator.cpp`
- **API Design**: Simple constructor-based approach with single `interpolate()` method
- **Algorithm**: Delaunay triangulation → Simplex location → Barycentric interpolation or nearest neighbor

### Dependencies
- **Qhull 2020.2**: Located in `third_party/qhull-2020.2/` for Delaunay triangulation
- **Eigen 3.4.0**: Located in `third_party/eigen-3.4.0/` for numerical computations (if used)

### Design Documents
The `.memo/` directory contains detailed design documentation:
- `01_project_overview.md`: Goals and implementation phases
- `02_theoretical_background.md`: Mathematical foundation of N-dimensional linear interpolation
- `03_qhull_api_analysis.md`: Qhull C++ API usage details
- `04_api_design.md`: Data structures and API specifications
- `05_implementation_plan.md`: Three-phase implementation strategy
- `06_implementation_notes.md`: Technical implementation details
- `07_scipy_compatibility_analysis.md`: SciPy compatibility analysis and test failure investigation
- `08_current_status_and_remaining_tasks.md`: Current implementation status and remaining tasks

## Planned API

```cpp
class LinearNdInterpolator {
public:
    LinearNdInterpolator(const std::vector<std::vector<double>>& points, 
                        const std::vector<double>& values);
    double interpolate(const std::vector<double>& point) const;
};
```

## Build Commands

```bash
# Generate Visual Studio project files (from project root)
cmake -B build -G "Visual Studio 17 2022"

# Build the project
cmake --build build

# Open generated solution file in Visual Studio
# build/Project.sln
```

## Implementation Status

### ✅ Phase 1: MVP (2D) - COMPLETED
- ✅ Basic 2D triangulation using Qhull
- ✅ Barycentric coordinate interpolation (SciPy-compatible)
- ✅ Correct vertex index mapping
- ✅ NaN handling for points outside convex hull (SciPy-compatible)

### ✅ Phase 2: N-Dimensional - COMPLETED  
- ✅ Generalized to arbitrary dimensions (2D-4D tested)
- ✅ N-dimensional simplex handling
- ✅ SciPy numerical compatibility achieved

### 🟡 Phase 3: Production Ready - IN PROGRESS
- ✅ SciPy compatibility tests (7/7 passing)
- 🟡 Performance optimization needed
- 🟡 Comprehensive error handling partially complete
- ❌ Independent unit test suite needs implementation

## Key Technical Details

- **Interpolation Method**: Linear interpolation using barycentric coordinates within simplices
- **Outside Hull Handling**: Returns NaN for points outside convex hull (SciPy-compatible)
- **Numerical Stability**: Uses appropriate tolerances for geometric computations
- **Vertex Mapping**: Coordinate-based matching between Qhull output and original points
- **SciPy Compatibility**: Passes all reference tests with exact numerical precision

## Test Commands

```bash
# Run SciPy compatibility tests
"build\bin\Debug\test_scipy_reference.exe"

# Run specific test case
"build\bin\Debug\test_scipy_reference.exe" --gtest_filter="SciPyReferenceTest.2D_Basic_Triangle"
```

## Current Status Summary

**✅ CORE FUNCTIONALITY COMPLETE**: The implementation successfully passes all SciPy compatibility tests and correctly handles N-dimensional linear interpolation. The critical vertex index mapping and barycentric coordinate calculation issues have been resolved.

**Next Priority**: Performance optimization and comprehensive unit test suite development.