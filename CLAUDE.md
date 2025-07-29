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

## Implementation Strategy

### Phase 1: MVP (2D)
- Basic 2D triangulation using Qhull
- Linear search for triangle location
- Barycentric coordinate interpolation
- Nearest neighbor fallback for points outside convex hull

### Phase 2: N-Dimensional
- Generalize to arbitrary dimensions
- N-dimensional simplex handling
- Robust numerical computation

### Phase 3: Production Ready
- Comprehensive error handling
- Edge case coverage
- Test suite completion

## Key Technical Details

- **Interpolation Method**: Linear interpolation using barycentric coordinates within simplices
- **Outside Hull Handling**: Returns value of nearest data point (not NaN or fill values)
- **Numerical Stability**: Uses appropriate tolerances for geometric computations
- **Memory Layout**: Efficient data structures for point storage and simplex representation