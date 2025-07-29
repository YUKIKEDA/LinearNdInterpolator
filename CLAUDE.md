# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

LinearNdInterpolator is a C++ library for n-dimensional linear interpolation using Qhull for convex hull computation. The project has a minimal structure with core interpolation functionality and uses Qhull as a third-party dependency for geometric computations.

## Architecture

- **Main library**: Located in `include/LinearNdInterpolator.h` and `src/LinearNdInterpolator.cpp`
- **Third-party dependency**: Qhull 2020.2 located in `third_party/qhull-2020.2/`
- **Build system**: CMake-based with Visual Studio project generation

The project follows a standard C++ library structure:
- Header files in `include/`
- Source files in `src/`
- Third-party dependencies in `third_party/`
- Build artifacts in `build/`

## Build Commands

The project uses CMake and is configured for Visual Studio on Windows:

```bash
# Generate Visual Studio project files (from project root)
cmake -B build -G "Visual Studio 17 2022"

# Build the project
cmake --build build

# Or open the generated solution file
# build/Project.sln in Visual Studio
```

### Qhull Dependency

Qhull is included as a third-party dependency with its own CMake configuration:
- Qhull provides convex hull, Delaunay triangulation, and Voronoi diagram computation
- The library is self-contained and builds with standard ANSI C/C++
- Qhull documentation available in `third_party/qhull-2020.2/html/index.htm`

## Development Notes

- The project appears to be in early development stage with minimal source files
- CMake configuration supports both static and shared library builds
- Qhull provides extensive geometric computation capabilities beyond basic convex hulls
- The build system is set up for Windows/Visual Studio but should be portable to other platforms