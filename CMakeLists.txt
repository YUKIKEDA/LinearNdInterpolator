cmake_minimum_required(VERSION 3.16)
project(LinearNdInterpolator VERSION 1.0.0 LANGUAGES CXX)

# C++標準設定
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 出力ディレクトリの設定
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# インクルードディレクトリ
include_directories(${CMAKE_SOURCE_DIR}/include)

# Qhullライブラリの設定
set(QHULL_DIR ${CMAKE_SOURCE_DIR}/third_party/qhull-2020.2)
add_subdirectory(${QHULL_DIR} EXCLUDE_FROM_ALL)

# メインライブラリの作成
set(LIBRARY_SOURCES
    src/LinearNdInterpolator.cpp
    src/Delaunay.cpp
    src/Qhull.cpp
)

set(LIBRARY_HEADERS
    include/LinearNdInterpolator.h
    include/Delaunay.h
    include/Qhull.h
)

# スタティックライブラリとして作成
add_library(LinearNdInterpolator STATIC
    ${LIBRARY_SOURCES}
    ${LIBRARY_HEADERS}
)

# Qhullライブラリとリンク
target_link_libraries(LinearNdInterpolator 
    PRIVATE qhullcpp qhullstatic_r
)

target_include_directories(LinearNdInterpolator 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${QHULL_DIR}/src
        ${QHULL_DIR}/src/libqhullcpp
        ${QHULL_DIR}/src/libqhull_r
)

# テストの設定
option(BUILD_TESTS "Build tests" ON)

if(BUILD_TESTS)
    # Google Testの設定
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    )
    
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    # テストファイルの追加
    enable_testing()

    # SciPy reference tests (auto-generated from JSON data)
    add_executable(test_scipy_reference
        tests/test_scipy_reference_auto_generated.cpp
    )
    
    target_link_libraries(test_scipy_reference
        LinearNdInterpolator
        gtest_main
    )
    
    include(GoogleTest)
    gtest_discover_tests(test_scipy_reference)
endif()


# インストール設定
install(TARGETS LinearNdInterpolator
    EXPORT LinearNdInterpolatorTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(FILES ${LIBRARY_HEADERS} DESTINATION include)