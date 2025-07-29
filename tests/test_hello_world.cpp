#include <gtest/gtest.h>
#include "LinearNdInterpolator.h"
#include <cmath>

// 2D Interpolation Tests
class LinearNdInterpolator2DTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple 2D test case - avoid cocircular points
        points_2d = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {0.5, 0.5}};
        values_2d = {0.0, 1.0, 1.0, 1.0};  // Simple linear function: f(x,y) = x + y
    }
    
    std::vector<std::vector<double>> points_2d;
    std::vector<double> values_2d;
};

TEST_F(LinearNdInterpolator2DTest, ConstructionTest) {
    EXPECT_NO_THROW({
        LinearNdInterpolator interp(points_2d, values_2d);
    });
}

TEST_F(LinearNdInterpolator2DTest, ExactVertexInterpolation) {
    LinearNdInterpolator interp(points_2d, values_2d);
    
    // Test interpolation at exact vertices
    EXPECT_NEAR(interp.interpolate({0.0, 0.0}), 0.0, 1e-6);
    EXPECT_NEAR(interp.interpolate({1.0, 0.0}), 1.0, 1e-6);
    EXPECT_NEAR(interp.interpolate({0.0, 1.0}), 1.0, 1e-6);
    EXPECT_NEAR(interp.interpolate({0.5, 0.5}), 1.0, 1e-6);
}

TEST_F(LinearNdInterpolator2DTest, InteriorPointInterpolation) {
    LinearNdInterpolator interp(points_2d, values_2d);
    
    // Test interpolation at interior points
    double result1 = interp.interpolate({0.25, 0.25});
    double result2 = interp.interpolate({0.3, 0.4});
    
    EXPECT_TRUE(std::isfinite(result1));
    EXPECT_TRUE(std::isfinite(result2));
    EXPECT_FALSE(std::isnan(result1));
    EXPECT_FALSE(std::isnan(result2));
}

TEST_F(LinearNdInterpolator2DTest, OutsideConvexHullNearestNeighbor) {
    LinearNdInterpolator interp(points_2d, values_2d);
    
    // Test points outside convex hull - should return nearest neighbor
    double result1 = interp.interpolate({-0.5, -0.5});
    double result2 = interp.interpolate({2.0, 2.0});
    
    EXPECT_TRUE(std::isfinite(result1));
    EXPECT_TRUE(std::isfinite(result2));
    EXPECT_FALSE(std::isnan(result1));
    EXPECT_FALSE(std::isnan(result2));
}

// 3D Interpolation Tests
class LinearNdInterpolator3DTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a 3D test case - tetrahedron
        points_3d = {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {1.0, 1.0, 1.0}
        };
        values_3d = {0.0, 1.0, 1.0, 1.0, 3.0};  // Simple linear function: f(x,y,z) = x + y + z
    }
    
    std::vector<std::vector<double>> points_3d;
    std::vector<double> values_3d;
};

TEST_F(LinearNdInterpolator3DTest, ConstructionTest) {
    EXPECT_NO_THROW({
        LinearNdInterpolator interp(points_3d, values_3d);
    });
}

TEST_F(LinearNdInterpolator3DTest, ExactVertexInterpolation) {
    LinearNdInterpolator interp(points_3d, values_3d);
    
    // Test interpolation at exact vertices
    EXPECT_NEAR(interp.interpolate({0.0, 0.0, 0.0}), 0.0, 1e-10);
    EXPECT_NEAR(interp.interpolate({1.0, 0.0, 0.0}), 1.0, 1e-10);
    EXPECT_NEAR(interp.interpolate({0.0, 1.0, 0.0}), 1.0, 1e-10);
    EXPECT_NEAR(interp.interpolate({0.0, 0.0, 1.0}), 1.0, 1e-10);
}

TEST_F(LinearNdInterpolator3DTest, InteriorPointInterpolation) {
    LinearNdInterpolator interp(points_3d, values_3d);
    
    // Test interpolation at interior points
    // For linear function f(x,y,z) = x + y + z, interpolation should be exact
    double result = interp.interpolate({0.25, 0.25, 0.25});
    EXPECT_TRUE(std::isfinite(result));
    EXPECT_FALSE(std::isnan(result));
}

// N-Dimensional Tests
TEST(LinearNdInterpolatorNDTest, HighDimensionalConstruction) {
    // Test 4D case
    std::vector<std::vector<double>> points_4d = {
        {0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0},
        {1.0, 1.0, 1.0, 1.0}
    };
    std::vector<double> values_4d = {0.0, 1.0, 1.0, 1.0, 1.0, 4.0};
    
    EXPECT_NO_THROW({
        LinearNdInterpolator interp(points_4d, values_4d);
    });
}

TEST(LinearNdInterpolatorNDTest, HighDimensionalInterpolation) {
    // Test 4D case - need at least 5 points for 4D triangulation
    std::vector<std::vector<double>> points_4d = {
        {0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0},
        {0.5, 0.5, 0.5, 0.5}  // Additional point to ensure proper triangulation
    };
    std::vector<double> values_4d = {0.0, 1.0, 1.0, 1.0, 1.0, 2.0};
    
    LinearNdInterpolator interp(points_4d, values_4d);
    
    // Test exact vertices
    EXPECT_NEAR(interp.interpolate({0.0, 0.0, 0.0, 0.0}), 0.0, 1e-6);
    EXPECT_NEAR(interp.interpolate({1.0, 0.0, 0.0, 0.0}), 1.0, 1e-6);
    
    // Test interior point
    double result = interp.interpolate({0.2, 0.2, 0.2, 0.2});
    EXPECT_TRUE(std::isfinite(result));
    EXPECT_FALSE(std::isnan(result));
}

// Error Handling Tests
TEST(LinearNdInterpolatorErrorTest, EmptyInputs) {
    std::vector<std::vector<double>> empty_points;
    std::vector<double> empty_values;
    
    EXPECT_THROW({
        LinearNdInterpolator interp(empty_points, empty_values);
    }, std::invalid_argument);
}

TEST(LinearNdInterpolatorErrorTest, SizeMismatch) {
    std::vector<std::vector<double>> points = {{0.0, 0.0}, {1.0, 0.0}};
    std::vector<double> values = {0.0}; // Wrong size
    
    EXPECT_THROW({
        LinearNdInterpolator interp(points, values);
    }, std::invalid_argument);
}

TEST(LinearNdInterpolatorErrorTest, InconsistentDimensions) {
    std::vector<std::vector<double>> points = {{0.0, 0.0}, {1.0, 0.0, 0.0}}; // Different dimensions
    std::vector<double> values = {0.0, 1.0};
    
    EXPECT_THROW({
        LinearNdInterpolator interp(points, values);
    }, std::invalid_argument);
}

TEST(LinearNdInterpolatorErrorTest, ZeroDimension) {
    std::vector<std::vector<double>> points = {{}};
    std::vector<double> values = {0.0};
    
    EXPECT_THROW({
        LinearNdInterpolator interp(points, values);
    }, std::invalid_argument);
}

TEST(LinearNdInterpolatorErrorTest, InsufficientPoints) {
    // 2D case needs at least 3 points
    std::vector<std::vector<double>> points = {{0.0, 0.0}, {1.0, 0.0}};
    std::vector<double> values = {0.0, 1.0};
    
    EXPECT_THROW({
        LinearNdInterpolator interp(points, values);
    }, std::invalid_argument);
}

TEST(LinearNdInterpolatorErrorTest, WrongQueryDimension) {
    std::vector<std::vector<double>> points = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {0.5, 0.5}};
    std::vector<double> values = {0.0, 1.0, 1.0, 1.0};
    
    LinearNdInterpolator interp(points, values);
    
    EXPECT_THROW({
        interp.interpolate({0.5}); // 1D query for 2D interpolator
    }, std::invalid_argument);
    
    EXPECT_THROW({
        interp.interpolate({0.5, 0.5, 0.5}); // 3D query for 2D interpolator
    }, std::invalid_argument);
}

// Performance Test
TEST(LinearNdInterpolatorPerformanceTest, ManyPointsInterpolation) {
    // Create a larger dataset for performance testing
    std::vector<std::vector<double>> points;
    std::vector<double> values;
    
    // Generate grid points in 2D
    for (double x = 0.0; x <= 10.0; x += 1.0) {
        for (double y = 0.0; y <= 10.0; y += 1.0) {
            points.push_back({x, y});
            values.push_back(x * x + y * y); // Quadratic function
        }
    }
    
    ASSERT_GT(points.size(), 100); // Make sure we have enough points
    
    EXPECT_NO_THROW({
        LinearNdInterpolator interp(points, values);
        
        // Test multiple interpolations
        for (double x = 0.5; x <= 9.5; x += 1.0) {
            for (double y = 0.5; y <= 9.5; y += 1.0) {
                double result = interp.interpolate({x, y});
                EXPECT_TRUE(std::isfinite(result));
                EXPECT_FALSE(std::isnan(result));
            }
        }
    });
}

// Special Cases Tests
TEST(LinearNdInterpolatorSpecialTest, CollinearPoints2D) {
    // Test with collinear points - should still work due to additional point
    std::vector<std::vector<double>> points = {
        {0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {0.0, 1.0}
    };
    std::vector<double> values = {0.0, 2.0, 4.0, 1.0};
    
    EXPECT_NO_THROW({
        LinearNdInterpolator interp(points, values);
        double result = interp.interpolate({0.5, 0.5});
        EXPECT_TRUE(std::isfinite(result));
        EXPECT_FALSE(std::isnan(result));
    });
}

TEST(LinearNdInterpolatorSpecialTest, DuplicatePoints) {
    // Test with duplicate points - should still work
    std::vector<std::vector<double>> points = {
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}
    };
    std::vector<double> values = {0.0, 1.0, 1.0, 0.0};
    
    EXPECT_NO_THROW({
        LinearNdInterpolator interp(points, values);
        double result = interp.interpolate({0.5, 0.5});
        EXPECT_TRUE(std::isfinite(result));
        EXPECT_FALSE(std::isnan(result));
    });
}