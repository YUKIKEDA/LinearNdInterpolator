#!/usr/bin/env python3
"""
SciPy LinearNDInterpolator テストデータ生成スクリプト

このスクリプトはSciPyのLinearNDInterpolatorを使用して、
C++実装との動作一致性を検証するためのテストデータを生成します。
"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator
import json
import os

def generate_2d_test_data():
    """2Dテストデータの生成"""
    print("=== 2D Test Data Generation ===")
    
    # テストケース1: 基本的な三角形
    points_2d_basic = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 0.5]
    ])
    values_2d_basic = np.array([0.0, 1.0, 1.0, 1.0])
    
    interp_2d_basic = LinearNDInterpolator(points_2d_basic, values_2d_basic)
    
    # クエリポイント
    query_points_2d = np.array([
        [0.0, 0.0],      # 頂点
        [1.0, 0.0],      # 頂点
        [0.0, 1.0],      # 頂点
        [0.5, 0.5],      # 頂点
        [0.25, 0.25],    # 内部点
        [0.3, 0.4],      # 内部点
        [0.75, 0.25],    # 内部点
        [-0.5, -0.5],    # 凸包外
        [2.0, 2.0],      # 凸包外
        [0.1, 0.1],      # 内部点（境界近く）
    ])
    
    results_2d_basic = []
    for point in query_points_2d:
        result = interp_2d_basic(point)
        results_2d_basic.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    test_case_2d_basic = {
        "name": "2D_Basic_Triangle",
        "description": "Basic 2D triangle interpolation",
        "input_points": points_2d_basic.tolist(),
        "input_values": values_2d_basic.tolist(),
        "queries": results_2d_basic
    }
    
    # テストケース2: より複雑な2D形状
    np.random.seed(42)
    points_2d_complex = np.random.rand(10, 2) * 2.0 - 1.0  # -1から1の範囲
    values_2d_complex = points_2d_complex[:, 0] + points_2d_complex[:, 1]  # x + y
    
    interp_2d_complex = LinearNDInterpolator(points_2d_complex, values_2d_complex)
    
    query_points_2d_complex = np.array([
        [0.0, 0.0],
        [0.5, 0.3],
        [-0.2, 0.4],
        [0.8, -0.1],
        [2.0, 2.0],  # 凸包外
        [-2.0, -2.0]  # 凸包外
    ])
    
    results_2d_complex = []
    for point in query_points_2d_complex:
        result = interp_2d_complex(point)
        results_2d_complex.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    test_case_2d_complex = {
        "name": "2D_Complex_Random",
        "description": "Complex 2D random points with f(x,y) = x + y",
        "input_points": points_2d_complex.tolist(),
        "input_values": values_2d_complex.tolist(),
        "queries": results_2d_complex
    }
    
    return [test_case_2d_basic, test_case_2d_complex]

def generate_3d_test_data():
    """3Dテストデータの生成"""
    print("=== 3D Test Data Generation ===")
    
    # テストケース1: 基本的な四面体
    points_3d_basic = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 1.0, 1.0]
    ])
    values_3d_basic = np.array([0.0, 1.0, 1.0, 1.0, 3.0])
    
    interp_3d_basic = LinearNDInterpolator(points_3d_basic, values_3d_basic)
    
    query_points_3d = np.array([
        [0.0, 0.0, 0.0],    # 頂点
        [1.0, 0.0, 0.0],    # 頂点
        [0.0, 1.0, 0.0],    # 頂点
        [0.0, 0.0, 1.0],    # 頂点
        [0.25, 0.25, 0.25], # 内部点
        [0.1, 0.2, 0.3],    # 内部点
        [2.0, 2.0, 2.0],    # 凸包外
        [-1.0, -1.0, -1.0]  # 凸包外
    ])
    
    results_3d_basic = []
    for point in query_points_3d:
        result = interp_3d_basic(point)
        results_3d_basic.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    test_case_3d_basic = {
        "name": "3D_Basic_Tetrahedron",
        "description": "Basic 3D tetrahedron interpolation",
        "input_points": points_3d_basic.tolist(),
        "input_values": values_3d_basic.tolist(),
        "queries": results_3d_basic
    }
    
    return [test_case_3d_basic]

def generate_4d_test_data():
    """4Dテストデータの生成"""
    print("=== 4D Test Data Generation ===")
    
    test_cases_4d = []
    
    # 4D基本テスト
    points_4d = np.array([
        [0.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.5, 0.5, 0.5, 0.5]
    ])
    values_4d = np.array([0.0, 1.0, 1.0, 1.0, 1.0, 2.0])
    
    interp_4d = LinearNDInterpolator(points_4d, values_4d)
    
    query_points_4d = np.array([
        [0.0, 0.0, 0.0, 0.0],    # 頂点
        [1.0, 0.0, 0.0, 0.0],    # 頂点
        [0.5, 0.5, 0.5, 0.5],    # 頂点
        [0.2, 0.2, 0.2, 0.2],    # 内部点
        [0.1, 0.3, 0.2, 0.4],    # 内部点
        [2.0, 2.0, 2.0, 2.0]     # 凸包外
    ])
    
    results_4d = []
    for point in query_points_4d:
        result = interp_4d(point)
        results_4d.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    test_cases_4d.append({
        "name": "4D_Basic_Simplex",
        "description": "Basic 4D simplex interpolation",
        "input_points": points_4d.tolist(),
        "input_values": values_4d.tolist(),
        "queries": results_4d
    })
    
    # 4Dベクトル値テスト
    values_4d_vector = np.array([
        [0.0, 0.0],
        [1.0, 2.0],
        [1.0, -1.0],
        [1.0, 0.5],
        [1.0, 1.5],
        [2.0, 1.0]
    ])
    
    interp_4d_vector = LinearNDInterpolator(points_4d, values_4d_vector)
    
    results_4d_vector = []
    for point in query_points_4d:
        result = interp_4d_vector(point)
        results_4d_vector.append({
            "query_point": point.tolist(),
            "result": result.flatten().tolist() if not np.any(np.isnan(result)) else None
        })
    
    test_cases_4d.append({
        "name": "4D_Vector_Values",
        "description": "4D interpolation with 2D vector values",
        "input_points": points_4d.tolist(),
        "input_values": values_4d_vector.tolist(),
        "queries": results_4d_vector
    })
    
    return test_cases_4d

def generate_5d_and_higher_test_data():
    """5次元以上のテストデータ生成（C++実装で追加オプション"Qx"が使用される）"""
    print("=== 5D and Higher Dimensional Test Data Generation ===")
    
    test_cases_high_dim = []
    
    # 5次元テスト（C++では"Qx"オプションが追加される）
    points_5d = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0],
        [0.2, 0.2, 0.2, 0.2, 0.2]
    ])
    values_5d = points_5d.sum(axis=1)  # 座標の合計値
    
    interp_5d = LinearNDInterpolator(points_5d, values_5d)
    
    query_points_5d = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0],    # 頂点
        [0.2, 0.2, 0.2, 0.2, 0.2],    # 頂点
        [0.1, 0.1, 0.1, 0.1, 0.1],    # 内部点
        [0.15, 0.15, 0.15, 0.15, 0.15], # 内部点
        [2.0, 2.0, 2.0, 2.0, 2.0]     # 凸包外
    ])
    
    results_5d = []
    for point in query_points_5d:
        result = interp_5d(point)
        results_5d.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    test_cases_high_dim.append({
        "name": "5D_Basic_Simplex",
        "description": "5D interpolation (triggers Qx option in C++)",
        "input_points": points_5d.tolist(),
        "input_values": values_5d.tolist(),
        "queries": results_5d
    })
    
    # 6次元テスト
    np.random.seed(789)
    points_6d = np.random.rand(15, 6) * 2.0 - 1.0  # -1から1の範囲
    values_6d = np.sum(points_6d**2, axis=1)  # 二乗和
    
    interp_6d = LinearNDInterpolator(points_6d, values_6d)
    
    query_points_6d = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.1, 0.2, -0.1, 0.3, -0.2, 0.4],
        [2.0, 2.0, 2.0, 2.0, 2.0, 2.0]  # 凸包外
    ])
    
    results_6d = []
    for point in query_points_6d:
        result = interp_6d(point)
        results_6d.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    test_cases_high_dim.append({
        "name": "6D_Random_Points",
        "description": "6D interpolation with random points",
        "input_points": points_6d.tolist(),
        "input_values": values_6d.tolist(),
        "queries": results_6d
    })
    
    return test_cases_high_dim

def generate_numerical_precision_test_data():
    """数値精度テストデータ生成（C++実装の1e-10許容範囲を検証）"""
    print("=== Numerical Precision Test Data Generation ===")
    
    test_cases_precision = []
    
    # 境界付近の数値精度テスト
    points_precision = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 0.5]
    ])
    values_precision = np.array([0.0, 1.0, 1.0, 1.0])
    
    interp_precision = LinearNDInterpolator(points_precision, values_precision)
    
    # C++実装の許容範囲1e-10周辺の点をテスト
    epsilon = 1e-10
    query_points_precision = np.array([
        [epsilon, epsilon],              # 1e-10の微小値
        [-epsilon, -epsilon],            # 負の1e-10
        [0.5 + epsilon, 0.5 - epsilon],  # 境界付近
        [0.5 - epsilon, 0.5 + epsilon],  # 境界付近
        [1.0 - epsilon, epsilon],        # 頂点近傍
        [epsilon, 1.0 - epsilon],        # 頂点近傍
        [0.333333333333, 0.333333333333], # 長い小数
        [1.0/3.0, 1.0/3.0],              # 1/3の表現
        [0.1 + 0.2, 0.3],                # 浮動小数点誤差が生じやすい計算
    ])
    
    results_precision = []
    for point in query_points_precision:
        result = interp_precision(point)
        results_precision.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    test_cases_precision.append({
        "name": "Numerical_Precision_Boundary",
        "description": "Numerical precision test around 1e-10 tolerance",
        "input_points": points_precision.tolist(),
        "input_values": values_precision.tolist(),
        "queries": results_precision
    })
    
    # 大きな値での精度テスト
    points_large = np.array([
        [1e6, 1e6],
        [1e6 + 1, 1e6],
        [1e6, 1e6 + 1],
        [1e6 + 0.5, 1e6 + 0.5]
    ])
    values_large = np.array([0.0, 1.0, 1.0, 1.0])
    
    interp_large = LinearNDInterpolator(points_large, values_large)
    
    query_points_large = np.array([
        [1e6 + 0.25, 1e6 + 0.25],       # 内部点
        [1e6 + 0.1, 1e6 + 0.1],         # 内部点
        [1e6 - 1, 1e6 - 1]              # 凸包外（スケールが大きい）
    ])
    
    results_large = []
    for point in query_points_large:
        result = interp_large(point)
        results_large.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    test_cases_precision.append({
        "name": "Large_Scale_Precision",
        "description": "Numerical precision test with large coordinate values",
        "input_points": points_large.tolist(),
        "input_values": values_large.tolist(),
        "queries": results_large
    })
    
    return test_cases_precision

def generate_boundary_test_data():
    """境界条件テストデータ生成"""
    print("=== Boundary Condition Test Data Generation ===")
    
    test_cases_boundary = []
    
    # 境界上の点のテスト
    points_boundary = np.array([
        [0.0, 0.0],
        [2.0, 0.0],
        [1.0, 2.0],
        [1.0, 1.0]
    ])
    values_boundary = np.array([1.0, 2.0, 3.0, 2.5])
    
    interp_boundary = LinearNDInterpolator(points_boundary, values_boundary)
    
    # 三角形の辺上の点をテスト
    query_points_boundary = np.array([
        [1.0, 0.0],     # 辺上
        [0.5, 1.0],     # 辺上
        [1.5, 1.0],     # 辺上
        [1.0, 0.5],     # 辺上
        [0.1, 0.05],    # 辺近傍
        [1.9, 0.05],    # 辺近傍
        [0.99, 1.99],   # 辺近傍
    ])
    
    results_boundary = []
    for point in query_points_boundary:
        result = interp_boundary(point)
        results_boundary.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    test_cases_boundary.append({
        "name": "Boundary_Edge_Points",
        "description": "Points on triangle edges and near boundaries",
        "input_points": points_boundary.tolist(),
        "input_values": values_boundary.tolist(),
        "queries": results_boundary
    })
    
    return test_cases_boundary

def generate_vector_value_test_data():
    """ベクトル値補間の追加テストデータ生成"""
    print("=== Vector Value Interpolation Test Data Generation ===")
    
    test_cases_vector = []
    
    # 2D→3Dベクトル値
    points_2d_to_3d = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0]
    ])
    values_2d_to_3d = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 1.0]
    ])
    
    interp_2d_to_3d = LinearNDInterpolator(points_2d_to_3d, values_2d_to_3d)
    
    query_points_2d_to_3d = np.array([
        [0.5, 0.5],     # 中心点
        [0.25, 0.25],   # 内部点
        [0.75, 0.25],   # 内部点
        [0.25, 0.75],   # 内部点
        [2.0, 2.0]      # 凸包外
    ])
    
    results_2d_to_3d = []
    for point in query_points_2d_to_3d:
        result = interp_2d_to_3d(point)
        results_2d_to_3d.append({
            "query_point": point.tolist(),
            "result": result.flatten().tolist() if not np.any(np.isnan(result)) else None
        })
    
    test_cases_vector.append({
        "name": "2D_to_3D_Vector",
        "description": "2D points to 3D vector values interpolation",
        "input_points": points_2d_to_3d.tolist(),
        "input_values": values_2d_to_3d.tolist(),
        "queries": results_2d_to_3d
    })
    
    # 3D→2Dベクトル値（逆変換的）
    points_3d_to_2d = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.5, 0.5, 0.5]
    ])
    values_3d_to_2d = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 0.5],
        [0.5, 0.25]
    ])
    
    interp_3d_to_2d = LinearNDInterpolator(points_3d_to_2d, values_3d_to_2d)
    
    query_points_3d_to_2d = np.array([
        [0.25, 0.25, 0.25], # 内部点
        [0.1, 0.2, 0.3],    # 内部点
        [2.0, 2.0, 2.0]     # 凸包外
    ])
    
    results_3d_to_2d = []
    for point in query_points_3d_to_2d:
        result = interp_3d_to_2d(point)
        results_3d_to_2d.append({
            "query_point": point.tolist(),
            "result": result.flatten().tolist() if not np.any(np.isnan(result)) else None
        })
    
    test_cases_vector.append({
        "name": "3D_to_2D_Vector",
        "description": "3D points to 2D vector values interpolation",
        "input_points": points_3d_to_2d.tolist(),
        "input_values": values_3d_to_2d.tolist(),
        "queries": results_3d_to_2d
    })
    
    return test_cases_vector

def generate_edge_case_test_data():
    """エッジケースのテストデータ生成"""
    print("=== Edge Case Test Data Generation ===")
    
    edge_cases = []
    
    # エッジケース1: 共線点（2D）
    points_collinear = np.array([
        [0.0, 0.0],
        [1.0, 1.0],
        [2.0, 2.0],
        [0.0, 1.0]  # 共線でない点を追加
    ])
    values_collinear = np.array([0.0, 2.0, 4.0, 1.0])
    
    interp_collinear = LinearNDInterpolator(points_collinear, values_collinear)
    
    query_collinear = np.array([
        [0.5, 0.5],
        [1.5, 1.5],
        [0.5, 0.0]
    ])
    
    results_collinear = []
    for point in query_collinear:
        result = interp_collinear(point)
        results_collinear.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    edge_cases.append({
        "name": "2D_Collinear_Points",
        "description": "2D interpolation with some collinear points",
        "input_points": points_collinear.tolist(),
        "input_values": values_collinear.tolist(),
        "queries": results_collinear
    })
    
    # エッジケース2: 重複点
    points_duplicate = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.0, 0.0]  # 重複点
    ])
    values_duplicate = np.array([0.0, 1.0, 1.0, 0.0])
    
    try:
        interp_duplicate = LinearNDInterpolator(points_duplicate, values_duplicate)
        
        query_duplicate = np.array([
            [0.0, 0.0],
            [0.5, 0.5]
        ])
        
        results_duplicate = []
        for point in query_duplicate:
            result = interp_duplicate(point)
            results_duplicate.append({
                "query_point": point.tolist(),
                "result": float(result) if not np.isnan(result) else None
            })
        
        edge_cases.append({
            "name": "2D_Duplicate_Points",
            "description": "2D interpolation with duplicate points",
            "input_points": points_duplicate.tolist(),
            "input_values": values_duplicate.tolist(),
            "queries": results_duplicate
        })
    except Exception as e:
        print(f"Warning: Duplicate points test failed: {e}")
    
    return edge_cases

def generate_performance_test_data():
    """パフォーマンステスト用のデータ生成"""
    print("=== Performance Test Data Generation ===")
    
    performance_cases = []
    
    # 2D大量ポイントテスト
    np.random.seed(123)
    n_points_2d = 100
    points_large_2d = np.random.rand(n_points_2d, 2) * 10.0  # 0-10の範囲
    values_large_2d = points_large_2d[:, 0]**2 + points_large_2d[:, 1]**2  # x^2 + y^2
    
    interp_large_2d = LinearNDInterpolator(points_large_2d, values_large_2d)
    
    # クエリポイント（少数）
    query_points_large_2d = np.random.rand(10, 2) * 10.0
    
    results_large_2d = []
    for point in query_points_large_2d:
        result = interp_large_2d(point)
        results_large_2d.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    performance_cases.append({
        "name": "2D_Large_Dataset",
        "description": f"Performance test with {n_points_2d} 2D points, f(x,y) = x^2 + y^2",
        "input_points": points_large_2d.tolist(),
        "input_values": values_large_2d.tolist(),
        "queries": results_large_2d
    })
    
    # 3D中規模ポイントテスト
    np.random.seed(456)
    n_points_3d = 50
    points_large_3d = np.random.rand(n_points_3d, 3) * 5.0  # 0-5の範囲
    values_large_3d = np.sqrt(np.sum(points_large_3d**2, axis=1))  # ユークリッド距離
    
    interp_large_3d = LinearNDInterpolator(points_large_3d, values_large_3d)
    
    query_points_large_3d = np.random.rand(8, 3) * 5.0
    
    results_large_3d = []
    for point in query_points_large_3d:
        result = interp_large_3d(point)
        results_large_3d.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    performance_cases.append({
        "name": "3D_Medium_Dataset",
        "description": f"Performance test with {n_points_3d} 3D points, f(x,y,z) = sqrt(x^2+y^2+z^2)",
        "input_points": points_large_3d.tolist(),
        "input_values": values_large_3d.tolist(),
        "queries": results_large_3d
    })
    
    # 高次元小規模テスト（計算複雑度確認）
    np.random.seed(789)
    n_points_high = 20
    n_dims_high = 7
    points_high_dim = np.random.rand(n_points_high, n_dims_high) * 2.0 - 1.0  # -1から1
    values_high_dim = np.sum(points_high_dim * np.arange(1, n_dims_high + 1), axis=1)  # 重み付き合計
    
    interp_high_dim = LinearNDInterpolator(points_high_dim, values_high_dim)
    
    query_points_high_dim = np.random.rand(5, n_dims_high) * 2.0 - 1.0
    
    results_high_dim = []
    for point in query_points_high_dim:
        result = interp_high_dim(point)
        results_high_dim.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    performance_cases.append({
        "name": f"{n_dims_high}D_High_Dimension",
        "description": f"Performance test with {n_points_high} {n_dims_high}D points",
        "input_points": points_high_dim.tolist(),
        "input_values": values_high_dim.tolist(),
        "queries": results_high_dim
    })
    
    return performance_cases

def save_test_data(test_cases, filename):
    """テストデータをJSONファイルに保存"""
    output_data = {
        "description": "SciPy LinearNDInterpolator reference test data for C++ implementation validation",
        "scipy_version": "1.x",
        "numpy_version": "1.x",
        "test_cases": test_cases
    }
    
    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, indent=2, ensure_ascii=False)
    
    print(f"Test data saved to: {filename}")

def main():
    """メイン処理"""
    print("SciPy LinearNDInterpolator テストデータ生成開始\n")
    
    # 出力ディレクトリの作成
    memo_dir = "out"
    if not os.path.exists(memo_dir):
        os.makedirs(memo_dir)
    
    # 全テストケースの生成
    all_test_cases = []
    
    # 基本次元テストデータ
    print("基本次元テストデータを生成中...")
    test_cases_2d = generate_2d_test_data()
    all_test_cases.extend(test_cases_2d)
    
    test_cases_3d = generate_3d_test_data()
    all_test_cases.extend(test_cases_3d)
    
    test_cases_4d = generate_4d_test_data()
    all_test_cases.extend(test_cases_4d)
    
    # 高次元テストデータ（5次元以上でC++の"Qx"オプションテスト）
    print("高次元テストデータを生成中...")
    high_dim_cases = generate_5d_and_higher_test_data()
    all_test_cases.extend(high_dim_cases)
    
    # 数値精度テストデータ
    print("数値精度テストデータを生成中...")
    precision_cases = generate_numerical_precision_test_data()
    all_test_cases.extend(precision_cases)
    
    # 境界条件テストデータ
    print("境界条件テストデータを生成中...")
    boundary_cases = generate_boundary_test_data()
    all_test_cases.extend(boundary_cases)
    
    # ベクトル値補間テストデータ
    print("ベクトル値補間テストデータを生成中...")
    vector_cases = generate_vector_value_test_data()
    all_test_cases.extend(vector_cases)
    
    # エッジケース
    print("エッジケーステストデータを生成中...")
    edge_cases = generate_edge_case_test_data()
    all_test_cases.extend(edge_cases)
    
    # パフォーマンステスト
    print("パフォーマンステストデータを生成中...")
    performance_cases = generate_performance_test_data()
    all_test_cases.extend(performance_cases)
    
    # JSONファイルに保存
    json_filename = os.path.join(memo_dir, "scipy_reference_data.json")
    save_test_data(all_test_cases, json_filename)
        
    print("\n=== 生成完了 ===")
    print(f"総テストケース数: {len(all_test_cases)}")
    print(f"JSONファイル: {json_filename}")
    
    # 詳細な統計情報
    total_queries = sum(len(tc['queries']) for tc in all_test_cases)
    print(f"総クエリ数: {total_queries}")
    
    dimensions = set()
    vector_dims = set()
    for tc in all_test_cases:
        if tc['input_points']:
            dimensions.add(len(tc['input_points'][0]))
        if tc['input_values'] and isinstance(tc['input_values'][0], list):
            vector_dims.add(len(tc['input_values'][0]))
    
    print(f"テスト空間次元: {sorted(dimensions)}")
    print(f"テスト値次元: {sorted(vector_dims)}")
    
    # テストカテゴリ別統計
    categories = {}
    for tc in all_test_cases:
        category = tc['name'].split('_')[0] + 'D' if tc['name'][0].isdigit() else tc['name'].split('_')[0]
        categories[category] = categories.get(category, 0) + 1
    
    print("\nテストカテゴリ別統計:")
    for category, count in sorted(categories.items()):
        print(f"  {category}: {count}ケース")
        
    print("\n拡充されたテスト項目:")
    print("  ✓ 高次元テスト（5D-6D、C++の'Qx'オプション検証）")
    print("  ✓ 数値精度テスト（1e-10許容範囲検証）")
    print("  ✓ 境界条件テスト（辺上の点）")
    print("  ✓ ベクトル値補間テスト（多次元出力値）")
    print("  ✓ 大規模データパフォーマンステスト")
    print("  ✓ 浮動小数点精度テスト")

if __name__ == "__main__":
    try:
        main()
    except ImportError as e:
        print(f"Error: Required package not found: {e}")
        print("Please install required packages: pip install numpy scipy")
    except Exception as e:
        print(f"Error: {e}")