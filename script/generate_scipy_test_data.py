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
    
    test_case_4d = {
        "name": "4D_Basic_Simplex",
        "description": "Basic 4D simplex interpolation",
        "input_points": points_4d.tolist(),
        "input_values": values_4d.tolist(),
        "queries": results_4d
    }
    
    return [test_case_4d]

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
    
    # 大量の2Dポイント
    np.random.seed(123)
    n_points = 100
    points_large = np.random.rand(n_points, 2) * 10.0  # 0-10の範囲
    values_large = points_large[:, 0]**2 + points_large[:, 1]**2  # x^2 + y^2
    
    interp_large = LinearNDInterpolator(points_large, values_large)
    
    # クエリポイント（少数）
    query_points_large = np.random.rand(10, 2) * 10.0
    
    results_large = []
    for point in query_points_large:
        result = interp_large(point)
        results_large.append({
            "query_point": point.tolist(),
            "result": float(result) if not np.isnan(result) else None
        })
    
    performance_case = {
        "name": "2D_Large_Dataset",
        "description": f"Performance test with {n_points} 2D points, f(x,y) = x^2 + y^2",
        "input_points": points_large.tolist(),
        "input_values": values_large.tolist(),
        "queries": results_large
    }
    
    return [performance_case]

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
    
    # 2Dテストデータ
    test_cases_2d = generate_2d_test_data()
    all_test_cases.extend(test_cases_2d)
    
    # 3Dテストデータ  
    test_cases_3d = generate_3d_test_data()
    all_test_cases.extend(test_cases_3d)
    
    # 4Dテストデータ
    test_cases_4d = generate_4d_test_data()
    all_test_cases.extend(test_cases_4d)
    
    # エッジケース
    edge_cases = generate_edge_case_test_data()
    all_test_cases.extend(edge_cases)
    
    # パフォーマンステスト
    performance_cases = generate_performance_test_data()
    all_test_cases.extend(performance_cases)
    
    # JSONファイルに保存
    json_filename = os.path.join(memo_dir, "scipy_reference_data.json")
    save_test_data(all_test_cases, json_filename)
        
    print("\n=== 生成完了 ===")
    print(f"総テストケース数: {len(all_test_cases)}")
    print(f"JSONファイル: {json_filename}")
    
    # 簡単な統計情報
    total_queries = sum(len(tc['queries']) for tc in all_test_cases)
    print(f"総クエリ数: {total_queries}")
    
    dimensions = set()
    for tc in all_test_cases:
        if tc['input_points']:
            dimensions.add(len(tc['input_points'][0]))
    print(f"テスト次元: {sorted(dimensions)}")

if __name__ == "__main__":
    try:
        main()
    except ImportError as e:
        print(f"Error: Required package not found: {e}")
        print("Please install required packages: pip install numpy scipy")
    except Exception as e:
        print(f"Error: {e}")