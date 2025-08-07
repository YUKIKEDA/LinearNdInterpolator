# LinearNdInterpolator

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](https://github.com/your-repo/LinearNdInterpolator)
[![SciPy Compatible](https://img.shields.io/badge/SciPy-compatible-blue)](https://scipy.org/)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue)](https://en.cppreference.com/w/cpp/17)

SciPyの`LinearNDInterpolator`と**完全な数値互換性**を実現する高性能C++ライブラリです。ドロネー三角分割を使用したN次元線形補間を提供し、科学計算アプリケーションに適した精密な補間結果を実現します。

## ✨ 特徴

- **🎯 SciPy互換性**: SciPyのLinearNDInterpolatorと100%の数値精度一致
- **📐 N次元サポート**: 2次元から4次元以上の補間をシームレスに処理
- **⚡ ドロネー三角分割**: Qhullライブラリを使用した堅牢な三角分割
- **🔢 重心座標補間**: 重心座標を使用した線形補間
- **🛡️ 堅牢なエラーハンドリング**: 包括的な入力検証とエッジケース管理
- **🔍 凸包検出**: 凸包外の点にはNaNを返す（SciPy互換）

## 🚀 クイックスタート

### 基本的な使用方法

```cpp
#include "LinearNdInterpolator.h"
#include <vector>

// 2次元の点とその値を定義
std::vector<std::vector<double>> points = {
    {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {0.5, 0.5}
};
std::vector<double> values = {0.0, 1.0, 1.0, 1.0};

// 補間器を作成
LinearNdInterpolator interp(points, values);

// 新しい点で補間
double result = interp.interpolate({0.25, 0.25});
std::cout << "補間値: " << result << std::endl;
```

### 複数点の補間

```cpp
// 複数の点を一度に補間
std::vector<std::vector<double>> query_points = {
    {0.1, 0.1}, {0.3, 0.4}, {0.7, 0.2}
};

auto results = interp.interpolate(query_points);
for (size_t i = 0; i < results.size(); ++i) {
    std::cout << "点 " << i << ": " << results[i][0] << std::endl;
}
```

## 🏗️ プロジェクトのビルド

### 前提条件

- **コンパイラ**: C++17対応コンパイラ（Visual Studio 2019+、GCC 7+、Clang 5+）
- **CMake**: バージョン3.15以上
- **依存関係**: Qhull（含まれています）、GoogleTest（自動ダウンロード）

### ビルド手順

```bash
# ビルドファイルを生成
cmake -B build -G "Visual Studio 17 2022"

# プロジェクトをビルド
cmake --build build

# テストを実行
"build\bin\Debug\test_scipy_reference.exe"
```

### 統合

```cmake
# CMakeLists.txtに追加
find_package(LinearNdInterpolator REQUIRED)
target_link_libraries(your_target LinearNdInterpolator::LinearNdInterpolator)
```

## 📊 パフォーマンスと互換性

### SciPy互換性テスト

✅ **全7つのテストスイートが通過:**
- 2次元基本三角形
- 2次元複雑ランダム点
- 3次元基本四面体
- 4次元基本単体
- 2次元共線点
- 2次元重複点
- 2次元大規模データセット（100点）

```bash
# 互換性テストを実行
"build\bin\Debug\test_scipy_reference.exe"
```

### 現在のパフォーマンス状況

| 機能 | 状況 | 備考 |
|------|------|------|
| **数値精度** | ✅ 完璧 | 100% SciPy互換性 |
| **基本機能** | ✅ 完了 | 全コア機能が動作 |
| **エラーハンドリング** | ✅ 良好 | 包括的な検証 |
| **変換行列** | ✅ 実装済み | SciPyスタイル最適化が有効 |
| **歩行アルゴリズム** | ✅ 完了 | SciPy互換の効率的な単体探索が実装済み |
| **パフォーマンス** | 🟡 部分的に最適化 | 変換行列が実装済み |
| **メモリ使用量** | 🟡 許容範囲 | 最適化の余地あり |

## 🛠️ APIリファレンス

### コンストラクタ

```cpp
LinearNdInterpolator(
    const std::vector<std::vector<double>>& points, 
    const std::vector<double>& values
);

LinearNdInterpolator(
    const std::vector<std::vector<double>>& points, 
    const std::vector<std::vector<double>>& values  // 多次元値
);
```

### メソッド

```cpp
// 単一点補間
double interpolate(const std::vector<double>& point) const;

// 複数点補間
std::vector<std::vector<double>> interpolate(
    const std::vector<std::vector<double>>& points
) const;
```

### エラーハンドリング

ライブラリは以下の場合に`std::invalid_argument`を投げます：
- 空の入力データ
- 点と値の数の不一致
- 次元の不整合
- 2次元未満の点
- NaNまたは無限大の入力値

## 📈 開発状況

### ✅ 完了済み（フェーズ1-2）

- **コアアルゴリズム**: N次元線形補間
- **SciPy互換性**: 完全な数値互換性を達成
- **変換行列**: SciPyスタイルの重心座標最適化を実装
- **歩行アルゴリズム**: SciPy互換の効率的な単体探索アルゴリズム
- **堅牢な実装**: 包括的なエラーハンドリングと検証
- **テストスイート**: 全17のSciPy参照テストが通過
- **クロスプラットフォームビルド**: CMakeベースのビルドシステム

### 🚧 現在の制限とロードマップ

#### 高優先度改善

1. **✅ 歩行アルゴリズム実装**（完了）
   - **達成**: SciPy互換の効率的な単体探索アルゴリズム
   - **実装**: findSimplex()、findSimplexDirected()、findSimplexBruteforce()
   - **状況**: 隣接関係データも完全実装済み

2. **✅ 変換行列実装**（完了）
   - **達成**: SciPyスタイルの事前計算変換行列
   - **影響**: 重心座標計算がO(d³) → O(d²)
   - **状況**: 全17のSciPy互換性テストが通過

#### 中優先度改善

3. **🟡 強化された単体テスト**（目標: 1週間）
   - 包括的なエッジケーステストの追加
   - パフォーマンスベンチマークスイート
   - メモリ使用量検証

4. **🟡 スレッド安全性**（目標: 3日）
   - スレッドセーフ操作の追加
   - 並行クエリ処理

5. **🟡 メモリ最適化**（目標: 1週間）
   - 頂点マッピングアルゴリズムの最適化
   - 大規模データセットのメモリフットプリント削減

#### 将来の拡張

- **増分点追加**: 動的ポイント挿入
- **カスタムフィル値**: 凸包外の点に対するNaN以外の選択肢
- **並列処理**: 大規模補間のためのOpenMPサポート
- **高度な補間**: 3次および高次メソッド

## 🔬 技術詳細

### アーキテクチャ

```
LinearNdInterpolator
├── Delaunay (Qhullベース三角分割)
│   ├── findSimplex()
│   ├── calculateBarycentricCoordinates()
│   └── getSimplices()
└── interpolate() (線形補間)
```

### アルゴリズムフロー

1. **三角分割**: Qhullを使用してドロネー三角分割を作成
2. **単体位置**: クエリ点を含む単体を検索
3. **重心座標計算**: 重心座標を計算
4. **線形補間**: 座標で頂点値を重み付け

### メモリレイアウト

- **点**: `std::vector<std::vector<double>>`として格納
- **値**: 多次元値サポート
- **三角分割**: RAII原則でQhullによって管理

## 📋 テスト

### テストカテゴリ

```bash
# SciPy互換性テスト（主要検証）
"build\bin\Debug\test_scipy_reference.exe"

# 特定のテストケース
"build\bin\Debug\test_scipy_reference.exe" --gtest_filter="SciPyReferenceTest.2D_Basic_Triangle"
"build\bin\Debug\test_scipy_reference.exe" --gtest_filter="SciPyReferenceTest.3D_Basic_Tetrahedron"
```

### テストカバレッジ

- ✅ **2D-4D補間**: 包括的な次元テスト
- ✅ **エッジケース**: 共線点、重複点、大規模データセット
- ✅ **エラー条件**: 無効な入力、凸包外の点
- ✅ **数値精度**: 正確なSciPy互換性検証

## 🤝 貢献

### 開発優先事項

1. **パフォーマンス最適化**: さらなる最適化とメモリ効率化
2. **テスト強化**: 独立した単体テストスイートの開発
3. **ドキュメント**: 使用例とベストプラクティス
4. **ベンチマーク**: パフォーマンス比較ツール

### 始め方

```bash
# クローンしてビルド
git clone https://github.com/your-repo/LinearNdInterpolator
cd LinearNdInterpolator
cmake -B build
cmake --build build

# テストを実行して検証
"build\bin\Debug\test_scipy_reference.exe"
```

## 📚 例

### 3次元表面補間

```cpp
// 3次元表面を定義: f(x,y,z) = x + y + z
std::vector<std::vector<double>> points = {
    {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}
};
std::vector<double> values = {0, 1, 1, 1, 3};

LinearNdInterpolator interp(points, values);

// 内部点で補間
double result = interp.interpolate({0.25, 0.25, 0.25});
// 期待値: 0.75 (線形結合)
```

### 大規模データセット処理

```cpp
// 大規模点群を処理
std::vector<std::vector<double>> large_points = loadPointCloud("data.csv");
std::vector<double> large_values = loadValues("values.csv");

LinearNdInterpolator interp(large_points, large_values);

// バッチ補間
std::vector<std::vector<double>> query_grid = generateGrid(100, 100);
auto results = interp.interpolate(query_grid);
```

## 📄 ライセンス

このプロジェクトはMITライセンスの下でライセンスされています - 詳細は[LICENSE](LICENSE)ファイルを参照してください。

## 🙏 謝辞

- **SciPyチーム**: 参照実装とテストケース
- **Qhull**: 堅牢な計算幾何学ライブラリ
- **GoogleTest**: テストフレームワーク

## 📞 サポート

- **Issues**: [GitHub Issues](https://github.com/your-repo/LinearNdInterpolator/issues)
- **ドキュメント**: 詳細な技術ドキュメントは`.memo/`ディレクトリを参照
- **パフォーマンス**: 最適化ロードマップは`10_interpolate_method_scipy_comparison.md`を参照

---

**状況**: ✅ **SciPy互換パフォーマンス最適化で本番準備完了**  
**最新の成果**: 🎯 **歩行アルゴリズム完全実装 - 効率的な単体探索**  
**次の焦点**: 🚀 **さらなるパフォーマンス最適化とメモリ効率化**