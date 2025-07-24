cmake_minimum_required(VERSION 3.10)
project(XFoil)

# (1) C++の標準指定（必要に応じて変更）
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# (2) 最適化フラグ -O2 の付与
#     Releaseビルドの場合は CMAKE_BUILD_TYPE=Release でも自動的に -O2 が付与されますが、
#     確実に指定したい場合は以下のように書くこともできます。
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# (4) 実行ファイルにするソースをまとめて指定
add_executable(testExe
    ./model/boundary_layer.cpp
    ./model/math_util.cpp
    ./model/spline.cpp
    ./XFoil.cpp
    ./main.cpp
)
target_link_libraries(testExe PRIVATE Eigen3::Eigen)

file(COPY sample/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/sample)

# (5) 出力ファイル名を "test.exe" にしたい場合の設定 (Linuxでも .exe にできます)
#     以下を追加すると、ビルド後に out/test.exe が生成されます。
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/out)
set_target_properties(testExe PROPERTIES OUTPUT_NAME "test")
set_target_properties(testExe PROPERTIES SUFFIX ".exe")
