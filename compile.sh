mkdir -p build
cd build
#cmake -DCMAKE_PREFIX_PATH=/Users/kristjan/lib/libtorch ..
#cmake -DCMAKE_PREFIX_PATH=/usr/local/opt/libtorch ..
cmake ..
make