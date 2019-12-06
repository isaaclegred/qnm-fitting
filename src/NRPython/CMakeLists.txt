add_library(mylib SHARED varpro.cpp)
include_directories(/usr/local/anaconda3/include/python3.7m/ ../includes
/usr/local/anaconda3/lib/python3.7/site-packages/numpy/core/include)
TARGET_LINK_LIBRARIES(
        mylib
        /usr/local/anaconda3/bin/python3.7m
        libnpymath.a
        )
