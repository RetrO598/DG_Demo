# add_executable(1D_Euler main.cpp)
# target_link_libraries(1D_Euler ${PYTHON3.12_LIB})

add_executable(1D_Euler main.cpp)
target_include_directories(1D_Euler PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(1D_Euler Python3::Python Python3::Module)