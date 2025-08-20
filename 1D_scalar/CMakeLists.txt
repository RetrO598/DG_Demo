add_executable(1D_scalar main.cpp)
target_include_directories(1D_scalar PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(1D_scalar Python3::Python Python3::Module)