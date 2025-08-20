add_executable(2D_Euler main.cpp)
target_include_directories(2D_Euler PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(2D_Euler Python3::Python Python3::Module)