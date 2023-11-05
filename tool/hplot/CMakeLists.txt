# include libraries for the hview binary
include_directories(tool/hplot/include)

# add hplot executable
add_executable(hplot
    tool/hplot/src/main.cpp
    tool/hplot/src/plot.cpp
)
