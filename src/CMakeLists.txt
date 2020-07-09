set(HEADER_FILES npm.h patch.h population.h visitors.h cmd_line.h individual.h rndutils.hpp)
add_executable(npm main.cpp npm.cpp patch.cpp population.cpp individual.cpp)
target_include_directories(npm PRIVATE "./")
install(TARGETS npm DESTINATION bin)
