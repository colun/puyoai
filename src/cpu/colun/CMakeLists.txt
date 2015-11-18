cmake_minimum_required(VERSION 2.8)

cpu_setup("colun")

function(colun_add_executable exe)
  cpu_add_executable(${exe} ${ARGN})
  cpu_target_link_common_libraries(${exe})
endfunction()

colun_add_executable(mc_simple mc_simple.cc)
colun_add_executable(mc_beam mc_beam.cc)
