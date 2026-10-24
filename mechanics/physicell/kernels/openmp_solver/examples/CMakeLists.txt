set(ALL_EXAMPLES benchmark)

foreach(example ${ALL_EXAMPLES})
  add_executable(mechanics.physicell.kernels.openmp_solver.${example}
                 ${example}.cpp)
  target_link_libraries(
    mechanics.physicell.kernels.openmp_solver.${example}
    PRIVATE mechanics.physicell.kernels.openmp_solver_internal_iface
            mechanics.physicell)
endforeach()
