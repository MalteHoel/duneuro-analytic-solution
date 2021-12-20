add_library("duneuroAnalyticSolutionPy" SHARED duneuro-analytic-solution.cc)
target_link_libraries(duneuroAnalyticSolutionPy ${PYTHON_LIBRARIES})
set_target_properties(duneuroAnalyticSolutionPy PROPERTIES PREFIX "")

