

if(CSPARSE_FOUND)
  add_subdirectory(csparse)
  add_subdirectory(slam2d_linear)
endif()

if(CHOLMOD_FOUND)
  add_subdirectory(cholmod)
endif(CHOLMOD_FOUND)

