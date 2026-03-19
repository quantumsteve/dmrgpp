add_executable (lanczos lanczos.cpp)
target_include_directories (lanczos PRIVATE
${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/Engine ${CMAKE_CURRENT_SOURCE_DIR}/Model)
target_link_libraries (lanczos PUBLIC psimaglite::psimaglite psimaglite)

add_executable (thermal thermal.cpp)
target_include_directories (thermal PRIVATE
${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/Engine)
target_link_libraries (thermal psimaglite::psimaglite)

add_executable (lorentzian lorentzian.cpp)



target_link_libraries (lorentzian
  psimaglite::psimaglite
  )
