#----------------------------------------------
set(LIB_NAME ${LC_PROJECT_NAME}WeightFunctions)
set(LIB_FOLDER engine/weight_functions)

# Library header files.
set(HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/weight_function.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/weight_factory.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cubic_spline.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cubic_spline.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/quartic_spline.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/quartic_spline.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gaussian_rbf.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gaussian_rbf.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/multiquadric_rbf.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/multiquadric_rbf.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/polyharmonic_rbf.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/polyharmonic_rbf.tpp 
)

#-------- Build library --------
add_library(${LIB_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${LIB_NAME} ALIAS ${LIB_NAME})
target_link_libraries(${LIB_NAME} INTERFACE IMP::IMP)

target_include_directories(${LIB_NAME} INTERFACE
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>/${PROJECT_NAME}
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
set_target_properties(${LIB_NAME}
    PROPERTIES
        PUBLIC_HEADER "${HEADERS};"
)
#--------------------------------------

install(TARGETS
        ${LIB_NAME}
    EXPORT
        ${LIB_NAME}Targets
    RUNTIME
        DESTINATION ${INSTALL_BINDIR}
        COMPONENT ${PROJECT_NAME}_Runtime
    LIBRARY
        DESTINATION ${INSTALL_LIBDIR}
        COMPONENT ${PROJECT_NAME}_Runtime
    ARCHIVE
        DESTINATION ${INSTALL_LIBDIR}
        COMPONENT ${PROJECT_NAME}_Development
    PUBLIC_HEADER
        DESTINATION ${INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${LIB_FOLDER}
        COMPONENT ${PROJECT_NAME}_Development
)

install(EXPORT ${LIB_NAME}Targets
    NAMESPACE "${PROJECT_NAME}::"
    DESTINATION ${INSTALL_CMAKEDIR}
    COMPONENT ${PROJECT_NAME}_Development
)