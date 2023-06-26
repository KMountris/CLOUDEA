#--------------------------------------------------------------
set(LIB_NAME ${PROJECT_NAME})
string(TOLOWER ${PROJECT_NAME} LC_PROJECT_NAME)

# Find IMP
find_package(IMP REQUIRED)
include_directories(${IMP_INCLUDE_DIRS})

# Find Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find Armadillo
find_package(Armadillo REQUIRED NO_MODULE)

# Find Boost
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

# Find TinyXML2
find_package(TinyXML2 REQUIRED)
include_directories(${TinyXML2_INCLUDE_DIRS})

# Find CGAL if it is available.
find_package(CGAL QUIET COMPONENTS Core)
include_directories(${CGAL_INCLUDE_DIRS})
# -------------------------

# The headers collection of the library.
set(HEADERS
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/CLOUDEA
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Approximants
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Conditions
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Elements
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Exporters
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Grid
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Integration
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Materials
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Mesh
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/MeshIO
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Models
    # ${CMAKE_SOURCE_DIR}/src/CLOUDEA/OptionsConfiguration
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Quadrature
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Sets
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Solvers
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/SupportDomain
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Utilities
    ${CMAKE_SOURCE_DIR}/src/CLOUDEA/Vectors
)

# Subdirectories of the library.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/approximants)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/conditions)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/elements)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/exporters)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/grid)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/integration)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/materials)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/mesh)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/mesh_io)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/models)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/monomial)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/options_configuration)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/quadrature)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/sets)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/solvers)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/support_domain)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/utilities)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/vectors)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA/engine/weight_functions)

set(${PROJECT_NAME}_LIBS
    ${LC_PROJECT_NAME}Approximants
    ${LC_PROJECT_NAME}Conditions
    ${LC_PROJECT_NAME}Elements
    ${LC_PROJECT_NAME}Exporters
    ${LC_PROJECT_NAME}Grid
    ${LC_PROJECT_NAME}Integration
    ${LC_PROJECT_NAME}Materials
    ${LC_PROJECT_NAME}Mesh
    ${LC_PROJECT_NAME}MeshIO
    ${LC_PROJECT_NAME}Models
    ${LC_PROJECT_NAME}Monomial
    # ${LC_PROJECT_NAME}OptionsConfiguration
    ${LC_PROJECT_NAME}Quadrature
    ${LC_PROJECT_NAME}Sets
    ${LC_PROJECT_NAME}Solvers
    ${LC_PROJECT_NAME}SupportDomain
    ${LC_PROJECT_NAME}Utilities
    ${LC_PROJECT_NAME}Vectors
    ${LC_PROJECT_NAME}WeightFunctions
)

#-------- Build library --------
add_library(${LIB_NAME} "")
target_sources(${LIB_NAME} PRIVATE ${HEADERS})
target_link_libraries(${LIB_NAME} PUBLIC ${${PROJECT_NAME}_LIBS}
    IMP::IMP
    Boost::filesystem
    Boost::program_options
    ${TinyXML2_LIBRARIES}
)
target_include_directories(${LIB_NAME} PUBLIC
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>/${PROJECT_NAME}
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    ${Boost_INCLUDE_DIR}
    ${TinyXML2_INCLUDE_DIRS}
    ${CGAL_INCLUDE_DIRS}
)
set_target_properties(${LIB_NAME}
    PROPERTIES
        LINKER_LANGUAGE CXX
        POSITION_INDEPENDENT_CODE 1
        DEBUG_POSTFIX "_d"
        PUBLIC_HEADER "${HEADERS}"
        $<$<BUILD_SHARED_LIBS>:CXX_VISIBILITY_PRESET hidden>
        $<$<BUILD_SHARED_LIBS>:VISIBILITY_INLINES_HIDDEN 1>
        $<$<BUILD_SHARED_LIBS>:SOVERSION ${PROJECT_VERSION}>
        $<$<BUILD_SHARED_LIBS>:OUTPUT_NAME ${PROJECT_NAME}>
        $<$<BUILD_SHARED_LIBS>:MACOSX_RPATH ON>
        $<$<NOT:$<BUILD_SHARED_LIBS>>:ARCHIVE_OUTPUT_NAME ${PROJECT_NAME}>
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
            DESTINATION ${INSTALL_INCLUDEDIR}/${PROJECT_NAME}
            COMPONENT ${PROJECT_NAME}_Development
)

install(EXPORT ${LIB_NAME}Targets
    NAMESPACE "${PROJECT_NAME}::"
    DESTINATION ${INSTALL_CMAKEDIR}
    COMPONENT ${PROJECT_NAME}_Development
)

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION
        ${INSTALL_CMAKEDIR}
    PATH_VARS
        INSTALL_INCLUDEDIR
        INSTALL_LIBDIR
)

write_basic_package_version_file(${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION
        ${PROJECT_VERSION}
    COMPATIBILITY
        SameMajorVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION
        ${INSTALL_CMAKEDIR}
)
