
# Check if Doxygen is installed in the OS.
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
    # Set input and output files.
    SET(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/CLOUDEA_doxyfile.in)
    SET(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/CLOUDEA_doxyfile)

    # Request to configure the Doxygen documentation file.
    CONFIGURE_FILE(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    
    # Inform user how to build the Doxygen documentation.
    MESSAGE(" ")
    MESSAGE("=============================================================================")
    MESSAGE("==   You have selected to build the doxygen documentation of CLOUDEA.      ==")
    MESSAGE("==   The required dependencies where found in your Operating System.       ==")
    MESSAGE("==                                                                         ==")
    MESSAGE("==   To build the documentation of CLOUDEA execute:                        ==")
    MESSAGE("==                                                                         ==")
    MESSAGE("==   make doc                                                              ==")
    MESSAGE("==                                                                         ==")
    MESSAGE("==   after completing the configuration and generation of the CMake build. ==")
    MESSAGE("=============================================================================")
    MESSAGE(" ")
    
    # Add custom target "doc" for documentation compilation.
    ADD_CUSTOM_TARGET(doc COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                          COMMENT "Generating the CLOUDEA Doxygen documentation..."
                          VERBATIM )
ELSE ()
  MESSAGE("Doxygen and Graphviz' dot need to be installed to generate the CLOUDEA Doxygen documentation.")
ENDIF ()
