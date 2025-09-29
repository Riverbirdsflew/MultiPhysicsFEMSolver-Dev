set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosThermDiffStructPhaseApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## ThermDiffStructPhase Core sources
file(GLOB_RECURSE KRATOS_THERM_DIFF_STRUCT_PHASE_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/therm_diff_struct_phase_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/therm_diff_struct_phase_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## ThermDiffStructPhase testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_THERM_DIFF_STRUCT_PHASE_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## ThermDiffStructPhase python interface sources
file(GLOB_RECURSE KRATOS_THERM_DIFF_STRUCT_PHASE_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosThermDiffStructPhaseCore SHARED ${KRATOS_THERM_DIFF_STRUCT_PHASE_APPLICATION_CORE} ${KRATOS_THERM_DIFF_STRUCT_PHASE_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosThermDiffStructPhaseCore PUBLIC KratosCore)
set_target_properties(KratosThermDiffStructPhaseCore PROPERTIES COMPILE_DEFINITIONS "THERM_DIFF_STRUCT_PHASE_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosThermDiffStructPhaseApplication MODULE THIN_LTO ${KRATOS_THERM_DIFF_STRUCT_PHASE_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosThermDiffStructPhaseApplication PRIVATE KratosThermDiffStructPhaseCore)
set_target_properties(KratosThermDiffStructPhaseApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosThermDiffStructPhaseApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosThermDiffStructPhaseApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ThermDiffStructPhaseApplication.py KratosMultiphysics/ThermDiffStructPhaseApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/ThermDiffStructPhaseApplication/python_registry_lists.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosThermDiffStructPhaseCore DESTINATION libs )
install(TARGETS KratosThermDiffStructPhaseApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosThermDiffStructPhaseCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosThermDiffStructPhaseApplication" PARENT_SCOPE)