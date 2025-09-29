set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosThermalConductionApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## ThermalConduction Core sources
file(GLOB_RECURSE KRATOS_THERMAL_CONDUCTION_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/thermal_conduction_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/thermal_conduction_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## ThermalConduction testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_THERMAL_CONDUCTION_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## ThermalConduction python interface sources
file(GLOB_RECURSE KRATOS_THERMAL_CONDUCTION_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosThermalConductionCore SHARED ${KRATOS_THERMAL_CONDUCTION_APPLICATION_CORE} ${KRATOS_THERMAL_CONDUCTION_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosThermalConductionCore PUBLIC KratosCore)
set_target_properties(KratosThermalConductionCore PROPERTIES COMPILE_DEFINITIONS "THERMAL_CONDUCTION_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosThermalConductionApplication MODULE THIN_LTO ${KRATOS_THERMAL_CONDUCTION_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosThermalConductionApplication PRIVATE KratosThermalConductionCore)
set_target_properties(KratosThermalConductionApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosThermalConductionApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosThermalConductionApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ThermalConductionApplication.py KratosMultiphysics/ThermalConductionApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/ThermalConductionApplication/python_registry_lists.py )

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
install(TARGETS KratosThermalConductionCore DESTINATION libs )
install(TARGETS KratosThermalConductionApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosThermalConductionCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosThermalConductionApplication" PARENT_SCOPE)