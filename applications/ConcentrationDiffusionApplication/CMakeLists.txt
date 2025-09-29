set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosConcentrationDiffusionApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## ConcentrationDiffusion Core sources
file(GLOB_RECURSE KRATOS_CONCENTRATION_DIFFUSION_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/concentration_diffusion_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/concentration_diffusion_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    
)

## ConcentrationDiffusion testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_CONCENTRATION_DIFFUSION_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## ConcentrationDiffusion python interface sources
file(GLOB_RECURSE KRATOS_CONCENTRATION_DIFFUSION_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

###############################################################
## ConcentrationDiffusionApplication core library (C++ parts)
add_library(KratosConcentrationDiffusionCore SHARED ${KRATOS_CONCENTRATION_DIFFUSION_APPLICATION_CORE} ${KRATOS_CONCENTRATION_DIFFUSION_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosConcentrationDiffusionCore PUBLIC KratosCore)
set_target_properties(KratosConcentrationDiffusionCore PROPERTIES COMPILE_DEFINITIONS "CONCENTRATION_DIFFUSION_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosConcentrationDiffusionApplication MODULE THIN_LTO ${KRATOS_CONCENTRATION_DIFFUSION_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosConcentrationDiffusionApplication PRIVATE KratosConcentrationDiffusionCore)
set_target_properties(KratosConcentrationDiffusionApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosConcentrationDiffusionApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosConcentrationDiffusionApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ConcentrationDiffusionApplication.py KratosMultiphysics/ConcentrationDiffusionApplication/__init__.py )
#kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/ConcentrationDiffusionApplication/python_registry_lists.py )

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
install(TARGETS KratosConcentrationDiffusionCore DESTINATION libs )
install(TARGETS KratosConcentrationDiffusionApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosConcentrationDiffusionCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosConcentrationDiffusionApplication" PARENT_SCOPE)
