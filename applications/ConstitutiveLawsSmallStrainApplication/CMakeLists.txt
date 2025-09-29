set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosConstitutiveLawsSmallStrainApplication ****")

################### PYBIND11
include(pybind11Tools)

# temp dependency until the remaining basic laws are moved to the Core
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralAnalysisApplication)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/StructuralAnalysisApplication )

## ConstitutiveLaws Core sources
file(GLOB_RECURSE KRATOS_CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws_small_strain_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws_small_strain_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strains/plasticity/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
)

## Exclude sources that might yield problems with the unity build
if(CMAKE_UNITY_BUILD MATCHES ON)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/advanced_constitutive_law_small_strain_utilities.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
endif(CMAKE_UNITY_BUILD MATCHES ON)

## ConstitutiveLaws testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    #file(GLOB_RECURSE KRATOS_CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## ConstitutiveLaws python interface sources
file(GLOB_RECURSE KRATOS_CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosConstitutiveLawsSmallStrainCore SHARED ${KRATOS_CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION_CORE} ${KRATOS_CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosConstitutiveLawsSmallStrainCore PUBLIC KratosCore KratosStructuralAnalysisCore)
set_target_properties(KratosConstitutiveLawsSmallStrainCore PROPERTIES COMPILE_DEFINITIONS "CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosConstitutiveLawsSmallStrainApplication MODULE THIN_LTO ${KRATOS_CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosConstitutiveLawsSmallStrainApplication PRIVATE KratosConstitutiveLawsSmallStrainCore)
set_target_properties(KratosConstitutiveLawsSmallStrainApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosConstitutiveLawsSmallStrainCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosConstitutiveLawsSmallStrainApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosConstitutiveLawsSmallStrainApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosConstitutiveLawsSmallStrainApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ConstitutiveLawsSmallStrainApplication.py KratosMultiphysics/ConstitutiveLawsSmallStrainApplication/__init__.py )

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
install(TARGETS KratosConstitutiveLawsSmallStrainCore DESTINATION libs )
install(TARGETS KratosConstitutiveLawsSmallStrainApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosConstitutiveLawsSmallStrainCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosConstitutiveLawsSmallStrainApplication" PARENT_SCOPE)
