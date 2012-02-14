# This project is designed to be built outside the Insight source tree.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(LandmarkIds)

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

IF( NOT ITK_USE_REVIEW )
  MESSAGE(FATAL_ERROR "ITK must be build with ITK_USE_REVIEW turned ON")
ENDIF( NOT ITK_USE_REVIEW )

ADD_EXECUTABLE(LandmarkIds LandmarkIds.cxx)

TARGET_LINK_LIBRARIES(LandmarkIds ITKQuadEdgeMesh ITKBasicFilters ITKCommon ITKIO ITKNumerics ITKAlgorithms)

INCLUDE(CTest)

SET(TEST_DATA_ROOT ${PROJECT_SOURCE_DIR}/Data)
SET(TEMP ${PROJECT_BINARY_DIR}/Testing/Temporary)

ADD_TEST( TEST_LandmarkIds ${CMAKE_CURRENT_BINARY_DIR}/LandmarkIds ${TEST_DATA_ROOT}/fixedMeshWithScalarsIC1.vtk ${TEST_DATA_ROOT}/movingMeshWithScalarsIC1.vtk ${TEST_DATA_ROOT}/fixedMeshWithScalarsIC1.vtk 100.0 ${TEST_DATA_ROOT}/LandmarkIDsOnFixedIC1.txt ${TEMP}/deformedFixedMeshWithScalarsIC1.vtk )
