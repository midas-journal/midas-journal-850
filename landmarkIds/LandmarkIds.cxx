/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: LandmarkIds.cxx,v $
  Language:  C++
  Date:      $Date: 2012-01-31 17:44:24 $
  Version:   $Revision: 1.0 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
/*
 * landmarFile: save pointIds for source and target landmarks
 * source landmark Ids are for inputFixedMesh
 * target landmark Ids are for inputMovingMesh
 */

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"

#include "itkVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "itkQuadEdgeMeshSphericalThinPlateSplineFilter.h"

#include <fstream>

int main( int argc, char * argv [] )
{
    if( argc < 5 )
    {
        std::cerr << "Missing arguments" << std::endl;
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << std::endl;
        std::cerr << "inputFixedMesh inputMovingMesh fixedMeshIC4 ";
        std::cerr << "radius landmarksFile ";
        std::cerr << "outputDeformedFixedMesh ";
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

    typedef float      MeshPixelType;
    const unsigned int Dimension = 3;

    typedef itk::QuadEdgeMesh< MeshPixelType, Dimension >   MeshType;

    typedef itk::VTKPolyDataReader< MeshType >     ReaderType;

    ReaderType::Pointer fixedReader = ReaderType::New();
    fixedReader->SetFileName( argv[1] );
    fixedReader->Update( );
    
    ReaderType::Pointer movingReader = ReaderType::New();
    movingReader->SetFileName( argv[2] );
    movingReader->Update( );
    
    ReaderType::Pointer fixedIC4Reader = ReaderType::New();
    fixedIC4Reader->SetFileName(argv[3]);
    fixedIC4Reader->Update();
    
    MeshType::Pointer  meshFixed  = fixedReader->GetOutput();
    MeshType::Pointer  meshMoving  = movingReader->GetOutput();
    
    typedef itk::QuadEdgeMeshSphericalThinPlateSplineFilter< MeshType, MeshType >    SphericalTPSType;

    SphericalTPSType::Pointer   tpsFilter  = SphericalTPSType::New();
    
    double radius = atof( argv[4] );
    tpsFilter->SetSphereRadius(radius);
    
    SphericalTPSType::PointType center;
    center.Fill( 0.0 );
    tpsFilter->SetSphereCenter(center);
    
    std::ifstream pointsFile;
    pointsFile.open( argv[5] );
    
    //Create source and target landmarks
    typedef MeshType::PointType LandmarkPointType;
    typedef SphericalTPSType::PointSetType LandmarkPointSetType;
    LandmarkPointSetType::Pointer sourceLandmarks = LandmarkPointSetType::New();
    LandmarkPointSetType::Pointer targetLandmarks = LandmarkPointSetType::New();

    MeshType::PointIdentifier sourceId;
    MeshType::PointIdentifier targetId;
    
    LandmarkPointType sourcePoint;
    LandmarkPointType targetPoint;
  
    unsigned int i = 0;
    
    pointsFile >> sourceId;
    pointsFile >> targetId;
    
    while( !pointsFile.fail() )
    {
        //get the point through point ID
        meshFixed->GetPoint(sourceId,&sourcePoint);
        sourceLandmarks->SetPoint( i, sourcePoint );

        meshMoving->GetPoint(targetId,&targetPoint);
        targetLandmarks->SetPoint( i, targetPoint );
        
        std::cout<<sourcePoint<<std::endl;
        std::cout<<targetPoint<<std::endl;
        i++;
        
        pointsFile >> sourceId;
        pointsFile >> targetId;
        
    }
    
    pointsFile.close();

    tpsFilter->SetFixedMesh(fixedIC4Reader->GetOutput());
    tpsFilter->SetSourceLandmarks(sourceLandmarks);
    tpsFilter->SetTargetLandmarks(targetLandmarks);
    tpsFilter->Update();
  
    typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
    WriterType::Pointer demonsWriter = WriterType::New();
    
    demonsWriter->SetInput( tpsFilter->GetOutput() );
    demonsWriter->SetFileName(argv[6]);
    demonsWriter->Update(); 

    return EXIT_SUCCESS;
}
