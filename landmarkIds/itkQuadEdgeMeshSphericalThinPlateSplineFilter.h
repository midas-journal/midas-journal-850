/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: ITKHeader.h,v $
 Language:  C++
 Date:      $Date: 2006-04-25 23:20:16 $
 Version:   $Revision: 1.1 $

 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkQuadEdgeMeshSphericalThinPlateSplineFilter_h
#define __itkQuadEdgeMeshSphericalThinPlateSplineFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
    /** \class QuadEdgeMeshSphericalThinPlateSplineFilter
     * This class defines the thin plate spline (TPS) transformation on spheres
     * It is implemented in as straightforward a manner as possible from
     * the MICCAI paper 2007;10(Pt 1):367-74.
     *
     * It takes source landmarks and target landmarks
     * Input mesh: fixed mesh (corresponding to the source landmarks)
     *
     * Output deformed fixed mesh
     *
     */
    template< class TFixedMesh, class TOutputMesh >
    class QuadEdgeMeshSphericalThinPlateSplineFilter :
    public QuadEdgeMeshToQuadEdgeMeshFilter< TFixedMesh, TOutputMesh >
    {
    public:
        typedef QuadEdgeMeshSphericalThinPlateSplineFilter                    Self;
        typedef SmartPointer< Self >                                          Pointer;
        typedef SmartPointer< const Self >                                    ConstPointer;
        typedef QuadEdgeMeshToQuadEdgeMeshFilter< TFixedMesh, TOutputMesh >   Superclass;

        /** Method that instantiates a new object */
        itkNewMacro( Self );

        /** Method that provides the name of the class as a string as well as the
         * name of the parent class. */
        itkTypeMacro( QuadEdgeMeshSphericalThinPlateSplineFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

        /** Input types. */
        typedef TFixedMesh                                        FixedMeshType;
        typedef typename  FixedMeshType::Pointer                  FixedMeshPointer;
        typedef typename  FixedMeshType::ConstPointer             FixedMeshConstPointer;
        typedef typename  FixedMeshType::PointType                PointType;
        typedef typename  FixedMeshType::PixelType                FixedPixelType;
        typedef typename  FixedMeshType::PointsContainer          FixedPointsContainer;
        typedef typename  FixedPointsContainer::Iterator          FixedPointsIterator;
        typedef typename  FixedPointsContainer::ConstIterator     FixedPointsConstIterator;

        typedef typename  PointType::VectorType                   VectorType;

        /** Output types. */
        typedef TOutputMesh                                       OutputMeshType;
        typedef typename  OutputMeshType::PointType               OutputPointType;
        typedef typename  OutputMeshType::Pointer                 OutputMeshPointer;
        typedef typename  OutputMeshType::PointsContainer         OutputPointsContainer;
        typedef typename  OutputPointsContainer::Iterator         OutputPointsIterator;

        /** Set/Get the Fixed mesh. */
        void SetFixedMesh( const FixedMeshType * fixedMesh );
        itkGetConstObjectMacro( FixedMesh, FixedMeshType );

        itkStaticConstMacro( PointDimension, unsigned int, FixedMeshType::PointDimension );

        typedef typename FixedMeshType::Traits                                FixedMeshTraits;
        typedef PointSet< FixedPixelType, PointDimension, FixedMeshTraits >   PointSetType;
        typedef typename PointSetType::Pointer                                PointSetPointer;
        typedef typename PointSetType::PointsContainer                        PointsContainer;
        typedef typename PointSetType::PointsContainerIterator                PointsIterator;
        typedef typename PointSetType::PointsContainerConstIterator           PointsConstIterator;
        typedef typename PointSetType::PointIdentifier                        PointIdentifier;

        /** Get the source landmarks. */
        itkGetObjectMacro( SourceLandmarks, PointSetType);

        /** Set the source landmarks. */
        itkSetObjectMacro( SourceLandmarks, PointSetType);

        /** Get the target landmarks. */
        itkGetObjectMacro( TargetLandmarks, PointSetType);

        /** Set the target landmarks. */
        itkSetObjectMacro( TargetLandmarks, PointSetType);

        /** Set Sphere Center.  The implementation of this filter assumes that the
         * Mesh surface has a spherical geometry (not only spherical topology). With
         * this method you can specify the coordinates of the center of the sphere
         * represented by the Mesh. This will be used in the coordinate transform
         * from Euclidean to Spherical and verse visa.
         */
        itkSetMacro( SphereCenter, PointType );
        itkGetConstMacro( SphereCenter, PointType );

        /** Set Sphere Radius.  The implementation of this filter assumes that the
         * Mesh surface has a spherical geometry (not only spherical topology). With
         * this method you can specify the radius of the sphere.
         * This will be used in the coordinate transform
         * from Euclidean to Spherical and verse visa.
         */
        itkSetMacro( SphereRadius, double );
        itkGetConstMacro( SphereRadius, double );

    protected:
        QuadEdgeMeshSphericalThinPlateSplineFilter();
        ~QuadEdgeMeshSphericalThinPlateSplineFilter();
        void PrintSelf(std::ostream& os, Indent indent) const;

        virtual void GenerateData( );

        /** The list of source landmarks, denoted 'p'. */
        PointSetPointer m_SourceLandmarks;

        /** The list of target landmarks, denoted 'q'. */
        PointSetPointer m_TargetLandmarks;

    private:
        QuadEdgeMeshSphericalThinPlateSplineFilter( const Self& );
        void operator = ( const Self& );

        /** matrix typedef. */

        /** 'K' matrix typedef. */
        typedef vnl_matrix<double> KMatrixType;

        /** 'd' matrix typedef. */
        typedef vnl_matrix<double> dMatrixType;

        /** 'c' matrix typedef. */
        typedef vnl_matrix<double> cMatrixType;

        /** 'T' matrix typedef. */
        typedef vnl_matrix<double> TMatrixType;

        /** 'z' matrix typedef. */
        typedef vnl_matrix<double> zMatrixType;

        /** 'TTInKTTTInK' matrix typedef. */
        typedef vnl_matrix<double> TTInKTTTInKMatrixType;

        /** Evaluate the value of K(X,Y) according to Eq.(2)
         * a = cos(angle(X,Y))
         */
        double EvaluateK (const double& a) const;

        /** Compute the n x n KMatrix with angles between
         *  each pair of landmarks in m_SourceLandMarks.
         *  Compute the TMatrix
         */
        void ComputeKAndT();

        /** Compute m_TTInKTTTInK=((T^T)(K^-1)T)^-1(T^T)(K^-1)
         *  the result is a row vector with n elements.
         *  n: number of landmarks
         *  which exits both in cMatrix and dMatrix
         *
         ***has to be called after ComputeT() and ComputeK()
         */
        void ComputeTTInKTTTInK();

        /** Compute dMatrix
         *  which is actually a vector with two elements (theta, phi)
         *  (m_TTInKTTTInK)z
         */
        void Compute_d();

        /** Compute cMatrix
         */
        void Compute_c();

        /** Compute zMatrix
         *  it's a n x 2 matrix to keep (theta, phi) for
         *  m_TargetLandmarks
         ***it must be called before Compute_c() and Compute_d()
         */
        void ComputeZ();

        FixedMeshConstPointer                 m_FixedMesh;

        PointType       m_SphereCenter;

        /** Radius of spherical mesh. We assume that both the source landmarks and
         * the target landmarks are from spheres that share the same
         * center and radius.
         */
        double          m_SphereRadius;

        /** KMatrix generated by the source landmarks and the target landmarks
         * with functionK defined in EvaluateK
         */
        KMatrixType m_KMatrix;

        cMatrixType m_cMatrix;

        dMatrixType m_dMatrix;

        /** TMatrix: a column vector with n ones
         *  n equals to the number of landmarks
         */
        TMatrixType m_TMatrix;

        zMatrixType m_zMatrix;

        /** ((T^T)(K^-1)T)^-1(T^T)(K^-1)
         *  1 x n row vector n: number of landmarks
         */
        TTInKTTTInKMatrixType m_TTInKTTTInK;

    };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadEdgeMeshSphericalThinPlateSplineFilter.txx"
#endif

#endif
