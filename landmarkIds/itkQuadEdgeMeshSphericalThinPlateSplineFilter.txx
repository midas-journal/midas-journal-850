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
#ifndef __itkQuadEdgeMeshSphericalThinPlateSplineFilter_txx
#define __itkQuadEdgeMeshSphericalThinPlateSplineFilter_txx

#include "itkQuadEdgeMeshSphericalThinPlateSplineFilter.h"

#include "itkVector.h"
#include "vnl/vnl_vector_fixed.h"
#include "itkMatrix.h"

#include "itkProgressReporter.h"
namespace itk
{

template< class TFixedMesh, class TOutputMesh >
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >
::QuadEdgeMeshSphericalThinPlateSplineFilter()
{
    this->SetNumberOfInputs( 1 );
    this->SetNumberOfOutputs( 1 );

    this->SetNumberOfRequiredInputs( 1 );
    this->SetNumberOfRequiredOutputs( 1 );

    this->SetNthOutput( 0, OutputMeshType::New() );

    this->m_SphereCenter.Fill( 0.0 );
    this->m_SphereRadius = 1.0;
}

template< class TFixedMesh, class TOutputMesh >
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >
::~QuadEdgeMeshSphericalThinPlateSplineFilter()
{
}

template< class TFixedMesh, class TOutputMesh >
void
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >
::SetFixedMesh( const FixedMeshType * fixedMesh )
{
    itkDebugMacro("setting Fixed Mesh to " << fixedMesh );

    if (this->m_FixedMesh.GetPointer() != fixedMesh )
    {
        this->m_FixedMesh = fixedMesh;

        // Process object is not const-correct so the const_cast is required here
        this->ProcessObject::SetNthInput(0, const_cast< FixedMeshType *>( fixedMesh ) );

        this->Modified();
    }
}


template< class TFixedMesh, class TOutputMesh >
void
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >
::GenerateData()
{
  if( this->m_TargetLandmarks.IsNull() )
    {
    itkGenericExceptionMacro( << "m_TargetLandmarks is null" );
    }

  if( this->m_SourceLandmarks.IsNull() )
    {
    itkGenericExceptionMacro( << "m_SourceLandmarks is null" );
    }

  // Prepare data
  this->CopyInputMeshToOutputMesh();

  //Prepare for matrices
  this->ComputeKAndT();
  this->ComputeZ();
  this->ComputeTTInKTTTInK();

  this->Compute_d();
  this->Compute_c();

  OutputPointsContainer * meshPoints = this->GetOutput()->GetPoints();

  OutputPointsIterator pointItr = meshPoints->Begin();
  OutputPointsIterator pointEnd = meshPoints->End();

  PointsContainer * sourcePoints = this->m_SourceLandmarks->GetPoints();

  PointsIterator sourceItr = sourcePoints->Begin();
  PointsIterator sourceEnd = sourcePoints->End();

  // Compute the displacement for each point on the output mesh
  // convert the displacement from (theta, phi) back to (x,y,z)
  PointType pi,destination_p;
  VectorType pToCenter, piToCenter, destination_pToCenter;
  double theta_p, phi_p, a, k_pTopi, new_theta_p, new_phi_p, radialDistance;

  typedef vnl_vector_fixed<double,2> PolarVectorType;
  //sum of ci*k_pTopi, i=1...n
  //u_p_polar: displacement in (theta,phi)
  PolarVectorType ci,sum_ck;

  while ( pointItr != pointEnd )
    {
    //transform to (theta,phi)
    pToCenter = pointItr.Value() - this->m_SphereCenter;

    theta_p = vcl_acos(pToCenter[2]/this->m_SphereRadius);
    phi_p = vnl_math::pi + vcl_atan2(-pToCenter[1],-pToCenter[0]);

    //we need m_SourceLandmarks to do the interpolation
    sum_ck[0] = 0.0;
    sum_ck[1] = 0.0;

    sourceItr = sourcePoints->Begin();
    int i=0;
    while ( sourceItr != sourceEnd )
    {
      piToCenter = sourceItr.Value() - this->m_SphereCenter;

      pToCenter.Normalize();
      piToCenter.Normalize();

      //a=cos( angle(p,pi) )
      a = pToCenter * piToCenter;
      k_pTopi = this->EvaluateK(a);

      //ci: i_th row of m_cMatrix
      ci[0] = this->m_cMatrix(i,0);
      ci[1] = this->m_cMatrix(i,1);

      sum_ck[0] += ci[0] * k_pTopi;
      sum_ck[1] += ci[1] * k_pTopi;

      ++sourceItr;
      ++i;
    }

    new_theta_p = theta_p + sum_ck[0] + this->m_dMatrix(0,0);
    new_phi_p = phi_p + sum_ck[1] + this->m_dMatrix(0,1);

    //project back to the sphere to get the new position for p
    destination_p[0] = this->m_SphereCenter[0] + this->m_SphereRadius * vcl_cos(new_phi_p) * vcl_sin(new_theta_p);
    destination_p[1] = this->m_SphereCenter[1] + this->m_SphereRadius * vcl_sin(new_phi_p) * vcl_sin(new_theta_p);
    destination_p[2] = this->m_SphereCenter[2] + this->m_SphereRadius * vcl_cos(new_theta_p);

    destination_pToCenter = destination_p - this->m_SphereCenter;
    radialDistance = destination_pToCenter.GetNorm();

    destination_pToCenter *= this->m_SphereRadius / radialDistance;

    destination_p = this->m_SphereCenter + destination_pToCenter;

    //assign the destination_p to the output mesh
    pointItr.Value() = destination_p;

    ++pointItr;
    }
}

/**
 * Evaluate the integral function from 0 to 1 using n points
 */
template< class TFixedMesh, class TOutputMesh >
double
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >
::EvaluateK(const double & a) const
{
  //****modified from itkFEMUtility
  double w[110] = {
        0.0, 1.0, 0.555555555555556, 0.888888888888889,
        0.347854845137454, 0.652145154862546, 0.236926885056189,
        0.478628670499366, 0.568888888888889, 0.171324492379170,
        0.360761573048139, 0.467913934572691, 0.129484966168870,
        0.279705391489277, 0.381830050505119, 0.417959183673469,
        0.101228536290376, 0.222381034453374, 0.313706645877887,
        0.362683783378362, 0.0812743883615744, 0.180648160694857,
        0.260610696402935, 0.312347077040003, 0.330239355001260,
        0.0666713443086881, 0.149451349150581, 0.219086362515982,
        0.269266719309996, 0.295524224714753, 0.0556685671161737,
        0.125580369464905, 0.186290210927734, 0.233193764591990,
        0.262804544510247, 0.272925086777901, 0.0471753363865118,
        0.106939325995318, 0.160078328543346, 0.203167426723066,
        0.233492536538355, 0.249147045813403, 0.0404840047653159,
        0.0921214998377284, 0.138873510219787, 0.178145980761946,
        0.207816047536889, 0.226283180262897, 0.232551553230874,
        0.0351194603317519,0.0801580871597602, 0.121518570687903,
        0.157203167158194, 0.185538397477938, 0.205198463721296,
        0.215263853463158, 0.0307532419961173,0.0703660474881081,
        0.107159220467172, 0.139570677926154, 0.166269205816994,
        0.186161000015562, 0.198431485327112, 0.202578241925561,
        0.0271524594117541,0.0622535239386478,0.0951585116824928,
        0.124628971255534, 0.149595988816577, 0.169156519395003,
        0.182603415044924, 0.189450610455068, 0.0241483028685479,
        0.0554595293739872, 0.0850361483171792,0.111883847193404,
        0.135136368468525, 0.154045761076810, 0.168004102156450,
        0.176562705366993, 0.179446470356207, 0.0216160135264833,
        0.0497145488949698, 0.0764257302548891,0.100942044106287,
        0.122555206711478, 0.140642914670651, 0.154684675126265,
        0.164276483745833, 0.169142382963144, 0.0194617882297265,
        0.0448142267656996,0.0690445427376412,0.0914900216224500,
        0.111566645547334, 0.128753962539336, 0.142606702173607,
        0.152766042065860, 0.158968843393954, 0.161054449848784,
        0.0176140071391521,0.0406014298003869,0.0626720483341091,
        0.0832767415767047, 0.101930119817240, 0.118194531961518,
        0.131688638449177, 0.142096109318382, 0.149172986472604,
        0.152753387130726
    };

  double z[110] = {
        0.0, 0.577350269189626, 0.774596669241483, 0.0,
        0.861136311594053, 0.339981043584856, 0.906179845938664,
        0.538469310105683, 0.0, 0.932469514203152,
        0.661209386466265, 0.238619186083197, 0.949107912342759,
        0.741531185599394, 0.405845151377397, 0.0,
        0.960289856497536, 0.796666477413627, 0.525532409916329,
        0.183434642495650, 0.968160239507626, 0.836031107326636,
        0.613371432700590, 0.324253423403809, 0.0,
        0.973906528517172, 0.865063366688985, 0.679409568299024,
        0.433395394129247, 0.148874338981631, 0.978228658146057,
        0.887062599768095, 0.730152005574049, 0.519096129206812,
        0.269543155952345, 0.0, 0.981560634246719,
        0.904117256370475, 0.769902674194305, 0.587317954286617,
        0.367831498998180, 0.125233408511469, 0.984183054718588,
        0.917598399222978, 0.801578090733310, 0.642349339440340,
        0.448492751036447, 0.230458315955135, 0.0,
        0.986283808696812, 0.928434883663574, 0.827201315069765,
        0.687292904811685, 0.515248636358154, 0.319112368927890,
        0.108054948707344, 0.987992518020485, 0.937273392400706,
        0.848206583410427, 0.724417731360170, 0.570972172608539,
        0.394151347077563, 0.201194093997435, 0.0,
        0.989400934991650, 0.944575023073233, 0.865631202387832,
        0.755404408355003, 0.617876244402644, 0.458016777657227,
        0.281603550779259, 0.0950125098376374,0.990575475314417,
        0.950675521768768, 0.880239153726986, 0.781514003896801,
        0.657671159216691, 0.512690537086477, 0.351231763453876,
        0.178484181495848, 0.0, 0.991565168420931,
        0.955823949571398, 0.892602466497556, 0.803704958972523,
        0.691687043060353, 0.559770831073948, 0.411751161462843,
        0.251886225691506, 0.0847750130417353,0.992406843843584,
        0.960208152134830, 0.903155903614818, 0.822714656537143,
        0.720966177335229, 0.600545304661681, 0.464570741375961,
        0.316564099963630, 0.160358645640225, 0.0,
        0.993128599185095, 0.963971927277914, 0.912234428251326,
        0.839116971822219, 0.746331906460151, 0.636053680726515,
        0.510867001950827, 0.373706088715420, 0.227785851141645,
        0.0765265211334973
   };

  double scale = (1.0 - 0.0)/2.0; // 0.5 ?
  int n = 20;
  int m, ibase,i;
  double sum,h,k_2,k_2_tl,k_2_tu,t,tl,tu;
  if ( (n&1) == 0 )
    {
    m = n/2;
    ibase = m*m;
    sum = 0.0;
    }
  else
    {
    m = (n - 1)/2;
    ibase = (n*n - 1)/4;
    h = (1.0 + 0.0)/2.0;
    k_2 = 0.0;
    if (h!=a)
      {
      k_2 = vcl_log10(h)*(1-1/h)*(1/vcl_sqrt(1-2*h*a+h*h)-1);
      }
    sum = w[ibase+m]*k_2;
    }

  for (i=1; i <= m; i++)
    {
    t = z[ibase + i - 1];
    tl = (1.0 - t)/2.0;
    tu = (1.0 + t)/2.0;
    k_2_tl = 0.0;
    k_2_tu = 0.0;

    if(tl!=0.0 && tl!=a)
      {
      k_2_tl = vcl_log10(tl)*(1-1/tl)*(1/vcl_sqrt(1-2*tl*a+tl*tl)-1);
      }
    if (tu!=0 && tu!=a)
      {
      k_2_tu = vcl_log10(tu)*(1-1/tu)*(1/vcl_sqrt(1-2*tu*a+tu*tu)-1);
      }
    sum = sum + w[ibase + i - 1]*(k_2_tl + k_2_tu);
    }

  //****modified from itkFEMUtility

  sum = scale*sum/(4.0*vnl_math::pi);
  return sum;
}

/**
* compute m_KMatrix with landmarks
*/
template< class TFixedMesh, class TOutputMesh >
void
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >::
ComputeKAndT()
{
  PointIdentifier numberOfLandmarks = this->m_SourceLandmarks->GetNumberOfPoints();

  //initiate m_KMatrix
  this->m_KMatrix.set_size(numberOfLandmarks,numberOfLandmarks);
  this->m_KMatrix.fill(0.0);

  PointType Pi;
  PointType Pj;
  double a = 0.0;
  double k_ij = 0.0;

  double k_ii = this->EvaluateK(1.0);
  for( PointIdentifier i = 0; i<numberOfLandmarks; i++ )
    {
    for( PointIdentifier j = i; j<numberOfLandmarks; j++)
      {
      if (i==j)
        {
        // cos(angle(pi,pj)) = 1.0
        this->m_KMatrix(i,i) = k_ii;
        }
      else if (this->m_KMatrix(i,j) == 0.0)
        {
        //(Pi,Pj) or (Pj,Pi) were never visited
        this->m_SourceLandmarks->GetPoint(i,&Pi);
        this->m_SourceLandmarks->GetPoint(j,&Pj);

        //compute cos(angle(Pi,Pj))
        VectorType PiToCenter = Pi - this->m_SphereCenter;
        VectorType PjToCenter = Pj - this->m_SphereCenter;

        PiToCenter.Normalize();
        PjToCenter.Normalize();

        //a=cos(angle(pi,pj))
        a = PiToCenter * PjToCenter;

        k_ij = this->EvaluateK(a);

        //save k_ij to KMatrix(i,j) and KMatrix(j,i)
        this->m_KMatrix(i,j) = k_ij;
        this->m_KMatrix(j,i) = k_ij;
        }
      }
    }

  this->m_TMatrix.set_size(numberOfLandmarks,1);
  this->m_TMatrix.fill(1.0);
}

/**
 * compute m_TTInKTTTInK
 * ((T^T)(K^-1)T)^-1(T^T)(K^-1)
 *  1 x n row vector n: number of landmarks
 *
 *  Has to be called after ComputeT and ComputeK
 */
template< class TFixedMesh, class TOutputMesh >
void
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >::
ComputeTTInKTTTInK()
{
  //m_TTInKTTTInK is a 1 x n row vector
  PointIdentifier numberOfLandmarks = this->m_SourceLandmarks->GetNumberOfPoints();

  this->m_TTInKTTTInK.set_size(1, numberOfLandmarks);
  this->m_TTInKTTTInK.fill(0.0);

  KMatrixType InK;
  InK.set_size(numberOfLandmarks, numberOfLandmarks);

  // not sure if it is the most efficient to inverse the matrix
  InK = vnl_matrix_inverse<double>(this->m_KMatrix);

  KMatrixType TTInKT;
  TTInKT.set_size(1,1);
  TTInKT = this->m_TMatrix.transpose() * InK * this->m_TMatrix;

  double TTInKT_value = TTInKT(0,0);

  //itkAssertInDebugAndIgnoreInReleaseMacro( TTInKT_value != 0. );

  if ( TTInKT_value != 0.0 ){
    this->m_TTInKTTTInK = 1.0 / TTInKT_value * this->m_TMatrix.transpose() * InK;
  }
}


/** Compute zMatrix
 *  it's a n x 2 matrix to keep (delta_theta, delta_phi)
 *  --displacement from m_SourceLandmarks to m_TargetLandmarks
 ***it must be called before Compute_c() and Compute_d()
 */
template< class TFixedMesh, class TOutputMesh >
void
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >::
ComputeZ()
{
  PointIdentifier numberOfLandmarks = this->m_SourceLandmarks->GetNumberOfPoints();

  this->m_zMatrix.set_size(numberOfLandmarks,2);
  this->m_zMatrix.fill(0.0);

  PointsContainer * sourcePoints = this->m_SourceLandmarks->GetPoints();

  PointsIterator sourceItr = sourcePoints->Begin();
  PointsIterator sourceEnd = sourcePoints->End();

  PointsContainer * targetPoints = this->m_TargetLandmarks->GetPoints();

  PointsIterator targetItr = targetPoints->Begin();
  PointsIterator targetEnd = targetPoints->End();

  int i=0;
  while ( sourceItr != sourceEnd )
  {
    //get one target landmark

    //transform it from (x,y,z) to (theta,phi)
    //using m_SphereCenter and m_SphereRadius
    VectorType siToCenter = sourceItr.Value() - this->m_SphereCenter;
    VectorType tiToCenter = targetItr.Value() - this->m_SphereCenter;

    //keep the displacement
    double theta = vcl_acos(tiToCenter[2]/this->m_SphereRadius) - vcl_acos(siToCenter[2]/this->m_SphereRadius);
    double phi = vcl_atan2(-tiToCenter[1],-tiToCenter[0])- vcl_atan2(-siToCenter[1],-siToCenter[0]);

    this->m_zMatrix(i,0) = theta;
    this->m_zMatrix(i,1) = phi;

    ++sourceItr;
    ++targetItr;
    ++i;
  }
}

/** Compute dMatrix
*  it's a 1 x 2 matrix to have (d_theta, d_phi)
***it must be called after ComputeTTInKTTTInK() and ComputeZ()
*/
template< class TFixedMesh, class TOutputMesh >
void
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >::
Compute_d()
{
    this->m_dMatrix.set_size(1,2);
    this->m_dMatrix.fill(0.0);

    //d=(m_TTInKTTTInK)z
    this->m_dMatrix = this->m_TTInKTTTInK * this->m_zMatrix;
}

/** Compute cMatrix
*  it's a n x 2 matrix
***it must be called after ComputeZ() and ComputeK()
*  and ComputeTTInKTTTInK()
*/
template< class TFixedMesh, class TOutputMesh >
void
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >::
Compute_c()
{
    //it's all about m_TargetLandmarks
    int numberOfLandmarks = this->m_SourceLandmarks->GetNumberOfPoints();

    this->m_cMatrix.set_size(numberOfLandmarks,2);
    this->m_cMatrix.fill(0.0);

    //inverse(m_KMatrix)
    KMatrixType InK;
    InK.set_size(numberOfLandmarks,numberOfLandmarks);

    InK = vnl_matrix_inverse<double>(this->m_KMatrix);

    //I: n x n identity matrix
    KMatrixType IMatrix;
    IMatrix.set_size(numberOfLandmarks,numberOfLandmarks);
    IMatrix.set_identity();

    this->m_cMatrix = InK * (IMatrix - this->m_TMatrix * this->m_TTInKTTTInK) * this->m_zMatrix;
}

template< class TFixedMesh, class TOutputMesh >
void
QuadEdgeMeshSphericalThinPlateSplineFilter< TFixedMesh, TOutputMesh >::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}
}

#endif
