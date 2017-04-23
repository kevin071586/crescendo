#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"

#include "Shards_CellTopology.hpp"              // For element topology

#include "Element.h"

using namespace Intrepid;
typedef FunctionSpaceTools fst;

// Constructor and destructor
// ============================================================================
Element::Element() {
}

Element::~Element() {
}


// Initialize cubature
// ============================================================================
void Element::initializeCubature()
{
  const int NUM_NODES_HEX8 = 8;
  const int CUBATURE_DEGREE = 2;
  const int SPATIAL_DIM = 3;

  // Define cell topology of the parent cell
  m_elemTopology = shards::getCellTopologyData<shards::Hexahedron<NUM_NODES_HEX8>>();

  // Create cubature
  Intrepid::DefaultCubatureFactory<double> cubFactory;
  Teuchos::RCP<Intrepid::Cubature<double>> hexCubature =
                                  cubFactory.create(m_elemTopology, CUBATURE_DEGREE);

  // Define basis
  typedef Intrepid::FieldContainer<double> FieldContainer;
  Intrepid::Basis_HGRAD_HEX_C1_FEM<double, FieldContainer> hexHGradBasis;
  m_numFields = hexHGradBasis.getCardinality();
  m_numCubPoints = hexCubature->getNumPoints();
  
  // Field containers for cubature
  m_cubDim = hexCubature->getDimension();
  m_cubWeights.resize(m_numCubPoints);
  m_cubPoints.resize(m_numCubPoints, m_cubDim);
  m_cubPointValues.resize(m_numFields, m_numCubPoints);
  m_cubPointGradient.resize(m_numFields, m_numCubPoints, SPATIAL_DIM); 

  // Get cubature points and weights from the parent element definition
  hexCubature->getCubature(m_cubPoints, m_cubWeights);
  hexHGradBasis.getValues(m_cubPointValues, m_cubPoints, Intrepid::OPERATOR_VALUE);

  // get gradient values for hex8 element (note: number of values here is
  // equal to cardinality of basis -- no need for # of cells).
  hexHGradBasis.getValues(m_cubPointGradient, m_cubPoints, Intrepid::OPERATOR_GRAD);

  return;
}


//-------------------------------------------------------------------------------- 
// Function to calculate stiffness matrix for a given shape function pair (a,b)
// such that
//
//    K_ab = int( N_a,j*Cijkl*N_b,l ) dV
//
// is an 8x8 matrix indexed by the cardinality (number of interpolation functions
// needed.  There are 9 such matrices when each node has 3 DOF corresponding to 
// Cijkl pairs with (i,k) indices fixed at (0,0), (1,0), (2,0), ... (2,2).
//
// Output: FieldContainer K_ab for the given (i,j) pair.
//
// This function wraps the FunctionSpaceTools::integrate() method.
//
//-------------------------------------------------------------------------------- 
void Element::integrateStiffnessMatrixAB(FieldContainer<double>& stiffnessMatrixAB, 
                                         FieldContainer<double> leftArg, 
                                         FieldContainer<double> rightArg, 
                                         FieldContainer<double> elasticModuli, int i, int k ) {

  // Get container dimensions
  const int num_elements = leftArg.dimension(0);
  const int num_shape_fields = leftArg.dimension(1);
  const int num_cub_points = leftArg.dimension(2);
  const int space_dim = leftArg.dimension(3);

  // Define some field containers
  FieldContainer<double> leftArgScaled(num_elements, num_shape_fields, 
    num_cub_points, space_dim);

  // Multiply shape function gradients (C,F,P,D1) by elastic tensor with a given (i,k) 
  // index pair held fixed, such that dimensions are (D1,D2) or (P,D1,D2) as needed.
  //
  // Output data should be (C,F,P,D1)
  //
  // Can't seem to find the right function to do this, so hacking one in myself.

  // Fill temp field container
  FieldContainer<double> C_temp(space_dim, space_dim);
  for( int j=0; j<space_dim; ++j ){
    for( int l=0; l<space_dim; ++l ){
      C_temp(j,l) = elasticModuli(i,j,k,l);  // fix i=0, k=0;
    }
  }

  // Inner product between field container and elasticity tensor
  for( int C=0; C < num_elements; ++C ){
    for( int F=0; F < num_shape_fields; ++F ){
      for( int P=0; P < num_cub_points; ++P ){
        for( int D=0; D < space_dim; ++D ){
          leftArgScaled(C,F,P,D) =
            leftArg(C,F,P,0)*C_temp(0,D) + 
            leftArg(C,F,P,1)*C_temp(1,D) + 
            leftArg(C,F,P,2)*C_temp(2,D); 
        }
      }
    }
  }

  // integrate stiffness matrix
  fst::integrate<double>(stiffnessMatrixAB, leftArgScaled, rightArg, Intrepid::COMP_CPP);

  return;
}

//-------------------------------------------------------------------------------- 
// 
// This function computes the complete stiffness matrix by evaluating 9 times, 
// once for each Cijkl pair (i,k).
//
// HEX8: The complete element stiffness matrix can be viewed as an 8x8 set of 
// submatrices, where each submatrix K_ik is a 3x3 corresponding to a single node.
// Displacements must then be in a vector u=[u_1, u_2, ... u_8] where u_1 is a
// sub-vector u_1 = [u_1x, u1y, u1z]
//
// The submatrices K_ik(a,b) are obtained by interlacing the K_ab submatrices:
//
// K = [ K_11(1,1), K_12(1,1), K_13(1,1), | K_11(1,2), ... K_13(1,8) ]
//     [ K_21(1,1), K_22(1,1), K_23(1,1), | K_21(1,2), ...           ]
//     [ K_31(1,1), K_32(1,1), K_33(1,1), | K_31(1,2), ...           ]
//     ---------------------------------------------------------------
//     [ K_11(2,1), K_12(2,1), K_13(2,1), | K_11(2,2), ...           ]
//     [                                                             ]
//     [                                                             ]
//     [ ...                                               K_33(8,8) ]
//
// Output: 
//  stiffnessMatrix - (C,24,24) field container (HEX8) 
//
//-------------------------------------------------------------------------------- 
void Element::integrateStiffnessMatrix(FieldContainer<double>& stiffnessMatrix, 
                                       FieldContainer<double> leftArg, 
                                       FieldContainer<double> rightArg, 
                                       FieldContainer<double> C) {

  const int num_elements = leftArg.dimension(0);
  const int space_dim = leftArg.dimension(3);
  const int num_shape_fields = leftArg.dimension(1);

  FieldContainer<double> stiffnessMatrixAB(num_elements, num_shape_fields, num_shape_fields);

  for( int i=0; i<space_dim; ++i ){
    for( int k=0; k<space_dim; ++k ){
      integrateStiffnessMatrixAB(stiffnessMatrixAB, leftArg, rightArg, C, i, k); 

      for( int a=0; a<num_shape_fields; ++a ){
        for( int b=0; b<num_shape_fields; ++b ){
          const int idx_row = a*space_dim + i;
          const int idx_col = b*space_dim + k;

          for( int idx_cell=0; idx_cell < num_elements; ++idx_cell ){
            stiffnessMatrix(idx_cell, idx_row, idx_col) = stiffnessMatrixAB(idx_cell,a,b);
          }
        }
      }
    }
  }

  return;
}


//-------------------------------------------------------------------------------- 
// Simple wrapper around FunctionSpaceTools::integrate for computing a 
// mass matrix.  Could also include mass-lumping functionality eventually.
// 
// integrate: mass matrix calculations
//
// Output:        2x8x8   (C,L,R)
// Left Fields:   2x8     (C,L,P)
// Right Fields:  2x8     (C,R,P)
//
// This calculation contracts, for all cells "C", point dimension "P" which
// denotes the number of integration points.
//-------------------------------------------------------------------------------- 
void Element::integrateMassMatrixAB(FieldContainer<double>& massMatrixAB, 
                                   FieldContainer<double> leftArg, 
                                   FieldContainer<double> rightArg,
                                   const double rho) {

  // scale left argument by density: assuming constant over the element.
  const int num_elements = leftArg.dimension(0);
  const int num_shape_fields = leftArg.dimension(1);
  const int num_cub_points = leftArg.dimension(2);

  FieldContainer<double> leftArgScaled(num_elements, num_shape_fields, num_cub_points);

  for(int C=0; C<num_elements; ++C) {
    for(int i=0; i<num_shape_fields; ++i){
      for(int j=0; j<num_cub_points; ++j){
        leftArgScaled(C,i,j) = rho*leftArg(C,i,j);
      }
    }
  }

  // integrate mass matrix
  fst::integrate<double>(massMatrixAB, leftArgScaled, rightArg, Intrepid::COMP_CPP);
  return;
}

//-------------------------------------------------------------------------------- 
// Compute full mass matrix for an element, in 3 dimensions.
//-------------------------------------------------------------------------------- 
void Element::integrateMassMatrix(FieldContainer<double>& massMatrix, 
                                  FieldContainer<double> leftArg, 
                                  FieldContainer<double> rightArg,
                                  const double rho,
                                  const unsigned space_dim) {

  const int num_elements = leftArg.dimension(0);
  const int num_shape_fields = leftArg.dimension(1);

  FieldContainer<double> massMatrixAB(num_elements, num_shape_fields, num_shape_fields);

  // Assemble total element mass matrix
  for(int idx_cell=0; idx_cell<num_elements; ++idx_cell) {
    integrateMassMatrixAB(massMatrixAB, leftArg, rightArg, rho); 

    for(int a=0; a<num_shape_fields; ++a) {
      for(int b=0; b<num_shape_fields; ++b) {
        for(int i=0; i<space_dim; ++i) {
          const int idx_row = a*space_dim + i;
          const int idx_col = b*space_dim + i;
          massMatrix(idx_cell, idx_row, idx_col) = massMatrixAB(idx_cell, a, b);
        }
      }
    }
  }

  return;
}

