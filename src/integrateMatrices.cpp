#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>
#include <integrateMatrices.hpp>

using namespace Intrepid;
typedef FunctionSpaceTools fst;

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
void integrateStiffnessMatrixAB(FieldContainer<double>& stiffnessMatrixAB, 
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

