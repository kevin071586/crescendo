#include <iostream>

//#include <gtest/gtest.h>                      // for AssertHelper, EXPECT_NE, etc
//#include <stddef.h>                           // for size_t, NULL
//#include <unistd.h>                           // for unlink
#include <string>                               // for string
#include <fstream>

#include "stk_mesh/base/Bucket.hpp"             // for Bucket
#include <stk_mesh/base/Entity.hpp>             // for EntityRank
#include <stk_mesh/base/GetEntities.hpp>        // for get_entities
#include "stk_mesh/base/Field.hpp"              // for Field
#include <stk_mesh/base/FieldBase.hpp>          // for field_data, FieldBase
#include <stk_mesh/base/MetaData.hpp>           // for MetaData
#include "stk_mesh/base/Types.hpp"              // for BucketVector
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian type, I think.
#include "stk_io/DatabasePurpose.hpp"           // for READ_MESH, WRITE_RESULTS, etc
#include <stk_io/StkMeshIoBroker.hpp>           // for StkMeshIoBroker
#include <stk_topology/topology.hpp>            // for topology::rank_t, etc.
#include <stk_util/parallel/Parallel.hpp>

#include <Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

// #include <Epetra_SerialDenseVector.h>
// #include <Epetra_Map.h>
// #include <Epetra_MpiComm.h>
// #include <Epetra_Vector.h>
// #include <Epetra_CrsMatrix.h>
 
#include <Intrepid_CellTools.hpp>
#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Shards_CellTopology.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

using namespace Intrepid;
typedef FunctionSpaceTools fst;

// Helper functions
FieldContainer<double> ElasticStiffnessTensor(double E, double nu);

void integrateStiffnessMatrixAB(FieldContainer<double>& stiffnessMatrix, 
                               FieldContainer<double> leftArg, 
                               FieldContainer<double> rightArg, 
                               FieldContainer<double> C, int i, int k);

void integrateMassMatrixAB(FieldContainer<double>& massMatrixAB, 
                               FieldContainer<double> leftArg, 
                               FieldContainer<double> rightArg,
                               const double rho);

void integrateStiffnessMatrix(FieldContainer<double>& stiffnessMatrix, 
                              FieldContainer<double> leftArg, 
                              FieldContainer<double> rightArg, 
                              FieldContainer<double> C);

void integrateMassMatrix(FieldContainer<double>& massMatrix, 
                              FieldContainer<double> leftArg,
                              FieldContainer<double> rightArg,
                              const double rho,
                              const unsigned space_dim); 

namespace utility {
  void writeElementMatrix(std::string fname, FieldContainer<double> inMatrix,
                          const int idx_elem);
}

//  
// Main Program
//
int main(int argc, char** argv) {
  // --------------------------------------------------------------------------
  //
  // Define model parameters
  //
  // --------------------------------------------------------------------------
  // Get elastic stiffness tensor
  double rho = 3.0;
  double E = 10.0;
  double nu = 0.0;

  // Define the input file
  std::string dbtype("exodusII");
  std::string in_filename("two_hex8_elements.g");
  //std::string in_filename("hex8_10x10x10.g");

  // --------------------------------------------------------------------------
  //
  // Initialize MPI, Read in STK Exodus Mesh 
  //
  // --------------------------------------------------------------------------
  // Set up a parallel communicator
  stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);
  stk::io::StkMeshIoBroker stkMeshIoBroker(comm);

  // Create the input mesh (read metadata only)
  stkMeshIoBroker.add_mesh_database(in_filename, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();

  // Populate/read in the bulk mesh data
  stkMeshIoBroker.populate_bulk_data(); 

  // Get bulk data object for the mesh
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  // Select local processor elements 
  stk::mesh::Selector localElementsSelector = stkMeshMetaData.locally_owned_part();
  const stk::mesh::BucketVector elementBuckets =
    stkMeshBulkData.get_buckets(stk::topology::ELEMENT_RANK, localElementsSelector);

  // Get the coordinates field  
  typedef stk::mesh::Field<double, stk::mesh::Cartesian> CoordinatesField_t;
  CoordinatesField_t const& coord_field = 
      static_cast<CoordinatesField_t const&>(stkMeshIoBroker.get_coordinate_field());

  // Loop over all elements and print out some information
  // TODO: get rid of these or replace with something better.
  const unsigned nodesPerHex = 8;
  const unsigned spatialDim = 3;

  // --------------------------------------------------------------------------
  //
  //  Cubature and Element Basis
  //
  // --------------------------------------------------------------------------
  // Define cell topology of the parent cell
  shards::CellTopology hexahedron_8( shards::getCellTopologyData<shards::Hexahedron<8>>() );

  // Create cubature
  Intrepid::DefaultCubatureFactory<double> cubFactory;
  Teuchos::RCP<Intrepid::Cubature<double>> hexCubature = cubFactory.create(hexahedron_8, 2);

  // Define basis
  Intrepid::Basis_HGRAD_HEX_C1_FEM<double, Intrepid::FieldContainer<double>> hexHGradBasis;
  int numFields = hexHGradBasis.getCardinality();
  int numCubPoints = hexCubature->getNumPoints();

  // Field containers for cubature
  int cubDim = hexCubature->getDimension();
  FieldContainer<double> cubWeights(numCubPoints);
  FieldContainer<double> cubPoints(numCubPoints, cubDim);
  FieldContainer<double> hexGVals(numFields, numCubPoints);
  FieldContainer<double> hexGGradient(numFields, numCubPoints, spatialDim ); 


  // Get cubature points and weights
  hexCubature->getCubature(cubPoints, cubWeights);
  hexHGradBasis.getValues(hexGVals, cubPoints, Intrepid::OPERATOR_VALUE);

  // get gradient values for hex8 element (note: number of values here is
  // equal to cardinality of basis -- no need for # of cells).
  hexHGradBasis.getValues(hexGGradient, cubPoints, Intrepid::OPERATOR_GRAD);


  // --------------------------------------------------------------------------
  //
  // Loop over element buckets
  //
  // --------------------------------------------------------------------------
  for (size_t bucketIndex = 0; bucketIndex < elementBuckets.size(); ++bucketIndex) {

    stk::mesh::Bucket &elemBucket = *elementBuckets[bucketIndex];

    // 
    // Loop over elements in each bucket
    //
    unsigned num_elements = elemBucket.size();
    std::cout << "num_elements: " << num_elements << std::endl;
    FieldContainer<double> hexNodes(num_elements, nodesPerHex, spatialDim);

    // TODO: make a getNodeCoordinates function for this loop.
    for (size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex) {
      stk::mesh::Entity elem = elemBucket[elemIndex];

      unsigned numNodes = stkMeshBulkData.num_nodes(elem);
      stk::mesh::Entity const* nodes = stkMeshBulkData.begin_nodes(elem);

      //
      // Loop over each node in the current element
      //
      for (unsigned inode = 0; inode < numNodes; ++inode) {
        double *coords = stk::mesh::field_data(coord_field, nodes[inode]);
        
        hexNodes(elemIndex, inode, 0) = coords[0];
        hexNodes(elemIndex, inode, 1) = coords[1];
        hexNodes(elemIndex, inode, 2) = coords[2];
      }
    }

    // Compute cell Jacobians
    FieldContainer<double> hexJacobian(num_elements , numCubPoints, spatialDim, spatialDim);
    FieldContainer<double> hexJacobianInv(num_elements, numCubPoints, spatialDim, spatialDim);
    FieldContainer<double> hexJacobianDet(num_elements, numCubPoints);
    FieldContainer<double> cellMeasure(num_elements, numCubPoints);
    FieldContainer<double> hexGValsTransformed(num_elements, numFields, numCubPoints);
    FieldContainer<double> hexGValsTransformedWeighted(num_elements, numFields, numCubPoints);

    CellTools<double>::setJacobian(hexJacobian, cubPoints, hexNodes, hexahedron_8);
    CellTools<double>::setJacobianInv(hexJacobianInv, hexJacobian);
    CellTools<double>::setJacobianDet(hexJacobianDet, hexJacobian);

    // simply replicates input hexGVals for all cells in the set.
    fst::HGRADtransformVALUE<double>(hexGValsTransformed, hexGVals);

    // compute and multiply cell measure
    fst::computeCellMeasure<double>(cellMeasure, hexJacobianDet, cubWeights );
    fst::multiplyMeasure<double>(hexGValsTransformedWeighted, cellMeasure, hexGValsTransformed);


    // -------------------------------------------------------------------------
    // Compute mass matrix
    // -------------------------------------------------------------------------
    const int num_elem_DOF = spatialDim*numFields;
    FieldContainer<double> massMatrix(num_elements, num_elem_DOF, num_elem_DOF);
    integrateMassMatrix(massMatrix, hexGValsTransformed, hexGValsTransformedWeighted,
        rho,spatialDim);

    // quick mass check
    double total_mass = 0.0;
    for( int i=0; i<massMatrix.dimension(1); ++i){
      for( int j=0; j<massMatrix.dimension(2); ++j){
        total_mass += massMatrix(0,i,j);
      }
    }
    std::cout << "Sum M(i,j) = " << total_mass << std::endl;


    // -------------------------------------------------------------------------
    // Compute stiffness matrix
    // -------------------------------------------------------------------------
    FieldContainer<double> hexGGradientTransformed(num_elements, numFields, 
        numCubPoints, spatialDim ); 
    FieldContainer<double> hexGGradientTransformedWeighted(num_elements, numFields, 
        numCubPoints, spatialDim ); 

    // transform to physical coordinates.  Each cell has its own physical coordinates,
    // so must have an array sized to accomodate each coordinate.
    fst::HGRADtransformGRAD<double>(hexGGradientTransformed, hexJacobianInv, hexGGradient);

    // multiply measure
    fst::multiplyMeasure<double>(hexGGradientTransformedWeighted, cellMeasure, 
        hexGGradientTransformed);

    // integrate element stiffness matrices
    FieldContainer<double> stiffnessMatrix(num_elements, num_elem_DOF, num_elem_DOF);
    FieldContainer<double> C = ElasticStiffnessTensor(E, nu);
    integrateStiffnessMatrix(stiffnessMatrix, hexGGradientTransformed,
      hexGGradientTransformedWeighted, C);

    // ------------------------------
    // Write out matrices
    // ------------------------------
    int world_rank = stk::parallel_machine_rank(comm);
    std::cout << "MPI rank: " << world_rank << std::endl;
    
    int idx_elem = 0;
    utility::writeElementMatrix("K.csv", stiffnessMatrix, idx_elem);
    utility::writeElementMatrix("M.csv", massMatrix, idx_elem);
  }

  // Call finalize for parallel (MPI prints an angry error message without this call!!)
  stk::parallel_machine_finalize();
  return 0;
}


//-------------------------------------------------------------------------------- 
// Write out matrices
//-------------------------------------------------------------------------------- 
void utility::writeElementMatrix(std::string fname, FieldContainer<double> inMatrix,
  const int idx_elem) {

  std::ofstream file_handle(fname);

  if( file_handle.is_open() ){
    for( int i=0; i<inMatrix.dimension(1); ++i ){
      for( int j=0; j<inMatrix.dimension(2); ++j ){
        file_handle << inMatrix(idx_elem,i,j);
        
        // new line and comma separator characters
        if( j==(inMatrix.dimension(2)-1) ){
          file_handle << "\n";
        }
        else{
          file_handle << ", ";
        }
      }
    }
    file_handle.close();
  }
  else {
    std::cout << "Unable to open matrix output file." << std::endl;
  }
  return;
}

//-------------------------------------------------------------------------------- 
// Function to return elastic stiffness tensor
//-------------------------------------------------------------------------------- 
FieldContainer<double> ElasticStiffnessTensor(double E, double nu) {
  unsigned spatial_dim = 3;
  FieldContainer<double> C(3,3,3,3);

  // Calculate lame parameters
  double lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
  double mu = E/(2.0*(1.0+nu));

  for( size_t i=0; i<spatial_dim; ++i ){
    for( size_t j=0; j<spatial_dim; ++j ){
      for( size_t k=0; k<spatial_dim; ++k ){
        for( size_t l=0; l<spatial_dim; ++l ){
          double dij = (i==j) ? 1.0 : 0.0;
          double dkl = (k==l) ? 1.0 : 0.0;
          double dik = (i==k) ? 1.0 : 0.0;
          double djl = (j==l) ? 1.0 : 0.0;
          double dil = (i==l) ? 1.0 : 0.0;
          double djk = (j==k) ? 1.0 : 0.0;

          // Definition of elastic stiffness tensor
          C(i,j,k,l) = lambda*dij*dkl + mu*(dik*djl + dil*djk);
        }
      }
    }
  }
  
  return C;
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
void integrateStiffnessMatrix(FieldContainer<double>& stiffnessMatrix, 
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
void integrateMassMatrixAB(FieldContainer<double>& massMatrixAB, 
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
void integrateMassMatrix(FieldContainer<double>& massMatrix, 
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
