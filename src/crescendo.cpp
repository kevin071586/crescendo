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

#include <integrateMatrices.hpp>
#include <ElasticMaterial.hpp>

using namespace Intrepid;
typedef FunctionSpaceTools fst;

// Helper functions
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

    ElasticMaterial elastic_mat;
    elastic_mat.setYoungsModulus(E);
    elastic_mat.setPoissonsRatio(nu);
    FieldContainer<double> C = elastic_mat.stiffnessTensor();

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

