#include <iostream>
#include "stk_mesh/base/Bucket.hpp"             // for Bucket
#include <stk_mesh/base/Entity.hpp>             // for EntityRank
#include <stk_mesh/base/GetEntities.hpp>        // for get_entities
// #include "stk_mesh/base/Field.hpp"              // for Field
// #include <stk_mesh/base/FieldBase.hpp>          // for field_data, FieldBase
#include <stk_mesh/base/MetaData.hpp>           // for MetaData
// #include "stk_mesh/base/Types.hpp"              // for BucketVector
// #include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian type, I think.
#include "stk_io/DatabasePurpose.hpp"           // for READ_MESH, WRITE_RESULTS, etc
#include <stk_io/StkMeshIoBroker.hpp>           // for StkMeshIoBroker
#include <stk_topology/topology.hpp>            // for topology::rank_t, etc.
#include <stk_util/parallel/Parallel.hpp>

// ============================================================================
// Objective: Demonstrate and experiment with a few different options for 
// stk_mesh selectors.  I am primarily looking at how one could (should?) 
// select element buckets for matrix assembly.  Therefore, I need to be able
// to obtain information about the 8 nodes of a HEX8 element so that I can
// correctly assemble local stiffness matrices into the correct global matrix
// position.
//
// This example is meant to be run on two processors.
// ============================================================================

void printElementNodes(stk::mesh::Selector allEntities, stk::mesh::BulkData &bulk_data);

int main(int argc, char** argv) {
  // Set up a parallel communicator.  Using Stk library for this to 
  // exercise initializing Epetra communicator using the Stk communicator.
  stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);
  stk::io::StkMeshIoBroker stkMeshIoBroker(comm);

  // Define the input file
  std::string dbtype("exodusII");
  std::string in_filename("test_2elem.g");
  
  // Create the input mesh (read metadata only)
  stkMeshIoBroker.add_mesh_database(in_filename, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();

  // Populate/read in the bulk mesh data
  stkMeshIoBroker.populate_bulk_data(); 
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::Selector allEntities;
  
  // ==========================================================================
  //
  // Below are several different selector options.  Comment or uncomment as 
  // many as desired to see the output from a given selector
  //
  // ==========================================================================

  // std::cout << "SELECTOR #1 -- Universal Part" << std::endl;
  // allEntities = stkMeshMetaData.universal_part();
  // printElementNodes(allEntities, stkMeshBulkData);

  // std::cout << "SELECTOR #2 -- Universal Part AND NOT Aura Part" << std::endl;
  // allEntities = stkMeshMetaData.universal_part() & !stkMeshMetaData.aura_part();
  // printElementNodes(allEntities, stkMeshBulkData);
 
  std::cout << "SELECTOR #3 -- Locally Owned Part" << std::endl;
  allEntities = stkMeshMetaData.locally_owned_part();
  printElementNodes(allEntities, stkMeshBulkData);

  // std::cout << "SELECTOR #4 -- Globally Shared Part" << std::endl;
  // allEntities = stkMeshMetaData.globally_shared_part();
  // printElementNodes(allEntities, stkMeshBulkData);

  // std::cout << "SELECTOR #5 -- Locally Owned OR Globally Shared Part" << std::endl;
  // allEntities = stkMeshMetaData.locally_owned_part() | stkMeshMetaData.globally_shared_part();
  // printElementNodes(allEntities, stkMeshBulkData);

  // std::cout << "SELECTOR #6 -- Locally Owned AND NOT Globally Shared Part" << std::endl;
  // allEntities = stkMeshMetaData.locally_owned_part() & !stkMeshMetaData.globally_shared_part();
  // printElementNodes(allEntities, stkMeshBulkData);

  // Call finalize for parallel.  MPI gets very angry without this.
  stk::parallel_machine_finalize();
  return 0;
}


// ==========================================================================
// HELPER FUNCTION: Print out nodes of the currently selected elements
// ==========================================================================
void printElementNodes(stk::mesh::Selector allEntities, stk::mesh::BulkData &bulk_data) {

  std::vector<unsigned> entityCounts;
  stk::mesh::count_entities(allEntities, bulk_data, entityCounts);
  int num_nodes = entityCounts[stk::topology::NODE_RANK];
  int num_elements = entityCounts[stk::topology::ELEMENT_RANK];
  std::cout << "num_selected_elements: " << num_elements << std::endl;
  std::cout << "num_selected_nodes:    " << num_nodes << std::endl;

  if( num_elements == 0 ){
    std::cout << "Proc [" << bulk_data.parallel_rank() << "]: ";
    std::cout << "No elements are selected." << std::endl;
  }

  const stk::mesh::BucketVector buckets =
    bulk_data.get_buckets(stk::topology::ELEMENT_RANK, allEntities);
  std::cout << "number of buckets: " << buckets.size() << std::endl;

  for (size_t bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    stk::mesh::Bucket &elemBucket = *buckets[bucketIndex];
    unsigned num_elements = elemBucket.size();

    for (size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex) {
      stk::mesh::Entity elem = elemBucket[elemIndex];
      unsigned num_nodes_in_elem = bulk_data.num_nodes(elem);
      std::cout << "Proc [" << bulk_data.parallel_rank() << "]: ";
      std::cout << "number of nodes in element bucket: " << num_nodes_in_elem << std::endl;

      std::cout << "Proc [" << bulk_data.parallel_rank() << "]: ";
      std::cout << "Global Node IDs: "; 
      stk::mesh::Entity const* nodes = bulk_data.begin_nodes(elem);
      for (unsigned inode = 0; inode < num_nodes_in_elem; ++inode) {
        int local_id = bulk_data.identifier(nodes[inode]);
        std::cout << local_id; 
        if( inode != num_nodes_in_elem-1 ){ std::cout << ", "; }
      }
      std::cout << std::endl;

    }
  }
  return;
}
