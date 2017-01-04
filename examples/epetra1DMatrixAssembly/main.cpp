#include <iostream>

#include <stk_util/parallel/Parallel.hpp>

#include <Epetra_SerialDenseVector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_FECrsMatrix.h>


// ============================================================================
// Objective: demonstrate finite-element matrix assembly using the FECrsMatrix
// to show:
//  (1) how global ID numbering can be used, and 
//  (2) how parallel matrix assembly and communication works
//
// The element configuration is:
//
//  (serial) local node ID:        0                 1                 2
//  (serial) local element ID:     o--------0--------o--------1--------o
//  global node ID:                10                20                30
//
// Assume that 
//  processor 0 owns nodes 10 and 20
//  processor 1 owns nodes 30
//
// The stiffness matrix for each element is assumed to be:
//
//  K_el = EA/L*[1, -1; -1, 1]
//
// For simplicity, assume EA/L=1.0
//
// This example is meant to be run on two processors.
// ============================================================================
int main(int argc, char** argv) {
    // Set up a parallel communicator.  Using Stk library for this to 
    // exercise initializing Epetra communicator using the Stk communicator.
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    // Wrap the MPI communicator so that Epetra can use it
    Epetra_MpiComm epetra_comm(comm);

    // Construct an Epetra_Map
    // (1) ... with num_elements and index base of 0
    // (2) Note that the map doesn't refer to actual row numbers.  It only
    //     refers to global ID's (10, 20, 30.. etc.) for its construction.
    //     Thus, it seems as though there is no concept of "first row" in
    //     whatever object (matrix, vector) is constructed from this map.
    //
    int num_local_nodes = 0;
    std::vector<int> global_node_ids;
    if( epetra_comm.MyPID() == 0 ){
      num_local_nodes = 2;
      global_node_ids.reserve(num_local_nodes);
      global_node_ids[0] = 10;
      global_node_ids[1] = 20;
    }
    else if( epetra_comm.MyPID() == 1 ){
      num_local_nodes = 1;
      global_node_ids.reserve(num_local_nodes);
      global_node_ids[0] = 30;
    }
    Epetra_Map row_map(-1, num_local_nodes, &global_node_ids[0], 0, epetra_comm);
    std::cout << row_map << std::endl;

    // Create matrix dimensions
    int num_nonzero_cols = 3;
    Epetra_FECrsMatrix stiffness_matrix(Copy, row_map, num_nonzero_cols);
    
    // Hard code the local element stiffness matrix
    Epetra_SerialDenseMatrix local_stiffness_matrix(2,2);
    local_stiffness_matrix(0,0) =  1.0;
    local_stiffness_matrix(0,1) = -1.0;
    local_stiffness_matrix(1,1) =  1.0;
    local_stiffness_matrix(1,0) = -1.0;

    std::cout << "Local stiffness matrix: " <<  local_stiffness_matrix << std::endl;

    // Insert local stiffness matrices into the global stiffness matrix
    // InsertGlobalValues(rowID, numEntries, values, colIDs)
    if( epetra_comm.MyPID() == 0 ){
      std::vector<int> row_id = {10, 20};
      std::vector<int> col_id = {10, 20};
      stiffness_matrix.InsertGlobalValues(row_id[0], 1, &local_stiffness_matrix(0,0), &col_id[0]);
      stiffness_matrix.InsertGlobalValues(row_id[0], 1, &local_stiffness_matrix(0,1), &col_id[1]);
      stiffness_matrix.InsertGlobalValues(row_id[1], 1, &local_stiffness_matrix(1,0), &col_id[0]);
      stiffness_matrix.InsertGlobalValues(row_id[1], 1, &local_stiffness_matrix(1,1), &col_id[1]);
    }
    else if( epetra_comm.MyPID() == 1 ){
      std::vector<int> row_id = {20, 30};
      std::vector<int> col_id = {20, 30};
      stiffness_matrix.InsertGlobalValues(row_id[0], 1, &local_stiffness_matrix(0,0), &col_id[0]);
      stiffness_matrix.InsertGlobalValues(row_id[0], 1, &local_stiffness_matrix(0,1), &col_id[1]);
      stiffness_matrix.InsertGlobalValues(row_id[1], 1, &local_stiffness_matrix(1,0), &col_id[0]);
      stiffness_matrix.InsertGlobalValues(row_id[1], 1, &local_stiffness_matrix(1,1), &col_id[1]);
    }

    int global_err = 0;
    global_err += stiffness_matrix.GlobalAssemble();

    std::cout << "Global stiffness matrix: " << stiffness_matrix << std::endl;
    std::cout << "Global error code: " << global_err << std::endl;

    // Call finalize for parallel.  MPI gets very angry without this.
    stk::parallel_machine_finalize();
    return 0;
}
