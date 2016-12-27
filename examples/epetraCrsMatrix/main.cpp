#include <iostream>

#include <stk_util/parallel/Parallel.hpp>

#include <Epetra_SerialDenseVector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_FECrsMatrix.h>

int main(int argc, char** argv) {
    // Set up a parallel communicator.  Using Stk library for this to 
    // exercise initializing Epetra communicator using the Stk communicator.
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    // Wrap the MPI communicator so that Epetra can use it
    Epetra_MpiComm epetra_comm(comm);

    // Basic Epetra serial vector example (because why not?)
    Epetra_SerialDenseVector my_vector(4);
    my_vector(0) = 10;
    my_vector(1) = 20;
    my_vector(2) = 30;
    my_vector(3) = 40;
    std::cout << my_vector;
    std::cout << "my_vector norm: " << my_vector.Norm2() << std::endl;

    // Construct a Map with num_elements and index base of 0
    int num_elements = 10;
    Epetra_Map Map(num_elements, 0, epetra_comm);

    // Create x and b vectors
    Epetra_Vector x(Map);
    Epetra_Vector b(Map);
    b.Random();
    x.Update(2.0, b, 0.0); // x = 2*b

    // ========================================================================
    //
    // Epetra matrix example: CrsMatrix vs. FECrsMatrix
    //
    // ========================================================================

    // Create map and set matrix dimensions.
    // NOTE: must create a second map if a rectangular matrix is needed, and 
    // use it with an alternative constructor.
    int num_rows = 3;
    int num_nonzero_cols = num_rows;
    Epetra_Map row_map(num_rows, 0, epetra_comm);
    std::cout << "Matrix map output:\n" << row_map << std::endl;

    // Demonstrating the difference between Epetra_CrsMatrix and Epetra_FECrsMatrix
    Epetra_CrsMatrix A(Copy, row_map, num_nonzero_cols);
    Epetra_FECrsMatrix B(Copy, row_map, num_nonzero_cols);
    
    // Fill matrix entries with value (i+1)*(j+1)
    for (int i = 0; i<num_rows; ++i ){
      for( int j = 0; j<num_nonzero_cols; ++j ){
        int num_entries = 1;
        double value = (i+1)*(j+1);

        // Epetra_CrsMatrix: values are not summed even if values inserted by 
        // multiple MPI processes.
        A.InsertGlobalValues(i, num_entries, &value, &j);

        // Epetra_FECrsMatrix: values get summed when inserted by multiple 
        // processes (even non-owning processes).
        B.InsertGlobalValues(i, num_entries, &value, &j);
      }
    } 

    // Demonstrate issues to understand when inserting additional entries when
    // using multiple MPI processes.  Only execute this section of code if the
    // calling process has ID == 1.
    if( epetra_comm.MyPID() == 1 ){
      int num_entries = 1;
      int i, j;
      double val;

      // CASE 1: Inserting an additional value into a row that is owned by the calling
      // process results in values being summed.
      i = 1;
      j = 1;
      val = 1000.0;
      A.InsertGlobalValues(i, num_entries, &val, &j);
      B.InsertGlobalValues(i, num_entries, &val, &j);

      // CASE 2: Inserting an additional value into a row that is NOT owned by the 
      // calling process does nothing for Epetra_CrsMatrix, but results in summing the
      // value for the Epetra_FECrsMatrix.
      i = 2;
      j = 1;
      val = 500.0;
      A.InsertGlobalValues(i, num_entries, &val, &j);
      B.InsertGlobalValues(i, num_entries, &val, &j);
    }
    
    // Tell the sparse matrix that we are done adding entries to it.
    int gblerr = 0;
    gblerr += A.FillComplete();
    gblerr += B.GlobalAssemble();   // calls FillComplete() by default, so no additional call needed.

    std::cout << A << std::endl;
    std::cout << B << std::endl;
    std::cout << "Global error code: " << gblerr << std::endl;

    // NOTE:
    //
    // 1) Final values should be (i+1)*(j+1) for CrsMatrix, except entry [1,1] which is modified.
    //      1, 2,    3
    //      2, 1004, 6
    //      3, 6,    9
    //
    // 2) Final values should be NumProcessors*(i+1)*(j+1) for FECrsMatrix since the matrix
    // filling code is called by every process.  Also, entries (1,1) and (2,1) will have been
    // modified by +1000 and +500, respectively.  Example output for 3 MPI processes:
    //      3, 6,     9
    //      6, 1012,  18
    //      9, 518,   27
    //

    // Call finalize for parallel.  MPI gets very angry without this.
    stk::parallel_machine_finalize();
    return 0;
}
