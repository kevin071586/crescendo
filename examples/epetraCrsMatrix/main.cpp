#include <iostream>

#include <stk_util/parallel/Parallel.hpp>

#include <Epetra_SerialDenseVector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

int main(int argc, char** argv) {
    // Set up a parallel communicator
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    // 
    // Basic Epetra serial vector example
    //
    Epetra_SerialDenseVector my_vector(4);
    my_vector(0) = 10;
    my_vector(1) = 20;
    my_vector(2) = 30;
    my_vector(3) = 40;
    std::cout << "\n" << my_vector;
    std::cout << my_vector.Norm2() << std::endl;

    // Wrap the MPI communicator so that Epetra can use it
    Epetra_MpiComm epetra_comm(comm);

    // Construct a Map with NumElements and index base of 0
    int NumElements = 10;
    Epetra_Map Map(NumElements, 0, epetra_comm);

    // Create x and b vectors
    Epetra_Vector x(Map);
    Epetra_Vector b(Map);
  
    b.Random();
    x.Update(2.0, b, 0.0); // x = 2*b
    std::cout << b << std::endl;
    std::cout << x << std::endl;

    // 
    // Basic Epetra matrix
    //
    int matrix_dim = 3;
    Epetra_Map MatrixMap( matrix_dim, 0, epetra_comm );
    Epetra_CrsMatrix A(Copy, MatrixMap, matrix_dim);
    
    double value = 10.0;
    int colIdx = 0;
    int colIdx2 = 2;
    A.PutScalar(0.0);

    A.InsertGlobalValues(0, 1, &value, &colIdx);
    A.InsertGlobalValues(2, 1, &value, &colIdx2); // need to figure out how to cast without
                                                  // defining an entire variable for it.
    // Tell the sparse matrix that we are done adding entries to it.
    int gblerr = 0;
    gblerr = A.FillComplete ();

    std::cout << A << std::endl;
    std::cout << "Matrix map output:\n" << MatrixMap << std::endl;

    // Call finalize for parallel
    stk::parallel_machine_finalize();
    return 0;
}
