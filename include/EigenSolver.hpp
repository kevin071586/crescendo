#ifndef EIGENSOLVER_HPP
#define EIGENSOLVER_HPP

#include <Epetra_FECrsMatrix.h>
#include <Epetra_MpiComm.h>

class EigenSolver {
  public:
    int Solve(Epetra_FECrsMatrix& K, Epetra_FECrsMatrix& M, Epetra_MpiComm comm);
};

#endif // EIGENSOLVER_HPP
