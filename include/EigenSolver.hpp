#ifndef EIGENSOLVER_HPP
#define EIGENSOLVER_HPP

#include <Epetra_FECrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

class EigenSolver {
  public:
    int Solve(Epetra_FECrsMatrix& K, Epetra_FECrsMatrix& M, Epetra_MpiComm& comm);
    int SolveIfpack(Epetra_FECrsMatrix& K, Epetra_FECrsMatrix& M, Epetra_MpiComm& comm);

    std::vector<double>* m_eigen_values;
    Teuchos::RCP<Epetra_MultiVector> m_eigen_vectors;
};

#endif // EIGENSOLVER_HPP
