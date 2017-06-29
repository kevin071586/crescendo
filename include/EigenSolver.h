#ifndef EIGENSOLVER_HPP
#define EIGENSOLVER_HPP

#include <Epetra_FECrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

#include "stk_util/parallel/Parallel.hpp"

#include "ParserCmdBlock.h"

class EigenSolver {
  public:
    // constructor & destructor
    EigenSolver(stk::ParallelMachine stkComm);
    ~EigenSolver() {};

    int Solve(Epetra_FECrsMatrix& K, Epetra_FECrsMatrix& M);

    int SolveIfpack(Epetra_FECrsMatrix& K, Epetra_FECrsMatrix& M);

    void setParams(const ParserCmdBlock& p) {m_solverParams = p;}

    std::vector<double>* m_eigen_values;
    Teuchos::RCP<Epetra_MultiVector> m_eigen_vectors;

  private:
    // Parser m_parserData; TODO: need to include header for this type
    stk::ParallelMachine m_stkComm;
    ParserCmdBlock m_solverParams;
};

#endif // EIGENSOLVER_HPP
