#ifndef EIGENSOLVER_HPP
#define EIGENSOLVER_HPP

#include "AnasaziBasicEigenproblem.hpp"

#include <Epetra_FECrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

#include "stk_util/environment/OutputLog.hpp"
#include "stk_util/parallel/Parallel.hpp"

#include "ParserCmdBlock.h"

class EigenSolver {
  public:
    // constructor & destructor
    EigenSolver(stk::ParallelMachine stkComm);
    ~EigenSolver() {};

    std::vector<double>* m_eigen_values;

    Teuchos::RCP<Epetra_MultiVector> m_eigen_vectors;

    int Solve(Epetra_FECrsMatrix& K, Epetra_FECrsMatrix& M);

    int SolveIfpack(Epetra_FECrsMatrix& K, Epetra_FECrsMatrix& M);

    void setParams(const ParserCmdBlock& p) {m_solverParams = p;}

  private:
    // Parser m_parserData; TODO: need to include header for this type
    stk::ParallelMachine m_stkComm;

    ParserCmdBlock m_solverParams;

    std::ostream& outputP0 = *(stk::get_log_ostream("output"));

    void computePrintResiduals(const Epetra_FECrsMatrix& K, 
            const Epetra_FECrsMatrix& M,
            Anasazi::Eigensolution<double, Epetra_MultiVector> sol,
            Anasazi::ReturnType returnCode,
            std::vector<double> compEvals);
};

#endif // EIGENSOLVER_HPP
