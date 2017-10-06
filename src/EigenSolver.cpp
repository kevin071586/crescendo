#include <iostream>

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Epetra_FECrsMatrix.h"
#include "Epetra_Map.h"
#include "EpetraExt_MatrixMatrix.h"

#include "stk_util/environment/Env.hpp"                 // for section_title

// Include header for Ifpack incomplete Cholesky preconditioner
#include "Ifpack.h"
#include "Epetra_InvOperator.h"

#include <EigenSolver.h>



// Constructor
EigenSolver::EigenSolver(stk::ParallelMachine stkComm) {
  m_stkComm = stkComm;
}

int EigenSolver::Solve(Epetra_FECrsMatrix& Kmat, Epetra_FECrsMatrix& Mmat) {
  outputP0 << sierra::Env::section_title("Eigenproblem Solution") << "\n" << std::endl;

  //
  //  Variables used for the Block Davidson Method
  //
  //  NOTE: See https://trilinos.org/pipermail/trilinos-users/2008-September/000822.html
  //  for some insight into choosing parameters.
  //
  const int    nev         = m_solverParams.getFieldInt("number of modes");
  const int    blockSize   = m_solverParams.getFieldInt("block size");
  const int    numBlocks   = m_solverParams.getFieldInt("number of blocks");
  const int    maxRestarts = m_solverParams.getFieldInt("maximum restarts");
  const double tol         = m_solverParams.getFieldDouble("target residual");
  const double shift       = m_solverParams.getFieldDouble("shift");

  // Create an Anasazi output manager
  Anasazi::BasicOutputManager<double> printer;
  printer.setOStream(Teuchos::rcp(&outputP0, false));

  // Set verbosity level
  bool verbose = true;    // Move this to input deck
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
    verbosity += Anasazi::IterationDetails;
    verbosity += Anasazi::OrthoDetails;
  }

  printer.stream(Anasazi::Errors) << Anasazi::Anasazi_Version() << std::endl << std::endl;

  Teuchos::RCP<Epetra_FECrsMatrix> K = Teuchos::rcp(const_cast<Epetra_FECrsMatrix*>(&Kmat), false);
  Teuchos::RCP<Epetra_FECrsMatrix> M = Teuchos::rcp(const_cast<Epetra_FECrsMatrix*>(&Mmat), false);

  // Create the shifted system K - sigma * M.
  double sigma = shift;
  Teuchos::RCP<Epetra_FECrsMatrix> Kshift = Teuchos::rcp( new Epetra_FECrsMatrix( *K ) );

  int addErr = EpetraExt::MatrixMatrix::Add( *M, false, -sigma, *Kshift, 1.0 );
  if (addErr != 0) {
    printer.print(Anasazi::Errors,"EpetraExt::MatrixMatrix::Add returned with error.\n");
    return -1;
  }

  // Preconditioner parameters
  bool usePrec = true;      
  double prec_dropTol = 1e-9;
  int prec_lofill = 1;

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;

  //************************************
  // Select the Preconditioner
  //************************************
  //
  Teuchos::RCP<Ifpack_Preconditioner> prec;
  Teuchos::RCP<Epetra_Operator> PrecOp;
  if (usePrec) {
    printer.stream(Anasazi::Errors) << "Constructing Incomplete Cholesky preconditioner..." << std::flush;
    Ifpack precFactory;
    // additive-Schwartz incomplete Cholesky with thresholding; see IFPACK documentation
    //std::string precType = "IC stand-alone";
    std::string precType = "ILU";

    int overlapLevel = 1;
    prec = Teuchos::rcp( precFactory.Create(precType,Kshift.get(),overlapLevel) );
    // parameters for preconditioner
    Teuchos::ParameterList precParams;
    precParams.set("fact: drop tolerance",prec_dropTol);
    precParams.set("fact: level-of-fill",prec_lofill);
    IFPACK_CHK_ERR(prec->SetParameters(precParams));
    IFPACK_CHK_ERR(prec->Initialize());
    IFPACK_CHK_ERR(prec->Compute());
    //
    printer.stream(Anasazi::Errors)
      << " done." << std::endl;
    // encapsulate this preconditioner into a IFPACKPrecOp class
    PrecOp = Teuchos::rcp( new Epetra_InvOperator(&*prec) );
  }

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();

  // Create the eigenproblem & define symmetry
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(Kshift, M, ivec) );

  MyProblem->setHermitian(true);

  MyProblem->setNEV( nev );

  if (usePrec) {
    MyProblem->setPrec(PrecOp);
  }

  // Done providing information to eigenproblem
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    printer.print(Anasazi::Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
    return -1;
  }

  // Create parameter list to pass into the solver manager
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", "SM" );  // "SM" or "LM" - smallest or largest eig vals.
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Num Restart Blocks", 1);
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );

  // Create the solver manager
  Anasazi::BlockDavidsonSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);

  // Solve the problem
  Anasazi::ReturnType returnCode = MySolverMan.solve();

  // print some precond info
  if (usePrec) {
    printer.stream(Anasazi::FinalSummary) << *prec << std::endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;

  // Undo shift transformation; computed eigenvalues are real
  int numev = sol.numVecs;
  std::vector<double> compEvals(numev);
  for (int i=0; i<numev; ++i) {
    compEvals[i] = evals[i].realpart + sigma;
  }

  // Assign result as public member variable
  m_eigen_values = &compEvals;
  m_eigen_vectors = sol.Evecs;

  computePrintResiduals(Kmat, Mmat, sol, returnCode, compEvals);

  std::ostringstream os;
  printer.print(Anasazi::Errors, os.str());

  return 0;
}


//-------------------------------------------------------------------
// Compute and print residuals for an eigen problem.
//
// TODO: Simplify this interface.. passing too many parameters in 
// for what this function does, some of this should be stored in the
// EigenSolver class itself as member variables or something.
//
void EigenSolver::computePrintResiduals(const Epetra_FECrsMatrix& K, 
              const Epetra_FECrsMatrix& M,
              Anasazi::Eigensolution<double, Epetra_MultiVector> sol,
              Anasazi::ReturnType returnCode,
              std::vector<double> compEvals) {

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;

  // Get eigenvectors & eigenvalues from solution 
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  
  // Compute the residuals
  std::vector<double> normR(sol.numVecs);

  if (sol.numVecs > 0) {
    Teuchos::SerialDenseMatrix<int,double> T(sol.numVecs, sol.numVecs);
    Epetra_MultiVector Kvec(K.OperatorDomainMap(), evecs->NumVectors());
    Epetra_MultiVector Mvec(M.OperatorDomainMap(), evecs->NumVectors());
    T.putScalar(0.0);
    for (int i=0; i<sol.numVecs; i++) {
      T(i,i) = evals[i].realpart;
    }
    K.Apply( *evecs, Kvec );
    M.Apply( *evecs, Mvec );
    MVT::MvTimesMatAddMv( -1.0, Mvec, T, 1.0, Kvec );
    MVT::MvNorm( Kvec, normR );
  }

  // Print the results only from processor rank 0
  // TODO: is there a nicer way to do this without the "if proc == 0" statement?
  if (stk::parallel_machine_rank(m_stkComm)==0) {

    std::ostringstream oss;
    std::string separator("------------------------------------------------------");

    oss.setf(std::ios_base::right, std::ios_base::adjustfield);

    oss << "Solver manager returned " << (returnCode == Anasazi::Converged ? "converged." : "unconverged.") << std::endl;
    oss << std::endl;
    oss << separator << std::endl;
    oss << std::setw(16)<<"Eigenvalue"
        << std::setw(18)<<"Direct Residual"
        << std::endl;
    oss << separator << std::endl;

    for (int i=0; i<sol.numVecs; i++) {
      oss << std::setw(16) << compEvals[i]
               << std::setw(18) << normR[i]/evals[i].realpart
               << std::endl;
    }
    oss << separator << std::endl;
    //outputP0 << "rank: " << stk::parallel_machine_rank(m_stkComm) << std::endl;
    outputP0 << oss.str();
  }

  return;
}