#include <iostream>

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "EpetraExt_MatrixMatrix.h"

#include <EigenSolver.h>

#include "stk_util/environment/OutputLog.hpp"

// Constructor
EigenSolver::EigenSolver(stk::ParallelMachine stkComm) {
  m_stkComm = stkComm;
}

int EigenSolver::Solve(Epetra_FECrsMatrix& Kmat, Epetra_FECrsMatrix& Mmat) {

  Epetra_MpiComm Comm(m_stkComm);

  std::ostream& outputP0 = *(stk::get_log_ostream("output"));
  outputP0 << "Solving Eigenproblem ..." << std::endl;

  // Create an Anasazi output manager
  Anasazi::BasicOutputManager<double> printer;
  printer.stream(Anasazi::Errors) << Anasazi::Anasazi_Version() << std::endl << std::endl;

  Teuchos::RCP<Epetra_FECrsMatrix> K = Teuchos::rcp(const_cast<Epetra_FECrsMatrix*>(&Kmat), false);
  Teuchos::RCP<Epetra_FECrsMatrix> M = Teuchos::rcp(const_cast<Epetra_FECrsMatrix*>(&Mmat), false);

  // Create the shifted system K - sigma * M.
  double sigma = 1.0;
  Teuchos::RCP<Epetra_CrsMatrix> Kshift = Teuchos::rcp( new Epetra_CrsMatrix( *K ) );

  int addErr = EpetraExt::MatrixMatrix::Add( *M, false, -sigma, *Kshift, 1.0 );
  if (addErr != 0) {
    printer.print(Anasazi::Errors,"EpetraExt::MatrixMatrix::Add returned with error.\n");
    return -1;
  }


  //************************************
  // Call the Block Davidson solver manager
  //***********************************
  //
  //  Variables used for the Block Davidson Method
  //
  //  NOTE: See https://trilinos.org/pipermail/trilinos-users/2008-September/000822.html
  //  for some insight into choosing parameters.
  //
  std::cout << "number of modes = " << m_solverParams.getFieldInt("number of modes") << std::endl;
  std::cout << "block size = " << m_solverParams.getFieldInt("block size") << std::endl;
  std::cout << "number of blocks = " << m_solverParams.getFieldInt("number of blocks") << std::endl;

  const int    nev         = m_solverParams.getFieldInt("number of modes");
  const int    blockSize   = m_solverParams.getFieldInt("block size");
  const int    numBlocks   = m_solverParams.getFieldInt("number of blocks");
  const int    maxRestarts = m_solverParams.getFieldInt("maximum restarts");
  const double tol         = m_solverParams.getFieldInt("target residual");

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();

  // Create the eigenproblem & define symmetry
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(Kshift, M, ivec) );

  MyProblem->setHermitian(true);

  MyProblem->setNEV( nev );

  // Done providing information to eigenproblem
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    printer.print(Anasazi::Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
    return -1;
  }

  // Create parameter list to pass into the solver manager
  Teuchos::ParameterList MyPL;
  MyPL.set( "Which", "SM" );  // "SM" or "LM" - smallest or largest eig vals.
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );

  // Create the solver manager
  Anasazi::BlockDavidsonSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);

  // Solve the problem
  Anasazi::ReturnType returnCode = MySolverMan.solve();

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

  // Compute residuals.
  std::vector<double> normR(sol.numVecs);
  if (sol.numVecs > 0) {
    Teuchos::SerialDenseMatrix<int,double> T(sol.numVecs, sol.numVecs);
    Epetra_MultiVector Kvec( K->OperatorDomainMap(), evecs->NumVectors() );
    Epetra_MultiVector Mvec( M->OperatorDomainMap(), evecs->NumVectors() );
    T.putScalar(0.0); 
    for (int i=0; i<sol.numVecs; i++) {
      T(i,i) = evals[i].realpart;
    }
    K->Apply( *evecs, Kvec );  
    M->Apply( *evecs, Mvec );  
    MVT::MvTimesMatAddMv( -1.0, Mvec, T, 1.0, Kvec );
    MVT::MvNorm( Kvec, normR );
  }

  // Print the results
  std::ostringstream os;
  os.setf(std::ios_base::right, std::ios_base::adjustfield);
  os<<"Solver manager returned " << (returnCode == Anasazi::Converged ? "converged." : "unconverged.") << std::endl;
  os<<std::endl;
  os<<"------------------------------------------------------"<<std::endl;
  os<<std::setw(16)<<"Eigenvalue"
    <<std::setw(18)<<"Direct Residual"
    <<std::endl;
  os<<"------------------------------------------------------"<<std::endl;
  for (int i=0; i<sol.numVecs; i++) {
    //os<<std::setw(16)<<evals[i].realpart
    os<<std::setw(16)<<compEvals[i]
      <<std::setw(18)<<normR[i]/evals[i].realpart
      <<std::endl;
  }
  os<<"------------------------------------------------------"<<std::endl;
  printer.print(Anasazi::Errors,os.str());

  return 0;
}

// Include header for Ifpack incomplete Cholesky preconditioner
#include "Ifpack.h"
#include "Epetra_InvOperator.h"

int EigenSolver::SolveIfpack(Epetra_FECrsMatrix& Kmat, Epetra_FECrsMatrix& Mmat) {
  using namespace Anasazi;

  Epetra_MpiComm Comm(m_stkComm);

  //************************************
  // Get the parameters from the command line
  //************************************
  int    nev       = 10;
  int    blockSize = 10;
  int    numBlocks   = 3; //4;
  int    maxRestarts = 100; 
  double tol       = 1.0e-8;
  bool verbose = false;    
  std::string which("SM");  
  bool usePrec = true;      
  double prec_dropTol = 1e-4;
  int prec_lofill = 0;

  //************************************
  // Create an Anasazi output manager
  //************************************
  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
  BasicOutputManager<double> printer(verbosity);
  printer.stream(Errors) << Anasazi_Version() << std::endl << std::endl;

  //************************************
  // Some useful typedefs
  //************************************
  //
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef MultiVecTraits<double, Epetra_MultiVector> MVT;

  // Mass and stiffness matrices
  Teuchos::RCP<Epetra_FECrsMatrix> K = Teuchos::rcp(const_cast<Epetra_FECrsMatrix*>(&Kmat), false);
  Teuchos::RCP<Epetra_FECrsMatrix> M = Teuchos::rcp(const_cast<Epetra_FECrsMatrix*>(&Mmat), false);

  // tell the user that we're done
  printer.stream(Errors) << " done." << std::endl;


  //************************************
  // Select the Preconditioner
  //************************************
  //
  Teuchos::RCP<Ifpack_Preconditioner> prec;
  Teuchos::RCP<Epetra_Operator> PrecOp;
  if (usePrec) {
    printer.stream(Errors) << "Constructing Incomplete Cholesky preconditioner..." << std::flush;
    Ifpack precFactory;
    // additive-Schwartz incomplete Cholesky with thresholding; see IFPACK documentation
    std::string precType = "IC stand-alone";
    int overlapLevel = 0;
    prec = Teuchos::rcp( precFactory.Create(precType,K.get(),overlapLevel) );
    // parameters for preconditioner
    Teuchos::ParameterList precParams;
    precParams.set("fact: drop tolerance",prec_dropTol);
    precParams.set("fact: level-of-fill",prec_lofill);
    IFPACK_CHK_ERR(prec->SetParameters(precParams));
    IFPACK_CHK_ERR(prec->Initialize());
    IFPACK_CHK_ERR(prec->Compute());
    //
    printer.stream(Errors)
      << " done." << std::endl;
    // encapsulate this preconditioner into a IFPACKPrecOp class
    PrecOp = Teuchos::rcp( new Epetra_InvOperator(&*prec) );
  }


  //************************************
  // Call the BlockDavidson solver manager
  //***********************************
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  //
  Teuchos::RCP<Epetra_MultiVector> ivec
    = Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();

  // Create the eigenproblem and define symmetry
  Teuchos::RCP<BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new BasicEigenproblem<double, MV, OP>(K, M, ivec) );

  MyProblem->setHermitian(true);

  if (usePrec) {
    MyProblem->setPrec(PrecOp);
  }

  MyProblem->setNEV( nev );

  // Done passing eigenproblem information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    printer.print(Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  // Parameter list for eigenproblem
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  
  // Create the solver manager
  BlockDavidsonSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);

  // Solve the problem
  printer.stream(Errors) << "Solving eigenvalue problem..." << std::endl;
  ReturnType returnCode = MySolverMan.solve();
  // print some precond info
  if (usePrec) {
    printer.stream(FinalSummary) << *prec << std::endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;

  // Compute the residuals
  std::vector<double> normR(sol.numVecs);
  if (sol.numVecs > 0) {
    Teuchos::SerialDenseMatrix<int,double> T(sol.numVecs, sol.numVecs);
    Epetra_MultiVector Kvec( K->OperatorDomainMap(), evecs->NumVectors() );
    Epetra_MultiVector Mvec( M->OperatorDomainMap(), evecs->NumVectors() );
    T.putScalar(0.0);
    for (int i=0; i<sol.numVecs; i++) {
      T(i,i) = evals[i].realpart;
    }
    K->Apply( *evecs, Kvec );
    M->Apply( *evecs, Mvec );
    MVT::MvTimesMatAddMv( -1.0, Mvec, T, 1.0, Kvec );
    MVT::MvNorm( Kvec, normR );
  }

  // Print the results
  std::ostringstream os;
  os.setf(std::ios_base::right, std::ios_base::adjustfield);
  os<<"Solver manager returned " << (returnCode == Converged ? "converged." : "unconverged.") << std::endl;
  os<<std::endl;
  os<<"------------------------------------------------------"<<std::endl;
  os<<std::setw(16)<<"Eigenvalue"
    <<std::setw(18)<<"Direct Residual"
    <<std::endl;
  os<<"------------------------------------------------------"<<std::endl;
  for (int i=0; i<sol.numVecs; i++) {
    os<<std::setw(16)<<evals[i].realpart
      <<std::setw(18)<<normR[i]/evals[i].realpart
      <<std::endl;
  }
  os<<"------------------------------------------------------"<<std::endl;
  printer.print(Errors,os.str());


  return 0;
}

