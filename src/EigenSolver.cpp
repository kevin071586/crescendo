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
#include "Teuchos_RefCountPtr.hpp"

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
  Teuchos::FancyOStream fancy_ostream(Teuchos::rcp(&outputP0, false));
  Anasazi::BasicOutputManager<double> printer(Anasazi::Errors, Teuchos::rcp(&fancy_ostream, false));
  printer.setFancyOStream(Teuchos::rcp(&fancy_ostream, false));

  // Set verbosity level
  bool verbose = true;    // Move this to input deck
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
    verbosity += Anasazi::IterationDetails;
    verbosity += Anasazi::OrthoDetails;
  }

  // printer.stream(Anasazi::Errors) << Anasazi::Anasazi_Version() << std::endl << std::endl;
  printer.print(Anasazi::Errors, Anasazi::Anasazi_Version());

  //Teuchos::RCP<Epetra_FECrsMatrix> K = Teuchos::rcp(const_cast<Epetra_FECrsMatrix*>(&Kmat), false);
  //Teuchos::RCP<Epetra_FECrsMatrix> M = Teuchos::rcp(const_cast<Epetra_FECrsMatrix*>(&Mmat), false);
  Teuchos::RCP<Epetra_FECrsMatrix> K = Teuchos::rcp(&Kmat, false);
  Teuchos::RCP<Epetra_FECrsMatrix> M = Teuchos::rcp(&Mmat, false);
  assert(K != Teuchos::null);
  assert(M != Teuchos::null);

  // Create the shifted system K - sigma * M.
  double sigma = shift;
  outputP0 << "\n\nShift set to sigma = " << shift << std::endl;
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
    outputP0 << "Constructing Preconditioner ...\n" << std::endl;
    //printer.stream(Anasazi::Errors) << "Constructing Incomplete Cholesky preconditioner..." << std::flush;
    printer.print(Anasazi::Errors, "Constructing Incomplete Cholesky preconditioner...\n");
    Ifpack precFactory;
    // additive-Schwartz incomplete Cholesky with thresholding; see IFPACK documentation
    //std::string precType = "IC stand-alone";
    std::string precType = "ILU";

    int overlapLevel = 0;
    outputP0 << "Calling Preconditioner Create() ...\n" << std::endl;
    // prec = Teuchos::rcp( precFactory.Create(precType,Kshift.get(),overlapLevel) );
    prec = Teuchos::rcp( precFactory.Create(precType,&*Kshift,overlapLevel) );
    assert(prec != Teuchos::null);
    // parameters for preconditioner
    Teuchos::ParameterList precParams;
    precParams.set("fact: drop tolerance",prec_dropTol);
    precParams.set("fact: level-of-fill",prec_lofill);
    precParams.set("schwarz: combine mode", "Add");
    IFPACK_CHK_ERR(prec->SetParameters(precParams));

    outputP0 << "Calling Preconditioner Initialize() ...\n" << std::endl;
    IFPACK_CHK_ERR(prec->Initialize());
    outputP0 << "Calling Preconditioner Compute() ...\n" << std::endl;
    IFPACK_CHK_ERR(prec->Compute());
    //
    // printer.stream(Anasazi::Errors)
    //  << " done." << std::endl;
    printer.print(Anasazi::Errors, "done.");
    // encapsulate this preconditioner into a IFPACKPrecOp class
    PrecOp = Teuchos::rcp( new Epetra_InvOperator(&*prec) );
  }

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->Map(), blockSize) );
  ivec->Random();

  // Create the eigenproblem & define symmetry
  outputP0 << "Create the eigenvalue problem ... \n";
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
  MyPL.set( "Use Locking", true );
  MyPL.set( "Locking Tolerance", tol/10. );
  MyPL.set( "Output Stream", Teuchos::rcp(&fancy_ostream, false));

  // Create the solver manager
  outputP0 << "Creating solver manager ..." << std::endl;
  Anasazi::BlockDavidsonSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);

  // Solve the problem
  outputP0 << "Solving problem ..." << std::endl;
  Anasazi::ReturnType returnCode = MySolverMan.solve();

  // print some precond info
  if (usePrec) {
    std::ostringstream os;
    os << *prec << std::endl;
    //printer.stream(Anasazi::FinalSummary) << *prec << std::endl;
    printer.print(Anasazi::FinalSummary, os.str());
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  std::cout << "getSolution()" << std::endl;
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::cout << "sol.Evals" << std::endl;
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  std::cout << "Teuchos::RCP" << std::endl;
  Teuchos::RCP<MV> evecs = sol.Evecs;

  // Undo shift transformation; computed eigenvalues are real
  int numev = sol.numVecs;
  std::cout << "Found numev eigenvectors: " << numev << std::endl;
  if (numev == 0) {
      std::cout << "Terminating without successful solution." << std::endl;
      return -1;
  }

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

#include "Epetra_LinearProblem.h"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AztecOO.h"
#include "AztecOO_Operator.h"

int EigenSolver::SolveBKS(Epetra_FECrsMatrix &Kmat, Epetra_FECrsMatrix &Mmat) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cerr;
  using std::cout;
  using std::endl;
  // Anasazi solvers have the following template parameters:
  //
  //   - Scalar: The type of dot product results.
  //   - MV: The type of (multi)vectors.
  //   - OP: The type of operators (functions from multivector to
  //     multivector).  A matrix (like Epetra_CrsMatrix) is an example
  //     of an operator; an Ifpack preconditioner is another example.
  //
  // Here, Scalar is double, MV is Epetra_MultiVector, and OP is
  // Epetra_Operator.
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, MV> MVT;

  const int MyPID = stk::parallel_machine_rank(m_stkComm);

  // Get the stiffness and mass matrices.
  //
  // rcp (T*, false) returns a nonowning (doesn't deallocate) RCP.
  RCP<Epetra_FECrsMatrix> K = Teuchos::rcp(&Kmat, false);
  RCP<Epetra_FECrsMatrix> M = Teuchos::rcp(&Mmat, false);

  //
  // Create linear solver for linear systems with K
  //
  // Anasazi uses shift and invert, with a "shift" of zero, to find
  // the eigenvalues of least magnitude.  In this example, we
  // implement the "invert" part of shift and invert by using an
  // AztecOO iterative linear solver with an Ifpack preconditioner.
  //

  //
  // Construct Ifpack preconditioner
  //

  // An Ifpack "factory" knows how to create Ifpack preconditioners.
  Ifpack Factory;

  // Use the Factory to create the preconditioner.  Factory.Create()
  // returns a raw Ifpack2_Preconditioner pointer.  The caller (that's
  // us!) is responsible for deallocation, so we use an owning RCP to
  // deallocate it automatically.
  //
  // The Factory needs a string name of the preconditioner, and the
  // overlap level.  (Almost) all Ifpack preconditioners are local to
  // each MPI process.  Ifpack automatically does an additive Schwarz
  // decomposition across processes.  The default overlap level for
  // additive Schwarz is zero, but you can use a different overlap
  // level here if you want.  The overlap level must be nonnegative,
  // and only matters if running with more than one MPI process.
  std::string PrecType = "ICT"; // incomplete Cholesky
  const int OverlapLevel = 0;
  // Create the preconditioner.
  RCP<Ifpack_Preconditioner> Prec =
      rcp (Factory.Create (PrecType, &*K, OverlapLevel));
  if (Prec.is_null ()) {
      throw std::logic_error ("Ifpack's factory returned a NULL preconditioner!");
  }

  //
  // Set Ifpack preconditioner parameters.
  //
  // Set parameters after creating the preconditioner.  Please refer
  // to Ifpack's documentation for a list of valid parameters.
  //
  Teuchos::ParameterList ifpackList;
  ifpackList.set ("fact: drop tolerance", 1e-4);
  ifpackList.set ("fact: ict level-of-fill", 0.0);
  // The "combine mode" describes how to combine contributions from
  // other MPI processes.  It only matters if the overlap level is
  // nonzero.  See Epetra_CombineMode.h for documentation of of the
  // various combine modes.
  ifpackList.set ("schwarz: combine mode", "Add");

  // Set the parameters.
  IFPACK_CHK_ERR(Prec->SetParameters(ifpackList));

  // Initialize the preconditioner.  This only looks at the structure
  // of the matrix, not the values.  Nevertheless, the matrix must
  // generally be fill complete at this point.  Compare Initialize()
  // to a symbolic factorization, and Compute() (see below) to a
  // numeric factorization.
  IFPACK_CHK_ERR(Prec->Initialize());

  // Compute the preconditioner.  This looks at both the structure and
  // the values of the matrix.  Compare Initialize() (see above) to a
  // symbolic factorization, and Compute() (see below) to a numeric
  // factorization.
  IFPACK_CHK_ERR(Prec->Compute());

  //
  // Set up AztecOO iterative solver for solving linear systems with K.
  //

  // Create Epetra linear problem class for solving linear systems
  // with K.  This implements the inverse operator for shift and
  // invert.
  Epetra_LinearProblem precProblem;
  // Tell the linear problem about the matrix K.  Epetra_LinearProblem
  // doesn't know about RCP, so we have to give it a raw pointer.
  precProblem.SetOperator (K.getRawPtr ());

  // Create AztecOO iterative solver for solving linear systems with K.
  AztecOO precSolver (precProblem);
//  // Tell the solver to use the Ifpack preconditioner we created above.
//  precSolver.SetPrecOperator (Prec.get ());
//  // Set AztecOO solver options.
//  precSolver.SetAztecOption (AZ_output, AZ_none); // Don't print output
//  precSolver.SetAztecOption (AZ_solver, AZ_cg); // Use CG
//
//  // Use the above AztecOO solver to create the AztecOO_Operator.
//  // This is the place where we tell the AztecOO solver the maximum
//  // number of iterations (here, we use the matrix dimension; in
//  // practice, you'll want a smaller number) and the convergence
//  // tolerance (here, 1e-12).
//  RCP<AztecOO_Operator> precOperator =
//      rcp (new AztecOO_Operator (&precSolver, K->NumGlobalRows (), 1e-12));

//  // Create an Operator that computes y = K^{-1} M x.
//  //
//  // This Operator object is the operator we give to Anasazi.  Thus,
//  // Anasazi just sees an operator that computes y = K^{-1} M x.  The
//  // matrix K got absorbed into precOperator via precProblem (the
//  // Epetra_LinearProblem object).  Later, when we set up the Anasazi
//  // eigensolver, we will need to tell it about M, so that it can
//  // orthogonalize basis vectors with respect to the inner product
//  // defined by M (since it is symmetric positive definite).
//  RCP<Anasazi::EpetraGenOp> Aop = rcp (new Anasazi::EpetraGenOp (precOperator, M));
//
//  //
//  // Set parameters for the block Krylov-Schur eigensolver
//  //
//
//  double tol = 1.0e-8;
//  int nev = 10;
//  int blockSize = 3;
//  int numBlocks = 3 * nev / blockSize;
//  int maxRestarts = 10;
//
//  // We're looking for the largest-magnitude eigenvalues of the
//  // _inverse_ operator, thus, the smallest-magnitude eigenvalues of
//  // the original operator.
//  std::string which = "LM";
//  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
//
//  // Create ParameterList to pass into eigensolver
//  Teuchos::ParameterList MyPL;
//  MyPL.set ("Verbosity", verbosity);
//  MyPL.set ("Which", which);
//  MyPL.set ("Block Size", blockSize);
//  MyPL.set ("Num Blocks", numBlocks);
//  MyPL.set ("Maximum Restarts", maxRestarts);
//  MyPL.set ("Convergence Tolerance", tol);
//
//  // Create an initial set of vectors to start the eigensolver.  Note:
//  // This needs to have the same number of columns as the block size.
//  RCP<MV> ivec = rcp (new MV (K->Map (), blockSize));
//  MVT::MvRandom (*ivec);
//
//  // This object holds all the stuff about your problem that Anasazi
//  // will see.
//  //
//  // Anasazi only needs M so that it can orthogonalize basis vectors
//  // with respect to the M inner product.  Wouldn't it be nice if
//  // Anasazi didn't require M in two different places?  Alas, this is
//  // not currently the case.
//  RCP<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem =
//      rcp (new Anasazi::BasicEigenproblem<double,MV,OP> (Aop, M, ivec));
//
//  // Tell the eigenproblem that the matrix pencil (K,M) is symmetric.
//  MyProblem->setHermitian (true);
//
//  // Set the number of eigenvalues requested
//  MyProblem->setNEV (nev);
//
//  // Tell the eigenproblem that you are finished passing it information.
//  bool boolret = MyProblem->setProblem();
//  if (boolret != true) {
//      if (MyPID == 0) {
//	  cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
//      }
//      return -1;
//  }
//
//  // Create the Block Krylov-Schur eigensolver.
//  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr (MyProblem, MyPL);
//
//  // Solve the eigenvalue problem.
//  //
//  // Note that creating the eigensolver is separate from solving it.
//  // After creating the eigensolver, you may call solve() multiple
//  // times with different parameters or initial vectors.  This lets
//  // you reuse intermediate state, like allocated basis vectors.
//  Anasazi::ReturnType returnCode = MySolverMgr.solve ();
//  if (returnCode != Anasazi::Converged && MyPID == 0) {
//      cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
//  }
//
//  // Get the eigenvalues and eigenvectors from the eigenproblem.
//  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution ();
//  // Anasazi returns eigenvalues as Anasazi::Value, so that if
//  // Anasazi's Scalar type is real-valued (as it is in this case), but
//  // some eigenvalues are complex, you can still access the
//  // eigenvalues correctly.  In this case, there are no complex
//  // eigenvalues, since the matrix pencil is symmetric.
//  std::vector<Anasazi::Value<double> > evals = sol.Evals;
//  Teuchos::RCP<MV> evecs = sol.Evecs;
//  int numev = sol.numVecs;
//
//  if (numev > 0) {
//      // Reconstruct the eigenvalues.  The ones that Anasazi gave back
//      // are the inverses of the original eigenvalues.  Reconstruct the
//      // eigenvectors too.
//      Teuchos::SerialDenseMatrix<int,double> dmatr(numev,numev);
//      MV tempvec (K->Map (), MVT::GetNumberVecs (*evecs));
//      K->Apply (*evecs, tempvec);
//      MVT::MvTransMv (1.0, tempvec, *evecs, dmatr);
//
//      if (MyPID == 0) {
//	  double compeval = 0.0;
//	  cout.setf (std::ios_base::right, std::ios_base::adjustfield);
//	  cout << "Actual Eigenvalues (obtained by Rayleigh quotient) : " << endl;
//	  cout << "------------------------------------------------------" << endl;
//	  cout << std::setw(16) << "Real Part"
//	      << std::setw(16) << "Rayleigh Error" << endl;
//	  cout << "------------------------------------------------------" << endl;
//	  for (int i = 0; i < numev; ++i) {
//	      compeval = dmatr(i,i);
//	      cout << std::setw(16) << compeval
//		  << std::setw(16)
//	      << std::fabs (compeval - 1.0/evals[i].realpart)
//	      << endl;
//	  }
//	  cout << "------------------------------------------------------" << endl;
//      }
//
//  }
  return -1;
//  return 0;
}

int EigenSolver::SolveNoOp(Epetra_FECrsMatrix &Kmat, Epetra_FECrsMatrix &Mmat) {
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  const int numEigVecs = 10;
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Kmat.Map(), numEigVecs) );
  ivec->Random();

  // Assign result as public member variable
  m_eigen_values->resize(numEigVecs);
  for (int i; i<m_eigen_values->size(); ++i) {
      (*m_eigen_values)[i] = 0;
  }
  m_eigen_vectors = ivec;

  std::cout << "done.. returning from SolveNoOp()" << std::endl;
  return -1;
}
