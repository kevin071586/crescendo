#include <iostream>
#include <string>
#include <stk_util/parallel/Parallel.hpp>
#include "stk_util/environment/Env.hpp"
#include "stk_util/environment/OutputLog.hpp"

// ============================================================================
// Objective: 
// Demonstrate basic usage of STK Output Logging.
//
// This example is meant to be run on one or two processors.
// ============================================================================


int main(int argc, char** argv) {
  // Set up a parallel communicator.  Using Stk library for this to 
  // exercise initializing Epetra communicator using the Stk communicator.
  stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);


  int myRank = stk::parallel_machine_rank(comm);
  int numProc = stk::parallel_machine_size(comm);
  std::cout << "Rank: " << myRank << std::endl;

  // Define some output streams
  stk::register_ostream(sierra::out(), "out");
  stk::register_ostream(sierra::pout(), "pout");

  // Define output rules for a single log file
  // -----------------------------------------------------------------
  if (myRank == 0) { 
    stk::bind_output_streams("logfile=\"output.log\" "    
                            "out>logfile+pout "          // Send output to the log file and to the per-processor stream
                            "pout>null");                // Throw per-processor output away
  }
  else {
    stk::bind_output_streams("out>null");
  }

  // -----------------------------------------------------------------
  // Define output rules for parallel log file output
  // -----------------------------------------------------------------
  std::string pLogFile = "parallel_output."+std::to_string(numProc)+"."+std::to_string(myRank)+".log";
  stk::bind_output_streams("plogfile=\""+pLogFile+"\" "    
                            "pout>plogfile");


  // Test the output streams
  sierra::out() << "[Proc " << myRank << "] " << "This is the sierra::out() stream" << std::endl;
  sierra::pout() << "[Proc " << myRank << "] " << "This is the sierra::pout() stream" << std::endl;

  // Use an alias for the output stream
  // std::ostream& myOut = sierra::out();
  // myOut << "[Proc " << myRank << "] " << "This is the myOut stream, aliased to the sierra::out() stream" << std::endl;

  // Get the output stream by key name
  // std::ostream& myOut2 = *(stk::get_ostream_ostream("out"));
  // myOut2 << "[Proc " << myRank << "] " << "This is the myOut2 stream, from stk::get_ostream_ostream()" << std::endl;


  // Call finalize for parallel.  MPI gets very angry without this.
  stk::parallel_machine_finalize();
  return 0;
}


