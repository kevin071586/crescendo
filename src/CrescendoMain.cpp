#include <iostream>

#include "boost/program_options.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/BroadcastArg.hpp"
#include "stk_util/environment/OutputLog.hpp"
#include "stk_util/environment/ProgramOptions.hpp"

#include "Parser.h"
#include "Simulation.h"

namespace po = boost::program_options;

namespace {
  const size_t SUCCESS = 0;
  const size_t ERROR = 1;
}

int main(int argc, char *argv[]) {
  // Get a parallel communicator
  stk::ParallelMachine stkComm = stk::parallel_machine_init(&argc, &argv);

  // Broadcast argc and argv to all processors
  stk::BroadcastArg b_arg(stkComm, argc, argv);

  // Program option variables
  std::string inputDeck; 
  std::string logFile("crescendo.log");

  // Populate program options
  po::options_description desc("Program options");
  desc.add_options()
    ("help,h", "Print help messages")
    ("input,i", po::value<std::string>(&inputDeck)->required(), "Provide input file")
    ("output,o", po::value<std::string>(&logFile), "Output log file");

  stk::get_options_description().add(desc);

  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // handle help option
    if (vm.count("help")) {
      std::cout << "Command line options:" << std::endl << desc << std::endl;
      return ERROR;
    }
    
    po::notify(vm);
  }
  catch(po::error& e) {
    std::cerr << "Error: " << e.what() << std::endl << std::endl;
    std::cerr << desc << std::endl;
    return ERROR;
  }

  // Open log file(s)
  // --------------------------------------------------------------------------
  stk::register_log_ostream(std::cout, "cout"); // Make standard out available
  stk::register_log_ostream(std::cerr, "cerr"); // Make standard error available
  stk::register_ostream(sierra::out(), "out");
  stk::register_ostream(sierra::pout(), "pout");
  stk::register_ostream(sierra::dout(), "dout");
  stk::register_ostream(sierra::tout(), "tout");

  //stk::create_log_file("logfile", logFile);              // Create main log file
  stk::bind_output_streams("logfile=\""+logFile+"\" "    
                            "out>logfile+pout "         // Send output to the log file and to the per-processor stream
                            "pout>null "                // Throw per-processor output away
                            "dout>out");                // Send diagnostic output to the regular output stream

  stk::bind_output_streams("out>pout "                  // Send output to the per-processor stream
                            "pout>null "                // Throw per-processor output away
                            "dout>out");                // Send diagnostic output to the regular output stream

  
  // std::ostream& outputP0 = *(stk::get_log_ostream("output"));
  std::ostream& outputP0 = *(stk::get_log_ostream("logfile"));
  outputP0 << "This is a test" << std::endl;
  //sierra::pout() << "Sierra parallel output stream" << stk::parallel_machine_rank(stkComm) << std::endl;

  // Parse the input deck
  // --------------------------------------------------------------------------
  Parser parser(inputDeck);
  parser.parse();

  // Create a simulation
  const int spatialDim = 3;
  Simulation simulation(parser, stkComm);

  try {
    simulation.setSpatialDim(spatialDim);
    simulation.Execute();
  }
  catch (const std::runtime_error& e) {
    stk::parallel_machine_finalize();
    return SUCCESS;
  }

  // Unregister log files
  stk::unregister_ostream(sierra::out());
  stk::unregister_ostream(sierra::pout());
  stk::unregister_ostream(sierra::dout());
  stk::unregister_ostream(sierra::tout());
  stk::unregister_log_ostream(std::cout);
  stk::unregister_log_ostream(std::cerr);

  // Call finalize for parallel (MPI prints an angry error message without this call!!)
  stk::parallel_machine_finalize();
  return SUCCESS;
}
