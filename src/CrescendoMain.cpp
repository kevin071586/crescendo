#include <iostream>

#include "boost/program_options.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/BroadcastArg.hpp"
#include "stk_util/environment/Env.hpp"
#include "stk_util/environment/OutputLog.hpp"
#include "stk_util/environment/ProgramOptions.hpp"
#include <stk_util/environment/ParsedOptions.hpp>
#include <stk_util/environment/ParseCommandLineArgs.hpp>

#include "Parser.h"
#include "Simulation.h"

namespace po = boost::program_options;

namespace {
  const size_t SUCCESS = 0;
  const size_t ERROR = 1;
}

int main(int argc, const char **argv) {
  char** m_argv = const_cast<char**>(argv);

  // Get a parallel communicator
  stk::ParallelMachine stkComm = stk::parallel_machine_init(&argc, &m_argv);

  // Broadcast argc and argv to all processors
  stk::BroadcastArg b_arg(stkComm, argc, m_argv);

  // Populate program options
  stk::OptionsSpecification optionsSpec;
  optionsSpec.add_options()
    ("help,h", false, false, "Print help messages")
    ("input,i", false, true, "Provide input file")
    ("output,o", "Output log file", stk::DefaultValue<std::string>("crescendo.log"));

  stk::get_options_specification().add(optionsSpec);
  stk::ParsedOptions parsedOptions;

  std::string inputDeck;
  std::string logFile;

  try {
    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);
    std::string inputDeck = parsedOptions["input"].as<std::string>();
    std::string logFile = parsedOptions["output"].as<std::string>();
  }
  catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl << std::endl;
    std::cerr << optionsSpec << std::endl;
    return ERROR;
  }

  // Open log file(s)
  // --------------------------------------------------------------------------
  // stk::register_log_ostream(std::cout, "cout"); // Make standard out available
  // stk::register_log_ostream(std::cerr, "cerr"); // Make standard error available
  stk::register_ostream(sierra::out(), "out");
  // stk::register_ostream(sierra::pout(), "pout");
  // stk::register_ostream(sierra::dout(), "dout");
  // stk::register_ostream(sierra::tout(), "tout");

  // stk::bind_output_streams("logfile=\""+logFile+"\" "    
  //                           "out>logfile+pout "         // Send output to the log file and to the per-processor stream
  //                           "pout>null "                // Throw per-processor output away
  //                           "dout>out");                // Send diagnostic output to the regular output stream

  // stk::bind_output_streams( "out>pout "                 // Send output to the per-processor stream
  //                         "pout>null "                // Throw per-processor output away
  //                         "dout>out");                // Send diagnostic output to the regular output stream

  // TODO: Testing to figure out how to control log output correctly ...
  // stk::bind_output_streams("pout>null "
  //                          "dout>null "
  //                          "tout>null ");
  if (stk::parallel_machine_rank(stkComm)==0) {
    stk::bind_output_streams("logfile=\""+logFile+"\" " 
                              "out>logfile ");
  }
  else {
    stk::bind_output_streams("out>null");
  }
  // *stk::get_ostream_ostream("out") << "Hello" << std::endl;

  // std::ostream& outputP0 = *(stk::get_log_ostream("logfile"));
  std::ostream& outputP0 = sierra::out();

  outputP0 << "CRESCENDO\n" << std::endl;
  // outputP0 << "Number of Processors: " << stk::parallel_machine_size(stkComm) << std::endl;
  // outputP0 << "Start Time: " << sierra::Env::startup_date() << std::endl;
  // outputP0 << sierra::Env::section_separator() << std::endl;

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
    std::cout << "CrescendoMain() :: Failed to execute.";
    return SUCCESS;
  }

  // Unregister log files
  stk::unregister_ostream(sierra::out());
  // stk::unregister_ostream(sierra::pout());
  // stk::unregister_ostream(sierra::dout());
  // stk::unregister_ostream(sierra::tout());
  // stk::unregister_log_ostream(std::cout);
  // stk::unregister_log_ostream(std::cerr);

  // Call finalize for parallel (MPI prints an angry error message without this call!!)
  stk::parallel_machine_finalize();
  return SUCCESS;
}
