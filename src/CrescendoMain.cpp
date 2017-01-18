#include <iostream>

#include "stk_util/parallel/Parallel.hpp"

#include "Parser.h"
#include "Simulation.h"

int main(int argc, char *argv[]) {
  // Get and process command line parameters
  // --------------------------------------------------------------------------
  if (argc == 1) {
    std::cout << "Error: Must supply an input file" << std::endl;
    return 1;
  }
  std::string inputDeck(argv[1]);
  std::cout << "Input file: " << inputDeck << std::endl;
  
  // Parse the input deck
  // --------------------------------------------------------------------------
  Parser parser(inputDeck);
  parser.parse();

  // Create a simulation
  const int spatialDim = 3;
  stk::ParallelMachine stkComm = stk::parallel_machine_init(&argc, &argv);
  Simulation simulation(parser, stkComm);

  try {
    simulation.setSpatialDim(spatialDim);
    simulation.Execute();
  }
  catch (const std::runtime_error& e) {
    stk::parallel_machine_finalize();
    return 0;
  }

  // Call finalize for parallel (MPI prints an angry error message without this call!!)
  stk::parallel_machine_finalize();
  return 0;
}
