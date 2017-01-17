#ifndef SIMULATION_H 
#define SIMULATION_H

#include "stk_util/parallel/Parallel.hpp"
#include <stk_io/StkMeshIoBroker.hpp>           // for StkMeshIoBroker

#include "Parser.h"

class Simulation
{
  public:
    Simulation(Parser& parserData, stk::ParallelMachine stkComm);
    ~Simulation();

    void Execute();

    // Set methods
    void setSpatialDim(int dim) { m_spatialDim = dim; return; }

    // Get methods
    int getSpatialDim() const {return m_spatialDim;}

  private:
    Parser m_parserData;

    int m_spatialDim;
    stk::ParallelMachine m_stkComm;
    stk::io::StkMeshIoBroker m_ioBroker;

    void initializeInputMesh();
    void setResultsOutput(size_t& resultOutputHandle);
};

#endif // SIMULATION_H
