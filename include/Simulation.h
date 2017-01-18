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

    // Initialize input mesh: read input mesh to create metadata 
    void initializeInputMesh();

    // Setup results output: define some global and nodal fields to be made
    // available for the output mesh, and then make call to populate_bulk_data().
    // The resulting output exodus mesh handle is returned in the argument
    // that is passed to it.  Eventually fields should be indicated elsewhere,
    // and assigned algorithmically (e.g., I probably don't want EigenFrequency
    // as a global variable in a statics analysis.
    void setResultsOutput(size_t& resultOutputHandle);

    // Initialize HEX8 cubature - make assumption that we have a HEX8 mesh, and
    // then populates the fields of the parent element.  No information about 
    // the mesh is required to do this, other than the topology of the element
    // itself.
    void initializeHex8Cubature();

};

#endif // SIMULATION_H
