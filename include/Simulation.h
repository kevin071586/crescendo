#ifndef SIMULATION_H 
#define SIMULATION_H

#include "stk_util/parallel/Parallel.hpp"
#include "stk_io/StkMeshIoBroker.hpp"           // for StkMeshIoBroker
#include "Epetra_Map.h"

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
    // Epetra_Map m_epetraRowMap;

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

    // Setup row map for distributed Epetra vectors/matrices/operators
    Epetra_Map setupEpetraRowMap();

    // Local-to-Global and Global-to-Local DOF Maps
    size_t localIdToLocalDof(size_t localId, int dofNum);
    size_t localDofToLocalId(size_t localId, int dofNum);
    size_t globalIdToGlobalDof(size_t globalId, int dofNum);
    size_t globalDofToGlobalId(size_t globalDof, int dofNum);

    // Element bucket loop

};

#endif // SIMULATION_H
