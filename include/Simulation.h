#ifndef SIMULATION_H 
#define SIMULATION_H

#include "stk_util/parallel/Parallel.hpp"
#include "stk_io/StkMeshIoBroker.hpp"           // for StkMeshIoBroker

#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"

#include "Intrepid_FieldContainer.hpp"

#include "Element.h"
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

    // Setup row map for distributed Epetra vectors/matrices/operators.  Also
    // providing an estimate of number of nonzero entries for constructing
    // FECrsMatrix efficiently
    Epetra_Map setupEpetraRowMap(int& numNonzeroEstimate);

    // Element Bucket Loop: iterate over element buckets and compute mass
    // and stiffness matrices, force vectors, jacobians, etc.  A lot of heavy
    // lifting is done in this loop.
    void elementBucketLoop(Element elementData, Epetra_FECrsMatrix& massMatrix,
        Epetra_FECrsMatrix& stiffMatrix);
    
    // Process all nodes on all elements in a bucket.  In particular, get
    // and return node coordinates.  Also define a global row ID for each
    // node degree of freedom.
    void processBucketNodes(const stk::mesh::Bucket& elemBucket,
                      Intrepid::FieldContainer<double>& nodeCoords,
                      Intrepid::FieldContainer<int>& nodeDofGlobalIds);

    // Local-to-Global and Global-to-Local DOF Maps (these are hard-coded
    // for 3 dimensions).  They are used to account for the fact that each 
    // node has multiple DOFs: u1 ==> (u1_x, u1_y, u1_z)
    size_t localIdToLocalDof(size_t localId, int dofNum);
    size_t localDofToLocalId(size_t localId, int dofNum);
    size_t globalIdToGlobalDof(size_t globalId, int dofNum);
    size_t globalDofToGlobalId(size_t globalDof, int dofNum);

};

#endif // SIMULATION_H
