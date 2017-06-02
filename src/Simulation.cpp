// ----------------------------------------------------------------------------
// TODO List:
//  1) Function to write log file to given stream, one proc output only.
//  2) Handle some constants better -- like m_spatialDim, NODES_PER_HEX, etc.
//     For example, NODES_PER_HEX can probably be evaluated based on element.
//  3) ElementLoop() needs to be broken up into smaller functions.  Some sections
//     can be moved into the Element calculation routine, like calculating the
//     jacobians.  Perhaps matrix calculation could be done there as well, provided
//     an input matrix to store data in?
//
//
// ----------------------------------------------------------------------------

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"

#include "stk_io/DatabasePurpose.hpp"           // for READ_MESH, WRITE_RESULTS, etc
#include "stk_io/StkMeshIoBroker.hpp"           // for StkMeshIoBroker
#include "stk_mesh/base/Field.hpp"              // for Field
#include <stk_mesh/base/GetEntities.hpp>        // for get_entities, count_entities
#include "stk_mesh/base/MetaData.hpp"           // for MetaData
#include "stk_mesh/base/CoordinateSystems.hpp"  // for Cartesian type

#include "Shards_CellTopology.hpp"              // For Hex8 topology
#include "Intrepid_CellTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"

#include "Teuchos_RCP.hpp"                      // For RCP

#include "Element.h"
#include "Simulation.h"
#include "ParserCmdBlock.h"
#include "ElasticMaterial.h"
#include "EigenSolver.h"


// Constructor and destructor
// ============================================================================
Simulation::Simulation(Parser& parserData, stk::ParallelMachine stkComm)
  : m_ioBroker(stkComm)
{
  m_parserData = parserData;
  m_stkComm = stkComm;
}

Simulation::~Simulation()
{
}

// Execute simulation 
// ============================================================================
void Simulation::Execute() 
{
  initializeInputMesh();

  // NOTE TO SELF: Currently working here.  I would like to store things like
  // number of nodes per hex, cubature degree, and element topology type information
  // in this class so that I can pass an object around with all that information,
  // and compute jacobians and things with that information.
  Element elementData;
  elementData.initializeCubature();

  size_t exoOutputMesh;
  setResultsOutput(exoOutputMesh);

  // Must come after populate_bulk_data is called
  int numNonzeroEstimate;
  Epetra_Map epetraRowMap = setupEpetraRowMap(numNonzeroEstimate);
  Epetra_FECrsMatrix stiffMatrix(Copy, epetraRowMap, numNonzeroEstimate);
  Epetra_FECrsMatrix massMatrix(Copy, epetraRowMap, numNonzeroEstimate);

  ParserCmdBlock cmdBlock = m_parserData.getCmdBlock("finite element model");
  const double rho = cmdBlock.getFieldDouble("density");
  const double E = cmdBlock.getFieldDouble("youngs modulus");
  const double nu = cmdBlock.getFieldDouble("poissons ratio");
  std::cout << "Density: " << rho << std::endl;
  std::cout << "Youngs Modulus: " << E << std::endl;
  std::cout << "Poissons Ratio: " << nu << std::endl;

  elementBucketLoop(elementData, massMatrix, stiffMatrix);

  // Solve eigen problem
  EigenSolver eigSolver;
  
  // Wrap the MPI communicator so Epetra can use it
  Epetra_MpiComm epetraComm(m_stkComm);
  eigSolver.Solve(stiffMatrix, massMatrix, epetraComm); // TODO: working here, not well tested yet.

  // Eigenvalue problem post-processing


  return; 
}

// Initialize input mesh metadata and bulk data
// ============================================================================
void Simulation::initializeInputMesh() 
{
  ParserCmdBlock cmdBlock = m_parserData.getCmdBlock("finite element model");
  std::string inputDatabase = cmdBlock.getFieldString("database name");

  // Read and create input deck
  try {
    m_ioBroker.add_mesh_database(inputDatabase, stk::io::READ_MESH);
    m_ioBroker.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = m_ioBroker.meta_data();
  }
  catch (const std::runtime_error& e) {
    if (stk::parallel_machine_rank(m_stkComm) == 0) {
      std::cout << "Error: Mesh '" << inputDatabase << "' not found.\n";
    }
    throw;
  }

  return; 
}

// Setup results output -- output variables names and mesh
// ============================================================================
void Simulation::setResultsOutput(size_t& resultOutputHandle)
{
  // Define some types for stk mesh output fields
  typedef stk::mesh::Field<double> ScalarField;
  typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> VectorField;

  // Get the metadata
  stk::mesh::MetaData& stkMeshMetaData = m_ioBroker.meta_data();

  // Get database name for output
  ParserCmdBlock cmdBlock = m_parserData.getCmdBlock("results output");
  std::string outputDatabase = cmdBlock.getFieldString("database name");

  // Define the output mesh
  std::string out_filename("output.e");
  resultOutputHandle = m_ioBroker.create_output_mesh(out_filename, stk::io::WRITE_RESULTS);

  // Define global fields
  std::string globalFieldFreq = "EigenFrequency";
  m_ioBroker.add_global(resultOutputHandle, globalFieldFreq, Ioss::Field::REAL);

  // Define node fields
  typedef stk::mesh::Field<double> ScalarField;
  typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> VectorField;

  VectorField& displacementsField =
    stkMeshMetaData.declare_field<VectorField>(stk::topology::NODE_RANK, "displacements");
  stk::mesh::put_field_on_entire_mesh(displacementsField);
  
  ScalarField& temperatureField = 
    stkMeshMetaData.declare_field<ScalarField>(stk::topology::NODE_RANK, "temperature");
  stk::mesh::put_field_on_entire_mesh(temperatureField);

  // Add defined fields to the "results output mesh"
  m_ioBroker.add_field(resultOutputHandle, displacementsField, "displ");
  m_ioBroker.add_field(resultOutputHandle, temperatureField);

  // Populate/read in the bulk mesh data
  // Note: MUST be called after defining fields
  m_ioBroker.populate_bulk_data();

  return;
}

// Setup Epetra Row map for distributed vectors/matrices/operators
// Note: must come AFTER populate_bulk_data() has been called
// ============================================================================
Epetra_Map Simulation::setupEpetraRowMap(int& numNonzeroEstimate)
{
  // Wrap the MPI communicator so Epetra can use it
  Epetra_MpiComm epetraComm(m_stkComm);

  // Get reference to bulk_data ... I'll use it a lot here.
  stk::mesh::BulkData& bulkData = m_ioBroker.bulk_data();

  // Select and count locally-owned elements
  std::vector<unsigned> entityCounts;
  stk::mesh::Selector localSelector = m_ioBroker.meta_data().locally_owned_part();
  stk::mesh::count_entities(localSelector, bulkData, entityCounts);

  // Calculate a rough estimate for number of non-zero entries in each row.  A
  // better estimate here may speed things up.
  const int numLocalNodes = entityCounts[stk::topology::NODE_RANK];
  const int numLocalDof = m_spatialDim*numLocalNodes;
  numNonzeroEstimate = m_spatialDim*numLocalDof;

  // A primary objective is to set up the global ID map for use in initializing
  // the Epetra_Map object
  std::vector<int> myGlobalIdMap(numLocalDof);

  // Bucket loop on local nodes to setup the Epetra Row Map 
  const stk::mesh::BucketVector nodeBuckets = 
    bulkData.get_buckets(stk::topology::NODE_RANK, localSelector);

  for (size_t bucketIdx = 0; bucketIdx < nodeBuckets.size(); ++bucketIdx) {
    stk::mesh::Bucket &nodeBucket = *nodeBuckets[bucketIdx];
    for (size_t nodeIdx = 0; nodeIdx < nodeBucket.size(); ++nodeIdx) {
      stk::mesh::Entity node = nodeBucket[nodeIdx];
      const int localId = bulkData.local_id(node);
      const int globalId = bulkData.identifier(node);

      // Assign each DOF to the map.  Note: using local Ids here to index the
      // map, but using the global_id (identifier) to derive a "global DOF
      // identifier" since each node has multiple DOFs.
      // NOTE: this is hard-coded for spatialDim = 3
      myGlobalIdMap[localIdToLocalDof(localId, 0)] = globalIdToGlobalDof(globalId, 0);
      myGlobalIdMap[localIdToLocalDof(localId, 1)] = globalIdToGlobalDof(globalId, 1);
      myGlobalIdMap[localIdToLocalDof(localId, 2)] = globalIdToGlobalDof(globalId, 2);
    }
  }

  // Finally, construct the Epetra_Map
  Epetra_Map epetraRowMap(-1, numLocalDof, &myGlobalIdMap[0], 0, epetraComm);

  return epetraRowMap;
}


// Element Bucket Loop 
// ============================================================================
void Simulation::elementBucketLoop(Element elementData, Epetra_FECrsMatrix& massMatrix,
    Epetra_FECrsMatrix& stiffMatrix) {
  using Intrepid::FieldContainer;
  using Intrepid::CellTools;
  typedef Intrepid::FunctionSpaceTools fst;

  // Get reference to bulk_data ... I'll use it a lot here.
  stk::mesh::BulkData& bulkData = m_ioBroker.bulk_data();

  // Define some important dimensions
  const int numCubPoints = elementData.m_numCubPoints;
  const int numFields = elementData.m_numFields;

  // Get element toplogy and cubature 
  const shards::CellTopology& elemTopology = elementData.getTopology();
  const FieldContainer<double>& cubPoints = elementData.getCubPoints();
  const FieldContainer<double>& cubWeights = elementData.getCubWeights();
  const FieldContainer<double>& cubPointValues = elementData.getCubPointValues();
  const FieldContainer<double>& cubPointGradient = elementData.getCubPointGradient();

  // Select local processor part (includes elements/nodes/etc.)
  stk::mesh::Selector localSelector = m_ioBroker.meta_data().locally_owned_part();
  const stk::mesh::BucketVector elementBuckets = 
    bulkData.get_buckets(stk::topology::ELEMENT_RANK, localSelector);

  // Loop over buckets (then later, loop over elements within a bucket)
  for (size_t bucketIndex = 0; bucketIndex < elementBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket& elemBucket = *elementBuckets[bucketIndex];

    unsigned numElements = elemBucket.size();
    std::cout << "Number of bucket elements: " << numElements << std::endl;

    // Arrays for node coordinates and global dof IDs.
    const int NUM_NODES_HEX8 = 8;
    Intrepid::FieldContainer<double> nodeCoords(numElements, NUM_NODES_HEX8, m_spatialDim);
    Intrepid::FieldContainer<int> nodeDofGlobalIds(numElements, NUM_NODES_HEX8*m_spatialDim);

    // Fill coordinate and global dof arrays
    processBucketNodes(elemBucket, nodeCoords, nodeDofGlobalIds);

    // Compute jacobian
    FieldContainer<double> hexJacobian(numElements, numCubPoints, m_spatialDim, m_spatialDim);
    FieldContainer<double> hexJacobianInv(numElements, numCubPoints, m_spatialDim, m_spatialDim);
    FieldContainer<double> hexJacobianDet(numElements, numCubPoints);
    FieldContainer<double> cellMeasure(numElements, numCubPoints);
    FieldContainer<double> hexGValsTransformed(numElements, numFields, numCubPoints);
    FieldContainer<double> hexGValsTransformedWeighted(numElements, numFields, numCubPoints);

    CellTools<double>::setJacobian(hexJacobian, cubPoints, nodeCoords, elemTopology);
    CellTools<double>::setJacobianInv(hexJacobianInv, hexJacobian);
    CellTools<double>::setJacobianDet(hexJacobianDet, hexJacobian);
    
    // TODO: 
    // simply replicates input cubPointValues for all cells in the set.
    fst::HGRADtransformVALUE<double>(hexGValsTransformed, cubPointValues);

    // compute and multiply cell measure
    fst::computeCellMeasure<double>(cellMeasure, hexJacobianDet, cubWeights);
    fst::multiplyMeasure<double>(hexGValsTransformedWeighted, cellMeasure, 
        hexGValsTransformed);


    // need # element DOF for both mass and stiffness matrices.  should
    // move this to element data class ...
    const int numElemDof = m_spatialDim*numFields;

   
    // TODO: ComputeElementMassMatrices(...)
    const int RHO_TEMP = 1.0; // TODO: Read density from input deck
    FieldContainer<double> elemMassMatrix(numElements, numElemDof, numElemDof);
    elementData.integrateMassMatrix(elemMassMatrix, hexGValsTransformed, 
        hexGValsTransformedWeighted, RHO_TEMP, m_spatialDim);
    
    // TODO: AssembleMassMatrix(...)
    for (size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex) {
      for (int i=0; i < elemMassMatrix.dimension(1); ++i) {
        for (int j=0; j < elemMassMatrix.dimension(2); ++j) {
          int rowIdx = nodeDofGlobalIds(elemIndex, i);
          int colIdx = nodeDofGlobalIds(elemIndex, j);
          double value = elemMassMatrix(elemIndex, i, j);
          massMatrix.InsertGlobalValues(rowIdx, 1, &value, &colIdx);
        }
      }
    }
    

    // TODO: ComputeElementStiffnessMatrices(...)
    FieldContainer<double> hexGGradientTransformed(numElements, numFields, 
        numCubPoints, m_spatialDim);
    FieldContainer<double> hexGGradientTransformedWeighted(numElements, numFields,
        numCubPoints, m_spatialDim);

    // transform to physical coordinates.  Each cell has its own physical coordinates,
    // so must have an array sized to accomodate each coordinate.
    fst::HGRADtransformGRAD<double>(hexGGradientTransformed, hexJacobianInv, cubPointGradient);

    // multiply measure
    fst::multiplyMeasure<double>(hexGGradientTransformedWeighted, cellMeasure, 
        hexGGradientTransformed);

    // integrate element stiffness matrices for work set
    FieldContainer<double> stiffnessMatrix(numElements, numElemDof, numElemDof);
    const double E_TEMP = 10.0; // TODO: Read density from input deck
    const double NU_TEMP = 0.0; // TODO: Read poissons ratio from input deck
    ElasticMaterial elasticMat;
    elasticMat.setYoungsModulus(E_TEMP);
    elasticMat.setPoissonsRatio(NU_TEMP);
    FieldContainer<double> C = elasticMat.stiffnessTensor();

    elementData.integrateStiffnessMatrix(stiffnessMatrix, hexGGradientTransformed,
      hexGGradientTransformedWeighted, C);

    // TODO: AsssembleStiffnessMatrix(...)
    for(size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex){
      for( int i=0; i<stiffnessMatrix.dimension(1); ++i ){
        for( int j=0; j<stiffnessMatrix.dimension(2); ++j ){
          int rowIdx = nodeDofGlobalIds(elemIndex, i);
          int colIdx = nodeDofGlobalIds(elemIndex, j);
          double value = stiffnessMatrix(elemIndex, i, j);
          stiffMatrix.InsertGlobalValues(rowIdx, 1, &value, &colIdx);
        }
      }
    }


  }

  // Call .GlobalAssemble() on mass and stiffness matrix
  massMatrix.GlobalAssemble();
  stiffMatrix.GlobalAssemble();
  
  return;
}

// Process Nodes (in an element bucket)
// ============================================================================
void Simulation::processBucketNodes(const stk::mesh::Bucket& elemBucket,
                              Intrepid::FieldContainer<double>& nodeCoords,
                              Intrepid::FieldContainer<int>& nodeDofGlobalIds) {
  // Get reference to bulk_data
  stk::mesh::BulkData& bulkData = m_ioBroker.bulk_data();

  // Get the coordinates field
  typedef stk::mesh::Field<double, stk::mesh::Cartesian> CoordinatesField_t;
  CoordinatesField_t const& coord_field = 
      static_cast<CoordinatesField_t const&>(m_ioBroker.get_coordinate_field());

  // Loop over each element in the bucket
  for (size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex) {
    stk::mesh::Entity elem = elemBucket[elemIndex];

    // Get nodes array belonging to this element
    unsigned numNodes = bulkData.num_nodes(elem);
    stk::mesh::Entity const* nodes = bulkData.begin_nodes(elem);

    // Loop over each node in the current element
    for (unsigned nodeIndex = 0; nodeIndex < numNodes; ++nodeIndex) {
      double* coords = stk::mesh::field_data(coord_field, nodes[nodeIndex]);

      // Assumption that spatial dimension is 3
      nodeCoords(elemIndex, nodeIndex, 0) = coords[0];
      nodeCoords(elemIndex, nodeIndex, 1) = coords[1];
      nodeCoords(elemIndex, nodeIndex, 2) = coords[2];

      // Assign 'spatial-dim' global IDs to each node for its DOFs
      int globalNodeId = bulkData.identifier(nodes[nodeIndex]);
      nodeDofGlobalIds(elemIndex, m_spatialDim*nodeIndex + 0) = m_spatialDim*globalNodeId + 0;
      nodeDofGlobalIds(elemIndex, m_spatialDim*nodeIndex + 1) = m_spatialDim*globalNodeId + 1;
      nodeDofGlobalIds(elemIndex, m_spatialDim*nodeIndex + 2) = m_spatialDim*globalNodeId + 2;
    }
  }

  return;
}



// Local-to-Global and Global-to-Local DOF Maps
size_t Simulation::localIdToLocalDof(size_t localId, int dofNum)
{ return m_spatialDim*localId + dofNum; }

size_t Simulation::localDofToLocalId(size_t localDof, int dofNum)
{ return (localDof - dofNum)/m_spatialDim; }

size_t Simulation::globalIdToGlobalDof(size_t globalId, int dofNum)
{ return m_spatialDim*globalId + dofNum; }

size_t Simulation::globalDofToGlobalId(size_t globalDof, int dofNum)
{ return (globalDof - dofNum)/m_spatialDim; }

