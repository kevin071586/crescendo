// ----------------------------------------------------------------------------
// TODO List:
//  1) Handle situation when input mesh doesnt exist
//  2) Function to write log file to given stream, one proc output only.
//
// ----------------------------------------------------------------------------

#include "stk_io/DatabasePurpose.hpp"           // for READ_MESH, WRITE_RESULTS, etc
#include "stk_io/StkMeshIoBroker.hpp"           // for StkMeshIoBroker
#include "stk_mesh/base/Field.hpp"              // for Field
#include "stk_mesh/base/MetaData.hpp"           // for MetaData
#include "stk_mesh/base/CoordinateSystems.hpp"  // for Cartesian type, I think.

#include "Simulation.h"
#include "ParserCmdBlock.h"

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

  size_t exoOutputMesh;
  setResultsOutput(exoOutputMesh);

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


  //double rho = cmdBlock.getFieldDouble("density"); 
  //double E = cmdBlock.getFieldDouble("youngs modulus"); 
  //double nu = cmdBlock.getFieldDouble("poissons ratio"); 

