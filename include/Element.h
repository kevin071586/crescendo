#ifndef ELEMENT_H 
#define ELEMENT_H

#include <Intrepid_FieldContainer.hpp>
#include "Shards_CellTopology.hpp"              // For element topology

using namespace Intrepid;

class Element
{
  public:
    // Constructor
    Element();

    // Destructor
    ~Element();

    // Initialize cubature - make assumption that we have a HEX8 mesh, and
    // then populates the fields of the parent element.  No information about 
    // the mesh is required to do this, other than the topology of the element
    // itself. 
    void initializeCubature();

    // Integrate stiffness matrix for shape function combination N_a, N_b
    void integrateStiffnessMatrixAB(FieldContainer<double>& stiffnessMatrix, 
                                    FieldContainer<double> leftArg, 
                                    FieldContainer<double> rightArg, 
                                    FieldContainer<double> C, int i, int k);

    // Integrate mass matrix for shape function combination N_a, N_b
    void integrateMassMatrixAB(FieldContainer<double>& massMatrixAB, 
                               FieldContainer<double> leftArg, 
                               FieldContainer<double> rightArg,
                               const double rho);

    // Integrate full stiffness matrix for a single element
    void integrateStiffnessMatrix(FieldContainer<double>& stiffnessMatrix, 
                                  FieldContainer<double> leftArg, 
                                  FieldContainer<double> rightArg, 
                                  FieldContainer<double> C);

    // Integrate full mass matrix for a single element
    void integrateMassMatrix(FieldContainer<double>& massMatrix, 
                             FieldContainer<double> leftArg,
                             FieldContainer<double> rightArg,
                             const double rho,
                             const unsigned space_dim); 
    
    // Get methods
    shards::CellTopology getTopology() const {return m_elemTopology;}
    FieldContainer<double> getCubPoints() const {return m_cubPoints;}
    FieldContainer<double> getCubWeights() const {return m_cubWeights;}
    FieldContainer<double> getCubPointValues() const {return m_cubPointValues;}
    FieldContainer<double> getCubPointGradient() const {return m_cubPointGradient;}

    int m_cubDim{-1};
    int m_numFields{-1};
    int m_numCubPoints{-1};

  private:
    shards::CellTopology m_elemTopology;

    Intrepid::FieldContainer<double> m_cubPoints;
    Intrepid::FieldContainer<double> m_cubWeights;
    Intrepid::FieldContainer<double> m_cubPointValues;
    Intrepid::FieldContainer<double> m_cubPointGradient;
};

#endif // ELEMENT_H
