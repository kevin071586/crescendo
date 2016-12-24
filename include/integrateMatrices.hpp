#ifndef INTEGRATEMATRICES_HPP
#define INTEGRATEMATRICES_HPP

#include <Intrepid_FieldContainer.hpp>

using namespace Intrepid;
void integrateStiffnessMatrixAB(FieldContainer<double>& stiffnessMatrix, 
                               FieldContainer<double> leftArg, 
                               FieldContainer<double> rightArg, 
                               FieldContainer<double> C, int i, int k);

// void integrateMassMatrixAB(FieldContainer<double>& massMatrixAB, 
//                                FieldContainer<double> leftArg, 
//                                FieldContainer<double> rightArg,
//                                const double rho);
// 
// void integrateStiffnessMatrix(FieldContainer<double>& stiffnessMatrix, 
//                               FieldContainer<double> leftArg, 
//                               FieldContainer<double> rightArg, 
//                               FieldContainer<double> C);
// 
// void integrateMassMatrix(FieldContainer<double>& massMatrix, 
//                               FieldContainer<double> leftArg,
//                               FieldContainer<double> rightArg,
//                               const double rho,
//                               const unsigned space_dim); 

#endif // INTEGRATEMATRICES_HPP
