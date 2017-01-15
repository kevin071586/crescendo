#include <Intrepid_FieldContainer.hpp>
#include <ElasticMaterial.h>

using namespace Intrepid;

FieldContainer<double> ElasticMaterial::stiffnessTensor() {
  double E = m_youngs_modulus;
  double nu = m_poissons_ratio;

  unsigned spatial_dim = 3;
  FieldContainer<double> C(3,3,3,3);

  // Calculate lame parameters
  double lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
  double mu = E/(2.0*(1.0+nu));

  for( size_t i=0; i<spatial_dim; ++i ){
    for( size_t j=0; j<spatial_dim; ++j ){
      for( size_t k=0; k<spatial_dim; ++k ){
        for( size_t l=0; l<spatial_dim; ++l ){
          double dij = (i==j) ? 1.0 : 0.0;
          double dkl = (k==l) ? 1.0 : 0.0;
          double dik = (i==k) ? 1.0 : 0.0;
          double djl = (j==l) ? 1.0 : 0.0;
          double dil = (i==l) ? 1.0 : 0.0;
          double djk = (j==k) ? 1.0 : 0.0;

          // Definition of elastic stiffness tensor
          C(i,j,k,l) = lambda*dij*dkl + mu*(dik*djl + dil*djk);
        }
      }
    }
  }
  
  return C;
}

void ElasticMaterial::setYoungsModulus(double E) {
  m_youngs_modulus = E;  
  return;
}

void ElasticMaterial::setPoissonsRatio(double nu) {
  m_poissons_ratio = nu;
  return;
}
