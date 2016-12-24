#ifndef ELASTICMATERIAL_HPP
#define ELASTICMATERIAL_HPP

#include <Intrepid_FieldContainer.hpp>

using namespace Intrepid;

class ElasticMaterial {
  private:
    double m_density;
    double m_youngs_modulus;
    double m_poissons_ratio;

  public:
    void setYoungsModulus(double E);
    void setPoissonsRatio(double nu);
    FieldContainer<double> stiffnessTensor();
};

#endif // ELASTICMATERIAL_HPP
