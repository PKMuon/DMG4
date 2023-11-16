#pragma once

#include <DMParticle.hh>

class DMParticleChiScalar : public DMParticle {
  private:
    static DMParticleChiScalar * theInstance;
    DMParticleChiScalar();
    ~DMParticleChiScalar();
  public:
    static DMParticleChiScalar * Definition();
};
