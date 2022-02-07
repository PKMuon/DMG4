#pragma once

#include <DMParticle.hh>

class DMParticleScalar : public DMParticle {
  private:
    static DMParticleScalar * theInstance;
    DMParticleScalar();
    ~DMParticleScalar();
  public:
    static DMParticleScalar * Definition();
};
