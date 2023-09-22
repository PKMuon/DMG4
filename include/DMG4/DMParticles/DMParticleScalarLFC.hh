#pragma once

#include <DMParticle.hh>

class DMParticleScalarLFC : public DMParticle {
  private:
    static DMParticleScalarLFC * theInstance;
    DMParticleScalarLFC();
    ~DMParticleScalarLFC();
  public:
    static DMParticleScalarLFC * Definition();
};
