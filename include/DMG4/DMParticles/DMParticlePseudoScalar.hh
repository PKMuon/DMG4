#pragma once

#include <DMParticle.hh>

class DMParticlePseudoScalar : public DMParticle {
  private:
    static DMParticlePseudoScalar * theInstance;
    DMParticlePseudoScalar();
    ~DMParticlePseudoScalar();
  public:
    static DMParticlePseudoScalar * Definition();
};
