#pragma once

#include <DMParticle.hh>

class DMParticleLFCScalar : public DMParticle {
  private:
    static DMParticleLFCScalar * theInstance;
    DMParticleLFCScalar();
    ~DMParticleLFCScalar();
  public:
    static DMParticleLFCScalar * Definition();
};
