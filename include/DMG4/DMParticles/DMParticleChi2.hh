#pragma once

#include <DMParticle.hh>

class DMParticleChi2 : public DMParticle {
  private:
    static DMParticleChi2 * theInstance;
    DMParticleChi2();
    ~DMParticleChi2();
  public:
    static DMParticleChi2 * Definition();
};
