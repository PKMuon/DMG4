#pragma once

#include <DMParticle.hh>

class DMParticleChi1 : public DMParticle {
  private:
    static DMParticleChi1 * theInstance;
    DMParticleChi1();
    ~DMParticleChi1();
  public:
    static DMParticleChi1 * Definition();
};
