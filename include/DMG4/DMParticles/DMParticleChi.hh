#pragma once

#include <DMParticle.hh>

class DMParticleChi : public DMParticle {
  private:
    static DMParticleChi * theInstance;
    DMParticleChi();
    ~DMParticleChi();
  public:
    static DMParticleChi * Definition();
};
