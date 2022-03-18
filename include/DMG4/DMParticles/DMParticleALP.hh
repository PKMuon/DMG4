#pragma once

#include <DMParticle.hh>

class DMParticleALP : public DMParticle {
  private:
    static DMParticleALP * theInstance;
    DMParticleALP();
    ~DMParticleALP();
  public:
    static DMParticleALP * Definition();
};
