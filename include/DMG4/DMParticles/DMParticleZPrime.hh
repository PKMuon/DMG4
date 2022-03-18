#pragma once

#include <DMParticle.hh>

class DMParticleZPrime : public DMParticle {
  private:
    static DMParticleZPrime * theInstance;
    DMParticleZPrime();
    ~DMParticleZPrime();
  public:
    static DMParticleZPrime * Definition();
};
