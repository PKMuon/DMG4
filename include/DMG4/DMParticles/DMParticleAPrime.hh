#pragma once

#include <DMParticle.hh>

class DMParticleAPrime : public DMParticle {
  private:
    static DMParticleAPrime * theInstance;
    DMParticleAPrime();
    ~DMParticleAPrime();
  public:
    static DMParticleAPrime * Definition();
};
