#pragma once

#include <DMParticle.hh>

class DMParticleAPrime : public DMParticle {
  private:
    static DMParticleAPrime * theInstance;
    DMParticleAPrime();
    ~DMParticleAPrime();
    double APrimeWidth(double, double, double);
  public:
    static DMParticleAPrime * Definition();
};
