#pragma once

#include <DMParticle.hh>

class DMParticleAxial : public DMParticle {
  private:
    static DMParticleAxial * theInstance;
    DMParticleAxial();
    ~DMParticleAxial();
  public:
    static DMParticleAxial * Definition();
};
