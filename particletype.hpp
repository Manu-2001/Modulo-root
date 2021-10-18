#ifndef PARTICLETYPE_HPP
#define PARTICLETYPE_HPP

#include <iostream>
#include <stdexcept>

class ParticleType {
 public:
  ParticleType() = delete;
  ParticleType(std::string const&, double, int);

  virtual ~ParticleType();
  virtual void Print() const;

  std::string const& GetParticleName() const;
  double GetMass() const;
  int GetCharge() const;

 private:
  std::string const particle_name;
  double const mass;
  int const charge;
};

#endif