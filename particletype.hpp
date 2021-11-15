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
  virtual double GetWidth() const;

  std::string const& GetParticleName() const;
  double GetMass() const;
  int GetCharge() const;

 private:
  std::string const fName;
  double const fMass;
  int const fCharge;
};

#endif