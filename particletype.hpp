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

  int GetCharge() const;
  double GetMass() const;
  std::string const& GetName() const;

 private:
  std::string const fName;
  double const fMass;
  int const fCharge;
};

#endif