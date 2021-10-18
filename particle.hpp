#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "impulse.hpp"
#include "particletype.hpp"
#include "resonancetype.hpp"

class Particle {
 public:
  static void AddParticleType(std::string const&, double, int, double);
  static void PrintParticleType();

  Particle(std::string const&, Impulse<double>);

  double Energy() const;
  double InvMass(Particle const&) const;
  void Print() const;

  Impulse<double> const& GetMomentum() const;
  unsigned int GetTypeIndex() const;

  void SetTypeIndex(unsigned int);
  void SetTypeIndex(std::string const&);
  void SetMomentum(Impulse<double> const&);

 private:
  static unsigned int FindParticle(std::string const&);

  using ParticleTypePtr = std::unique_ptr<ParticleType>;

  static inline std::vector<ParticleTypePtr> fParticleType{};

  unsigned int typeIndex;
  Impulse<double> momentum;
};

#endif