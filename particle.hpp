#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "particletype.hpp"
#include "point.hpp"
#include "resonancetype.hpp"

class Particle {
 public:
  static void AddParticleType(std::string const&, double, int, double = 0.);
  static void PrintParticleTypes();

  Particle();
  Particle(std::string const&, Point<double> const& = Point<double>{});

  int Decay2body(Particle&, Particle&) const;
  double Energy() const;
  double InvMass(Particle const&) const;
  void Print() const;

  int GetCharge() const;
  double GetMass() const;
  Point<double> const& GetMomentum() const;
  std::string const& GetName() const;
  unsigned int GetIndex() const;
  double GetWidth() const;

  void SetMomentum(Point<double> const&);
  void SetIndex(std::string const&);
  void SetIndex(unsigned int);

 private:
  static unsigned int FindParticle(std::string const&);

  void Boost(Point<double> const&);

  using ParticleTypePtr = std::unique_ptr<ParticleType>;

  static std::vector<ParticleTypePtr> fParticleType;

  unsigned int fIndex;
  Point<double> fMomentum;
};

#endif