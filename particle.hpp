#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "point.hpp"
#include "particletype.hpp"
#include "resonancetype.hpp"

class Particle {
 public:
  static void AddParticleType(std::string const&, double, int,
                              double = double{});
  static void PrintParticleType();

  Particle(std::string const&, Point<double> const& = Point<double>{});
  Particle();

  double Energy() const;
  double InvMass(Particle const&) const;
  int Decay2body(Particle &,Particle &) const;
  void Print() const;

  Point<double> const& GetMomentum() const;
  unsigned int GetTypeIndex() const;
  double GetMass() const;
  int GetCharge() const;

  void SetTypeIndex(unsigned int);
  void SetTypeIndex(std::string const&);
  void SetMomentum(Point<double> const&);

 private:
  static unsigned int FindParticle(std::string const&);

  void Boost(Point<double> const&);

  using ParticleTypePtr = std::unique_ptr<ParticleType>;

  static std::vector<ParticleTypePtr> fParticleType;

  unsigned int fTypeIndex;
  Point<double> fMomentum;
};

#endif