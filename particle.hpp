#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "particletype.hpp"
#include "point.hpp"
#include "resonancetype.hpp"

class Particle {
 public:
  static void AddParticleType(std::string const&, double, int,
                              double = double{});
  static void PrintParticleType();

  Particle();
  Particle(std::string const&, Point<double> const& = Point<double>{});

  int Decay2body(Particle&, Particle&) const;
  double Energy() const;
  double InvMass(Particle const&) const;
  void Print() const;

  int GetCharge() const;
  double GetMass() const;
  Point<double> const& GetMomentum() const;
  unsigned int GetTypeIndex() const;

  void SetMomentum(Point<double> const&);
  void SetTypeIndex(std::string const&);
  void SetTypeIndex(unsigned int);

 private:
  static unsigned int FindParticle(std::string const&);

  void Boost(Point<double> const&);

  using ParticleTypePtr = std::unique_ptr<ParticleType>;

  static std::vector<ParticleTypePtr> fParticleType;

  unsigned int fTypeIndex;
  Point<double> fMomentum;
};

#endif