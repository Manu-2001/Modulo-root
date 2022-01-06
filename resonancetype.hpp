#ifndef RESONANCETYPE_HPP
#define RESONANCETYPE_HPP

#include <iostream>
#include <stdexcept>

#include "particletype.hpp"

class ResonanceType : public ParticleType {
 public:
  ResonanceType() = delete;
  ResonanceType(std::string const&, double, int, double);

  ~ResonanceType();

  void Print() const override;

  double GetWidth() const override;

 private:
  double const fWidth;
};

#endif
