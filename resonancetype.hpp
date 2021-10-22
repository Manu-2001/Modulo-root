#ifndef RESOUNCETYPE_HPP
#define RESOUNCETYPE_HPP

#include <iostream>
#include <stdexcept>

#include "particletype.hpp"

class ResonanceType : public ParticleType {
 public:
  ResonanceType() = delete;
  ResonanceType(std::string const&, double, int, double);

  ~ResonanceType();

  void Print() const override;

  double GetWidth() const;

 private:
  double const fWidth;
};

#endif