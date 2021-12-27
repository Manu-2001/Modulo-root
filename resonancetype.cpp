#include "resonancetype.hpp"

// costructor
ResonanceType::ResonanceType(std::string const& name, double mass, int charge,
                             double width)
    : ParticleType(name, mass, charge), fWidth{width} {
  if (fWidth <= 0.) {
    throw std::runtime_error{
        "Resonacetype::Resonancetype(std::string&, double, int, double) : "
        "Invalid inizialization"};
  }
}

// distructor
ResonanceType::~ResonanceType() {}

// public methods
void ResonanceType::Print() const {
  std::cout << "Resonance ";
  ParticleType::Print();
  std::cout << "  Resonance width: " << this->fWidth << " Hz\n";
}

// get methods
double ResonanceType::GetWidth() const { return this->fWidth; }