#include "resonancetype.hpp"

// costruttori
ResonanceType::ResonanceType(std::string const& name, double mass, int charge,
                             double width)
    : ParticleType(name, mass, charge), fWidth{width} {
  if (fWidth <= 0.) {
    throw std::runtime_error{
        "Resonacetype::Resonancetype(std::string&, double, int, double) : "
        "Invalid inizialization"};
  }
}

// distruttore
ResonanceType::~ResonanceType() {}

// metodi pubblici
void ResonanceType::Print() const {
  std::cout << "Resonance ";
  ParticleType::Print();
  std::cout << "  Resonance width: " << this->fWidth << " Hz\n";
}

// metodi get
double ResonanceType::GetWidth() const { return this->fWidth; }