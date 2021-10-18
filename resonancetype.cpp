#include "resonancetype.hpp"

// costruttori
ResonanceType::ResonanceType(std::string const& name_, double mass_,
                             int charge_, double fWidth_)
    : ParticleType(name_, mass_, charge_), fWidth{fWidth_} {
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
double ResonanceType::GetFWidth() const { return this->fWidth; }