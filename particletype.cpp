#include "particletype.hpp"

// costructor
ParticleType::ParticleType(std::string const& name, double mass, int charge)
    : fName{name}, fMass{mass}, fCharge{charge} {
  if (fMass < 0.) {
    throw std::runtime_error{
        "Particletype::Particletype(std::string&, double, int) : Invalid "
        "inizialization"};
  }
}

// distructor
ParticleType::~ParticleType() {}

// public methods
void ParticleType::Print() const {
  std::cout << "Particle " << fName << "\n  Mass: " << fMass << " GeV/c^2"
            << "\n  Charge: " << fCharge << " e\n";
}

// get methods
double ParticleType::GetWidth() const { return double{}; }

int ParticleType::GetCharge() const { return fCharge; }

double ParticleType::GetMass() const { return fMass; }

std::string const& ParticleType::GetParticleName() const { return fName; }