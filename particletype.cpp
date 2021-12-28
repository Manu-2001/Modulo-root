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
  std::cout << "Particle " << this->fName 
            << "\n  Mass: " << this->fMass << " GeV/c^2"
            << "\n  Charge: " << this->fCharge << " e\n";
}

// get methods
double ParticleType::GetWidth() const { return double{}; }

int ParticleType::GetCharge() const { return this->fCharge; }

double ParticleType::GetMass() const { return this->fMass; }

std::string const& ParticleType::GetParticleName() const { return this->fName; }