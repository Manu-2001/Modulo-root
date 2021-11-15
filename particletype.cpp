#include "particletype.hpp"

// costruttore
ParticleType::ParticleType(std::string const& name, double mass, int charge)
    : fName{name}, fMass{mass}, fCharge{charge} {
  if (fMass < 0.) {
    throw std::runtime_error{
        "Particletype::Particletype(std::string&, double, int) : Invalid "
        "inizialization"};
  }
}

// distruttore
ParticleType::~ParticleType() {}

// metodi pubblici
void ParticleType::Print() const {
  std::cout << "Particle " << this->fName << "\n  Mass: " << this->fMass
            << " GeV/c^2"
            << "\n  Charge: " << this->fCharge << " e\n";
}

double ParticleType::GetWidth() const { return double{}; }

// metodi get
std::string const& ParticleType::GetParticleName() const { return this->fName; }

double ParticleType::GetMass() const { return this->fMass; }

int ParticleType::GetCharge() const { return this->fCharge; }