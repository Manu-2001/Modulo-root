#include "particletype.hpp"

// costruttore
ParticleType::ParticleType(std::string const& particle_name_, double mass_,
                           int charge_)
    : particle_name{particle_name_}, mass{mass_}, charge{charge_} {
  if (mass <= 0.) {
    throw std::runtime_error{
        "Particletype::Particletype(std::string&, double, int) : Invalid "
        "inizialization"};
  }
}

// distruttore
ParticleType::~ParticleType() {}

// metodi pubblici
void ParticleType::Print() const {
  std::cout << "Particle " << this->particle_name << "\n  Mass: " << this->mass
            << " GeV/c^2"
            << "\n  Charge: " << this->charge << " e\n";
}

// metodi get
std::string const& ParticleType::GetParticleName() const {
  return this->particle_name;
}

double ParticleType::GetMass() const { return this->mass; }

int ParticleType::GetCharge() const { return this->charge; }