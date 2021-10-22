#include "particle.hpp"

// costruttore
Particle::Particle(std::string const& name, Impulse<double> const& momentum)
    : fTypeIndex{}, fMomentum{momentum} {
  auto const index = FindParticle(name);

  if (index == fParticleType.size()) {
    std::cerr
        << "\nParticle::Particle(std::string&, Impulse<double>&) : Particle "
           "name not found\n";
    return;
  }

  this->fTypeIndex = index;
}

// metodi pubblici
double Particle::Energy() const {
  auto const mass = fParticleType[this->fTypeIndex]->GetMass();
  return std::hypot(mass, fMomentum.Norm());
}

double Particle::InvMass(Particle const& particle) const {
  return sqrt(pow(this->Energy() + particle.Energy(), 2) -
              pow((this->fMomentum).Norm() + particle.fMomentum.Norm(), 2));
}

void Particle::Print() const {
  std::cout << "\nIndex: " << fTypeIndex
            << "\n  Type: " << fParticleType[fTypeIndex]->GetParticleName()
            << "\n  Momentum: P(" << fMomentum.x << ", " << fMomentum.y << ", "
            << fMomentum.z << ") GeV/c\n";
}

// metodi get
Impulse<double> const& Particle::GetMomentum() const { return this->fMomentum; }

unsigned int Particle::GetTypeIndex() const { return this->fTypeIndex; }

// metodi set
void Particle::SetTypeIndex(unsigned int const typeIndex) {
  if (typeIndex >= fParticleType.size()) {
    throw std::runtime_error{
        "Particle::SetTypeIndex(unsigned int) : Invalid typeIndex"};
  }

  this->fTypeIndex = typeIndex;
}

void Particle::SetTypeIndex(std::string const& name) {
  auto const index = FindParticle(name);

  if (index == fParticleType.size()) {
    throw std::runtime_error{
        "Particle::SetTypeIndex(std::string&) : Particle name not found"};
  }

  this->fTypeIndex = index;
}

void Particle::SetMomentum(Impulse<double> const& momentum) {
  this->fMomentum = momentum;
}

// metodi statici
void Particle::AddParticleType(std::string const& name, double const mass,
                               int const charge, double const width) {
  if (FindParticle(name) != fParticleType.size()) {
    return;
  }

  if (width == double{}) {
    fParticleType.push_back(
        std::make_unique<ParticleType>(ParticleType{name, mass, charge}));
    return;
  }

  fParticleType.push_back(std::make_unique<ResonanceType>(
      ResonanceType{name, mass, charge, width}));
}

void Particle::PrintParticleType() {
  for (auto const& ParticlePointer : fParticleType) {
    std::cout << '\n';
    ParticlePointer->Print();
  }
}

unsigned int Particle::FindParticle(std::string const& name) {
  auto const typePosition =
      std::find_if(fParticleType.begin(), fParticleType.end(),
                   [&name](ParticleTypePtr const& pointer) -> bool {
                     return pointer.get()->GetParticleName() == name;
                   });

  return typePosition - fParticleType.begin();
}