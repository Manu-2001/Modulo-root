#include "particle.hpp"

// costruttore
Particle::Particle(std::string const& name_, Impulse<double> const& momentum_)
    : typeIndex{}, momentum{momentum_} {
  auto const index = FindParticle(name_);

  if (index == fParticleType.size()) {
    std::cerr
        << "\nParticle::Particle(std::string&, Impulse<double>&) : Particle "
           "name not found\n";
    return;
  }

  this->typeIndex = index;
}

// metodi pubblici
double Particle::Energy() const {
  auto const mass = fParticleType[this->typeIndex]->GetMass();
  return std::hypot(mass, momentum.Norm());
}

double Particle::InvMass(Particle const& particle) const {
  return sqrt(pow(this->Energy() + particle.Energy(), 2) -
              pow((this->momentum).Norm() + particle.momentum.Norm(), 2));
}

void Particle::Print() const {
  std::cout << "\nIndex: " << typeIndex
            << "\n  Type: " << fParticleType[typeIndex]->GetParticleName()
            << "\n  Momentum: P(" << momentum.x << ", " << momentum.y << ", "
            << momentum.z << ") GeV/c\n";
}

// metodi get
Impulse<double> const& Particle::GetMomentum() const { return this->momentum; }

unsigned int Particle::GetTypeIndex() const { return this->typeIndex; }

// metodi set
void Particle::SetTypeIndex(unsigned int const typeIndex_) {
  if (typeIndex_ >= fParticleType.size()) {
    throw std::runtime_error{
        "Particle::SetTypeIndex(unsigned int) : Invalid typeIndex"};
  }

  this->typeIndex = typeIndex_;
}

void Particle::SetTypeIndex(std::string const& name_) {
  auto const index = FindParticle(name_);

  if (index == fParticleType.size()) {
    throw std::runtime_error{
        "Particle::SetTypeIndex(std::string&) : Particle name not found"};
  }

  this->typeIndex = index;
}

void Particle::SetMomentum(Impulse<double> const& momentum_) {
  this->momentum = momentum_;
}

// metodi statici
void Particle::AddParticleType(std::string const& particleName_,
                               double const mass_, int const charge_,
                               double const width_) {
  if (FindParticle(particleName_) != fParticleType.size()) {
    return;
  }

  if (width_ == double{}) {
    fParticleType.push_back(std::make_unique<ParticleType>(
        ParticleType{particleName_, mass_, charge_}));
    return;
  }

  fParticleType.push_back(std::make_unique<ResonanceType>(
      ResonanceType{particleName_, mass_, charge_, width_}));
}

void Particle::PrintParticleType() {
  for (auto const& ParticlePointe : fParticleType) {
    std::cout << '\n';
    ParticlePointe->Print();
  }
}

unsigned int Particle::FindParticle(std::string const& particleName_) {
  auto const typePosition =
      std::find_if(fParticleType.begin(), fParticleType.end(),
                   [&particleName_](ParticleTypePtr const& pointee) -> bool {
                     return pointee.get()->GetParticleName() == particleName_;
                   });

  return typePosition - fParticleType.begin();
}