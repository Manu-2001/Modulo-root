#include "particle.hpp"

// costruttore
Particle::Particle(std::string const& name, Point<double> const& momentum)
    : fTypeIndex{}, fMomentum{momentum} {
  auto const index = FindParticle(name);

  if (index == fParticleType.size()) {
    std::cerr
        << "\nParticle::Particle(std::string&, Point<double>&) : Particle "
           "name not found\n";
    return;
  }

  this->fTypeIndex = index;
}

Particle::Particle() : fTypeIndex{}, fMomentum{} {}

// metodi pubblici
double Particle::Energy() const {
  auto const mass = fParticleType[this->fTypeIndex]->GetMass();
  return std::hypot(mass, fMomentum.Norm());
}

double Particle::InvMass(Particle const& particle) const {
  return sqrt(pow(this->Energy() + particle.Energy(), 2) -
              pow((this->fMomentum).Norm() + particle.fMomentum.Norm(), 2));
}

int Particle::Decay2body(Particle& dau1, Particle& dau2) const {
  if (this->GetMass() == 0.0) {
    std::cerr << "Decayment cannot be preformed if mass is zero\n";
    return 1;
  }

  double massMot = this->GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if (this->fTypeIndex != fParticleType.size()) {  // add width effect

    // gaussian random numbers

    float x1, x2, w, y1 /*, y2 (?)*/;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    // y2 = x2 * w;     (?)

    massMot += fParticleType[this->fTypeIndex]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    std::cerr
        << "Decayment cannot be preformed because mass is too low in this "
           "channel\n";
    return 2;
  }

  double pout =
      sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;

  Point<double> Momentum{};

  Momentum.x = pout * sin(theta) * cos(phi);
  Momentum.y = pout * sin(theta) * sin(phi);
  Momentum.z = pout * cos(theta);

  dau1.SetMomentum(Momentum);

  Momentum.x = -Momentum.x;
  Momentum.y = -Momentum.y;
  Momentum.z = -Momentum.z;

  dau2.SetMomentum(Momentum);

  double energy = sqrt(pow(this->fMomentum.Norm(), 2) + massMot * massMot);

  Point<double> b{};
  b.x = this->fMomentum.x / energy;
  b.y = this->fMomentum.y / energy;
  b.z = this->fMomentum.z / energy;

  dau1.Boost(b);
  dau2.Boost(b);

  return 0;
}

void Particle::Print() const {
  std::cout << "\nIndex: " << fTypeIndex
            << "\n  Type: " << fParticleType[fTypeIndex]->GetParticleName()
            << "\n  Momentum: P(" << fMomentum.x << ", " << fMomentum.y << ", "
            << fMomentum.z << ") GeV/c\n";
}

// metodi get
Point<double> const& Particle::GetMomentum() const { return this->fMomentum; }

unsigned int Particle::GetTypeIndex() const { return this->fTypeIndex; }

double Particle::GetMass() const {
  return fParticleType[this->fTypeIndex]->GetMass();
}

int Particle::GetCharge() const {
  return fParticleType[this->fTypeIndex]->GetCharge();
}

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

void Particle::SetMomentum(Point<double> const& momentum) {
  this->fMomentum = momentum;
}

// metodi privati
void Particle::Boost(Point<double> const& b) {
  double const energy = this->Energy();

  // Boost this Lorentz vector
  double b2 = b * b;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = fMomentum.x * b.x + fMomentum.y * b.y + fMomentum.z * b.z;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  this->fMomentum.x = b.x * (gamma2 * bp + gamma * energy);
  this->fMomentum.y = b.y * (gamma2 * bp + gamma * energy);
  this->fMomentum.z = b.z * (gamma2 * bp + gamma * energy);
}

// metodi statici
void Particle::AddParticleType(std::string const& name, double const mass,
                               int const charge, double const width) {
  if (FindParticle(name) != fParticleType.size()) {
    return;
  }

  if (width == double{}) {
    fParticleType.push_back(
        std::unique_ptr<ParticleType>(new ParticleType{name, mass, charge}));
    return;
  }

  fParticleType.push_back(std::unique_ptr<ResonanceType>(
      new ResonanceType{name, mass, charge, width}));
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

std::vector<std::unique_ptr<ParticleType>> Particle::fParticleType{};