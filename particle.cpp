#include "particle.hpp"

// public static methods
void Particle::AddParticleType(std::string const& name, double mass, int charge,
                               double width) {
  if (FindParticle(name) != fParticleType.size()) {
    return;
  }

  if (width == 0.) {
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

// costructor
Particle::Particle() : fTypeIndex{}, fMomentum{} {}

Particle::Particle(std::string const& name, Point<double> const& momentum)
    : fTypeIndex{}, fMomentum{momentum} {
  auto const index = FindParticle(name);

  if (index == fParticleType.size()) {
    std::cerr << "\nParticle::Particle(std::string&, Point<double>&) : "
                 "Particle name not found\n";
    return;
  }

  fTypeIndex = index;
}

// public methods
int Particle::Decay2body(Particle& dau1, Particle& dau2) const {
  if (GetMass() == 0.) {
    std::cerr << "Decayment cannot be preformed if mass is zero\n";
    return 1;
  }

  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if (fTypeIndex != fParticleType.size()) {
    double x1{}, x2{}, w{}, y1{}, y2{};
    double invnum = 1. / RAND_MAX;

    do {
      x1 = 2. * rand() * invnum - 1.;
      x2 = 2. * rand() * invnum - 1.;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.);

    w = sqrt((-2. * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;
    massMot += fParticleType[fTypeIndex]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    std::cerr << "Decayment cannot be performed because mass is too low in "
                 "this channel\n";
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

  Point<double> P{pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi),
                  pout * cos(theta)};

  dau1.SetMomentum(P);
  dau2.SetMomentum(-P);

  double energy = sqrt(fMomentum.Norm2() + massMot * massMot);

  P = fMomentum / energy;

  dau1.Boost(P);
  dau2.Boost(P);

  return 0;
}

double Particle::Energy() const {
  return std::hypot(GetMass(), fMomentum.Norm());
}

double Particle::InvMass(Particle const& particle) const {
  return sqrt(pow(Energy() + particle.Energy(), 2) -
              (fMomentum + particle.fMomentum).Norm2());
}

void Particle::Print() const {
  std::cout << "\nIndex: " << fTypeIndex << "\n  Type: " << GetName()
            << "\n  Momentum: P(" << fMomentum.x << ", " << fMomentum.y << ", "
            << fMomentum.z << ") GeV/c\n";
}

// get methods
int Particle::GetCharge() const {
  return fParticleType[fTypeIndex]->GetCharge();
}

double Particle::GetMass() const {
  return fParticleType[fTypeIndex]->GetMass();
}

Point<double> const& Particle::GetMomentum() const { return fMomentum; }

std::string const& Particle::GetName() const {
  return fParticleType[fTypeIndex]->GetName();
};

unsigned int Particle::GetTypeIndex() const { return fTypeIndex; }

double Particle::GetWidth() const {
  return fParticleType[fTypeIndex]->GetWidth();
}

// set methods
void Particle::SetMomentum(Point<double> const& momentum) {
  fMomentum = momentum;
}

void Particle::SetTypeIndex(std::string const& name) {
  auto const index = FindParticle(name);

  if (index == fParticleType.size()) {
    std::cerr
        << "Particle::SetTypeIndex(std::string&) : Particle name not found";

    return;
  }

  fTypeIndex = index;
}

void Particle::SetTypeIndex(unsigned int const typeIndex) {
  if (typeIndex >= fParticleType.size()) {
    std::cerr << "Particle::SetTypeIndex(unsigned int) : Invalid typeIndex";

    return;
  }

  fTypeIndex = typeIndex;
}

// private static methods
void Particle::Boost(Point<double> const& b) {
  double b2 = b.Norm2();
  double gamma = 1. / sqrt(1. - b2);
  double bp = fMomentum * b;
  double gamma2 = (b2 > 0.) ? (gamma - 1.) / b2 : 0.;

  fMomentum = b * (gamma2 * bp + gamma * Energy());
}

// private methods
unsigned int Particle::FindParticle(std::string const& name) {
  auto const typePosition =
      std::find_if(fParticleType.begin(), fParticleType.end(),
                   [&name](ParticleTypePtr const& pointer) {
                     return pointer.get()->GetName() == name;
                   });

  return typePosition - fParticleType.begin();
}

// private static attributes
std::vector<std::unique_ptr<ParticleType>> Particle::fParticleType{};