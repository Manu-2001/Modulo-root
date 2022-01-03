#include "particle.hpp"

// public static methods
void Particle::AddParticleType(std::string const &name, double mass, int charge,
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

void Particle::PrintParticleTypes() {
  for (auto const &ParticlePointer : fParticleType) {
    ParticlePointer->Print();
  }
}

// costructor
Particle::Particle() : fIndex{}, fMomentum{} {}

Particle::Particle(std::string const &name, Point<double> const &momentum)
    : fIndex{}, fMomentum{momentum} {
  auto const index = FindParticle(name);

  if (index == fParticleType.size()) {
    throw std::runtime_error{
        "Particle::Particle(std::string&, Point<double>&) : Particle name not "
        "found"};

    return;
  }

  fIndex = index;
}

// public method
int Particle::Decay2body(Particle &dau1, Particle &dau2) const {
  if (GetMass() == 0.) {
    std::cerr << "Decayment cannot be preformed if mass is zero\n";
    return 1;
  }

  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if (fIndex != fParticleType.size()) {
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
    massMot += fParticleType[fIndex]->GetWidth() * y1;
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
  return sqrt(GetMass() * GetMass() + fMomentum.Norm2());
}

double Particle::InvMass(Particle const &particle) const {
  return sqrt(pow(Energy() + particle.Energy(), 2) -
              (fMomentum + particle.fMomentum).Norm2());
}

void Particle::Print() const {
  std::cout << "\nIndex: " << fIndex << "\n  Type: " << GetName()
            << "\n  Momentum: P(" << fMomentum.x << ", " << fMomentum.y << ", "
            << fMomentum.z << ") GeV/c\n";
}

// get methods
int Particle::GetCharge() const { return fParticleType[fIndex]->GetCharge(); }

double Particle::GetMass() const { return fParticleType[fIndex]->GetMass(); }

Point<double> const &Particle::GetMomentum() const { return fMomentum; }

std::string const &Particle::GetName() const {
  return fParticleType[fIndex]->GetName();
};

unsigned int Particle::GetIndex() const { return fIndex; }

double Particle::GetWidth() const { return fParticleType[fIndex]->GetWidth(); }

// set methods
void Particle::SetMomentum(Point<double> const &momentum) {
  fMomentum = momentum;
}

void Particle::SetIndex(std::string const &name) {
  auto const index = FindParticle(name);

  if (index == fParticleType.size()) {
    throw std::runtime_error{
        "Particle::SetIndex(std::string&) : Particle name not found"};

    return;
  }

  fIndex = index;
}

void Particle::SetIndex(unsigned int const typeIndex) {
  if (typeIndex >= fParticleType.size()) {
    throw std::runtime_error{
        "Particle::SetIndex(unsignd int) : Invalid typeIndex"};

    return;
  }

  fIndex = typeIndex;
}

// private static methods
unsigned int Particle::FindParticle(std::string const &name) {
  auto const typePosition =
      std::find_if(fParticleType.begin(), fParticleType.end(),
                   [&name](ParticleTypePtr const &pointer) {
                     return pointer.get()->GetName() == name;
                   });

  return typePosition - fParticleType.begin();
}

// private method
void Particle::Boost(Point<double> const &b) {
  double b2 = b.Norm2();
  double gamma = 1. / sqrt(1. - b2);
  double bp = b * fMomentum;
  double gamma2 = (b2 > 0) ? (gamma - 1.) / b2 : 0.;

  fMomentum += b * (gamma2 * bp + gamma * Energy());
}

// private static attributes
std::vector<std::unique_ptr<ParticleType>> Particle::fParticleType{};
