#include <array>
#include <cmath>
#include <stdexcept>

#include "particle.hpp"
#include "particletype.hpp"
#include "point.hpp"
#include "resonancetype.hpp"

int main() {
  try {
    R__LOAD_LIBRARY(point_cpp.so);
    R__LOAD_LIBRARY(particletype_cpp.so);
    R__LOAD_LIBRARY(resonancetype_cpp.so);
    R__LOAD_LIBRARY(particle_cpp.so);

    // tipi di partielle
    Particle::AddParticleType("pion+", 0.13957, +1);
    Particle::AddParticleType("pion-", 0.13957, -1);
    Particle::AddParticleType("kaon+", 0.49367, +1);
    Particle::AddParticleType("kaon-", 0.49367, -1);
    Particle::AddParticleType("proton", 0.93827, +1);
    Particle::AddParticleType("antiproton", 0.93827, -1);
    Particle::AddParticleType("kaon*", 0.89166, 0, 0.050);

    // canvas e istogrammi
    TCanvas* myCanvas = new TCanvas("Canvas", "Particle");
    myCanvas->Divide(3, 3);

    TH1F* hParticleType = new TH1F("Type", "Particle Types", 7, 0, 7);
    TH1F* hTheta = new TH1F("Theta", "Theta", 1000, 0., M_PI);
    TH1F* hPhi = new TH1F("Phi", "Phi", 1000, 0., 2 * M_PI);
    TH1F* hImpulse = new TH1F("Impulse", "Momentum distribution", 1000, 0, 7);
    TH1F* hTrasversalImpulse = new TH1F(
        "TrasversalImpulse", "Trasversla momentum distribution", 1000, 0, 7);
    TH1F* hEnergy = new TH1F("Energy", "Particle energy", 1000, 0, 7);
    TH1F* hInvMass = new TH1F("InvMass", "Invariant Mass", 1000, 0, 5);

    TH1F* hInvMassOppCharge = new TH1F(
        "InvMassOppCharge", "Invariant Mass of opposite charge", 1000, 0, 5);
    TH1F* hInvMassSameCharge = new TH1F(
        "InvMassSameCharge", "Invariant Mass of same charge", 1000, 0, 5);
    TH1F* hInvMassPpKp =
        new TH1F("InvMassPpK", "Invariant Mass pione+/kaone+", 1000, 0, 5);
    TH1F* hInvMassPpKm =
        new TH1F("hInvMassPpKm", "Invariant Mass pione+/kaone-", 1000, 0, 5);
    TH1F* hInvMassPmKp =
        new TH1F("hInvMassPmKp", "Invariant Mass pione-/kaone+", 1000, 0, 5);
    TH1F* hInvMassPmKm =
        new TH1F("hInvMassPmKm", "Invariant Mass pione-/kaone-", 1000, 0, 5);
    TH1F* hInvMassDecay =
        new TH1F("hInvMassDecay", "Invariant Mass Decay particle", 1000, 0, 5);

    // array delle particelle 100 + 20
    std::array<Particle, 120> myParticleArray({});

    // ultimo elemento dell'array è il kaone* usato per i decaimenti
    auto const kaon = myParticleArray.end() - 1;
    kaon->SetTypeIndex("kaon*");

    // variabili usate per generare i valori random
    double phi{};
    double theta{};
    double pNorm{};
    double result{};

    Point<double> P{};

    // variabili usate negli ultimi cicli for annidati
    // per tenere salvati indici e cariche
    unsigned int pTypeIndex{};
    unsigned int itTypeIndex{};

    int pCharge{};

    // variabile per segnalare errori dovuti al decadimento
    int typeError{};

    // iteratore all'ultima particella da inserire/inserita
    auto lastParticle = myParticleArray.end() - 20;

    for (int event{}; event < 10E5; ++event) {
      lastParticle = myParticleArray.end() - 20;

      for (auto particle = myParticleArray.begin(); particle != lastParticle;
           ++particle) {
        phi = gRandom->Uniform(0., 2 * M_PI);
        theta = gRandom->Uniform(0., M_PI);
        pNorm = gRandom->Exp(1.);

        P.x = sin(theta) * cos(phi) * pNorm;
        P.y = sin(theta) * sin(phi) * pNorm;
        P.z = cos(theta) * pNorm;

        particle->SetMomentum(P);

        result = gRandom->Rndm();

        if (result <= 0.4) {
          // ho deciso di fare la fill di questi due
          // instogrammi qui per via dei kaoni* e per
          // risparmiare tempo sul primo
          particle->SetTypeIndex("pion+");
          hParticleType->Fill(0);
          hEnergy->Fill(particle->Energy());
        } else if (result <= 0.8) {
          particle->SetTypeIndex("pion-");
          hParticleType->Fill(1);
          hEnergy->Fill(particle->Energy());
        } else if (result <= 0.85) {
          particle->SetTypeIndex("kaon+");
          hParticleType->Fill(2);
          hEnergy->Fill(particle->Energy());
        } else if (result <= 0.9) {
          particle->SetTypeIndex("kaon-");
          hParticleType->Fill(3);
          hEnergy->Fill(particle->Energy());
        } else if (result <= 0.945) {
          particle->SetTypeIndex("proton");
          hParticleType->Fill(4);
          hEnergy->Fill(particle->Energy());
        } else if (result <= 0.99) {
          particle->SetTypeIndex("antiproton");
          hParticleType->Fill(5);
          hEnergy->Fill(particle->Energy());
        } else {
          // la particlella puntata da particle e la successiva
          // saranno gli esiti del decadimento, 'setto' due
          // particelle e non una, quindi incrementa lastParticle
          ++lastParticle;
          if (lastParticle == kaon) {
            throw std::runtime_error{
                "error in main() : myParticleArray is full"};
          }

          kaon->SetMomentum(P);

          hParticleType->Fill(6);
          hEnergy->Fill(kaon->Energy());

          result = gRandom->Rndm();

          if (result <= 0.5) {
            particle->SetTypeIndex("pion+");
            ++particle;
            particle->SetTypeIndex("kaon-");
          } else {
            particle->SetTypeIndex("pion-");
            ++particle;
            particle->SetTypeIndex("kaon+");
          }

          typeError = kaon->Decay2body(*particle, *(particle - 1));

          if (typeError) {
            std::cerr << "\ntype error: " << typeError << ".\n";
          }

          hInvMassDecay->Fill(particle->InvMass(*(particle - 1)));
        }

        // fill degli istogrammi, indipendenti da particle
        hTheta->Fill(theta);
        hPhi->Fill(phi);
        hImpulse->Fill(pNorm);
        hTrasversalImpulse->Fill(sqrt(P.x * P.x + P.y * P.y));
      }  // fine ciclo for

      // doppio ciclo for per gli ultimi istogrammi, non avendo kaoni*
      // non devo ogni volta effettuare il controllo
      for (auto particle = myParticleArray.begin(); particle != lastParticle;
           ++particle) {
        pTypeIndex = particle->GetTypeIndex();
        pCharge = particle->GetCharge();

        for (auto it = particle + 1; it != lastParticle; ++it) {
          itTypeIndex = it->GetTypeIndex();

          hInvMass->Fill(particle->InvMass(*it));

          (it->GetCharge() * pCharge > 0)
              ? hInvMassSameCharge->Fill(particle->InvMass(*it))
              : hInvMassOppCharge->Fill(particle->InvMass(*it));

          if (pTypeIndex == 0) {
            if (itTypeIndex == 2) {
              hInvMassPpKp->Fill(particle->InvMass(*it));
            } else if (itTypeIndex == 3) {
              hInvMassPpKm->Fill(particle->InvMass(*it));
            }
          } else if (pTypeIndex == 1) {
            if (itTypeIndex == 2) {
              hInvMassPmKp->Fill(particle->InvMass(*it));
            } else if (itTypeIndex == 3) {
              hInvMassPmKm->Fill(particle->InvMass(*it));
            }
          }

        }  // fine ciclo for
      }    // fine ciclo for
    }      // fine ciclo for

  } catch (std::exception const& Exception) {
    std::cerr << Exception.what() << '\n';
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "an unknown excwption was caught";
    return EXIT_FAILURE;
  }
}