#include <array>
#include <cmath>
#include <stdexcept>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
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

    using iterator = std::array<Particle, 120>::iterator;

    // tipi di partielle
    Particle::AddParticleType("pion+", 0.13957, +1);
    Particle::AddParticleType("pion-", 0.13957, -1);
    Particle::AddParticleType("kaon+", 0.49367, +1);
    Particle::AddParticleType("kaon-", 0.49367, -1);
    Particle::AddParticleType("proton", 0.93827, +1);
    Particle::AddParticleType("antiproton", 0.93827, -1);
    Particle::AddParticleType("kaon*", 0.89166, 0, 0.050);

    // canvas e istogrammi
    TCanvas* ParCanvas = new TCanvas("ParticleCanvas", "Particle");
    TCanvas* InvCanvas = new TCanvas("InvMassCanvas", "Massa invariante");
    TCanvas* InvPKCanvas =
        new TCanvas("InvPioneKaoneCanvas", "Massa invariante pione kaone");

    TH1D* hParticleType = new TH1D("Type", "Particle Types", 7, 0, 7);
    TH1D* hTheta = new TH1D("Theta", "Theta", 1000, 0., M_PI);
    TH1D* hPhi = new TH1D("Phi", "Phi", 1000, 0., 2 * M_PI);
    TH1D* hImpulse = new TH1D("Impulse", "Momentum distribution", 1000, 0, 5);
    TH1D* hTrasversalImpulse = new TH1D(
        "TrasversalImpulse", "Trasversla momentum distribution", 1000, 0, 5);
    TH1D* hEnergy = new TH1D("Energy", "Particle energy", 1000, 0.1, 5);
    TH1D* hInvMass = new TH1D("InvMass", "Invariant Mass", 1000, 0.225, 5);

    TH1D* hInvMassOppCharge = new TH1D(
        "InvMassOppCharge", "Invariant Mass of opposite charge", 1000, 0.225, 5);
    TH1D* hInvMassSameCharge = new TH1D(
        "InvMassSameCharge", "Invariant Mass of same charge", 1000, 0.225, 5);
    TH1D* hInvMassPpKp =
        new TH1D("InvMassPpK", "Invariant Mass pione+/kaone+", 1000, 0.62, 3);
    TH1D* hInvMassPpKm =
        new TH1D("hInvMassPpKm", "Invariant Mass pione+/kaone-", 1000, 0.62, 3);
    TH1D* hInvMassPmKp =
        new TH1D("hInvMassPmKp", "Invariant Mass pione-/kaone+", 1000, 0.62, 3);
    TH1D* hInvMassPmKm =
        new TH1D("hInvMassPmKm", "Invariant Mass pione-/kaone-", 1000, 0.62, 3);
    TH1D* hInvMassDecay =
        new TH1D("hInvMassDecay", "Invariant Mass Decay particle", 1000, 0.62, 0.88);

    // array delle particelle 100 + 20
    std::array<Particle, 120> myParticleArray({});

    // variabili per generare i valori random
    double phi{};
    double theta{};
    double pNorm{};
    double result{};

    Point<double> P{};

    // variabili di riferimento
    unsigned int pTypeIndex{};
    unsigned int itTypeIndex{};

    int pCharge{};
    int typeError{};

    double invMass{};

    // iteratori utili
    iterator lastParticle{};
    iterator particle{};
    iterator it{};
    iterator const first = myParticleArray.begin();
    iterator const last = myParticleArray.end() - 20;
    iterator const kaon = myParticleArray.end() - 1;

    kaon->SetTypeIndex("kaon*");

    for (int event{}; event != 10E5; ++event) {
      lastParticle = last;

      for (particle = first; particle != lastParticle; ++particle) {
        phi = gRandom->Uniform(0., 2 * M_PI);
        theta = gRandom->Uniform(0., M_PI);
        pNorm = gRandom->Exp(1.);

        P.x = sin(theta) * cos(phi) * pNorm;
        P.y = sin(theta) * sin(phi) * pNorm;
        P.z = cos(theta) * pNorm;

        particle->SetMomentum(P);

        result = gRandom->Rndm();

        if (result <= 0.4) {
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

          if (result < 0.5) {
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

        // fill istogrammi
        hTheta->Fill(theta);
        hPhi->Fill(phi);
        hImpulse->Fill(pNorm);
        hTrasversalImpulse->Fill(sqrt(P.x * P.x + P.y * P.y));
      }  // fine ciclo for

      for (particle = first; particle != lastParticle; ++particle) {
        pTypeIndex = particle->GetTypeIndex();
        pCharge = particle->GetCharge();

        for (it = particle + 1; it != lastParticle; ++it) {
          itTypeIndex = it->GetTypeIndex();
          invMass = particle->InvMass(*it);

          hInvMass->Fill(invMass);

          (it->GetCharge() * pCharge > 0) ? hInvMassSameCharge->Fill(invMass)
                                          : hInvMassOppCharge->Fill(invMass);

          if (particle->GetMass() != it->GetMass() && pTypeIndex < 4 &&
              itTypeIndex < 4) {
            if ((pTypeIndex == 0 && itTypeIndex == 2) ||
                (itTypeIndex == 0 && pTypeIndex == 2)) {
              hInvMassPpKp->Fill(invMass);
            } else if ((pTypeIndex == 0 && itTypeIndex == 3) ||
                       (itTypeIndex == 0 && pTypeIndex == 3)) {
              hInvMassPpKm->Fill(invMass);
            } else if ((pTypeIndex == 1 && itTypeIndex == 2) ||
                       (itTypeIndex == 1 && pTypeIndex == 2)) {
              hInvMassPmKp->Fill(invMass);
            } else if ((pTypeIndex == 1 && itTypeIndex == 3) ||
                       (itTypeIndex == 1 && pTypeIndex == 3)) {
              hInvMassPmKm->Fill(invMass);
            }
          }

        }  // fine ciclo for
      }    // fine ciclo for
    }      // fine ciclo for

    ParCanvas->Divide(3, 2);
    InvCanvas->Divide(2, 2);
    InvPKCanvas->Divide(2, 2);

    gStyle->SetOptStat(112210);

    ParCanvas->cd(1);
    hParticleType->Draw();
    ParCanvas->cd(2);
    hTheta->Draw();
    ParCanvas->cd(3);
    hPhi->Draw();
    ParCanvas->cd(4);
    hImpulse->Draw();
    ParCanvas->cd(5);
    hTrasversalImpulse->Draw();
    ParCanvas->cd(6);
    hEnergy->Draw();

    InvCanvas->cd(1);
    hInvMass->Draw();
    InvCanvas->cd(2);
    hInvMassOppCharge->Draw();
    InvCanvas->cd(3);
    hInvMassSameCharge->Draw();
    InvCanvas->cd(4);
    hInvMassDecay->Draw();

    InvPKCanvas->cd(1);
    hInvMassPpKp->Draw();
    InvPKCanvas->cd(2);
    hInvMassPpKm->Draw();
    InvPKCanvas->cd(3);
    hInvMassPmKp->Draw();
    InvPKCanvas->cd(4);
    hInvMassPmKm->Draw();

    ParCanvas->Draw();
    InvCanvas->Draw();
    InvPKCanvas->Draw();

    ParCanvas->Print("Particle.pdf");
    InvCanvas->Print("InvMass.pdf");
    InvPKCanvas->Print("InvPKMass.pdf");

    ParCanvas->Print("Particle.root");
    InvCanvas->Print("InvMass.root");
    InvPKCanvas->Print("InvPKMass.root");

    ParCanvas->Print("Particle.C");
    InvCanvas->Print("InvMass.C");
    InvPKCanvas->Print("InvPKMass.C");

  } catch (std::exception const& Exception) {
    std::cerr << Exception.what() << '\n';
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "an unknown excwption was caught";
    return EXIT_FAILURE;
  }
}