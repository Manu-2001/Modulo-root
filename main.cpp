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
  // bloco try-catch
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

    /*  INIZIALIZZAZIONE VARIABILI/OGGETTI  */

    // canvas
    TCanvas* ParCanvas = new TCanvas("ParticleCanvas", "Particle");
    TCanvas* InvCanvas = new TCanvas("InvMassCanvas", "Massa invariante");
    TCanvas* InvPKCanvas =
        new TCanvas("InvPioneKaoneCanvas", "Massa invariante pione kaone");

    // istogrammi
    TH1D* hParticleType = new TH1D("Type", "Particle Types", 7, 0, 7);
    TH1D* hTheta = new TH1D("Theta", "Theta", 1000, 0., M_PI);
    TH1D* hPhi = new TH1D("Phi", "Phi", 1000, 0., 2 * M_PI);
    TH1D* hImpulse = new TH1D("Impulse", "Momentum distribution", 1000, 0, 5);
    TH1D* hTrasversalImpulse = new TH1D(
        "TrasversalImpulse", "Trasversla momentum distribution", 1000, 0, 5);
    TH1D* hEnergy = new TH1D("Energy", "Particle energy", 1000, 0.1, 5);
    TH1D* hInvMass = new TH1D("InvMass", "Invariant Mass", 1000, 0.225, 5);

    TH1D* hInvMassOppCharge =
        new TH1D("InvMassOppCharge", "Invariant Mass of opposite charge", 1000,
                 0.225, 5);
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
    TH1D* hInvMassDecay = new TH1D(
        "hInvMassDecay", "Invariant Mass Decay particle", 1000, 0.62, 0.88);

    // array delle particelle
    std::array<Particle, 120> myParticleArray({});

    // iteratori
    iterator lastParticle{};
    iterator particle{};
    iterator next{};
    iterator const first = myParticleArray.begin();
    iterator const last = myParticleArray.end() - 20;
    iterator const kaon = myParticleArray.end() - 1;

    // indice di particle e next
    unsigned int pIndex{};
    unsigned int npIndex{};
    int pCharge{};
    int typeError{};

    // variabili
    double phi{};
    double theta{};
    double pNorm{};
    double result{};
    double invMass{};
    Point<double> P{};

    /*  FINE INIZIALIZZAZIONE VARIABILI/OGGETTI, ISTRUZIONI PROGRAMMA */

    // indice kaone*
    kaon->SetTypeIndex("kaon*");

    // inizio generazioni eventi
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

        // indice, fill istogrammi energia e tipo
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
          // particle e particle+1 saranno gli esiti del decadimento
          // incremento last particle per generare sempre 100 particelle
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

          // fill istogramma massa invariante particelle decadute
          hInvMassDecay->Fill(particle->InvMass(*(particle - 1)));
        }

        // fill istogrammi proprietà geometriche e impulso
        hTheta->Fill(theta);
        hPhi->Fill(phi);
        hImpulse->Fill(pNorm);
        hTrasversalImpulse->Fill(sqrt(P.x * P.x + P.y * P.y));
      }  // fine ciclo for, riempimento array completato

      for (particle = first; particle != lastParticle; ++particle) {
        pIndex = particle->GetTypeIndex();
        pCharge = particle->GetCharge();

        for (next = particle + 1; next != lastParticle; ++next) {
          npIndex = next->GetTypeIndex();
          invMass = particle->InvMass(*next);

          hInvMass->Fill(invMass);

          // fill istogrammi carica discorde e concorde
          (next->GetCharge() * pCharge > 0) ? hInvMassSameCharge->Fill(invMass)
                                            : hInvMassOppCharge->Fill(invMass);

          // fill istogrammi invMass pione/kaone
          if (particle->GetMass() != next->GetMass() && pIndex < 4 &&
              npIndex < 4) {
            if ((pIndex == 0 && npIndex == 2) ||
                (npIndex == 0 && pIndex == 2)) {
              hInvMassPpKp->Fill(invMass);
            } else if ((pIndex == 0 && npIndex == 3) ||
                       (npIndex == 0 && pIndex == 3)) {
              hInvMassPpKm->Fill(invMass);
            } else if ((pIndex == 1 && npIndex == 2) ||
                       (npIndex == 1 && pIndex == 2)) {
              hInvMassPmKp->Fill(invMass);
            } else if ((pIndex == 1 && npIndex == 3) ||
                       (npIndex == 1 && pIndex == 3)) {
              hInvMassPmKm->Fill(invMass);
            }
          }

        }  // fine ciclo for
      }    // fine ciclo for

    }  // fine ciclo for eventi

    /*  STAMPA E SALVATAGGIO DEI DATI */

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

    /*  GESTIONE ERRORI E FINE PROGRAMMA  */
  } catch (std::exception const& Exception) {
    std::cerr << Exception.what() << '\n';
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "an unknown excwption was caught";
    return EXIT_FAILURE;
  }  // fine blocco try-catch
}