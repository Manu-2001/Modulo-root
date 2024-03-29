#include <array>
#include <cmath>
#include <stdexcept>

#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TRandom.h"
#include "particle.hpp"

R__LOAD_LIBRARY(point_cpp.so);
R__LOAD_LIBRARY(particletype_cpp.so);
R__LOAD_LIBRARY(resonancetype_cpp.so);
R__LOAD_LIBRARY(particle_cpp.so);

int main() {
  try {
    /*  INITIALIZATION OF VARIABLES AND OBJECT  */

    using iterator = std::array<Particle, 120>::iterator;

    // type of particle
    Particle::AddParticleType("pion+", 0.13957, +1);        // pion+  Index = 0
    Particle::AddParticleType("pion-", 0.13957, -1);        // pion-  Index = 1
    Particle::AddParticleType("kaon+", 0.49367, +1);        // kaon+  Index = 2
    Particle::AddParticleType("kaon-", 0.49367, -1);        // kaon-  Index = 3
    Particle::AddParticleType("proton", 0.93827, +1);       // proton Index = 4
    Particle::AddParticleType("antiproton", 0.93827, -1);   // antip  Index = 5
    Particle::AddParticleType("kaon*", 0.89166, 0, 0.050);  // kaon*  Index = 6

    // file
    TFile* myFile = new TFile("histogram.root", "RECREATE");

    // histograms
    TH1D* hParticleType = new TH1D("ParticleType", "Particle Types", 7, 0, 7);
    TH1D* hTheta = new TH1D("Theta", "Theta", 1000, 0., M_PI);
    TH1D* hPhi = new TH1D("Phi", "Phi", 1000, 0., 2 * M_PI);
    TH1D* hImpulse = new TH1D("Impulse", "Momentum distribution", 1000, 0, 5);
    TH1D* hTrasversalImpulse = new TH1D(
        "TrasversalImpulse", "Trasversal momentum distribution", 1000, 0, 5);
    TH1D* hEnergy = new TH1D("Energy", "Particle energy", 1000, 0, 5);
    TH1D* hInvMass = new TH1D("InvMass", "Invariant Mass", 1000, 0, 5);
    TH1D* hInvMassOppCharge = new TH1D(
        "InvMassOppCharge", "Invariant Mass of opposite charge", 1000, 0, 5);
    TH1D* hInvMassSameCharge = new TH1D(
        "InvMassSameCharge", "Invariant Mass of same charge", 1000, 0, 5);
    TH1D* hInvMassDecay = new TH1D(
        "InvMassDecay", "Invariant Mass Decay particle", 1000, 0.3, 1.5);
    TH1D* hInvMassOppChargePK =
        new TH1D("InvMassOppChargePK",
                 "Invariant Mass opposite charge pione/kaone", 1000, 0, 5);
    TH1D* hInvMassSameChargePK =
        new TH1D("InvMassSameChargePK",
                 "Invariant Mass same charge pione/kaone", 1000, 0, 5);

    // particle array
    std::array<Particle, 120> myParticleArray({});

    // iterator
    iterator lastParticle{};
    iterator particle{};
    iterator next{};
    iterator const first = myParticleArray.begin();
    iterator const last = myParticleArray.end() - 20;
    iterator const kaon = myParticleArray.end() - 1;

    // particle index charge and mass
    unsigned int pIndex{};
    int pCharge{};
    double pMass{};

    // variables
    double phi{};
    double theta{};
    double pNorm{};
    double invMass{};
    double randomNumber{};

    // momentum
    Point<double> P{};

    /*  END OF INITIALIZATION OF VARIABLES / OBJECTS, PROGRAM INSTRUCTIONS */

    // sum in quadrature of the errors
    hInvMassOppCharge->Sumw2();
    hInvMassSameCharge->Sumw2();
    hInvMassSameChargePK->Sumw2();
    hInvMassOppChargePK->Sumw2();

    Particle::PrintParticleTypes();

    // kaone* index
    kaon->SetIndex(6);

    // Seed
    gRandom->SetSeed();

    // start of event generation
    for (int event{}; event != 1E5; ++event) {
      lastParticle = last;

      // array filling
      for (particle = first; particle != lastParticle; ++particle) {
        phi = gRandom->Uniform(0., 2 * M_PI);
        theta = gRandom->Uniform(0., M_PI);
        pNorm = gRandom->Exp(1.);

        P.x = sin(theta) * cos(phi) * pNorm;
        P.y = sin(theta) * sin(phi) * pNorm;
        P.z = cos(theta) * pNorm;

        particle->SetMomentum(P);

        randomNumber = gRandom->Rndm();

        // index, fill energy and type histograms
        if (randomNumber < 0.4) {
          particle->SetIndex(0);  // pion+

          (void)hParticleType->Fill(0);
          (void)hEnergy->Fill(particle->Energy());
        } else if (randomNumber < 0.8) {
          particle->SetIndex(1);  // pion-

          (void)hParticleType->Fill(1);
          (void)hEnergy->Fill(particle->Energy());
        } else if (randomNumber < 0.85) {
          particle->SetIndex(2);  // kaon+

          (void)hParticleType->Fill(2);
          (void)hEnergy->Fill(particle->Energy());
        } else if (randomNumber < 0.9) {
          particle->SetIndex(3);  // kaon-

          (void)hParticleType->Fill(3);
          (void)hEnergy->Fill(particle->Energy());
        } else if (randomNumber < 0.945) {
          particle->SetIndex(4);  //  proton

          (void)hParticleType->Fill(4);
          (void)hEnergy->Fill(particle->Energy());
        } else if (randomNumber < 0.99) {
          particle->SetIndex(5);  //  antiproton

          (void)hParticleType->Fill(5);
          (void)hEnergy->Fill(particle->Energy());
        } else {  //  kaone*
          // particle and particle + 1 will be the result of the decay
          // increment lastParticle to always generate 100 particles
          ++lastParticle;

          if (lastParticle == kaon) {
            throw std::runtime_error{
                "error in main() : myParticleArray is full"};
          }

          kaon->SetMomentum(P);

          (void)hParticleType->Fill(6);
          (void)hEnergy->Fill(kaon->Energy());

          if (gRandom->Rndm() < 0.5) {
            particle->SetIndex(0);  //  pion+
            ++particle;
            particle->SetIndex(3);  //  kaon-
          } else {
            particle->SetIndex(1);  //  pion-
            ++particle;
            particle->SetIndex(2);  //  kaon+
          }

          (void)kaon->Decay2body(*particle, *(particle - 1));

          // fill histogram invariant mass decayed particles
          (void)hInvMassDecay->Fill(particle->InvMass(*(particle - 1)));
        }

        // fill histograms of geometric properties and impulse
        (void)hTheta->Fill(theta);
        (void)hPhi->Fill(phi);
        (void)hImpulse->Fill(pNorm);
        (void)hTrasversalImpulse->Fill(sqrt(P.x * P.x + P.y * P.y));
      }  // end of for loop, array filling completed

      for (particle = first; particle != lastParticle; ++particle) {
        pIndex = particle->GetIndex();
        pCharge = particle->GetCharge();
        pMass = particle->GetMass();

        for (next = particle + 1; next != lastParticle; ++next) {
          invMass = particle->InvMass(*next);

          (void)hInvMass->Fill(invMass);

          // fill histograms of opposite/same charge
          (next->GetCharge() == pCharge)
              ? (void)hInvMassSameCharge->Fill(invMass)
              : (void)hInvMassOppCharge->Fill(invMass);

          // fill histograms of pion/kaon opposite/same charge
          if (pMass != next->GetMass() && pIndex < 4 && next->GetIndex() < 4) {
            (next->GetCharge() == pCharge)
                ? (void)hInvMassSameChargePK->Fill(invMass)
                : (void)hInvMassOppChargePK->Fill(invMass);
          }

        }  // end of for loop
      }    // end of for loop

    }  // end of event for loop

    /*  SAVE HISTOGRAMS */
    (void)myFile->Write();

    myFile->Close();

    /*  ERROR HANDLING, END OF PROGRAM  */
  } catch (std::exception const& Exception) {
    std::cerr << Exception.what() << '\n';
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "an unknown excwption was caught";
    return EXIT_FAILURE;
  }
}