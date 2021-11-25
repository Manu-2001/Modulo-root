#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TStyle.h"

void analyze() {
  TFile* myFile = new TFile("histogram.root", "READ");

  TH1D* ParticleType = (TH1D*)myFile->Get("hParticleType");
  TH1D* Theta = (TH1D*)myFile->Get("hTheta");
  TH1D* Phi = (TH1D*)myFile->Get("hPhi");
  TH1D* Impulse = (TH1D*)myFile->Get("hImpulse");
  TH1D* TrasversalImpulse = (TH1D*)myFile->Get("hTrasversalImpulse");
  TH1D* Energy = (TH1D*)myFile->Get("hEnergy");
  TH1D* InvMass = (TH1D*)myFile->Get("hInvMass");
  TH1D* InvMassOppCharge = (TH1D*)myFile->Get("hInvMassOppCharge");
  TH1D* InvMassSameCharge = (TH1D*)myFile->Get("hInvMassSameCharge");
  TH1D* InvMassPpKp = (TH1D*)myFile->Get("hInvMassPpKp");
  TH1D* InvMassPpKm = (TH1D*)myFile->Get("hInvMassPpKm");
  TH1D* InvMassPmKp = (TH1D*)myFile->Get("hInvMassPmKp");
  TH1D* InvMassPmKm = (TH1D*)myFile->Get("hInvMassPmKm");
  TH1D* InvMassDecay = (TH1D*)myFile->Get("hInvMassDecay");

  // canvas
  TCanvas* ParCanvas = new TCanvas("ParticleCanvas", "Particle");
  TCanvas* InvCanvas = new TCanvas("InvMassCanvas", "Massa invariante");
  TCanvas* InvPKCanvas =
      new TCanvas("InvPioneKaoneCanvas", "Massa invariante pione kaone");

  ParCanvas->Divide(3, 2);
  InvCanvas->Divide(2, 2);
  InvPKCanvas->Divide(2, 2);

  gStyle->SetOptStat(112210);

  ParCanvas->cd(1);
  ParticleType->DrawCopy();
  ParCanvas->cd(2);
  Theta->DrawCopy();
  ParCanvas->cd(3);
  Phi->DrawCopy();
  ParCanvas->cd(4);
  Impulse->DrawCopy();
  ParCanvas->cd(5);
  TrasversalImpulse->DrawCopy();
  ParCanvas->cd(6);
  Energy->DrawCopy();

  InvCanvas->cd(1);
  InvMass->DrawCopy();
  InvCanvas->cd(2);
  InvMassOppCharge->DrawCopy();
  InvCanvas->cd(3);
  InvMassSameCharge->DrawCopy();
  InvCanvas->cd(4);
  InvMassDecay->DrawCopy();

  InvPKCanvas->cd(1);
  InvMassPpKp->DrawCopy();
  InvPKCanvas->cd(2);
  InvMassPpKm->DrawCopy();
  InvPKCanvas->cd(3);
  InvMassPmKp->DrawCopy();
  InvPKCanvas->cd(4);
  InvMassPmKm->DrawCopy();

  ParCanvas->Print("Particle.pdf");
  InvCanvas->Print("InvMass.pdf");
  InvPKCanvas->Print("InvPKMass.pdf");
}