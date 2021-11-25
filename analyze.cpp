#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TStyle.h"

void analyze() {
  TFile* hFile = new TFile("histogram.root", "WRITE");

  TH1D* ParticleType = (TH1D*)hFile->Get("hParticleType");
  TH1D* Theta = (TH1D*)hFile->Get("hTheta");
  TH1D* Phi = (TH1D*)hFile->Get("hPhi");
  TH1D* Impulse = (TH1D*)hFile->Get("hImpulse");
  TH1D* TrasversalImpulse = (TH1D*)hFile->Get("hTrasversalImpulse");
  TH1D* Energy = (TH1D*)hFile->Get("hEnergy");
  TH1D* InvMass = (TH1D*)hFile->Get("hInvMass");
  TH1D* InvMassOppCharge = (TH1D*)hFile->Get("hInvMassOppCharge");
  TH1D* InvMassSameCharge = (TH1D*)hFile->Get("hInvMassSameCharge");
  TH1D* InvMassPpKp = (TH1D*)hFile->Get("hInvMassPpKp");
  TH1D* InvMassPpKm = (TH1D*)hFile->Get("hInvMassPpKm");
  TH1D* InvMassPmKp = (TH1D*)hFile->Get("hInvMassPmKp");
  TH1D* InvMassPmKm = (TH1D*)hFile->Get("hInvMassPmKm");
  TH1D* InvMassDecay = (TH1D*)hFile->Get("hInvMassDecay");

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

  ParCanvas->Draw();
  InvCanvas->Draw();
  InvPKCanvas->Draw();

  ParCanvas->Print("Particle.pdf");
  InvCanvas->Print("InvMass.pdf");
  InvPKCanvas->Print("InvPKMass.pdf");

  hFile->Close();
}