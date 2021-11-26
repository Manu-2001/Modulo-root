#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TStyle.h"

void analyze() {
  // apertura file
  TFile* hFile = new TFile("histogram.root", "READ");

  // funzioni
  TF1* uniform = new TF1("Uniform", "[0]", 0., 2 * M_PI);
  TF1* exp = new TF1("Exp", "[0]*exp([1] * x + [2]) + [3]", 0., 8.);

  // istogrammi da file salvato
  TH1D* particleType = (TH1D*)hFile->Get("ParticleType");
  TH1D* theta = (TH1D*)hFile->Get("Theta");
  TH1D* phi = (TH1D*)hFile->Get("Phi");
  TH1D* impulse = (TH1D*)hFile->Get("Impulse");
  TH1D* trasversalImpulse = (TH1D*)hFile->Get("TrasversalImpulse");
  TH1D* energy = (TH1D*)hFile->Get("Energy");
  TH1D* invMass = (TH1D*)hFile->Get("InvMass");
  TH1D* invMassOppCharge = (TH1D*)hFile->Get("InvMassOppCharge");
  TH1D* invMassSameCharge = (TH1D*)hFile->Get("InvMassSameCharge");
  TH1D* invMassPpKp = (TH1D*)hFile->Get("InvMassPpKp");
  TH1D* invMassPpKm = (TH1D*)hFile->Get("InvMassPpKm");
  TH1D* invMassPmKp = (TH1D*)hFile->Get("InvMassPmKp");
  TH1D* invMassPmKm = (TH1D*)hFile->Get("InvMassPmKm");
  TH1D* invMassDecay = (TH1D*)hFile->Get("InvMassDecay");

  // canvas
  TCanvas* ParCanvas = new TCanvas("ParticleCanvas", "Particle");
  TCanvas* InvCanvas = new TCanvas("InvMassCanvas", "Massa invariante");
  TCanvas* InvPKCanvas =
      new TCanvas("InvPioneKaoneCanvas", "Massa invariante pione kaone");

  double binContent{};
  double binError{};
  double entries{};

  /*  PARTE 9 */

  std::cout << "\n\tHistogram: ParticleType"
            << "\nbin\t\tproportion\n";

  entries = particleType->GetEntries();

  for (int bin{1}; bin != 8; ++bin) {
    binContent = particleType->GetBinContent(bin);
    binError = particleType->GetBinError(bin);

    std::cout << bin << '\t' << binContent / entries << " +- "
              << binError / entries << '\n';
  }

  std::cout << "\n\n\tHistogram: Theta\n";

  theta->Fit("Uniform");

  std::cout << "\nTheta fit result: y = " << uniform->GetParameter(0) << " +- "
            << uniform->GetParError(0) << ".\n"
            << "Chisquare/NDF: " << uniform->GetChisquare() / uniform->GetNDF()
            << ".\n";

  std::cout << "\n\n\tHistogram: Phi\n";

  phi->Fit("Uniform");

  std::cout << "\nPhi fit result: y = " << uniform->GetParameter(0) << " +- "
            << uniform->GetParError(0) << ".\n"
            << "Chisquare/NDF: " << uniform->GetChisquare() / uniform->GetNDF()
            << ".\n";

  std::cout << "\n\n\tHistogram: Momentum Norm"
            << "\n function: par[0] * exp(par[1] * x + par[2]) + par[3]\n\n";

  exp->SetParameters(800E3, -1, 0, 0);

  impulse->Fit("Exp");

  TF1* fitResult = impulse->GetFunction("Exp");

  const Double_t* par = fitResult->GetParameters();
  const Double_t* dpar = fitResult->GetParErrors();

  for (Int_t i{}; i != 4; i++) {
    std::cout << "\n par[" << i << "] = " << par[i] << " dpar[" << i
              << "] = " << dpar[i];
  }

  std::cout << "\n\nChisquare/NDF ="
            << fitResult->GetChisquare() / fitResult->GetNDF() << '\n';

  /*  PARTE 11  */

  // stampa canvas
  ParCanvas->Divide(3, 2);
  InvCanvas->Divide(2, 2);
  InvPKCanvas->Divide(2, 2);

  gStyle->SetOptStat(112210);

  ParCanvas->cd(1);
  particleType->DrawCopy();
  ParCanvas->cd(2);
  theta->DrawCopy();
  ParCanvas->cd(3);
  phi->DrawCopy();
  ParCanvas->cd(4);
  impulse->DrawCopy();
  ParCanvas->cd(5);
  trasversalImpulse->DrawCopy();
  ParCanvas->cd(6);
  energy->DrawCopy();

  InvCanvas->cd(1);
  invMass->DrawCopy();
  InvCanvas->cd(2);
  invMassOppCharge->DrawCopy();
  InvCanvas->cd(3);
  invMassSameCharge->DrawCopy();
  InvCanvas->cd(4);
  invMassDecay->DrawCopy();

  InvPKCanvas->cd(1);
  invMassPpKp->DrawCopy();
  InvPKCanvas->cd(2);
  invMassPpKm->DrawCopy();
  InvPKCanvas->cd(3);
  invMassPmKp->DrawCopy();
  InvPKCanvas->cd(4);
  invMassPmKm->DrawCopy();

  ParCanvas->Draw();
  InvCanvas->Draw();
  InvPKCanvas->Draw();

  ParCanvas->Print("Particle.pdf");
  InvCanvas->Print("InvMass.pdf");
  InvPKCanvas->Print("InvPKMass.pdf");

  hFile->Close();
}