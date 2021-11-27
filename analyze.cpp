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
  TF1* exp = new TF1("Exp", "[0] * exp([1] * x + [2]) + [3]", 0., 8.);

  // istogrammi da file salvato
  TH1D* particleType = (TH1D*)hFile->Get("ParticleType");
  TH1D* theta = (TH1D*)hFile->Get("Theta");
  TH1D* phi = (TH1D*)hFile->Get("Phi");
  TH1D* impulse = (TH1D*)hFile->Get("Impulse");

  TH1D* trasversalImpulse = (TH1D*)hFile->Get("TrasversalImpulse");
  TH1D* energy = (TH1D*)hFile->Get("Energy");
  TH1D* invMass = (TH1D*)hFile->Get("InvMass");
  TH1D* invMassDecay = (TH1D*)hFile->Get("InvMassDecay");

  TH1D* invMassSameCharge = (TH1D*)hFile->Get("InvMassSameCharge");
  TH1D* invMassOppCharge = (TH1D*)hFile->Get("InvMassOppCharge");
  TH1D* invMassSameChargePK = (TH1D*)hFile->Get("InvMassSameChargePK");
  TH1D* invMassOppChargePK = (TH1D*)hFile->Get("InvMassOppChargePK");

  // istogrammi differenza
  TH1D* hDiffSameOppPK =
      new TH1D("DiffSameOppPK", "Pione/Kaone: Opposite Charge - Same Charge",
               1000, 0.62, 3.2);
  TH1D* hDiffSameOppCharge = new TH1D(
      "DiffSameOppCharge", "Opposite Charge - Same Charge", 1000, 0.225, 5);

  // canvas
  TCanvas* RandomValueC =
      new TCanvas("RandomValueCanvas", "Proprietà generate");
  TCanvas* InvMassEnergyC =
      new TCanvas("InvMassEnergyCanvas", "Energia e massa invariante");
  TCanvas* InvSameOppC = new TCanvas(
      "InvOppositeSameCanvas", "Massa invariante carica concorde e discorde");
  TCanvas* InvSubCanvas =
      new TCanvas("InvMassSubCanvas", "Invariant Mass Subtraction");

  double binContent{};
  double binError{};
  double entries{};

  Bool_t typeError{};

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

  /*  PARTE 10 */
  hDiffSameOppPK->Sumw2();
  hDiffSameOppCharge->Sumw2();

  typeError =
      hDiffSameOppPK->Add(invMassOppChargePK, invMassSameChargePK, 1, -1);
  typeError =
      hDiffSameOppCharge->Add(invMassSameCharge, invMassOppCharge, 1, -1);

  /*  PARTE 11  */

  // stampa canvas
  gStyle->SetOptStat(112210);
 
  RandomValueC->Divide(2, 2);
  InvMassEnergyC->Divide(2, 2);
  InvSameOppC->Divide(2, 2);
  InvSubCanvas->Divide(2);

  RandomValueC->cd(1);
  particleType->DrawCopy();
  RandomValueC->cd(2);
  theta->DrawCopy();
  RandomValueC->cd(3);
  phi->DrawCopy();
  RandomValueC->cd(4);
  impulse->DrawCopy();

  InvMassEnergyC->cd(1);
  trasversalImpulse->DrawCopy();
  InvMassEnergyC->cd(2);
  energy->DrawCopy();
  InvMassEnergyC->cd(3);
  invMass->DrawCopy();
  InvMassEnergyC->cd(4);
  invMassDecay->DrawCopy();

  InvSameOppC->cd(1);
  invMassSameCharge->DrawCopy();
  InvSameOppC->cd(2);
  invMassOppCharge->DrawCopy();
  InvSameOppC->cd(3);
  invMassSameChargePK->DrawCopy();
  InvSameOppC->cd(4);
  invMassOppChargePK->DrawCopy();

  InvSubCanvas->cd(1);
  hDiffSameOppPK->DrawCopy();
  InvSubCanvas->cd(2);
  hDiffSameOppCharge->DrawCopy();

  RandomValueC->Print("RandomValue.pdf");
  InvMassEnergyC->Print("EnergiaInvMass.pdf");
  InvSameOppC->Print("InvMasSameOppoCharge.pdf");
  InvSubCanvas->Print("InvMassSubtraction.pdf");

  hFile->Close();
}