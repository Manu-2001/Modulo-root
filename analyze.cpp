#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

void analyze() {
  // apertura file
  TFile* hFile = new TFile("histogram.root", "READ");

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
               1000, 0, 5);
  TH1D* hDiffSameOppCharge = new TH1D(
      "DiffSameOppCharge", "Opposite Charge - Same Charge", 1000, 0, 5);

  double binContent{};
  double binError{};
  double entries{};

  double mean{};
  double meanError{};

  /*  PARTE 9 */

  gStyle->SetOptFit(111);/*
  gStyle->SetOptStat(112210);*/

  // istogramma Tipi di particelle
  std::cout<<"\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ParticleType  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << "\n\n bin\t\tproportion\n";

  entries = particleType->GetEntries();

  for (int bin = 1; bin != 8; ++bin) {
    binContent = particleType->GetBinContent(bin);
    binError = particleType->GetBinError(bin);

    std::cout << ' ' << bin << '\t' << binContent / entries << " ± "
              << binError / entries << '\n';
  }

  std::cout<< "\n\n bin\t\tentries\n";

  for (int bin = 1; bin != 8; ++bin) {
    binContent = particleType->GetBinContent(bin);
    binError = particleType->GetBinError(bin);

    std::cout << ' ' << bin << '\t' << binContent << " ± " << binError << '\n';
  }
  std::cout<< "\n                    ______________   .   ______________\n\n";

  // istogramma theta
  theta->Fit("pol0", "Q");
  TF1* FitTheta = theta->GetFunction("pol0");
  FitTheta->SetLineColor(kOrange + 10);
  FitTheta->SetLineWidth(2);
  std::cout<<"\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Theta  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << "\n\n\tHistogram: Theta\n\n"
            << " Theta fit result: y = " << FitTheta->GetParameter(0) << " ± "
            << FitTheta->GetParError(0) << ".\n"
            << " Chisquare: " << FitTheta->GetChisquare()
            << " \n NDF: " << FitTheta->GetNDF()
            << " \n Chisquare/NDF: "
            << FitTheta->GetChisquare() / FitTheta->GetNDF() << ".\n"
            << "\n                    ______________   .   ______________\n\n";

  // istogramma phi
  phi->Fit("pol0", "Q");
  TF1* fitPhi = phi->GetFunction("pol0");
  fitPhi->SetLineColor(kOrange + 10);
  fitPhi->SetLineWidth(2);

  std::cout<<"\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Phi  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << " \nPhi fit result: y = " << fitPhi->GetParameter(0) << " ± "
            << fitPhi->GetParError(0) << ".\n"
            << " Chisquare: " << fitPhi->GetChisquare()
            << " \n NDF: " << fitPhi->GetNDF()
            << " \n Chisquare/NDF: "
            << fitPhi->GetChisquare() / fitPhi->GetNDF() << ".\n"
            << "\n                    ______________   .   ______________\n\n";

  // istogramma dell'impulso
  impulse->Fit("expo");
  TF1* fitResultExp = impulse->GetFunction("expo");
  fitResultExp->SetLineColor(kOrange + 10);
  fitResultExp->SetLineWidth(2);
  mean = 1 / (-fitResultExp->GetParameter(1));
  meanError =
      fitResultExp->GetParError(1) / (pow(fitResultExp->GetParameter(1), 2));

  std::cout<<"\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Momentum  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << "\n\n\tHistogram: Momentum Norm"
            << "\n function: A * exp(B * x)\n\n"
            << " fit mean: " << mean << " ± " << meanError << '\n'
            << " histo mean: " << impulse->GetMean() << " ± "
            << impulse->GetMeanError() << '\n'
            << " \n\n fit\n Chisquare: " << fitResultExp->GetChisquare()
            << " \n NDF: " << fitResultExp->GetNDF()
            << "\n Chisquare/NDF : "
            << fitResultExp->GetChisquare() / fitResultExp->GetNDF()
            << "\n                    ______________   .   ______________\n\n";

  /*  PARTE 10 */

  hDiffSameOppPK->Sumw2();
  hDiffSameOppCharge->Sumw2();

  (void)hDiffSameOppPK->Add(invMassOppChargePK, invMassSameChargePK, 1, -1);
  (void)hDiffSameOppCharge->Add(invMassOppCharge, invMassSameCharge, 1, -1);

  TF1* fGaus = new TF1("myGauss", "gaus", 0.6, 1.25);

  invMassDecay->Fit("myGauss","","",0.6, 1.25);
  TF1* DecayFit = invMassDecay->GetFunction("myGauss");
  double const decayMean = DecayFit->GetParameter(1);
  double const decayMeanError = DecayFit->GetParError(1);
  double const decayStd = DecayFit->GetParameter(2);
  double const decayStdError = DecayFit->GetParError(2);
  double const decayChi = DecayFit->GetChisquare();
  double const decayNDF =  DecayFit->GetNDF();

  hDiffSameOppPK->Fit("myGauss","","",0.6, 1.25);
  TF1* PKfit = hDiffSameOppPK->GetFunction("myGauss");
  double const PKMean = PKfit->GetParameter(1);
  double const PKMeanError = PKfit->GetParError(1);
  double const PKStd = PKfit->GetParameter(2);
  double const PKStdError = PKfit->GetParError(2);
  double const PKChi = PKfit->GetChisquare();
  double const PKNDF =  PKfit->GetNDF();

  hDiffSameOppCharge->Fit("myGauss","","",0.6, 1.25);
  TF1* DiffSOFit = hDiffSameOppCharge->GetFunction("myGauss");
  double const DiffSOMean = DiffSOFit->GetParameter(1);
  double const DiffSOMeanError = DiffSOFit->GetParError(1);
  double const DiffSOStd = DiffSOFit->GetParameter(2);
  double const DiffSOStdError = DiffSOFit->GetParError(2);
  double const DiffSOSChi = DiffSOFit->GetChisquare();
  double const DiffSOSNDF =  DiffSOFit->GetNDF();

  std::cout<<"\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Kaon*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
           << "\n\nform decay particles ( Chi/NDF: " << decayChi/decayNDF << " )"
           << "\n chi: " << decayChi
           << "\n NDF: " << decayNDF
           <<" \n mean: " << decayMean <<" ± " << decayMeanError
           <<".\tstd: " << decayStd << " ± " << decayStdError
           << "\nform difference same opposite pion/kaon charge histos"
           << "( Chi/NDF: " << PKChi/PKNDF << " )"
           << "\n chi: " << PKChi
           << "\n NDF: " << PKNDF
           <<"\n mean: " << PKMean <<" ± " << PKMeanError
           <<".\tstd: " << PKStd << " ± " << PKStdError
           << "\nform difference same opposite charge histos"
           << "( Chi/NDF: " << DiffSOSChi/DiffSOSNDF << " )"
           << "\n chi: " << DiffSOSChi
           << "\n NDF: " << DiffSOSNDF
           <<"\n mean: " << DiffSOMean <<" ± " << DiffSOMeanError
           <<".\tstd: " << DiffSOStd << " ± " << DiffSOStdError <<"\n"
           << "                    ______________   .   ______________\n\n";
           
  /*  PARTE 11  */

  // particleType
  particleType->GetXaxis()->SetBinLabel(1, "pion+");
  particleType->GetXaxis()->SetBinLabel(2, "pion-");
  particleType->GetXaxis()->SetBinLabel(3, "kaon+");
  particleType->GetXaxis()->SetBinLabel(4, "kaon-");
  particleType->GetXaxis()->SetBinLabel(5, "proton");
  particleType->GetXaxis()->SetBinLabel(6, "antiproton");
  particleType->GetXaxis()->SetBinLabel(7, "kaon*");
  particleType->GetXaxis()->SetTitle("Types of particles");
  particleType->GetYaxis()->SetTitle("Occurrences");
  particleType->SetFillColor(kGreen + 1);
  particleType->SetLineColor(kGreen + 1);

  // phi distribution
  phi->GetXaxis()->SetTitle("Phi distribution (rad)");
  phi->GetYaxis()->SetTitle("Occurrences");
  phi->SetLineColor(kGreen + 1);
  phi->SetFillColor(kGreen + 1);

  // theta distribution
  theta->GetXaxis()->SetTitle("Theta distribution (rad)");
  theta->GetYaxis()->SetTitle("Occurrences");
  theta->SetLineColor(kGreen + 1);
  theta->SetFillColor(kGreen + 1);

  // momentum distribution
  impulse->GetXaxis()->SetTitle("Momentum distribution (kg*m/s)");
  impulse->GetYaxis()->SetTitle("Occurrences");
  impulse->SetLineColor(kGreen + 1);
  impulse->SetFillColor(kGreen + 1);

  // trasversal momentum distribution
  trasversalImpulse->GetXaxis()->SetTitle(
      "Transverse momentum distribution kg*m/s)");
  trasversalImpulse->GetYaxis()->SetTitle("Occurrences");
  trasversalImpulse->SetLineColor(kGreen + 1);
  trasversalImpulse->SetFillColor(kGreen + 1);

  // energy distribution
  energy->GetXaxis()->SetTitle("Energy distribution (J)");
  energy->GetYaxis()->SetTitle("Occurrences");
  energy->SetLineColor(kGreen + 1);
  energy->SetFillColor(kGreen + 1);

  // ivariant mass
  invMass->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
  invMass->GetYaxis()->SetTitle("Occurrences");
  invMass->SetLineColor(kGreen + 1);
  invMass->SetFillColor(kGreen + 1);

  // invariant mass Same Charge
  invMassSameCharge->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
  invMassSameCharge->GetYaxis()->SetTitle("Occurrences");
  invMassSameCharge->SetLineColor(kGreen + 1);
  invMassSameCharge->SetFillColor(kGreen + 1);

  // invariant mass Opposite Charge
  invMassOppCharge->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
  invMassOppCharge->GetYaxis()->SetTitle("Occurrences");
  invMassOppCharge->SetLineColor(kGreen + 1);
  invMassOppCharge->SetFillColor(kGreen + 1);

  // invariant mass Same Charge Pione Kaone
  invMassSameChargePK->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
  invMassSameChargePK->GetYaxis()->SetTitle("Occurrences");
  invMassSameChargePK->SetLineColor(kGreen + 1);
  invMassSameChargePK->SetFillColor(kGreen + 1);

  // invariant mass Opposite Charge Pione Kaone
  invMassOppChargePK->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
  invMassOppChargePK->GetYaxis()->SetTitle("Occurrences");
  invMassOppChargePK->SetLineColor(kGreen + 1);
  invMassOppChargePK->SetFillColor(kGreen + 1);

  // invariant mass decay particle
  invMassDecay->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
  invMassDecay->GetYaxis()->SetTitle("Occurrences");
  invMassDecay->SetLineColor(kGreen + 1);
  invMassDecay->SetFillColor(kGreen + 1);

  // histo differenza pk
  hDiffSameOppPK->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
  hDiffSameOppPK->GetYaxis()->SetTitle("Occurrences");
  hDiffSameOppPK->SetLineColor(kGreen + 1);

  // histo differenza
  hDiffSameOppCharge->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
  hDiffSameOppCharge->GetYaxis()->SetTitle("Occurrences");
  hDiffSameOppCharge->SetLineColor(kGreen + 1);

  // canvas con angoli, impulso e distribuzione particelle
  TCanvas* RandomValueC =
      new TCanvas("RandomValueCanvas", "Proprietà generate");

  RandomValueC->Divide(2, 2);

  RandomValueC->cd(1);
  particleType->DrawCopy();
  RandomValueC->cd(2);
  theta->DrawCopy();
  RandomValueC->cd(3);
  phi->DrawCopy();
  RandomValueC->cd(4);
  impulse->DrawCopy();

  // canvas con momento trasversale, energia massa invariante
  TCanvas* InvMassEnergyC =
      new TCanvas("InvMassEnergyCanvas", "Energia e massa invariante");

  InvMassEnergyC->Divide(2, 2);

  InvMassEnergyC->cd(1);
  trasversalImpulse->DrawCopy();
  InvMassEnergyC->cd(2);
  energy->DrawCopy();
  InvMassEnergyC->cd(3);
  invMass->DrawCopy();
  InvMassEnergyC->cd(4);
  invMassDecay->DrawCopy();

  // canvas massa invariante segno
  TCanvas* InvSameOppC = new TCanvas(
      "InvOppositeSameCanvas", "Massa invariante carica concorde e discorde");

  InvSameOppC->Divide(2, 2);

  InvSameOppC->cd(1);
  invMassSameCharge->DrawCopy("histo");
  InvSameOppC->cd(2);
  invMassOppCharge->DrawCopy("histo");
  InvSameOppC->cd(3);
  invMassSameChargePK->DrawCopy("histo");
  InvSameOppC->cd(4);
  invMassOppChargePK->DrawCopy("histo", "E");

  // canvas differenza
  TCanvas* InvKaonKanvas =
      new TCanvas("InvMassSubCanvas", "Invariant Mass Subtraction");

  InvKaonKanvas->Divide(1, 3);

  InvKaonKanvas->cd(1);
  invMassDecay->DrawCopy();
  InvKaonKanvas->cd(2);
  hDiffSameOppCharge->DrawCopy();
  InvKaonKanvas->cd(3);
  hDiffSameOppPK->DrawCopy();

  RandomValueC->Print("RandomValue.pdf");
  InvMassEnergyC->Print("EnergiaInvMass.pdf");
  InvSameOppC->Print("InvMasSameOppoCharge.pdf");
  InvKaonKanvas->Print("InvMassSubtraction.pdf");

  RandomValueC->Print("RandomValue.jpg");
  InvMassEnergyC->Print("EnergiaInvMass.jpg");
  InvSameOppC->Print("InvMasSameOppoCharge.jpg");
  InvKaonKanvas->Print("InvMassSubtraction.jpg");

  hFile->Close();
}