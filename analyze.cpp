#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

void printData(double, double, double = 0.);

void analyze(bool const print = false) {
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
  double value{};

  /*  PARTE 9 */

  gStyle->SetOptFit(111);

  // istogramma ParticleType
  std::cout<<"\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ParticleType  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << "\n\nbin\t\tproportion\n";

  entries = particleType->GetEntries();

  for (int bin = 1; bin != 8; ++bin) {
    binContent = particleType->GetBinContent(bin);
    binError = particleType->GetBinError(bin);

    if (bin == 1 || bin == 2){
        value = 0.4;
    } else if (bin == 3 || bin == 4){
        value = 0.05;
    } else if (bin == 5 || bin == 6){
        value = 0.045;
    } else {
        value = 0.01;
    }

    std::cout << ' ' << bin << '\t';
    printData(binContent / entries, binError / entries, value);
    std::cout<<'\n';
  }

  std::cout<< "\n\nbin\t\tentries\n";

  for (int bin = 1; bin != 8; ++bin) {
    binContent = particleType->GetBinContent(bin);
    binError = particleType->GetBinError(bin);

    if (bin == 1 || bin == 2){
        value = 0.4E7;
    } else if (bin == 3 || bin == 4){
        value = 0.05E7;
    } else if (bin == 5 || bin == 6){
        value = 0.045E7;
    } else {
        value = 0.01E7;
    }

    std::cout << ' ' << bin << '\t';
    printData(binContent, binError, value);
    std::cout<<'\n';
  }
  std::cout<< "\n                    ______________   .   ______________\n\n";

  // istogramma theta
  theta->Fit("pol0", "Q");
  TF1* FitTheta = theta->GetFunction("pol0");
  FitTheta->SetLineColor(kOrange + 10);
  FitTheta->SetLineWidth(2);
  std::cout<<"\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Theta  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << "\n\n fit result: y = ";
            printData(FitTheta->GetParameter(0), FitTheta->GetParError(0), 10000);
  std::cout<< "\n Chisquare: " << FitTheta->GetChisquare()
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
            << "\n\n fit result: y = ";
            printData(fitPhi->GetParameter(0), fitPhi->GetParError(0), 10000);
  std::cout<< "\n Chisquare: " << fitPhi->GetChisquare()
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
            << "\n\n function: A * exp(B * x)\n\n"
            << " fit mean: ";
            printData(mean, meanError, 1);
   std::cout<< "\n histo mean: "; 
            printData(impulse->GetMean(), impulse->GetMeanError(), 1);
   std::cout<< "\n\n\n fit\n Chisquare: " << fitResultExp->GetChisquare()
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
  double const deacyAmp = DecayFit->GetParameter(0);
  double const decayAmpError = DecayFit->GetParError(0);
  double const decayMean = DecayFit->GetParameter(1);
  double const decayMeanError = DecayFit->GetParError(1);
  double const decayStd = DecayFit->GetParameter(2);
  double const decayStdError = DecayFit->GetParError(2);
  double const decayChi = DecayFit->GetChisquare();
  double const decayNDF =  DecayFit->GetNDF();

  hDiffSameOppPK->Fit("myGauss","","",0.6, 1.25);
  TF1* PKfit = hDiffSameOppPK->GetFunction("myGauss");
  double const PKAmp = PKfit->GetParameter(0);
  double const PKAmpError = PKfit->GetParError(0);
  double const PKMean = PKfit->GetParameter(1);
  double const PKMeanError = PKfit->GetParError(1);
  double const PKStd = PKfit->GetParameter(2);
  double const PKStdError = PKfit->GetParError(2);
  double const PKChi = PKfit->GetChisquare();
  double const PKNDF =  PKfit->GetNDF();

  hDiffSameOppCharge->Fit("myGauss","","",0.6, 1.25);
  TF1* DiffSOFit = hDiffSameOppCharge->GetFunction("myGauss");
  double const DiffSOAmp = DiffSOFit->GetParameter(0);
  double const DiffSOAmpError = DiffSOFit->GetParError(0);
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
           << "\n Amplitude: ";
           printData(deacyAmp, decayAmpError);
  std::cout<< " \n mean: ";
           printData(decayMean, decayMeanError, 0.89166);
  std::cout<<"\n std: ";
           printData(decayStd, decayStdError, 0.050);
  std::cout<< "\nform difference same opposite pion/kaon charge histos"
           << "( Chi/NDF: " << PKChi/PKNDF << " )"
           << "\n chi: " << PKChi
           << "\n NDF: " << PKNDF
           << "\n Amplitude: ";
           printData(PKAmp, PKAmpError);
  std::cout<< " \n mean: ";
           printData(PKMean, PKMeanError, 0.89166);
  std::cout<<"\n std: ";
           printData(PKStd, PKStdError, 0.050);
  std::cout<< "\nform difference same opposite charge histos" 
           << "( Chi/NDF: " << DiffSOSChi/DiffSOSNDF << " )"
           << "\n chi: " << DiffSOSChi
           << "\n NDF: " << DiffSOSNDF
           << "\n Amplitude: ";
           printData(DiffSOAmp, DiffSOAmpError);
  std::cout<< " \n mean: ";
           printData(DiffSOMean, DiffSOMeanError, 0.89166);
  std::cout<<"\n std: ";
           printData(DiffSOStd, DiffSOStdError, 0.050);
  std::cout<< "\n                    ______________   .   ______________\n\n";
           
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
  impulse->GetXaxis()->SetTitle("Momentum distribution (Gev/c)");
  impulse->GetYaxis()->SetTitle("Occurrences");
  impulse->SetLineColor(kGreen + 1);
  impulse->SetFillColor(kGreen + 1);

  // trasversal momentum distribution
  trasversalImpulse->GetXaxis()->SetTitle(
      "Transverse momentum distribution Gev/c)");
  trasversalImpulse->GetYaxis()->SetTitle("Occurrences");
  trasversalImpulse->SetLineColor(kGreen + 1);
  trasversalImpulse->SetFillColor(kGreen + 1);

  // energy distribution
  energy->GetXaxis()->SetTitle("Energy distribution (GeV)");
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
  TLegend* ThetaLeg = new TLegend (.1, .7, .5, .9, "Legend");
  ThetaLeg->SetFillColor(0);
  ThetaLeg->AddEntry(theta, "Theta distribution");
  ThetaLeg->AddEntry(FitTheta, "Fit: y = p0");
  theta->DrawCopy();
  ThetaLeg->Draw("SAME");
  
  RandomValueC->cd(3);
  TLegend* PhiLeg = new TLegend (.1, .7, .5, .9, "Legend");
  PhiLeg->SetFillColor(0);
  PhiLeg->AddEntry(phi, "Phi distribution");
  PhiLeg->AddEntry(fitPhi, "Fit: y = p0");
  phi->DrawCopy();
  ThetaLeg->Draw("SAME");
  
  RandomValueC->cd(4);
  TLegend* ImpLegend = new TLegend (.1, .7, .5, .9, "Legend");
  ImpLegend->SetFillColor(0);
  ImpLegend->AddEntry(impulse, "impulse distribution");
  ImpLegend->AddEntry(fitResultExp, "Exponential fit");
  impulse->DrawCopy();
  ImpLegend->Draw("SAME");

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
  TLegend* DecayLegend = new TLegend(.1, .7, .3, .9, "Decay K*");
  DecayLegend->SetFillColor(0);
  DecayLegend->AddEntry(invMassDecay,"Invariant mass decay particle");
  DecayLegend->AddEntry(DecayFit, "Gauss fit");
  invMassDecay->DrawCopy();
  DecayLegend->Draw("SAME");

  InvKaonKanvas->cd(2);
  TLegend* DiffSOLegend = new TLegend(.1, .7, .3, .9, "Difference opposite-equal charge");
  DiffSOLegend->SetFillColor(0);
  DiffSOLegend->AddEntry(hDiffSameOppCharge,"Invariant mass");
  DiffSOLegend->AddEntry(DiffSOFit, "Gauss fit");
  hDiffSameOppCharge->DrawCopy();
  DiffSOLegend->Draw("SAME");

  InvKaonKanvas->cd(3);
  TLegend* PKLegend = new TLegend(.1, .7, .3, .9, "Difference p/k charge of opposite - equale sign");
  PKLegend->SetFillColor(0);
  PKLegend->AddEntry(hDiffSameOppCharge,"Invariant mass");
  PKLegend->AddEntry(DiffSOFit, "Gauss fit");
  hDiffSameOppPK->DrawCopy();
  PKLegend->Draw("SAME");

  if (print){
    RandomValueC->Print("RandomValue.pdf");
    InvMassEnergyC->Print("EnergiaInvMass.pdf");
    InvSameOppC->Print("InvMasSameOppoCharge.pdf");
    InvKaonKanvas->Print("InvMassSubtraction.pdf");
  }
  
  hFile->Close();
  
}

void printData(double data, double dataError, double value){
  if (value == 0.){
    std::cout<< data << " ± " << dataError;
  } else if (value >= data - dataError && value <= data + dataError){
    std::cout<<"\033[32m" << data << " ± " << dataError << "\033[0m";
  } else {
    std::cout<<"\033[35m" << data << " ± " << dataError
             <<"   dis.: " << (value - data) / dataError << " error\033[0m";
  }
}
