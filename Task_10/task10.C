#include "TLegend.h"
#include <TCanvas.h>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <ROOT/TSeq.hxx>
#include <TVirtualFitter.h>
#include <TFitter.h>
#include <TMinuit.h>

#include <iostream>
#include <fstream>

enum ParIndex_t 
{
    Bkg = 0,
    SigScale,
    SigSigma,
    SigMean,
    N_PAR
};

const std::map<ParIndex_t, TString> parNames 
{
    {Bkg,      "Background"},
    {SigScale, "Gauss scale"},
    {SigSigma,    "Gauss #sigma"},
    {SigMean,     "Gauss #mu"}
};

// Background function
Double_t Background(Double_t *x, Double_t *par)
{
    return par[Bkg];
}

// Gauss Peak function
Double_t Signal(Double_t *x, Double_t *par)
{
    return par[SigScale]*TMath::Gaus(x[0], par[SigMean], par[SigSigma], true);
}

// Sum of Background and peak funciton 
Double_t FitFunction(Double_t *x, Double_t *par)
{
    return Background(x, par) + Signal(x, par);
}

void FitFcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    TH1D *histOne = nullptr; TH1D *histTwo = nullptr;

    gROOT->GetObject("histOne", histOne);
    gROOT->GetObject("histTwo", histTwo);

    const int nBins1 = histOne->GetNbinsX();
    const int nBins2 = histTwo->GetNbinsX();

    Double_t chisq = 0;
    Double_t delta, n, x;
    for (auto i : ROOT::TSeqI(nBins1))
    {
        n = histOne->GetBinContent(i);
        if (n != 0)
        {
            x = histOne->GetBinCenter(i);

            delta = (n - FitFunction(&x, par)) / TMath::Sqrt(n);
            chisq += delta*delta;
        }
    }
    
    for (auto i : ROOT::TSeqI(nBins2))
    {
        n = histTwo->GetBinContent(i);
        if (n != 0)
        {
            x = histTwo->GetBinCenter(i);

            delta = (n - Background(&x, par)) / TMath::Sqrt(n);
            chisq += delta*delta;
            //std::cout << "chisq = " << n << std::endl;
        }
    }

    f = chisq;
}

void FillHist(TH1D *hist, TString filename)
{
    std::ifstream ifile(filename, std::ios::in); 
    if (ifile.is_open())
    {
        Double_t n = 0;
        Double_t x;
        while (ifile >> x)
        {
            hist->Fill(x);
            n++;
        }
    } else 
    {
        std::cerr << "File " << filename << " is not open!!!";
        exit(-1);
    }
}

void task10()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    const Int_t nBins1 = 100;
    const Int_t nBins2 = 100;
    TH1D *histOne = new TH1D("histOne", "Distribution data from data_1.dat;Energy, [MeV];quantity",
                             nBins1, 500, 600);
    TH1D *histTwo = new TH1D("histTwo", "Distribution data from data_2.dat;Energy, [MeV];quantity",
                             nBins2, 500, 600);

    FillHist(histOne, "data_1.dat");
    FillHist(histTwo, "data_2.dat");

    TVirtualFitter *minuit = TVirtualFitter::Fitter(0, 4);
 
    Double_t vstart[N_PAR] = {7,   1,   10, 550};
    Double_t step[N_PAR]   = {0.1, 0.1, 0.01, 0.1};
    for (auto i : {Bkg, SigScale, SigSigma, SigMean})
    {
        minuit->SetParameter(i, (*parNames.find(i)).second,
                             vstart[i], step[i], 0, 0);
        std::cout << minuit->GetParameter(i) << std::endl;;
    }

    minuit->SetFCN(FitFcn);

    // Set Print level 
    Double_t arglist[100];
    arglist[0] = 0;
    minuit->ExecuteCommand("SET PRINT", arglist, 1);

    // Start Fiting 
    arglist[0] = 5000;   // number of function calls 
    arglist[1] = 1;  // tolerance 
    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    Double_t minParams[N_PAR];
    Double_t parErrors[N_PAR];

    std::cout << std::endl << std::endl;
    std::cout << "------------------------------" << std::endl;
    
    for (auto i : {Bkg, SigScale, SigSigma, SigMean})
    {
        minParams[i] = minuit->GetParameter(i);
        parErrors[i] = minuit->GetParError(i);
    }

    std::cout << "Number of Background events: " 
        << (Int_t)(minParams[Bkg] * nBins1) << std::endl;

    //std::cout << "Chi2 equal " << minuit->Chisquare(N_PAR, minParams) << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << std::endl << std::endl;

    TF1 *fit1 = new TF1("fit1", FitFunction, 500, 600, 4);
    TF1 *fit2 = new TF1("fit2", FitFunction, 500, 600, 1);

    for (auto i : {Bkg, SigScale, SigSigma, SigMean})
        fit1->SetParameter(i, minParams[i]);
    fit2->SetParameter(Bkg, minParams[Bkg]);

    fit2->SetParName(0, "Bkg");

    TFile *file = new TFile("newfile.root", "recreate");
    TCanvas *canvas = new TCanvas("Canvas", "Fitting two Histogram", 10, 10, 700, 500);
    canvas->Divide(2, 1);

    canvas->cd(1);
    histOne->Draw("E");
    fit1->Draw("Same");
    histOne->Write();
    fit1->Write();

    canvas->cd(2);
    histTwo->Draw("E");
    fit2->Draw("Same");
    histTwo->Write();
    fit2->Write();

    canvas->SaveAs("FittingHistogram.png");
}
