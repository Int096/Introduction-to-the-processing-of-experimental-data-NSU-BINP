#include <TH1D.h>
#include <TROOT.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>

enum ParIndex_t 
{
    PScale = 0,
    Pmean,
    BWMean,
    BWWidth,
    N_PAR 
};

const std::map<ParIndex_t, TString> parNames 
{
    {PScale,  "Poisson scale"},
    {Pmean,   "Poisson mean"},
    {BWMean,  "Breit Wigner mean"},
    {BWWidth, "Breit Wigner width"},
};

void FillHistogram(TH1D *histogram, TString filename) 
{
    std::ifstream ifile(filename, std::ios::in);
    if (ifile.is_open())
    {
        Double_t x; 
        while (ifile >> x)
            histogram->Fill(x);
    } else {
        std::cerr << "File " << filename << " is not open!\n";
        exit(-1);
    }
}

Double_t bww(Double_t *x, Double_t *par)
{
    return TMath::BreitWigner(x[0], par[0], par[1]);
}

TF1 *bw = new TF1("bw", bww, -10, 10, 2);

Double_t fitFunc(Double_t *x, Double_t *par)
{
    TH1D *hist;
    gROOT->GetObject("hist", hist);

    bw->SetParameter(0, par[BWMean]);
    bw->SetParameter(1, par[BWWidth]);
    
    return 
    par[PScale] * TMath::Poisson(x[0], par[Pmean]) +
    (150 - par[PScale])* TMath::BreitWigner(x[0], par[BWMean], par[BWWidth]);
}

void task12()
{   
    gStyle->SetOptFit(11111);
    
    TH1D *histogram = new TH1D("hist", 
                               "Distribution data from task10Nov.dat;mm;yields / 0.1mm",
                               100, 0, 10);

    FillHistogram(histogram, "task10Nov.dat");

    TF1 *fit = new TF1("fit", fitFunc, 0, 10, 4);

    // Set parameters
    fit->SetParameter(PScale, 1000);
    fit->SetParameter(Pmean, 0.2);
    fit->SetParameter(BWMean, 5);
    fit->SetParameter(BWWidth, 10);

    fit->SetParLimits(PScale, 0, histogram->GetSum());
    fit->SetParLimits(Pmean, 0, 100);
    fit->SetParLimits(BWMean, 1, 10);
    fit->SetParLimits(BWWidth, 0, 100);

    for (auto i : {PScale, Pmean, BWMean, BWWidth}) {
        fit->SetParName(i, (*parNames.find(i)).second);
    }

    TCanvas *c1 = new TCanvas("c1", "Distribution data", 1200, 600);
    histogram->Draw("E");
    TFitResultPtr fitRes = histogram->Fit(fit, "SMPE", "L", 0, 10);

}
