#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>

enum ParIndex_t 
{
    PScale = 0,
    BWScale,
    BWMean,
    BWWidth,
    N_PAR 
};

const std::map<ParIndex_t, TString> parNames 
{
    {PScale,  "Poisson scale"},
    {BWScale, "Breit Wigner scale"},
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

Double_t fitFunc(Double_t *x, Double_t *par)
{
    return par[PScale]  * TMath::Poisson(x[0], 0.2) 
         + par[BWScale] * TMath::BreitWigner(x[0], par[BWMean], par[BWWidth]);
}

void task12()
{   
    gStyle->SetOptFit(11111);
    TH1D *histogram = new TH1D("hist", 
                               "Distribution data from task10Nov.dat;mm;yields / 0.1mm",
                               100, 0, 10);

    FillHistogram(histogram, "task10Nov.dat");

    TF1 *fit = new TF1("fit", fitFunc, 0, 10, 4);
    for (auto i : {PScale, BWScale, BWMean, BWWidth}) {
        fit->SetParameter(i, 1);
        fit->SetParLimits(i, 0, 1000);
        fit->SetParName(i, (*parNames.find(i)).second);
    }

    TCanvas *c1 = new TCanvas("c1", "Distribution data", 1200, 600);
    histogram->Draw("E");
    TFitResultPtr fitRes = histogram->Fit(fit, "SMPE", "L", 0, 10);

}
