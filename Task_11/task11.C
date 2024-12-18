#include "TMinuit.h"
#include "TVirtualFitter.h"
#include <TF1.h>
#include <TGraph.h>
#include <TDataType.h>
#include <TROOT.h>
#include <TFile.h>
#include <ROOT/TSeq.hxx>
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>

#include <map>
#include <fstream>
#include <iostream>

// Для избавления от магических чисел
// N_PAR при такой индексации имеет смысл
// количества параметров
enum ParIndex_t 
{
    Bkg = 0,
    SigScale,
    SigSigma,
    SigMean,
    N_PAR
};

// Названия параметров для удобства визуализации
const std::map<ParIndex_t, TString> parNames 
{
    { Bkg,      "Background"   },
    { SigScale, "Gauss scale"  },
    { SigSigma, "Gauss #sigma" },
    { SigMean,  "Gauss #mu"    }
};

// Функция фона (константа)
Double_t Background(Double_t *x, Double_t *par)
{ return par[Bkg]; }

// Функция сигнала (гаусс)
Double_t Signal(Double_t *x, Double_t *par)
{ return par[SigScale]*TMath::Gaus(x[0], par[SigMean], par[SigSigma], true); }

// Фитируящая функция (гаусс + подложка)
Double_t FitFunction(Double_t *x, Double_t *par) 
{ return Background(x, par) + Signal(x, par); }

// Заполнение гистограммы данными из файла
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

void FitFcn(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t)
{
    TH1D *histOne = nullptr; TH1D *histTwo = nullptr;
    gROOT->GetObject("histGauss", histOne);
    gROOT->GetObject("histBkgrd", histTwo);

    const Int_t nBins = histOne->GetNbinsX();

    Double_t chisq = 0;
    Double_t delta, n, x;

    for (auto i : ROOT::TSeqI(nBins))
    {
        n = histOne->GetBinContent(i);
        x = histOne->GetBinCenter(i);

        delta = TMath::Poisson(n, FitFunction(&x, par));
        chisq += -2 * TMath::Log(delta);
    }
    for (auto i : ROOT::TSeqI(nBins))
    {
        n = histTwo->GetBinContent(i);
        x = histTwo->GetBinCenter(i);

        delta = TMath::Poisson(n, Background(&x, par));
        chisq += -2 * TMath::Log(delta);
    }
    f = chisq;
}

void GetParametersOfMinuit(Double_t *params, Double_t *errpar)
{
    TVirtualFitter *minuit = TVirtualFitter::Fitter(0, N_PAR);

    Double_t vstart[N_PAR] = {7.0, 1.0, 10.0, 550};
    Double_t step[N_PAR]   = {0.1, 0.1, 0.1,  0.1};

    for (auto i : {Bkg, SigScale, SigSigma, SigMean})
         minuit->SetParameter(i,
                              (*parNames.find(i)).second,
                              vstart[i],
                              step[i],
                              0,
                              1500);

    minuit->SetFCN(FitFcn);

    Double_t arglist[100];
    arglist[0] = 10000;
    arglist[1] = 0.1;
    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    for (auto i : {Bkg, SigScale, SigSigma, SigMean})
    {
        params[i] = minuit->GetParameter(i);
        errpar[i] = minuit->GetParError(i);
    }
}

// Первая часть задания, фит двух гистограмм
void FittingTwoHistogram(TH1D *histOne, TH1D *histTwo)
{
    Double_t minParams[N_PAR];
    Double_t errParams[N_PAR];
    GetParametersOfMinuit(minParams, errParams);

    Double_t bkg = (Int_t)(minParams[Bkg] * histOne->GetNbinsX());
    Double_t sig = histOne->GetEntries() - bkg;
    
    // Write result to file  
    std::ofstream ofile("result/result.txt", std::ios::out);
    if (ofile.is_open())
    {
        ofile << "\tFor nBins = 100" << std::endl;
        ofile << "-------------------------------------------" << std::endl;
        ofile << "Number of Background events: " << bkg << std::endl;
        ofile << "Number of Signal events: " << sig << std::endl;
        ofile << "-------------------------------------------" << std::endl;
        ofile << std::endl;
        ofile << "\tFit parameters:" << std::endl;
        ofile << "-------------------------------------------" << std::endl;
        ofile << "Par. Names\t | Value\t | Error" << std::endl;
        ofile << "-------------------------------------------" << std::endl;
        for (auto i : {Bkg, SigScale, SigSigma, SigMean})
            ofile << (*parNames.find(i)).second << "\t | " 
                << minParams[i] << "\t | "
                << errParams[i] << std::endl; 
        ofile.close();
    }

    TF1 *fit1 = new TF1("fit1", FitFunction, 500, 600, 4);
    TF1 *fit2 = new TF1("fit2", Background, 500, 600, 1);

    for (auto i : {Bkg, SigScale, SigSigma, SigMean})
    {
        fit1->SetParameter(i, minParams[i]);
        fit2->SetParameter(i, minParams[i]);
    }

    TFile *file = new TFile("result/result.root", "recreate");
    TCanvas *c1 = new TCanvas("Canvas", "Fitting two Histogram", 10, 10, 1200, 600);
    c1->Divide(2, 1);

    c1->cd(1);
    histOne->Draw("E");
    fit1->Draw("SAME");

    c1->cd(2);
    histTwo->Draw("E");
    fit2->Draw("SAME");

    c1->SaveAs("result/Histograms.png");
    histOne->Write();
    histTwo->Write();
    fit1->Write();
    fit2->Write();
}

void BinsDependency()
{
    TH1D *histOne = nullptr; TH1D *histTwo = nullptr;
    gROOT->GetObject("histGauss", histOne);
    gROOT->GetObject("histBkgrd", histTwo);

    const Int_t n = 10;
    Double_t min = 100, step = 100;
    Double_t bin[n], x[n];
    
    std::ofstream ofile("test.txt", std::ios::out);
    for (Int_t i = 0; i < n; i++)
    {
        bin[i] = min + step*i;

        histOne->SetBins(bin[i], 500, 600);
        histTwo->SetBins(bin[i], 500, 600);

        Double_t minParams[N_PAR]; Double_t errParams[N_PAR];
        GetParametersOfMinuit(minParams, errParams);

        Double_t bkg = (Int_t)(minParams[Bkg] * bin[i]);
        Double_t sig = histOne->GetEntries() - bkg;
        ofile << "For " << min + step*i 
            << " nBins sig = " << sig 
            << " bkg = " << bkg << std::endl;
        
        x[i] = sig;
    }
    ofile.close();
    
    TCanvas *c2 = new TCanvas("Bins", "BinsDependency", 10, 10, 700, 500);
    TGraph *gr = new TGraph(n, bin, x);
    gr->SetTitle("N_{signal} dependency on nBins;nBins;N_{signal}");
    gr->SetMarkerColor(4);
    gr->SetMarkerSize(1.5);
    gr->SetMarkerStyle(21);
    gr->Draw("AP");

    c2->SaveAs("result/NSignalNbins.png");

    histOne->SetBins(100, 500, 600);
    histTwo->SetBins(100, 500, 600);
}

void Countour()
{
    TCanvas *c3 = new TCanvas("c3", "Contours", 10, 10, 600, 800);
    TGraph *gr12 = (TGraph*)gMinuit->Contour(40, 1, 2);

    gMinuit->SetErrorDef(2.25);

    TGraph *gr2 = (TGraph*)gMinuit->Contour(40, 1, 2);
    gr2->SetTitle("Errors Contour");
    gr2->SetFillColor(42);
    gr2->Draw("Alf");
    gr12->Draw("C");
}

void task11()
{
    const Int_t nBins = 100;
    TH1D *histGauss = new TH1D("histGauss", 
                               "Distribution data from data_1.dat;Energy, [MeV];quantity",
                               nBins, 500, 600);
    TH1D *histBkgrd = new TH1D("histBkgrd", 
                               "Distribution data from data_2.dat;Energy, [MeV];quantity",
                               nBins, 500, 600);
    
    FillHistogram(histGauss, "data/data_1.dat");
    FillHistogram(histBkgrd, "data/data_2.dat");

    FittingTwoHistogram(histGauss, histBkgrd);

    Countour();
    
    BinsDependency();

}
