#include <TH1D.h>
#include <TGraph.h>
#include <TMinuit.h>
#include <TFile.h>
#include <TStyle.h>
#include <TVirtualFitter.h>
#include <ROOT/TSeq.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TROOT.h>
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
    {SigSigma, "Gauss #sigma"},
    {SigMean,  "Gauss #mu"}
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

    const Int_t nBins = histOne->GetNbinsX();
    
    Double_t chisq = 0;
    Double_t delta, n, x;
    for (auto i : ROOT::TSeqI(nBins))
    {
        n = histOne->GetBinContent(i);
        x = histOne->GetBinCenter(i);
        delta = TMath::Poisson(n, FitFunction(&x, par));
        chisq += -2 *TMath::Log(delta);
    }
    for (auto i : ROOT::TSeqI(nBins))
    {
        n = histTwo->GetBinContent(i);
        x = histTwo->GetBinCenter(i);
            
        delta = TMath::Poisson(n, Background(&x, par));
        chisq += -2 *log(delta);
    }
    f = chisq;
}

void FillHist(TH1D *hist, TString filename)
{
    std::ifstream ifile(filename, std::ios::in);
    if (ifile.is_open())
    {
        Double_t x; 
        while (ifile >> x)
        {
            hist->Fill(x);
        }
    } else {
        std::cerr << "File " << filename << " is not open!";
        exit(-1);
    }
}

void GetParametersOfMinuit(Double_t *params, Double_t *errpar, Int_t npar)
{
    TVirtualFitter *minuit = TVirtualFitter::Fitter(0, 4);
 
    Double_t vstart[N_PAR] = {7.0, 1.0, 10.0, 550};
    Double_t step[N_PAR]   = {0.1, 0.1, 0.01, 0.1};

    for (auto i : {Bkg, SigScale, SigSigma, SigMean})
    {
        minuit->SetParameter(i, (*parNames.find(i)).second,
                             vstart[i], step[i], 0, 2000);
        std::cout << minuit->GetParameter(i) << std::endl;;
    }

    minuit->SetFCN(FitFcn);

    // Set Print level 
    Double_t arglist[100];
    arglist[0] = 0;
    minuit->ExecuteCommand("SET PRINT", arglist, 1);

    // Start Fiting 
    arglist[0] = 50000;   // number of function calls 
    arglist[1] = 0.1;  // tolerance 
    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    
    for (auto i : {Bkg, SigScale, SigSigma, SigMean})
    {
        params[i] = minuit->GetParameter(i);
        errpar[i] = minuit->GetParError(i);
    }

    
    TCanvas *canvas2 = new TCanvas("c2", "contours", 10, 10, 600, 800);
    TGraph *gr12 = (TGraph*)gMinuit->Contour(40, 0, 1);

    gMinuit->SetErrorDef(2.25);

    TGraph *gr2 = (TGraph*)gMinuit->Contour(40, 0, 1);
    gr2->SetFillColor(42);
    gr2->Draw("Alf");
    gr12->Draw("C");
}

void BinsDependency()
{
    TH1D *histOne = nullptr; TH1D *histTwo = nullptr;
     
    gROOT->GetObject("histOne", histOne);
    gROOT->GetObject("histTwo", histTwo);
    
    TGraph *graph = new TGraph();

    const Int_t n = 10;
    Double_t min = 100, max = 1000, step = max / n;
    Double_t bin[n], x[n];

    for (Int_t i = 100, j = 0; i < 1000; i += 100, j++)
    {
        histOne->SetBins(i, 500, 600);
        histTwo->SetBins(i, 500, 600);

        Double_t minParams[N_PAR]; Double_t errParams[N_PAR];
        GetParametersOfMinuit(minParams, errParams, N_PAR);

        Double_t bkg = (Int_t)(minParams[Bkg] * i);
        Double_t sig = histOne->Integral() - bkg;

        bin[j] = i;
        x[j] = sig;
    }

    TGraph *gr = new TGraph(n, bin, x);
    gr->Draw("AC*");

    histOne->SetBins(100, 500, 600);
    histTwo->SetBins(100, 500, 600);
}

void task11() 
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    const Int_t nBins = 100;
    TH1D *histOne = new TH1D("histOne", "Distribution data from data_1.dat",
                             nBins, 500, 600);
    TH1D *histTwo = new TH1D("histTwo", "Distribution data from data_1.dat",
                             nBins, 500, 600);

    FillHist(histOne, "data_1.dat");
    FillHist(histTwo, "data_2.dat");

    Double_t minParams[N_PAR];
    Double_t errParams[N_PAR];

    GetParametersOfMinuit(minParams, errParams, N_PAR);

    Double_t bkg = (Int_t)(minParams[Bkg] * nBins);
    Double_t sig = histOne->Integral() - bkg;
    
    // Write result to file  
    std::ofstream ofile("result.txt", std::ios::out);
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
    TF1 *fitG = new TF1("fitG", Signal, 500, 600, 4);
    TF1 *fitB = new TF1("fitB", Background, 500, 600, 4);

    TF1 *fit2 = new TF1("fit2", FitFunction, 500, 600, 1);

    for (auto i : {Bkg, SigScale, SigSigma, SigMean})
    {
        fit1->SetParameter(i, minParams[i]);
        fitG->SetParameter(i, minParams[i]);
        fitB->SetParameter(i, minParams[i]);
        fit2->SetParameter(i, minParams[i]);
    }

    fitG->SetLineColor(kBlue);
    fitB->SetLineColor(kGreen);

    BinsDependency();


    // Draw Histograms and fit
    TFile *file = new TFile("newfile.root", "recreate");
    TCanvas *canvas = new TCanvas("Canvas", "Fitting two Histogram", 10, 10, 700, 500);
    canvas->Divide(2, 1);

    canvas->cd(1);
    histOne->Draw("E");
    fit1->Draw("SAME");
    fitG->Draw("SAME");
    fitB->Draw("SAME");
    histOne->Write();
    fit1->Write();

    canvas->cd(2);
    histTwo->Draw("E");
    fit2->Draw("SAME");
    histTwo->Write();
    fit2->Write();

    canvas->SaveAs("Histograms.png");
}
