#include "TF1Convolution.h"
#include <TStyle.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TMath.h>
#include <fstream>

enum Index_t {
    LineB = 0, 
    LineK, 
    BWScale,
    GSigma1,
    GSigma2, 
    GMean2,
    N_PAR
};

std::map<Index_t, TString> IndexNames {
    {LineB,   "Line B"},
    {LineK,   "Line K"},
    {BWScale, "Breit Wigner scale"},
    {GSigma1, "Gauss #sigma_{1}"},
    {GSigma2, "Gauss #sigma_{2}"},
    {GMean2,  "Gauss #mu_{2}"}
};

std::map<Index_t, Double_t> IndexInit {
    {LineB,   -10.},
    {LineK,   1.},
    {BWScale, 11.},
    {GSigma1, 1},
    {GSigma2, 0.01},
    {GMean2,  0.1}
};

Double_t Response(Double_t x, Double_t GSigma1, Double_t GSigma2, Double_t GMean2)
{
    using namespace TMath;
    Double_t mu1 = 0;
    Double_t G1 = Exp(- 1./2 * ((x-mu1)/GSigma1) * ((x-mu1)/GSigma1)) 
                / (Sqrt(2 * Pi() * GSigma1));
    Double_t G2 = Exp(- 1./2 * ((x-GMean2)/GSigma2) * ((x-GMean2)/GSigma2))
                / (Sqrt(2 * Pi() * GSigma2));

    Double_t g1 = TMath::Gaus(x, 0, GSigma1, kTRUE);
    Double_t g2 = TMath::Gaus(x, GMean2, GSigma2, kTRUE);

    return G1 + G2;
}

Double_t Signal(Double_t x)
{
    return TMath::BreitWigner(x, 3.0969, 0.000093);
}

Double_t BG(Double_t x, Double_t LineB, Double_t LineK)
{
    return LineB + LineK * x;
}

Double_t Convolution(Double_t *x, Double_t *par)
{
    Double_t TMin = 2;
    Double_t TMax = 3.2;
    Int_t NStep = 400;
    Double_t Steps = (TMax - TMin) / NStep;

    Double_t Sum = 0.;
    for (Int_t i = 0; i < NStep; i++) {
        Double_t t = TMin + (i + 0.5) * Steps;
        Double_t BW = Signal(t);
        Double_t RP = Response(x[0] - t, par[GSigma1], par[GMean2], par[GSigma2]);
        Sum += BW * RP;
    }
    Sum *= Steps;
    return par[BWScale] * Sum + BG(x[0], par[LineB], par[LineK]);
}

Int_t convolution()
{
    gStyle->SetOptStat(111);
    gStyle->SetOptFit(1);

    // Create hist for data
    Int_t N = 400;
    TH1 *hist = new TH1D("hist", 
                         "Invariant mass of 3#pi, Gev/c^{2};m_{3#pi}, [GeV/c^{2}]",
                         N, 3, 3.2);

    // Read data from file
    std::ifstream ifile("m3piJPSI_cut.dat", std::ios::in);
    if (ifile.is_open()) {
        Double_t x;
        while (ifile >> x) 
            hist->Fill(x);
    } else {
        std::cerr << "File m3piJPSI_cut.dat is not open!\n";
        return -1;
    }

    TF1 *fit = new TF1("fit", Convolution, 3.0, 3.2, 6);
    for (auto i : {LineK, LineB, BWScale, GSigma2, GSigma1, GMean2}) {
        fit->SetParameter(i, IndexInit[i]);
        fit->SetParName(i, IndexNames[i]);
    }
    fit->SetParLimits(GSigma1, 1e-5, 0.01);
    fit->SetParLimits(GSigma2, 1e-5, 0.01);

    
    TCanvas *c = new TCanvas("c", "Mass of 3pi", 800, 600);

    c->SetLogy();
    hist->Draw("E");
    hist->Fit(fit); 
    
    TF1 *fBG = new TF1("fBG", "[0] + [1]*x", 3, 3.2);
    fBG->SetParameters(fit->GetParameter(LineB), fit->GetParameter(LineK));
    fBG->SetLineColor(kBlue);
    fBG->SetLineStyle(5);
    fBG->Draw("SAME");

    std::cout << std::endl << "Parameters of Gauss:" << std::endl;
    for (auto i : {GSigma1, GMean2, GSigma2}) 
        std::cout << IndexNames[i] << " = " << fit->GetParameter(i) << std::endl;

    return 0;
}

