#include <TH1.h>
#include <TMath.h>
#include <TF1.h>
#include <ROOT/TSeq.hxx>
#include <TVirtualFFT.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>

#define DEBUG 1
#define CD3P 0


void FillHist(TH1D *hist, TString filename)
{
    std::ifstream ifile(filename, std::ios::in); 
    if (ifile.is_open())
    {
        Double_t n = 0;
        Double_t x;
        Int_t i = 1;
        while (ifile >> x)
        {
            hist->SetBinContent(i, x); 
            i++;
        }
    } else 
    {
        std::cerr << "File " << filename << " is not open!!!";
        exit(-1);
    }
}

void fft()
{
    Int_t N = 500;
    TH1D *signalH = new TH1D("sighanH", "Signal;t, [ns];U, [mV]", N+1, 0, 500);
    FillHist(signalH, "dataFFT.dat");
    
    TCanvas *fftc = new TCanvas("fftc", "Fast Fourier Transform", 1200, 1000);
    fftc->Divide(2, 2);

    // Исходная гистограмма
    fftc->cd(1);
    signalH->Draw();

    fftc->cd(2);
    TH1 *hm = nullptr;
    TVirtualFFT::SetTransform(nullptr);
    hm = signalH->FFT(hm, "MAG");
    hm->SetTitle("Magnitude of the 1st transform");
    hm->Scale(1./N);
    hm->Draw();
    
    fftc->cd(3);
    TH1 *hp = nullptr;
    hp = signalH->FFT(hp, "PH");
    hp->SetTitle("Phase of the 1st transform");
    hp->Draw();
    hp->SetStats(kFALSE);

    fftc->cd(4);
    Double_t re, im; 
    TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
    Double_t *re_full = new Double_t[N];
    Double_t *im_full = new Double_t[N];
    fft->GetPointsComplex(re_full, im_full);
    TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2R M K");
    
    // Filter
    Int_t p = 8;
    for (Int_t i = 1; i <= N; i++) {
        #ifdef DEBUG
        std::cout <<
            TString::Format("j = %d, Mag = %lf", i, TMath::Sqrt(re_full[i]*re_full[i] + im_full[i]*im_full[i])) << 
        std::endl;
        #endif
        if (i > p) {
            re_full[i] = 0;
            im_full[i] = 0;
        }
    }
    
    fft_back->SetPointsComplex(re_full, im_full);
    fft_back->Transform();

    TH1 *hb = nullptr;
    hb = TH1::TransformHisto(fft_back, hb, "Re");
    hb->SetTitle("The backward transform result;t, [ns];U, [mV]");
    hb->Scale(1./N);
    hb->Draw("HIST SAME C");
    hb->SetStats(kFALSE);

    Double_t maxSignal = hb->GetMaximum();
    Double_t startSignalT = 0, startSignalU = 0;
    std::cout << "Maximum of signal equal " << maxSignal << " [mV]" << std::endl;
    for (Int_t i = 1; i <= N; i++) {
        Double_t content = hb->GetBinContent(i);
        if (content >= 0.2*maxSignal) {
            startSignalT = i;
            startSignalU = content;
            break;
        }
    }
    std::cout << "Start of signal: T = " << startSignalT
        << " ns; U = " << startSignalU << " mV" << std::endl;

    auto f = [=](Double_t *x, Double_t *par) { return startSignalU; };
    TF1 *sigSt = new TF1("sigSt", f, 0, 500, 0);
    sigSt->Draw("SAME");

    delete fft_back;
    fft_back = nullptr;

    delete [] re_full;
    delete [] im_full;
}
