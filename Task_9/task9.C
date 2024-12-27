#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TMath.h>
#include <iostream>
#include <vector>

void task9() {
    TFile *resultsFile = new TFile("results.root", "UPDATE");

    TTree *candTree = (TTree *)resultsFile->Get("Pi0CandidatesTree");

    Float_t theta_pi0, phi_pi0;
    candTree->SetBranchAddress("PiTheta", &theta_pi0);
    candTree->SetBranchAddress("PiPhi", &phi_pi0);

    const Int_t nbins_theta = 18; // 0-180 градусов, шаг 10 градусов
    const Int_t nbins_phi = 36;   // 0-360 градусов, шаг 10 градусов
    std::vector<Double_t> thetaCenters(nbins_theta), phiCenters(nbins_phi);
    std::vector<Double_t> thetaCounts(nbins_theta, 0), phiCounts(nbins_phi, 0);
    std::vector<Double_t> thetaErrors(nbins_theta, 0), phiErrors(nbins_phi, 0);

    for (Int_t i = 0; i < nbins_theta; ++i)
        thetaCenters[i] = 10.0 * i + 5.0;
    for (Int_t i = 0; i < nbins_phi; ++i)
        phiCenters[i] = 10.0 * i + 5.0;  

    Long64_t nEntries = candTree->GetEntries();
    Long64_t Events = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        candTree->GetEntry(i);
        if (phi_pi0 < 0) phi_pi0 += 2 * TMath::Pi();

        Int_t thetaBin = (Int_t)(theta_pi0 * 180.0 / TMath::Pi() / 10);
        Int_t phiBin   = (Int_t)(phi_pi0 * 180.0 / TMath::Pi() / 10);

        if (thetaBin >= 0 && thetaBin < nbins_theta) thetaCounts[thetaBin]++;
        if (phiBin >= 0 && phiBin < nbins_phi) phiCounts[phiBin]++;
    }

    for (Int_t i = 0; i < nbins_theta; ++i) Events += thetaCounts[i];

    for (Int_t i = 0; i < nbins_theta; ++i)
        thetaErrors[i] = TMath::Sqrt(thetaCounts[i]);
    for (Int_t i = 0; i < nbins_phi; ++i)
        phiErrors[i] = TMath::Sqrt(phiCounts[i]);

    TCanvas *canvas = new TCanvas("c", "Task 9 Results", 1200, 600);
    canvas->Divide(2, 1);

    TGraphErrors *thetaGraph = new TGraphErrors(nbins_theta, thetaCenters.data(), thetaCounts.data(), nullptr, thetaErrors.data());
    thetaGraph->SetTitle("#pi^{0} Candidates from #theta;#theta [deg];N of candidates");
    thetaGraph->SetMarkerStyle(21);
    thetaGraph->SetMarkerColor(kBlue);

    TGraphErrors *phiGraph = new TGraphErrors(nbins_phi, phiCenters.data(), phiCounts.data(), nullptr, phiErrors.data());
    phiGraph->SetTitle("#pi^{0} Candidates from #phi;#phi [deg];N of candidates");
    phiGraph->SetMarkerStyle(21);
    phiGraph->SetMarkerColor(kRed);

    canvas->cd(1);
    thetaGraph->Draw("AP");
    auto *legend1 = new TLegend(0.1, 0.7, 0.48, 0.9);
    legend1->AddEntry(thetaGraph, "#pi_{0} number", "lep");
    legend1->Draw();

    canvas->cd(2);
    phiGraph->Draw("AP");

    TDirectory *task9Dir = (TDirectory *)resultsFile->Get("Task9");
    if (!task9Dir) {
        task9Dir = resultsFile->mkdir("Task9");
    }
    task9Dir->cd();
    thetaGraph->Write("Graph_theta");
    phiGraph->Write("Graph_phi");
    canvas->Write("Task9_Canvas");

    std::cout << "Events: " << Events << std::endl;

    resultsFile->Close();
}
