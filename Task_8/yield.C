#define yield_cxx
#include <TTree.h>
#include "TFile.h"
#include "yield.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <ROOT/TSeq.hxx>
#include <iostream> 
#include <Math/Vector3D.h>

#define N 0
#define DEBUG 0

int main() 
{
    TChain *chain = new TChain("MyTree");
    chain->Add("newtree.root");
    yield as(chain);
    as.Loop();
    return 0;
}

void yield::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    
    TFile *newfile = new TFile("pi0tree.root", "recreate");

    TH1D *histM = new TH1D("histM", 
                           "Invariant mass of candidates for #pi^{0};M_{#gamma#gamma}, [MeV/c^{2}];",
                           100, 0.1 * 1000, 0.2 * 1000);
    TH1D *histTheta = new TH1D("histTheta",
                               "Angle between #gamma_{1} and #gamma_{2} for all pairs;#eta, [rad]",
                               100, 0, TMath::Pi());

    Int_t countCandidates = 0;
    Int_t candidatesNumber = 0;

    TTree *Pi0CandTree = new TTree("Pi0CandidatesTree", "Tree of pi0 candidates");
    Double_t PiMass[300], PiTheta[300], PiPhi[300];
    Int_t PiN;
    Pi0CandTree->Branch("Pi0Mass", PiMass, "Pi0Mass[300]/D");
    Pi0CandTree->Branch("Pi0Theta", PiTheta, "Pi0Theta[300]/D");
    Pi0CandTree->Branch("Pi0Phi", PiPhi, "Pi0Phi[300]/D");
    Pi0CandTree->Branch("Pi0N", &PiN, "Pi0N/I");
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        candidatesNumber = 0;    
        //////
        /// My code

        Double_t s[300], theta[300];
        Double_t invariantMass = 0;

        TLorentzVector ph1(0,0,0,0), ph2(0,0,0,0);

        PiN = 0;
        for (Int_t i = 0; i < (nph - 1); i++) {
            for (Int_t j = i + 1; j < nph; j++) {
                using namespace TMath;

                ph1.SetPx(eph[i]*Sin(thetaph[i])*Cos(phiph[i]));
                ph1.SetPy(eph[i]*Sin(thetaph[i])*Sin(phiph[i]));
                ph1.SetPz(eph[i]*Cos(thetaph[i]));
                ph1.SetE(eph[i]);

                ph2.SetPx(eph[j]*Sin(thetaph[j])*Cos(phiph[j]));
                ph2.SetPy(eph[j]*Sin(thetaph[j])*Sin(phiph[j]));
                ph2.SetPz(eph[j]*Cos(thetaph[j]));
                ph2.SetE(eph[j]);

                Double_t eta = (ph1).Angle(ph2.Vect());
                invariantMass = (ph1 + ph2).M(); 

                histTheta->Fill((eta));

                if ( (invariantMass >= 0.1) && ( invariantMass <= 0.2 ) ) {
                    auto Pi0Direction = ph1.Vect() + ph2.Vect();

                    s[candidatesNumber] = invariantMass;
                    //theta[candidatesNumber] = eta;

                    PiMass[candidatesNumber] = invariantMass;
                    PiTheta[candidatesNumber] = Pi0Direction.Theta();
                    PiPhi[candidatesNumber] = Pi0Direction.Phi();
                    candidatesNumber++;
                    }
                }
            }
   
        PiN = candidatesNumber;
        Pi0CandTree->Fill();
   
        if ( 2 == candidatesNumber ) {
            histM->Fill(s[0]*1000); histM->Fill(s[1]*1000);
            //histTheta->Fill(theta[0]); histTheta->Fill(theta[1]);
            countCandidates++;
            }
        }

    Pi0CandTree->Write();

    std::cout << "Count candidates: " << countCandidates << std::endl;
    
    TCanvas *c1 = new TCanvas("c1", "#pi^{0} candidates info", 1200, 600);
    c1->Divide(2, 1);

    c1->cd(1);
    histM->Draw();
    histM->Write();

    c1->cd(2);
    histTheta->Draw();
    histTheta->Write();
}
