#include "TBranch.h"
#include <TTree.h>
#include <TFile.h>
#include <ROOT/TSeq.hxx>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>

void task7()
{
    TString filename = "m3pimc.root";
    TString newfilename  = "newtree.root";
 
    TFile *oldfile = TFile::Open(filename);
    TTree *oldtree = (TTree*)oldfile->Get("h10");

    TFile *newfile = new TFile(newfilename, "recreate");
    TTree *temptree = oldtree->CopyTree("Isrfilter==1&&chi2_3p<30");
    temptree->SetBranchStatus("*", 0);
    for (auto activeBranchName : {"nph", "eph", "phiph", "thetaph"})
         temptree->SetBranchStatus(activeBranchName, 1);

    TTree *newtree = temptree->CloneTree();
    newtree->SetName("MyTree");
    
    newtree->Print();
    newtree->Write();

    TH1F *energyHistogram = new TH1F("energyHistogram", "Distibution of phtons energy;energy, [MeV];quintity", 100, 0, 9);
    TF1 *fit = new TF1("fit", "[0]/x^[1]", 0, 9);
    fit->SetParameter(0, 3);
    fit->SetParameter(1, 3);

    TCanvas *canvas = new TCanvas();
    newtree->Draw("eph[]>>energyHistogram", "", "", newtree->GetEntries(), 0);
    energyHistogram->Fit(fit);
    energyHistogram->Write();

    Float_t energyMax = newtree->GetMaximum("eph");
    Float_t energyMin = newtree->GetMinimum("eph");

    std::cout << std::endl
        << "Maximum energy of photons is " << energyMax << " MeV" << std::endl
        << "Minimum energy of photons is " << energyMin << " MeV" << std::endl
        << std::endl;
}   
