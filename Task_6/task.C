#include <ROOT/TSeq.hxx>
#include <TF1.h>
#include <TLorentzRotation.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <iostream>

namespace T6 {
    // Mass of particles (MeV) from PDG
    constexpr Double_t FULL_ENERGY = 1020;
    constexpr Double_t MASS_OF_PI = 139.57039;
    constexpr Double_t MASS_OF_K = 497.611;

    // Geometry constants (cm)
    constexpr Double_t RADIUS_OF_DETECTOR = 30;
    constexpr Double_t LEHGTH_OF_DETECTOR = 50;
    const Double_t MAX_LENGTH = TMath::Sqrt(RADIUS_OF_DETECTOR*RADIUS_OF_DETECTOR + LEHGTH_OF_DETECTOR*LEHGTH_OF_DETECTOR/4);

    // Other constants
    constexpr Double_t MIN_ENERGY_FOR_REGISTER_OF_PI = 40;
    constexpr Double_t TWOPI = TMath::Pi();
    constexpr Double_t AVERAGE_MILEAGE = 2.6;
    constexpr Double_t C_TAU = 2.6844;
}

void task()
{
    using namespace TMath;

    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed(10);

    // Create histograms
    TH1D *histKsThetaLab = new TH1D("histKsThetaLab", "", 100, 0, TMath::Pi());
    TH1D *histKsPhiLab   = new TH1D("histKsPhiLab", "", 100, 0, TMath::Pi());

    TH1D *histPiPhiLab = new TH1D("histPiPhiLab", "", 100, 0, Pi());
    TH1D *histPiThetaLab = new TH1D("histPiPhiLab", "", 100, 0, Pi());

    TH1D *histKsLength = new TH1D("histKslength", "", 100, 0, 20);
    TH1D *histAngle = new TH1D("histAngle", "", 100, 0, Pi());

    // For Ks theta distribution
    TF1 *KsThetaDist = new TF1("KsTheta", "sin(x)*sin(x)*sin(x)", 0, Pi());

    UInt_t NTotal = 10000000;
    UInt_t NReg   = 0;
    for (UInt_t i = 0; i < NTotal; i++) {
        // Generate angle Ks in lab system
        Double_t KsTheta = KsThetaDist->GetRandom();
        Double_t KsPhi   = 2 * Pi() * rnd->Rndm() - Pi();

        Double_t KsCosTheta = Cos(KsTheta);
        Double_t KsSinTheta = Sin(KsTheta);
            
        Double_t KsEnergy = (T6::FULL_ENERGY - 2*T6::MASS_OF_PI) * rnd->Rndm() + T6::MASS_OF_PI;
            Double_t KsMomentum = Sqrt(KsEnergy*KsEnergy - T6::MASS_OF_K*T6::MASS_OF_K);

            // Ks in lab system
            TLorentzVector Ks(0, 0, 0, T6::MASS_OF_K);

            // For boost and rotate
            Double_t KsBx = KsMomentum*KsSinTheta*Cos(KsPhi) / KsEnergy;
            Double_t KsBy = KsMomentum*KsSinTheta*Sin(KsPhi) / KsEnergy;
            Double_t KsBz = KsMomentum*KsCosTheta            / KsEnergy;

        // Generate mileage Ks in detector
        Double_t Length = - Log(rnd->Rndm()) * T6::C_TAU * KsMomentum/KsEnergy; 
        TVector3 Decay(Length*KsSinTheta*Cos(KsPhi), 
                       Length*KsSinTheta*Sin(KsPhi),
                       Length*KsCosTheta);

        if (Abs(Decay.Z()) < T6::LEHGTH_OF_DETECTOR && Abs(Decay.X()) < T6::RADIUS_OF_DETECTOR) {
            // Generate Pi+ and Pi-
            Double_t PiCosTheta = 2 * rnd->Rndm() - 1;
            Double_t PiSinTheta = Sqrt(1 - PiCosTheta*PiCosTheta);

            Double_t PiTheta = ACos(PiCosTheta);
            Double_t PiPhi   = 2 * Pi() * rnd->Rndm() - Pi();

            Double_t PiEnergy = (T6::MASS_OF_K - 2*T6::MASS_OF_PI) * rnd->Rndm() + T6::MASS_OF_PI;
            Double_t PiMomentum = Sqrt(PiEnergy*PiEnergy - T6::MASS_OF_PI*T6::MASS_OF_PI);

            Double_t PiPx = PiMomentum*PiSinTheta*Cos(PiPhi);
            Double_t PiPy = PiMomentum*PiSinTheta*Sin(PiPhi);
            Double_t PiPz = PiMomentum*PiCosTheta;

            TLorentzVector PiPlus(PiPx, PiPy, PiPz, PiEnergy);
            TLorentzVector PiMinus(-PiPx, -PiPy, -PiPz, PiEnergy);

            TLorentzRotation T;
            T.Boost(KsBx, KsBy, KsBz);

            Ks = T.VectorMultiplication(Ks);
            PiPlus = T.VectorMultiplication(PiPlus);
            PiMinus = T.VectorMultiplication(PiMinus);

            if (PiPlus.P() < T6::MIN_ENERGY_FOR_REGISTER_OF_PI || PiMinus.P() < T6::MIN_ENERGY_FOR_REGISTER_OF_PI) continue;
           
            NReg++;
            histKsThetaLab->Fill(KsTheta);
            histKsPhiLab->Fill(KsPhi);

            histPiThetaLab->Fill(PiMinus.Theta());
            histPiThetaLab->Fill(PiMinus.Theta());
            histPiPhiLab->Fill(PiPlus.Phi());
            histPiPhiLab->Fill(PiMinus.Phi());

            histKsLength->Fill(Length);
            histAngle->Fill(Ks.Angle(PiPlus.Vect()));
            histAngle->Fill(Ks.Angle(PiMinus.Vect()));
        }
    }

    std::cout << "NReg = " << NReg << std::endl;

    TCanvas *c = new TCanvas("c", "Task 6", 1200, 800);
    c->Divide(2, 3);

    c->cd(1);
    histKsPhiLab->SetTitle("Ks phi angle distribution in lab;#phi_{Ks}, [rad]");
    histKsPhiLab->Draw();

    c->cd(2);
    histKsThetaLab->SetTitle("Ks theta angle distribution in lab;#theta_{Ks}, [rad]");
    histKsThetaLab->Draw();

    c->cd(3);
    histKsLength->SetTitle("Ks mileage distribution;length, [cm]");
    histKsLength->Draw();

    c->cd(4);
    histAngle->SetTitle("Angle between Pi+/Pi- and Ks;#eta_{Pi_{#pm}-Ks}, [rad]");
    histAngle->Draw();
   
    c->cd(5);
    histPiPhiLab->SetTitle("Pi phi angle distribution in lab;#phi_{Pi_{#pm}}, [rad]");
    histPiPhiLab->Draw();

    c->cd(6);
    histPiThetaLab->SetTitle("Pi theta angle distribution in lab;#theta_{Pi_{#pm}}, [rad]");
    histPiThetaLab->Draw();
}
