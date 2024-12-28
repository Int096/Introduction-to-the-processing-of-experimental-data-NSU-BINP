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
    constexpr Double_t AVERAGE_MILEAGE = 6;
}

void task6(UInt_t N = 100000)
{
    using namespace TMath;

    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed(0);

    TH1D *histKsLength = new TH1D("histKsLength", "histKsLength", 100, 0,     T6::MAX_LENGTH);
    TH1D *histKsTheta  = new TH1D("histKsTheta",  "histKsTheta",  100, 0,     Pi()/2);
    TH1D *histKsPhi    = new TH1D("histKsPhi",    "histKsPhi",    100, 0,     Pi());
    TH1D *histBTW      = new TH1D("histBTW",      "histBTW",      100, -Pi(), Pi());

    TF1 *ThetaDistribution = new TF1("ThetaDistribution", "sin(x)*sin(x)", 0, TMath::Pi());

    // Generate one event
    for (auto i : ROOT::TSeqI(N))
    {
        // Generate Ks meson in moving frame
        Double_t KsEnergy = (T6::FULL_ENERGY - T6::MASS_OF_K) * rnd->Rndm() - T6::MASS_OF_K;
        TLorentzVector *Ks = new TLorentzVector(0., 0., 0., T6::MASS_OF_K);

        // Create Rotate Matrix and Boost 
        TLorentzRotation T;

        Double_t KsMomentum = Sqrt(KsEnergy*KsEnergy - T6::MASS_OF_K*T6::MASS_OF_K);
        Double_t Phi        = 2 * Pi() * rnd->Rndm();
        Double_t CosTheta   = 2 * rnd->Rndm() - 1;
        Double_t SinTheta   = Sqrt( 1 - CosTheta*CosTheta ); 

        Double_t Length = -Log(rnd->Rndm()) * T6::AVERAGE_MILEAGE;
        histKsLength->Fill(Length);
       
        if (kTRUE) {
            Double_t bx = KsMomentum*SinTheta*TMath::Cos(Phi);
            Double_t by = KsMomentum*SinTheta*TMath::Sin(Phi);
            Double_t bz = KsMomentum*CosTheta;

            T.Boost(bx, by, bz);
            
            // Pi -
            Double_t PiPlusEnergy    = (T6::MASS_OF_K - T6::MASS_OF_PI) * rnd->Rndm() - T6::MASS_OF_PI;
            Double_t PiPlusMoment    = (Sqrt(PiPlusEnergy*PiPlusEnergy - T6::MASS_OF_PI*T6::MASS_OF_PI));
 
            Double_t PiPlusPhi       = 2*Pi()*rnd->Rndm();
            Double_t PiPlusCosTheta  = 2 * rnd->Rndm() - 1;
            Double_t PiPlusSinTheta  = Sqrt( 1 - PiPlusCosTheta*PiPlusCosTheta );

            // Pi +
            Double_t PiMinusEnergy   = T6::MASS_OF_PI - PiPlusEnergy;
            Double_t PiMinusMoment   = Sqrt(PiMinusEnergy*PiMinusEnergy - T6::MASS_OF_PI*T6::MASS_OF_PI);

            Double_t PiMinusPhi      = 2*Pi() - PiPlusPhi;
            Double_t PiMinusCosTheta = Cos(Pi()/2 + ACos(PiPlusCosTheta));
            Double_t PiMinusSinTheta = Sqrt( 1 - PiMinusCosTheta*PiMinusCosTheta );

            TLorentzVector *PiPlus = new TLorentzVector(
                                                        PiPlusMoment*PiPlusSinTheta*TMath::Cos(PiPlusPhi),
                                                        PiPlusMoment*PiPlusSinTheta*TMath::Sin(PiPlusPhi),
                                                        PiPlusMoment*PiPlusCosTheta,
                                                        PiPlusEnergy);
            TLorentzVector *PiMinus = new TLorentzVector(
                                                        PiMinusMoment*PiMinusSinTheta*Cos(PiMinusPhi),
                                                        PiMinusMoment*PiMinusSinTheta*Sin(PiMinusPhi),
                                                        PiMinusMoment*PiMinusCosTheta,
                                                        PiMinusEnergy);

            *PiPlus = T.VectorMultiplication(*PiPlus);
            *PiMinus = T.VectorMultiplication(*PiMinus);

            histKsPhi->Fill(Phi);
            histKsLength->Fill(Length);
            histKsTheta->Fill(ACos(CosTheta));
            histBTW->Fill((Ks->Vect()).Angle(PiPlus->Vect()));
            histBTW->Fill((Ks->Vect()).Angle(PiMinus->Vect()));
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Histograms", 1500, 800);
    c1->Divide(2, 3);

    c1->cd(1);
    histKsLength->Draw();

    c1->cd(2);
    histKsTheta->Draw();

    c1->cd(3);
    histKsPhi->Draw();

    c1->cd(4);
    histBTW->Draw();
    
}
