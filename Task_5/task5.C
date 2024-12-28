#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <ROOT/TSeq.hxx>
#include <iostream>

void GeneratorOne(TRandom3*);
Double_t functionTwo(Double_t);
Double_t EmpErrors(Double_t, Double_t);

Double_t RejectionMethod(TRandom3*, UInt_t);
Double_t AverageMethod(TRandom3*, UInt_t);
Double_t ImportanceSampling(TRandom3*, UInt_t);

void PrintHistograms(TRandom3*, UInt_t, UInt_t, Double_t);

void task5(UInt_t N = 100000, UInt_t count = 1000)
{
    TH1::AddDirectory(false);
    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed(time(NULL));

    GeneratorOne(rnd);

    constexpr Double_t exactMeaning = 6 - 16/TMath::E();
    Double_t neumann   = RejectionMethod(rnd, N);
    Double_t isampling = ImportanceSampling(rnd, N);
    Double_t average   = AverageMethod(rnd, N);

    std::cout << std::fixed;
    std::cout.precision(7);

    std::cout << "\t\tFor N = " << N << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Methods name\t|  Result\t|  EmpError"  << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Exactly     \t|  " << exactMeaning << "\t|  "
        << "---" << std::endl;
    std::cout << "Rejection   \t|  " << neumann << "\t|  " 
        << EmpErrors(exactMeaning, neumann) <<  std::endl;
    std::cout << "Average     \t|  " << average << "\t|  " 
        << EmpErrors(exactMeaning, average) << std::endl;
    std::cout << "Importance  \t|  " << isampling << "\t|  " 
        << EmpErrors(exactMeaning, isampling) << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    PrintHistograms(rnd, N, count, exactMeaning);
}

// Генерация и вывод первой гистограммы
void GeneratorOne(TRandom3 *rnd)
{
    TCanvas *canvas = new TCanvas("c1", "Histogram for one formula", 800, 600);
    TH1D *histogramOne = new TH1D("histogramOne", "Distribution", 100, 0, 2.5);

    Int_t N = 100000;
    for (auto i : ROOT::TSeqI(N))
        histogramOne->Fill(TMath::Sqrt(-TMath::Log(rnd->Rndm())));

    TString Nstr = ""; Nstr.Form("%d", N);
    TString title = "Distribution 2#upoint x#upoint exp(-x^{2}) for N = " + Nstr;
    histogramOne->SetTitle(title);

    histogramOne->Draw("E");
    canvas->SaveAs("One.png");
}

Double_t functionTwo(Double_t x) { return TMath::Exp(-x) * x*x*x; }
Double_t EmpErrors(Double_t trueValue, Double_t calculatedValue) 
{ return TMath::Abs(trueValue - calculatedValue); };

// Оценка интеграла методом выбраковки Неймана
Double_t RejectionMethod(TRandom3 *rnd, const UInt_t N)
{
    Double_t a = 0, b = 1, C = 0.9;
    Double_t ksi, eta; Double_t n = 0;
    for (auto i : ROOT::TSeqI(N))
    {
        ksi = (b - a) * rnd->Rndm() + a; 
        eta = C * rnd->Rndm();
        if (eta <= functionTwo(ksi))
            n++;
    }
    return n/N * (b-a)*C;
}

// Подсчет методом среднего
Double_t AverageMethod(TRandom3 *rnd, const UInt_t N)
{
    Double_t a = 0, b = 1, R = 0;
    for (auto i : ROOT::TSeqI(N))
        R += (b - a) * functionTwo((b - a) * rnd->Rndm() + a);

    return R/N;
}

// Оценка методом существенной выборки (выборка по значимости)
// g(x) = 4 * x^3
Double_t ImportanceSampling(TRandom3 *rnd, const UInt_t N)
{
    Double_t a = 0, b = 1, R = 0;
    Double_t r;
    for (auto i : ROOT::TSeqI(N))
    {
        r = TMath::Sqrt(TMath::Sqrt(rnd->Rndm()));
        R += functionTwo(r) / 4 / r/r/r;
    }
    return R/N;
}

void PrintHistograms(TRandom3 *rnd, UInt_t N, UInt_t count, Double_t exactMeaning)
{
    TCanvas *canvas2 = new TCanvas("c2", "Errors Histograms", 800, 600);
    canvas2->Divide(2, 2);
   
    Double_t min = exactMeaning - 1/TMath::Sqrt(N);
    Double_t max = exactMeaning + 1/TMath::Sqrt(N);

    TH1D *RejectionHistogram  = new TH1D("rejectionHist", "Rejection Method;result value;quantity", 100, min, max);
    TH1D *AverageHistogram    = new TH1D("averageHist", "Average Value Method;result value;quantity", 100, min, max);
    TH1D *ImportanceHistogram = new TH1D("importanceHist", "Importance Sampling Method;result value;quantity", 100, min, max);

    for (auto i : ROOT::TSeqI(count))
    {
        RejectionHistogram->Fill(RejectionMethod(rnd, N));   
        AverageHistogram->Fill(AverageMethod(rnd, N));   
        ImportanceHistogram->Fill(ImportanceSampling(rnd, N));   
    }
    canvas2->cd(1);
    RejectionHistogram->Draw("E");

    canvas2->cd(2);
    AverageHistogram->Draw("E");

    canvas2->cd(3);
    ImportanceHistogram->Draw("E");

    canvas2->SaveAs("Two.png");
}
