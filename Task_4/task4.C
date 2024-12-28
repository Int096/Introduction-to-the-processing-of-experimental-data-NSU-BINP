#include <TMath.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TGraph2D.h>

#include <iostream>
#include <ctime>
#include <numeric>


Double_t CalculatePi(TRandom3 *rnd, Double_t L, Double_t l, UInt_t N);

Double_t SearchOpt(TRandom3 *rnd, UInt_t N)
{   
    Double_t Lmin = 1, Lmax = 10,
             lmin = 1, lmax = 10;
    Double_t Lstep = 0.01, lstep = 0.01;

    TGraph2D *graph = new TGraph2D((Int_t)( (lmax - lmin) / lstep * (Lmax - Lmin) / Lstep ));

    Double_t Lopt = 0, lopt = 0, piopt = 0, pierr = 100;

    Double_t L = Lmin, l = lmin; 
    Int_t i = 0;
    while (L <= Lmax)
    {
        l = lmin;
        while (l <= lmax)
        {
            Double_t sum = 0; UInt_t K = 100;
            for (int i = 0; i < K; i++)
            {
                std::cout << "";
                sum += CalculatePi(rnd, L, l, N);
            }

            Double_t mean = sum / K;

            graph->SetPoint(i, l, L, mean);
            
            if (TMath::Abs(mean - TMath::Pi()) < pierr)
            {
                pierr = TMath::Abs(mean - TMath::Pi());
                piopt = mean; 
                Lopt = L; 
                lopt = l;
            }
            std::cout << "help: " << L << " " << l << std::endl;
            l += lstep;
            i++;
        }
        L += Lstep;
        i++;
    }

    std::cout
        << "L optimal = " << Lopt
        << "; l optimal = " << lopt 
        << "; pi = " << piopt 
        << std::endl;

    graph->SetTitle("Pi(L, l);l;L;Pi"); 
    graph->Draw();

    return 0;
}

int task4(UInt_t N = 1000) 
{
    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed(time(NULL));
   
    SearchOpt(rnd, N);
    
    Double_t L, l;

    std::cout << "For generate Pi you need to insert N, L and l: ";
    std::cin >> N >> L >> l;
    std::cout << std::endl;

    std::cout 
        << "For L = "  << L 
        << ", l = "    << l 
        << " and N = " << N 
        << " Pi = "    << CalculatePi(rnd, L, l, N)
        << std::endl;

    return 0;
}


Double_t CalculatePi(TRandom3 *rnd, Double_t L, Double_t l, UInt_t N)
{
    Double_t x, phi;
    Double_t n;

    for (Int_t i = 0; i < N; i++)
    {
        x   = L *rnd->Rndm();
        phi = TMath::Pi() * rnd->Rndm();

        if ((x + l/2*TMath::Sin(phi)) > L ) n++;
        if ((x - l/2*TMath::Sin(phi)) < 0.) n++;
    }
    return 2 * N * l / n / L;
}
