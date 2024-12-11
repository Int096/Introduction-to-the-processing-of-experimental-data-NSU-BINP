#include <TMath.h> 
#include <TRandom3.h>

#include <iostream>
#include <ctime>

#define DEBUG 0

Double_t CalculatePi(TRandom3*, Double_t, Double_t, UInt_t);
void SearchOpt(TRandom3*, UInt_t, Double_t, Double_t);

void buffonsMethod(UInt_t K = 1000) 
{
    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed(std::time(NULL));

    Double_t N, L, l;
    std::cout << "For generate Pi you need to insert N, L and l: ";
    std::cin >> N >> L >> l;
    std::cout << std::endl;

    std::cout 
        << "For L = "   << L
        << ", l = "     << l
        << " and N = "  << N
        << " Pi = "     << CalculatePi(rnd, L, l, N)
        << std::endl;

    std::cout << std::endl << std::endl;
    SearchOpt(rnd, K, 0.5, 6);
    std::cout << std::endl << std::endl;
}

void SearchOpt(TRandom3 *rnd, UInt_t N, Double_t step, Double_t L)
{
    Double_t lopt, Lopt, err = 1e10, pi, l;
    UInt_t K = 10000;

    
    for (l = 1; l <= 1.5*L; l += step) {
        Double_t sum = 0;
        for (Int_t i = 0; i < K; i++) {
            std::cout << "";
            sum += CalculatePi(rnd, L, l, N);
        }

        Double_t mean = sum / K;

        #ifdef DEBUG
        std::cout 
            << "mean = " << mean
            << " L = " << L
            << " l = " << l
            << std::endl;
        #endif

        if ( TMath::Abs(mean - TMath::Pi()) < err ) {
            err = TMath::Abs(mean - TMath::Pi());

            pi = mean;
            Lopt = L;
            lopt = l;
        }
    }
    TString result 
        = TString::Format("\nOptimal parameters: L = %.2lf, l = %.2lf, N = %u and pi = %.4lf\n",
                          L, l, N, pi);
    std::cout << result << std::endl;
}

Double_t CalculatePi(TRandom3 *rnd, Double_t L, Double_t l, UInt_t N)
{
    Double_t x, phi, n;
    for (Int_t i = 0; i < N; i++) {
        x   = L * rnd->Rndm(); 
        phi = TMath::Pi() * rnd->Rndm();

        if ((x + l/2*TMath::Sin(phi)) > L) n++;
        if ((x - l/2*TMath::Sin(phi)) < 0) n++;
    }
    return 2 * N * l / n / L;
}
