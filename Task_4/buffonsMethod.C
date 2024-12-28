#include <TRandom3.h>
#include <TMath.h>

#include <iostream>

// Function for Calculate Pi
//
Double_t Calculate(UInt_t N, Double_t L, Double_t l, TRandom3 *rnd)
{
    Double_t x = 0, phi = 0;
    UInt_t n = 0;

    for (UInt_t i = 0; i < N; i++)
    {
        x   = L/2.           * rnd->Rndm();
        phi = TMath::Pi()/2. * rnd->Rndm();
        if (x <= l/2. * TMath::Sin(phi)) n++;
    }

    std::cout << "N = " << N
        << " n = " << n << std::endl;

    if (l <= L)
        return 2. * l/L * N/n;
    else
        return 2.*N/n*((l-TMath::Sqrt(l*l-L*L))/L + TMath::ACos((Double_t)L/l));
}

void buffonsMethod()
{
    using namespace TMath;

    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed(0);
    
    std::cout << ""
}
