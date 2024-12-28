#include <TMath.h>
#include <TRandom3.h>

#include <iostream>
#include <ctime>
#include <numeric>

Double_t CalculatePi(TRandom3 *rnd, Double_t L, Double_t l, Double_t N);
std::pair<Double_t, Double_t> Parametric(TRandom3 *rnd,
                Double_t L_min, Double_t L_max, Double_t L_delta,
                Double_t l_min, Double_t l_delta, Int_t N, Int_t depth);
Double_t mean(std::vector<Double_t> &data);


int main() 
{
    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed(time(NULL));

    Int_t N = 100000;
    Int_t depth = 10;

    std::cout 
        << "For historical data (L = 4 cm, l = 3 cm, N = 530) Pi equals "
        << CalculatePi(rnd, 4, 3, 530)
        << std::endl;

    std::pair<Double_t, Double_t> opt = Parametric(rnd, 1, 1.5, 0.1, 0.5, 0.1, N, depth);

    std::cout 
        << "\nFor optimal parameters (L = " << opt.first << ", l = "
        << opt.second << " and N = " << N << ") Pi = "
        << CalculatePi(rnd, opt.first, opt.second, N)
        << std::endl;
    
    std::cout 
        << "\nFor optimal parameters and historical N (L = " 
        << opt.first  << " cm, l = " 
        << opt.second << " cm, N = 530) Pi equals "
        << CalculatePi(rnd, 4, 3, 530)
        << std::endl;
}

std::pair<Double_t, Double_t> Parametric(TRandom3 *rnd,
                Double_t L_min, Double_t L_max, Double_t L_delta,
                Double_t l_min, Double_t l_delta, Int_t N, Int_t depth)
{
    Double_t L_opt = -100,
             l_opt = -100;
    Double_t eps_opt = 100;

    Double_t L = L_min;
    Double_t l = l_min;
    Double_t pi, eps;

    std::vector<Double_t> data(depth);
    for (Int_t i = 0; i < depth; i++)
        data[i] = 0;

    while(L < L_max)
    {
        l = l_min;
        while(l < L)
        {
            for (Int_t i = 0; i < depth; i++)
            {
                data[i] = CalculatePi(rnd, L, l, N);
                std::cout << "";
            }

            pi = mean(data); 
            eps = TMath::Abs(TMath::Pi() - pi);
            std::cout << pi << std::endl;
            if (eps < eps_opt)
            {
                L_opt = L;
                l_opt = l;
                eps_opt = eps;
            }
            l += l_delta;
        }
        std::cout << "For L = " << L << " optimal parameters equals L_opt = " 
            << L_opt << ", l_opt = " << l_opt << std::endl;
        L += L_delta;
    }

    return std::make_pair(L_opt, l_opt);
}

Double_t CalculatePi(TRandom3 *rnd, Double_t L, Double_t l, Double_t N)
{
    Double_t x, phi;
    Double_t n;

    for (Int_t i = 0; i < N; i++)
    {
        x = L *rnd->Rndm();
        phi = TMath::Pi() * rnd->Rndm() - TMath::Pi() / 2;

        if (x <= l * TMath::Cos(phi))
            n++;
    }
    return 2 * N * l / n / L;
}

Double_t mean(std::vector<Double_t> &data)
{
    Int_t count = data.size();
    Double_t sum = 0;
    for (Int_t i = 0; i < count; i++)
    {
        sum += data[i];
    }

    return sum / count;
}
