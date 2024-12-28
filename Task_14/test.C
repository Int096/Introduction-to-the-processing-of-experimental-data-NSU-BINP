#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>

void test()
{
    TString filename = "test.dat";
    std::ofstream ofile(filename, std::ios::out);
    for (int i = 0; i < 500; i++) {
        Double_t x =  TMath::Gaus(i, 250, 50);
        ofile << x << std::endl; 
        std::cout << x << std::endl;
    }
    ofile.close();
}
