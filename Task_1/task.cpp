#include <cmath>
#include <iostream>

double Func(double x, double arg) { return -2 + x + std::exp(arg * x); }

double FuncDerivative(double x, double arg) { return 1 + arg * std::exp(arg * x); }

int main() 
{
    double arg,
           epsilon,
           firstX0,
           secondX0;

    std::cout << "Please, enter the argument, the precision and both start points: ";
    std::cin >> arg >> epsilon >> firstX0 >> secondX0;


    std::cout << "Arg: " << arg << ", Eps: " << epsilon
              << ", FX0: " << firstX0 << ", SX0: " << secondX0 << std::endl;

    double prevX;
    do 
    {
        prevX = firstX0;
        firstX0 = firstX0 - Func(firstX0, arg) / FuncDerivative(firstX0, arg);
        
    } while(std::abs(firstX0 - prevX) > epsilon);

    std::cout << "First solution: " << firstX0 << std::endl;
    
    if (arg < 0) 
    {
        do 
        {
            prevX = secondX0;
            secondX0 = secondX0 - Func(secondX0, arg) / FuncDerivative(secondX0, arg);
        
        } while(std::abs(secondX0 - prevX) > epsilon);        
        std::cout << "Second solution: " << secondX0 << std::endl;
    }

}
