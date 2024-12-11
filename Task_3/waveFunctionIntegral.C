#include <TF1.h>
#include <TCanvas.h>

Double_t psif(Double_t *x, Double_t *par)
{
    using namespace TMath;

    Double_t xx = x[0];
    Double_t a  = par[0];
 
    return Power(2./Pi(), 0.25) * 1/Sqrt(a) * Exp(-xx*xx/a/a);
}

Double_t integrandf(Double_t *x, Double_t *par)
{
    using namespace TMath;

    Double_t xx = x[0];
    Double_t a  = par[0];

    TF1 psi("psi", psif, -10*a, 10*a, 1);
    psi.SetParameter(0, a);

    auto U = [](Double_t x) -> Double_t 
        { return TMath::Abs(x) <= 10 ? -0.5 : 0; };
    
    return psi.Eval(xx) * ( -3.8 * psi.Derivative2(xx) + U(xx) * psi.Eval(xx) );
}

Double_t integral(Double_t *x, Double_t *par)
{
    TF1 integrand("integrand", integrandf, -10*x[0], 10*x[0], 1);
    integrand.SetParameter(0, x[0]);

    return integrand.Integral(-10*x[0], 10*x[0]);
}

void waveFunctionIntegral()
{
    TF1 *intgrl = new TF1("intgrl", integral, 5, 45);
    TF1 *function = new TF1("func", psif, -25, 25, 1);
    
    Double_t a_opt = intgrl->GetMinimumX();
    function->SetParameter(0, a_opt);

    TCanvas *c1 = new TCanvas("c1", "Energy minimazer", 1200, 600);
    c1->Divide(2, 1);
    
    c1->cd(1);
    intgrl->SetTitle("<#psi|H|#psi>;a, [A];<#psi|H|#psi>, [eV]");
    intgrl->Draw();

    c1->cd(2);
    TString sum = ""; sum.Form("%.2f", a_opt);
    TString title = "#psi(x) for a = " + sum + ";x, [A]; #psi(x), [A^{-1/2}]";
    function->SetTitle(title);
    function->Draw();
}
