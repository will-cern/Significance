#include "Significance.h"
#include "RooStats/NumberCountingUtils.h"
#include "TMath.h"

#include <stdexcept>
#include <math.h>




Significance::Significance() : TObject() { }

double Significance::applyPrescription(double n, double b, double Ze, double Zd, Prescription p) {

    switch(p) {
        case PRESCRIPTION1:
            //Return excess significance if n>=b, otherwise return deficit significance
            return ((n>=b) ? Ze : Zd);
        case PRESCRIPTION2:
            //Like Prescription1, but zeros unintuitive case of negative/positive significance for excesses/deficits
            if (n>=b && Ze > 0) return Ze;
            if (n<b && Zd < 0) return Zd;
            return 0;
        case PRESCRIPTION3:
            //Allows for negative significances for small excesses, because technically they have pValues > 0.5
            if (Ze > 0 && Zd > 0) return Ze;
            if (Zd < 0 && Zd < 0) return Zd;
            return 0;
    }

    throw std::runtime_error("Unknown prescription");


}


double Significance::GaussApproxNoUncert(double n, double b) {
    return (n-b)/sqrt(b);
}

double Significance::GaussApproxWithUncert(double n, double b, double sigma) {
    return (n-b)/sqrt(b + sigma*sigma);
}

double Significance::PoissonNoUncert(double n, double b, Prescription p) {
    if(n==0 && b==0) return 0; //special case

    double Ze = -999;
    if (n>0) Ze = RooStats::PValueToSignificance( ROOT::Math::poisson_cdf_c(n-1,b) ); //tail integral for x >= n
    double Zd = RooStats::PValueToSignificance( ROOT::Math::poisson_cdf_c(n,b) );     //tail integral for x > n

    return applyPrescription( n, b, Ze, Zd, p );


}

double Significance::AsymptoticPoissonPoissonModel(double n, double b, double sigma) {
    double t0 = 0;
    if(sigma<=0.) {
        //use simplified expression ...
        t0 = 2.*( ((n==0)?0:n*log(n/b)) - (n-b) );
    } else {

        double tau = b / (sigma*sigma);
        double m = b*tau;

        double b_hathat = (n+m)/(1.+tau);
        //double s_hat = n - m/tau;
        //double b_hat = m/tau;

        t0 = 2.*( ((n==0)?0:n*log(n/b_hathat)) +
                  ((m==0)?0:m*log(m/(tau*b_hathat))));
    }

    return (n>=b) ? sqrt(t0) : -sqrt(t0);

}

double Significance::AsymptoticPoissonGaussianModel(double n, double b, double sigma) {
    double t0 = 0;

    if(sigma<=0.) {
        //use simplified expression ...
        t0 = 2.*( ((n==0)?0:n*log(n/b)) - (n-b) );
    } else {

        double sigma2 = sigma*sigma;
        double b_hathat = 0.5*( b - sigma2 + sqrt( pow(b-sigma2,2) + 4*n*sigma2 ) );
        //double s_hat = n - m;
        //double b_hat = m;

        t0 = 2.*( ((n==0)?0:n*log(n/b_hathat)) + b_hathat - n + pow(b-b_hathat,2)/(2.*sigma2) );
    }

    return (n>=b) ? sqrt(t0) : -sqrt(t0);

}

double Significance::BinomialModel(double n, double b, double sigma) {
    return (n>0) ? RooStats::NumberCountingUtils::BinomialObsZ(n,b,sigma/b) : -999.; // BinomialObsZ return +inf for n=0, not -inf
}

double Significance::ModifiedBinomialModel(double n, double b, double sigma, Prescription p) {
    double tau = b/(sigma*sigma);
    double rho = 1./(1.+tau); //probability of event being in the aux region
    double m = b*tau; //inferred events in the aux region

    double Ze = (n <= 0) ? -999 : RooStats::PValueToSignificance( TMath::BetaIncomplete(rho,n,m+1) ); //this is exactly what BinomialObsZ method uses
    double Zd = RooStats::PValueToSignificance( TMath::BetaIncomplete(rho,n+1,m) );

    return applyPrescription( n, b, Ze, Zd, p );
}

double Significance::PoissonGammaModel(double n, double b, double sigma, Prescription p) {
    if (n<=0 && b<=0) return 0.;

    double B = b/(sigma*sigma);
    double A = b*B;

    // suitable for excess
    double Ze = 0.;
    if (n == 0) {
        Ze = -999.;
    } else {
        if (A>100) { // need to use logarithms

            /// NB: must work in log-scale otherwise troubles!
            double logProb = A*log(B/(1+B));
            double sum = exp(logProb); // P(n=0)
            for (unsigned u=1; u<=n-1; ++u) { // NOTE u <= n-1
                logProb += log((A+u-1)/(u*(1+B)));
                sum += exp(logProb);
            }
            Ze = RooStats::PValueToSignificance(1. - sum);

        } else {

            // Recursive formula: P(n;A,B) = P(n-1;A,B) (A+n-1)/(n*(1+B))
            double p0 = pow(B/(1+B),A); // P(0;A,B)
            //double nExp = A/B;
            double pLast = p0;
            double sum = p0;
            for (unsigned k=1; k<=n-1; ++k) { // NOTE k <= n-1
                double p = pLast * (A+k-1) / (k*(1+B));
                sum += p;
                pLast = p;
            }
            Ze = RooStats::PValueToSignificance(1. - sum);
        }
    }

    // suitable for deficit
    double Zd = 0.;
    if (A>100) { // need to use logarithms

        /// NB: must work in log-scale otherwise troubles!
        double logProb = A*log(B/(1+B));
        double sum = exp(logProb); // P(n=0)
        for (unsigned u=1; u<=n; ++u) { // NOTE u <= n
            logProb += log((A+u-1)/(u*(1+B)));
            sum += exp(logProb);
        }
        Zd = RooStats::PValueToSignificance(1. - sum);

    } else {

        // Recursive formula: P(n;A,B) = P(n-1;A,B) (A+n-1)/(n*(1+B))
        double p0 = pow(B/(1+B),A); // P(0;A,B)
        //double nExp = A/B;
        double pLast = p0;
        double sum = p0;
        for (unsigned k=1; k<=n; ++k) { // NOTE k <= n
            double p = pLast * (A+k-1) / (k*(1+B));
            sum += p;
            pLast = p;
        }
        Zd = RooStats::PValueToSignificance(1. - sum);
    }

    return applyPrescription( n, b, Ze, Zd, p );

}


