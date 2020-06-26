#pragma once

#include "TObject.h"

#include "RooStats/RooStatsUtils.h"


class Significance : public TObject {


public:
    enum Prescription {
        PRESCRIPTION1 = 1,
        PRESCRIPTION2 = 2,
        PRESCRIPTION3 = 3
    };

    static const Prescription RECOMMENDED = PRESCRIPTION3;

    Significance();

    /**
     * Compute the recommended significance
     */
    static double Recommended(double n, double b, double sigma) { return AsymptoticPoissonPoissonModel(n,b,sigma); }


    static double GaussApproxNoUncert(double n, double b);

    static double GaussApproxWithUncert(double n, double b, double sigma);

    static double PoissonNoUncert(double n, double b, Prescription p = RECOMMENDED);

    static double AsymptoticPoissonPoissonModel(double n, double b, double sigma);

    static double AsymptoticPoissonGaussianModel(double n, double b, double sigma);

    static double BinomialModel(double n, double b, double sigma);

    static double ModifiedBinomialModel(double n, double b, double sigma, Prescription p = RECOMMENDED);

    static double PoissonGammaModel(double n, double b, double sigma, Prescription p = RECOMMENDED);


    static double applyPrescription(double n, double b, double Ze, double Zd, Prescription p);


    ClassDef(Significance,1);



};



