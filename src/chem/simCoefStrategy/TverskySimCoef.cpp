
#include "TverskySimCoef.hpp"

TverskySimCoef::TverskySimCoef() :
    mA(1.0), mB(1.0)
{
    // no-op
}

TverskySimCoef::TverskySimCoef(double a, double b) :
    mA(a), mB(b)
{
    if (mA < 0.0) {
        mA = 0.0;
    } else if (mA > 1.0) {
        mA = 1.0;
    }

    if (mB < 0.0) {
        mB = 0.0;
    } else if (mB > 1.0) {
        mB = 1.0;
    }
}

double TverskySimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return TverskySimilarity(*fp1, *fp2, mA, mB);
}
double TverskySimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}