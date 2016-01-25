
#include "OnBitSimCoef.hpp"

double OnBitSimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return OnBitSimilarity(*fp1, *fp2);
}
double OnBitSimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}