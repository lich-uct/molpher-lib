
#include "SokalSimCoef.hpp"

double SokalSimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return SokalSimilarity(*fp1, *fp2);
}
double SokalSimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}