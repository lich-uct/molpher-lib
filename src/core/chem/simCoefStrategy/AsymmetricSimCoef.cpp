
#include "AsymmetricSimCoef.hpp"

double AsymmetricSimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return AsymmetricSimilarity(*fp1, *fp2);
}
double AsymmetricSimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}