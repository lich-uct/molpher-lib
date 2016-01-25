
#include "KulczynskiSimCoef.hpp"

double KulczynskiSimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return KulczynskiSimilarity(*fp1, *fp2);
}
double KulczynskiSimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}