
#include "AllBitSimCoef.hpp"

double AllBitSimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return AllBitSimilarity(*fp1, *fp2);
}
double AllBitSimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}