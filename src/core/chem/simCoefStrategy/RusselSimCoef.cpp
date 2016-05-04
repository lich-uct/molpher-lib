
#include "RusselSimCoef.hpp"

double RusselSimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return RusselSimilarity(*fp1, *fp2);
}
double RusselSimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}
