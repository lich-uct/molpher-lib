
#include "TanimotoSimCoef.hpp"

double TanimotoSimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return TanimotoSimilarity(*fp1, *fp2);
}
double TanimotoSimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}