
#include "McConnaugheySimCoef.hpp"

double McConnaugheySimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return McConnaugheySimilarity(*fp1, *fp2);
}
double McConnaugheySimCoef::ConvertToDistance(double coef)
{
    return 1 - (coef + 1) / 2;
}