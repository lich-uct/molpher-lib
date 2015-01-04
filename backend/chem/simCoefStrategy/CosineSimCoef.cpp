
#include "CosineSimCoef.hpp"

double CosineSimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return CosineSimilarity(*fp1, *fp2);
}
double CosineSimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}