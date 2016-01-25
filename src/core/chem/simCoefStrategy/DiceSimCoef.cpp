
#include "DiceSimCoef.hpp"

double DiceSimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return DiceSimilarity(*fp1, *fp2);
}
double DiceSimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}