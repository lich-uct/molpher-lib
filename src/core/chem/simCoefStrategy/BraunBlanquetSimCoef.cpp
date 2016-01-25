
#include "BraunBlanquetSimCoef.hpp"

double BraunBlanquetSimCoef::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return BraunBlanquetSimilarity(*fp1, *fp2);
}
double BraunBlanquetSimCoef::ConvertToDistance(double coef)
{
    return 1 - coef;
}