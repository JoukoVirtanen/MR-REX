#ifndef _EmpiricalMaximumLikelihood_include_
#define _EmpiricalMaximumLikelihood_include_

struct ErrorHistogramStruct
{
	Matrix histogram;
	Real qBinSize, ratioBinSize, minRatio;
};

void readErrorHistogram(string errorHistogramFile, ErrorHistogramStruct &errorHistogram);
void addScatteringToDistribution(XRayStruct &xray, Matrix &distribution);
void randomPlacement(ProteinStruct &Protein);
void calcDistributionForEveryPoint(XRayStruct xray, ProteinStruct Protein, XRayParamStruct &params, Matrix &distribution);
Real calcEmpiricalMaximumLikelihood(Real calc, Real correctedObs, Vector &distribution);
Real calcEmpiricalMaximumLikelihood(Real calc, Real obs, ErrorHistogramStruct &errorHistogram, int qBin, Vector &distribution);
Real calcEmpiricalMaximumLikelihood(XRayStruct &xray, XRayStruct &expXRay, ErrorHistogramStruct &errorHistogram, Matrix &distribution);
Real calcEmpiricalMaximumLikelihood(XRayStruct &xray, XRayStruct &expXRay, ProteinStruct &Protein, XRayParamStruct &params);
#endif
