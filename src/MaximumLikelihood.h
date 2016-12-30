#ifndef _MaximumLikelihood_include_
#define _MaximumLikelihood_include_
Real calcWoolfson(Real fObs, Real fCalc, Real fVar);
Real calcRice(Real fObs, Real fCalc, Real fVar);
Real calcMLRF0_NoScale(vector<Real> &fObs, vector<Real> &fCalc, vector<Real> &fError, vector<bool> &centric);
Real calcMLRF0(vector<Real> &fObs, vector<Real> &fCalc, vector<Real> &fError, vector<bool> &centric);
Real calcMLRF0(XRayStruct &xray, XRayStruct &expXRay);
Real calcMLTF(XRayStruct &xray, XRayStruct &expXRay, int nmove);
void calcRootMeanIntensity(XRayStruct &xray, Vector &rootMeanI);
void calcIntensityHistogram(XRayStruct &xray, Vector &histogram);
void calcIntensityHistogram(ProteinStruct &Protein, XRayStruct &xray, Vector &histogram);
void calcAverageAndStdDevIntensity(XRayStruct &xray, Vector &rootMeanIntensity, Vector &stdDevIntensity);
Real calcRotationPossibleScore(XRayStruct &xray, XRayStruct &expXRayStruct);
Real calcRotationPossibleScore(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRayStruct);
#endif
