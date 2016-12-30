#ifndef _CalcXRay_include_
#define _CalcXRay_include_
void calcScattering(ProteinStruct &Protein, vector<CubeStruct> &Cubes, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &scores);
void calcMultipleScattering(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
void calcScattering(string inputPdbFile, XRayParamStruct &params);
int main(int argc, char *argv[]);
#endif
