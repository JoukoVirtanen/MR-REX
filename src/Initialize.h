#ifndef _Initialize_include_
#define _Initialize_include_
# include "LatticeStruct.h"
# include "XRayStruct.h"
# include "Params.h"
void initializeLattice(LatticeStruct &lattice, XRayParamStruct &params);
void initializeSegments(ProteinStruct &Protein, XRayParamStruct &params);
void initializePatterson(PattersonStruct &patterson);
void initializeDWF(Matrix &dwf, vector<PosStruct> &q, ProteinStruct &Protein);
void initializeFA(Matrix &f, Matrix &dwf, Matrix &fa, ProteinStruct &Protein);
void setExcludedVolumeRadii(vector<Real> &excludedVolumeRadii, XRayParamStruct &params);
void initializeScatteringFactors(XRayStruct &xray, ProteinStruct &Protein, XRayParamStruct &params);
void initializeCentric(XRayStruct &xray, vector<SymmetryOperationStruct> &symmetryOperators);
void initializeCentric(XRayStruct &xray, ProteinStruct &Protein);
void readCorrectionFile(string correctionFile, XRayStruct &xray);
void calcUnitCellVectors(ProteinStruct &Protein, XRayParamStruct &params);
void initializeContinuousMiller(XRayStruct &xray, XRayParamStruct &params);
void initializeContinuous(XRayStruct &xray, XRayParamStruct &params);
void initializeMiller(XRayStruct &xray, XRayParamStruct &params);
void initializeXRay(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
void initializeXRay(XRayStruct &xray, XRayStruct &expXRay, ProteinStruct &Protein, XRayParamStruct &params);
void initializeRand(int randomSeed);
void initialize(LatticeStruct &lattice, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
#endif
