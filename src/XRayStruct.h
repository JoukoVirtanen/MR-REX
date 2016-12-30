#ifndef _XRayStruct_included_
#define _XRayStruct_included_

# include "../../LibraryFiles/MathStructures.h"
# include "LatticeStruct.h"
# include "Params.h"

const int CLASH_CC=0;
const int CLASH_CS=1;
const int CLASH_SS=2;
const int RFACTOR=3;
const int SQR_DIFF=4;
const int PEARSON=5;
const int PRODUCT=6;
const int ROTATION=7;
const int LOG_SQR_DIFF=8;
const int PHASE=9;
const int DENSITY_CC=10;
const int MLTF=11;

const int NUM_SCORE_TYPES=12;

struct PattersonStruct
{
	vector<PosStruct> u;
	vector<Real> real, imag;
};

struct OccupancyCorrectionStruct
{
	Real occupancyBinSize;
	Real correctionBinSize;
	Matrix correction;
};

struct XRayStruct
{
	//bool boolVerbose;
	vector<bool> centric, free;
	vector< vector< vector<int> > > updatedContinuous;
	//vector< vector< vector<bool> > > updatedContinuous;
	//vector< vector< vector< vector<bool> > > > updatedContinuous;
	Vector i, imag, real, waterReal, waterImag, phase, error, qMag;
	Vector excludedReal, excludedImag;
	Vector F, Ferror, Ffix, Fmove;
	Vector realBlock, imagBlock;
	Vector correction;
	Vector cosx1d, cosy1d, cosz1d;
	Vector sinx1d, siny1d, sinz1d;
	Vector score, stdDev;
	vector<PosStruct> q, miller;
	Real a, b, c, alpha, beta, gamma, correctionBinSize;
	Real hBinSize, kBinSize, lBinSize;  //For continuous Miller indexes
	Real clashScore;
	int maxH, maxK, maxL;
	int minH, minK, minL;
	int maxContinuousH, maxContinuousK, maxContinuousL;
	vector<int> scoreIndex;
	Array3D realContinuous, imagContinuous;	//Continuous Miller indexes are used for fast rotation and translation
	Array4D realContinuousSegment, imagContinuousSegment;
	//vector< vector< vector<PosStruct> > > continuousMiller, qContinuous;
	//vector< vector< vector<PosStruct> > > realContinuousDer;
	//vector< vector< vector<PosStruct> > > imagContinuousDer;
	//vector< vector< vector< vector< vector<Real> > > > > realContinuousHes;
	//vector< vector< vector< vector< vector<Real> > > > > imagContinuousHes;
	Real ***i3d, ***imag3d, ***real3d, ***phase3d;
	vector< vector< vector<PosStruct> > > q3d, miller3d;
	Matrix dwf, f, fa, complexReal, complexImag, symReal, symImag;	//Debye-Waller Factor, form factor, dwf*f
	Matrix realBlockTranslated, imagBlockTranslated;
	Matrix symTempReal, symTempImag;
	Matrix cosx, cosy, cosz;
	Matrix sinx, siny, sinz;
	Matrix dwfx, dwfy, dwfz;
	Matrix matTrig1, matTrig2;
	Matrix distribution;
	Array3D complexSymReal, complexSymImag;
	Array3D complexSymTempReal, complexSymTempImag;
	PattersonStruct patterson;
	LookUpStruct lookUp;
	vector<ProteinStruct> Proteins;
	OccupancyCorrectionStruct OccupancyCorrection;
	XRayStruct& operator = (const XRayStruct &xray);
};

struct MRStruct
{
        ProteinStruct Protein;
	vector<ProteinStruct> Proteins;
        XRayStruct xray, expXRay;
        LatticeStruct lattice;
        XRayParamStruct params;
};

struct EstimateScoreStruct
{
	ProteinStruct Protein, referenceProtein;
	Matrix cartGrad; //Cartesian gradient. How much does the X-ray score
			 //change when an atom is moved.
};

struct SolventStruct
{
	XRayStruct xray, expXRay, xraySol, xrayCorrected;
	XRayParamStruct params;
};

#endif
