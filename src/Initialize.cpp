# include <iostream>
# include <vector>

using namespace std;

# include "../../LibraryFiles/TypeDef.h"
# include "../../LibraryFiles/Structures.h"
# include "../../LibraryFiles/AtomUtils.h"
# include "../../LibraryFiles/VectorManip.h"
# include "../../LibraryFiles/MathUtils.h"
# include "../../LibraryFiles/PrintPdb.h"
# include "LatticeStruct.h"
# include "Params.h"
# include "XRay.h"
# include "FlexibleEditing.h"

void initializeLattice(LatticeStruct &lattice, XRayParamStruct &params)
{
	cout <<"In initializeLattice "<<params.XCubeLength<<endl;
	lattice.xCubeLength=params.a/params.NumQxValues;
	lattice.yCubeLength=params.b/params.NumQyValues;
	lattice.zCubeLength=params.c/params.NumQzValues;
	lattice.a=params.a;
	lattice.b=params.b;
	lattice.c=params.c;
	cout <<"a= "<<lattice.a<<" params.XCubeLength= "<<params.XCubeLength
	<<" params.XCubes= "<<params.XCubes<<endl;
	lattice.firstPos.Pos.resize(3);
	//lattice.firstPos.Pos[X]=-params.a*0.5;
	//lattice.firstPos.Pos[Y]=-params.b*0.5;
	//lattice.firstPos.Pos[Z]=-params.c*0.5;
	lattice.firstPos.Pos[X]=0;
	lattice.firstPos.Pos[Y]=0;
	lattice.firstPos.Pos[Z]=0;
	if (lattice.density.size()==0) 
	{
		Safe3DAlloc(lattice.density, params.NumQxValues, params.NumQyValues, params.NumQzValues, "lattice.density");
	}
	cout <<"NumQxValues= "<<params.NumQxValues<<endl;
	cout <<"NumQyValues= "<<params.NumQyValues<<endl;
	cout <<"NumQzValues= "<<params.NumQzValues<<endl;
	Safe3DAlloc(lattice.real, params.NumQxValues, params.NumQyValues, params.NumQzValues, "lattice.real");
	Safe3DAlloc(lattice.imag, params.NumQxValues, params.NumQyValues, params.NumQzValues, "lattice.imag");
	cout <<"Leaving initializeLattice"<<endl;
}

void initializePatterson(PattersonStruct &patterson)
{
	int npoint=16;
	Real xmin, ymin, zmin;
	Real xmax, ymax, zmax;
	Real xinc, yinc, zinc;
	PosStruct tempPos;

	xmin=-10.0;
	ymin=-10.0;
	zmin=-10.0;
	xmax=-xmin;
	ymax=-ymin;
	zmax=-zmin;

	xinc=(xmax-xmin)/Real(npoint);
	yinc=(ymax-ymin)/Real(npoint);
	zinc=(zmax-zmin)/Real(npoint);
	tempPos.Pos.resize(3);
	for (Real xgrid=xmin;xgrid<=xmax;xgrid+=xinc)
	{
		tempPos.Pos[X]=xgrid;	
		for (Real ygrid=ymin;ygrid<=ymax;ygrid+=yinc)
		{
			tempPos.Pos[Y]=ygrid;	
			for (Real zgrid=zmin;zgrid<=zmax;zgrid+=zinc)
			{
				tempPos.Pos[Z]=zgrid;	
				SafePushBack(patterson.u, tempPos, "patterson.u");
			}
		}
	}
}


void initializeDWF(Matrix &dwf, vector<PosStruct> &q, ProteinStruct &Protein)
{
	int npoint=q.size();
	int natom=Protein.Atoms.size();
	cout <<"In initializeDWF"<<endl;
	printSymmetryOperations(Protein.symmetryOperations);
	cout <<"npoint= "<<npoint<<" natom= "<<natom<<endl;	
	Safe2DAlloc(dwf, npoint, natom, "dwf");
	cout <<"After SafeAlloc"<<endl;
	for (int i=0;i<npoint;i++)
	{
		for (int j=0;j<natom;j++)
		{
			Protein.Atoms[j].BFactor=1.0;
			dwf[i][j]=calcDWF(q[i], Protein.Atoms[j]);
		}
	}
}

void initializeFA(Matrix &f, Matrix &dwf, Matrix &fa, ProteinStruct &Protein)
{
	int npoint=f.size();
	int natom=Protein.Atoms.size(), atomID;

	Safe2DAlloc(fa, npoint, natom, "fa");
	print2DVectorSize(fa, "fa");
	cout <<"In initializeFZ"<<endl;
	printSymmetryOperations(Protein.symmetryOperations);
	
	for (int i=0;i<npoint;i++)
	{
		for (int j=0;j<natom;j++)
		{
			atomID=Protein.Atoms[j].atomid;
			fa[i][j]=f[atomID][i]*dwf[i][j];
		}
	}
}

void setExcludedVolumeRadii(vector<Real> &excludedVolumeRadii, XRayParamStruct &params)
{
	SafeAlloc(excludedVolumeRadii, NumAtomTypes, "excludedVolumeRadii");
	excludedVolumeRadii[HYDROGEN]=params.ExcludedVolumeRadiusHydrogen;
	excludedVolumeRadii[CARBON]=params.ExcludedVolumeRadiusCarbon;
	excludedVolumeRadii[NITROGEN]=params.ExcludedVolumeRadiusNitrogen;
	excludedVolumeRadii[OXYGEN]=params.ExcludedVolumeRadiusOxygen;
	excludedVolumeRadii[SULFUR]=params.ExcludedVolumeRadiusSulfur;
}

void initializeScatteringFactors(XRayStruct &xray, ProteinStruct &Protein, XRayParamStruct &params)
{
	cout <<"In initializeScatteringFactors"<<endl;
	solventCorrectedScatteringFactor(xray.f, xray.q, params);
	cout <<"Before initializeDWF"<<endl;
	cout <<"UseExpLookUp= "<<params.UseExpLookUp<<endl;
	initializeCosLookUp(xray.lookUp);
	initializeSinLookUp(xray.lookUp);
	if (params.UseExpLookUp)
	{
		cout <<"Initialize ExpLookUp"<<endl;
		//initializeExpLookUp(xray.lookUp);
		initializeNegativeExpLookUp(xray.lookUp);
		cout <<"After initializeExpLookUp"<<endl;
		//initializeCosLookUp(xray.lookUp);
		//initializeSinLookUp(xray.lookUp);
	}
	else
	{
		initializeDWF(xray.dwf, xray.q, Protein);
		cout <<"Before initializeFA"<<endl;
		initializeFA(xray.f, xray.dwf, xray.fa, Protein);
		cout <<"Leaving initializeScatteringFactors"<<endl;
	}
	cout <<"Leaving initializeScatteringFactors"<<endl;
}

void initializeCentric(XRayStruct &xray, vector<SymmetryOperationStruct> &symmetryOperators)
{
	int npoint=xray.miller.size();
	int noperations=symmetryOperators.size();
	Real h, k, l;

	SafeAlloc(xray.centric, false, npoint);
	for (int i=0;i<npoint;i++)
	{
		xray.centric[i]=false;
		for (int j=0;j<noperations;j++)
		{
			h=xray.miller[i].Pos[X]*symmetryOperators[j].rotationMat[X][X];
			k=xray.miller[i].Pos[Y]*symmetryOperators[j].rotationMat[Y][Y];
			l=xray.miller[i].Pos[Z]*symmetryOperators[j].rotationMat[Z][Z];
			if (h==-xray.miller[i].Pos[X] && k==-xray.miller[i].Pos[Y] && l==-xray.miller[i].Pos[Z])
			{
				xray.centric[i]=true;
			}
		}
	}
}

void initializeCentric(XRayStruct &xray, ProteinStruct &Protein)
{
	initializeCentric(xray, Protein.symmetryOperations);
}

void readCorrectionFile(string correctionFile, Real &correctionBinSize, Vector &correction)
{
	vector<string> lines;
	int nline;
	Real temp;

	readLines(correctionFile, lines, "correctionFile");
	nline=lines.size();
	
	correctionBinSize=toReal(GetWord2(lines[0], 2));
	for (int i=1;i<nline;i++)
	{
		temp=toReal(GetWord2(lines[i], 3));
		SafePushBack(correction, temp, "correction");
	}

}

void readCorrectionFile(string correctionFile, XRayStruct &xray)
{
	readCorrectionFile(correctionFile, xray.correctionBinSize, xray.correction);
}

void readOccupancyCorrectionFile(string correctionFile, XRayStruct &xray)
{
	string file;
	vector<string> lines;
	int nline;
	Vector temp;
	readLines(correctionFile, lines, "correctionFile");
	nline=lines.size();

	xray.OccupancyCorrection.occupancyBinSize=toReal(GetWord2(lines[0], 2));
	printLocation(__LINE__, __FILE__);
	cout <<"occupancyBinSize= "<<xray.OccupancyCorrection.occupancyBinSize<<endl;
	for (int i=1;i<nline;i++)
	{
		file=lines[i];
		readCorrectionFile(file, xray.OccupancyCorrection.correctionBinSize, temp);
		SafePushBack(xray.OccupancyCorrection.correction, temp, "temp");
	}
}

void initializeSegments(ProteinStruct &Protein, XRayParamStruct &params)
{
	if (params.SegmentsFile!="")
	{
		readSegmentsFile(params.SegmentsFile, Protein.segmentStart, Protein.segmentEnd);	
		residueSegmentsToAtomSegments(Protein.Atoms, Protein.segmentStart, Protein.segmentEnd);
	}
}

void calcUnitCellVectors(ProteinStruct &Protein, XRayParamStruct &params)
{
	//if (Protein.XBoxLength==0) Protein.XBoxLength=params.a;
	//if (Protein.YBoxLength==0) Protein.YBoxLength=params.b;
	//if (Protein.ZBoxLength==0) Protein.ZBoxLength=params.c;
	//if (Protein.alpha==0) Protein.alpha=params.alpha;
	//if (Protein.beta==0) Protein.beta=params.beta;
	//if (Protein.gamma==0) Protein.gamma=params.gamma;
	if (params.a!=0) Protein.XBoxLength=params.a;
	if (params.b!=0) Protein.YBoxLength=params.b;
	if (params.c!=0) Protein.ZBoxLength=params.c;
	if (params.alpha!=0) Protein.alpha=params.alpha;
	if (params.beta!=0) Protein.beta=params.beta;
	if (params.gamma!=0) Protein.gamma=params.gamma;
	calcUnitCellVectors(Protein);
}

void initializeContinuousMiller(XRayStruct &xray, XRayParamStruct &params)
{
	//For fast rotation translation function
	int hBins, kBins, lBins;
	PosStruct tempPos;

	tempPos.Pos.resize(3);
	hBins=params.ContinuousHValues;
	kBins=params.ContinuousKValues;
	lBins=params.ContinuousLValues;
	findMaxMillerForFastRotation(xray, params);
	xray.hBinSize=Real(xray.maxContinuousH)/Real(hBins);
	xray.kBinSize=Real(xray.maxContinuousK)/Real(kBins);
	xray.lBinSize=Real(xray.maxContinuousL)/Real(lBins);
/*
	Safe3DAlloc(xray.continuousMiller, tempPos, hBins, kBins*2, lBins*2, "continuousMiller");

	for (int h=0;h<hBins;h++)
	{
		//for (int k=-kBins+1;k<kBins;k++)
		for (int k=-kBins;k<kBins;k++)
		{
			//for (int l=-lBins+1;l<lBins;l++)
			for (int l=-lBins;l<lBins;l++)
			{
				xray.continuousMiller[h][k+kBins][l+lBins].Pos[X]=Real(h)*xray.hBinSize;
				xray.continuousMiller[h][k+kBins][l+lBins].Pos[Y]=Real(k)*xray.kBinSize;
				xray.continuousMiller[h][k+kBins][l+lBins].Pos[Z]=Real(l)*xray.lBinSize;
			}
		}
	}
*/
}

void initializeContinuous(XRayStruct &xray, XRayParamStruct &params)
{
	int hBins, kBins, lBins;
	PosStruct tempPos;
	Matrix Mat;
	
	tempPos.Pos.resize(3);
	hBins=params.ContinuousHValues;
	kBins=params.ContinuousKValues*2;
	lBins=params.ContinuousLValues*2;
	cout <<"Before initializeContinuousMiller"<<endl;	
	initializeContinuousMiller(xray, params);
	cout <<"After initializeContinuousMiller"<<endl;	
	//continuousMillerToContinuousQ(xray, xlength, ylength, zlength);
	Safe2DAlloc(Mat, 3, 3, "Mat");
	Safe3DAlloc(xray.realContinuous, hBins, kBins, lBins);
	cout <<"After Safe3DAlloc"<<endl;
	Safe3DAlloc(xray.imagContinuous, hBins, kBins, lBins);
/*
	if (params.FastTranslationRotationInterpolation=="Hessian")
	{
		Safe3DAlloc(xray.realContinuousDer, tempPos, hBins, kBins, lBins, "Der");
		Safe3DAlloc(xray.imagContinuousDer, tempPos, hBins, kBins, lBins, "Der");
		Safe3DAlloc(xray.realContinuousHes, Mat, hBins, kBins, lBins, "Mat");
		Safe3DAlloc(xray.imagContinuousHes, Mat, hBins, kBins, lBins, "Mat");
	}
*/
	cout <<"After Safe3DAlloc"<<endl;

	cout <<"maxH= "<<xray.maxH<<" maxK= "<<xray.maxK<<" maxL= "<<xray.maxL<<endl;
	cout <<"hBins= "<<hBins<<" kBins= "<<kBins<<" lBins= "<<lBins<<endl;
	cout <<"hBinSize= "<<xray.hBinSize<<" kBinSize= "<<xray.kBinSize<<" lBinSize= "<<xray.lBinSize<<endl;
}

void initializeMiller(XRayStruct &xray, XRayParamStruct &params)
{
	PosStruct tempPos;
	tempPos.Pos.resize(3);
	cout <<"In initializeMiller"<<endl;
	for (int i=0;i<=params.NumQxValues;i++)
	{
		tempPos.Pos[X]=Real(i);
		for (int j=0;j<=params.NumQyValues;j++)
		{
			tempPos.Pos[Y]=Real(j);
			for (int k=0;k<=params.NumQzValues;k++)
			{
				tempPos.Pos[Z]=Real(k);
				SafePushBack(xray.miller, tempPos, "xray.miller");
			}
		}
	}
}

void initializeXRay(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	int npoint;
	//Real xlength, ylength, zlength;
	cout <<"In initializeXRay"<<endl;

	xray.maxH=params.NumQxValues;
	xray.maxK=params.NumQyValues;
	xray.maxL=params.NumQzValues;
	cout <<"maxH= "<<xray.maxH<<" maxK= "<<xray.maxK<<" maxL= "<<xray.maxL<<endl;
	cout <<"NumQxValues= "<<params.NumQxValues<<endl;
	cout <<"NumQyValues= "<<params.NumQyValues<<endl;
	cout <<"NumQzValues= "<<params.NumQzValues<<endl;
	if (params.UseMillerIndexesFromExperiment) 
	{
		cout <<"xray.miller=expXRay.miller"<<endl;
		xray.miller=expXRay.miller;
		cout <<"xray.miller.size()= "<<xray.miller.size()<<endl;
	}
	else initializeMiller(xray, params);
	if (params.FilterOutHighResolution) 
	{
		filterOutHighResolution(xray, params.NumQxValues, params.NumQyValues, params.NumQzValues);
		cout <<"Before removeNonoverlapping"<<endl;
		removeNonoverlapping(xray, expXRay);
		cout <<"After removeNonoverlapping"<<endl;
	}
	//xlength=params.XCubeLength*params.XCubes;
	//ylength=params.YCubeLength*params.YCubes;
	//zlength=params.ZCubeLength*params.ZCubes;
	cout <<"Before millerToQ"<<endl;
	cout <<"xray.miller.size()= "<<xray.miller.size()<<endl;
	//millerToQ(xray, params.a, params.b, params.c);
	millerToQ(xray, params);
	cout <<"After millerToQ"<<endl;
	cout <<"xray.miller.size()= "<<xray.miller.size()<<endl;
	cout <<"xray.q.size()= "<<xray.q.size()<<endl;
//	initializePatterson(xray.patterson);
	cout <<"After initializePatterson"<<endl;
	npoint=xray.q.size();
	cout <<"npoint= "<<npoint<<endl;
	SafeAlloc(xray.real, npoint, "xray.real");
	SafeAlloc(xray.imag, npoint, "xray.imag");
	SafeAlloc(xray.waterReal, npoint, "xray.waterReal");
	SafeAlloc(xray.waterImag, npoint, "xray.waterImag");
	cout <<"After SafeAlloc"<<endl;
	cout <<"xray.q.size()= "<<xray.q.size()<<endl;
	SafeAlloc(xray.patterson.real, npoint, "xray.patterson.real");
	SafeAlloc(xray.patterson.imag, npoint, "xray.patterson.imag");
	cout <<"After SafeAlloc"<<endl;
	cout <<"xray.q.size()= "<<xray.q.size()<<endl;
	solventCorrectedScatteringFactor(xray.f, xray.q, params);
	initializeCosLookUp(xray.lookUp);
	initializeSinLookUp(xray.lookUp);
	cout <<"Before initializeExpLookUp"<<endl;
	initializeExpLookUp(xray.lookUp);
	initializeNegativeExpLookUp(xray.lookUp);
	cout <<"After solventCorrectedScatteringFactor"<<endl;
	cout <<"xray.q.size()= "<<xray.q.size()<<endl;
	cout <<"UseFastTranslationRotation= "<<params.UseFastTranslationRotation<<endl;
	if (params.UseFastTranslationRotation)
	{
		initializeContinuous(xray, params);
	}
	cout <<"xray.q.size()= "<<xray.q.size()<<endl;
	if (params.CorrectionFile!="") 
	{
		readCorrectionFile(params.CorrectionFile, xray);
	}
	if (params.OccupancyCorrectionFile!="")
	{
		cout <<"OccupancyCorrectionFile= "<<params.OccupancyCorrectionFile<<endl;
		readOccupancyCorrectionFile(params.OccupancyCorrectionFile, xray);
	}
}

void initializeXRay(XRayStruct &xray, XRayStruct &expXRay, ProteinStruct &Protein, XRayParamStruct &params)
{
	initializeXRay(xray, expXRay, params);
	initializeCentric(xray, Protein);
	calcUnitCellVectors(Protein, params);
	findCoreAtoms(Protein.Atoms, params.CalphaRadius, params.AsaCutoff);
	initializeScatteringFactors(xray, Protein, params);
}

void initializeRand(int randomSeed)
{
	srand(randomSeed);
}

void initializeScores(XRayStruct &xray, XRayParamStruct &params)
{
	int count=0;

	cout <<"In inializeScores"<<endl;
	cout <<"NUM_SCORE_TYPES= "<<NUM_SCORE_TYPES<<endl;
	SafeAlloc(xray.scoreIndex, NUM_SCORE_TYPES, "scoreIndex");
	SafeAlloc(xray.score, NUM_SCORE_TYPES, "scores");
	SafeAlloc(xray.stdDev, Real(1.0), NUM_SCORE_TYPES, "stdDev");
	cout <<"xray.scoreIndex.size()= "<<xray.scoreIndex.size()<<endl;
	if (params.RFactorWeight!=0)
	{
		xray.scoreIndex[RFACTOR]=count;
		count++;
	}
	if (params.DiffOfSquaresWeight!=0)
	{
		xray.scoreIndex[SQR_DIFF]=count;
		count++;
	}
	if (params.PearsonWeight!=0)
	{
		xray.scoreIndex[PEARSON]=count;
		count++;
	}
	if (params.ProductWeight!=0)
	{
		xray.scoreIndex[PRODUCT]=count;
		count++;
	}
	if (params.ClashWeightCoreCore!=0)
	{
		xray.scoreIndex[CLASH_CC]=count;
		count++;
	}
	if (params.ClashWeightCoreSurface!=0)
	{
		xray.scoreIndex[CLASH_CS]=count;
		count++;
	}
	if (params.ClashWeightSurfaceSurface!=0)
	{
		xray.scoreIndex[CLASH_SS]=count;
		count++;
	}
}

void initialize(LatticeStruct &lattice, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	//initializeLattice(lattice, params);
	initializeXRay(xray, expXRay, params);
	initializeRand(params.RandomSeed);
	initializeScores(xray, params);
}

