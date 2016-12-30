# include <iostream>
# include <vector>
# include <sys/time.h>
# include <algorithm>

using namespace std;


# include "../../LibraryFiles/TypeDef.h"
# include "../../LibraryFiles/Structures.h"
# include "../../LibraryFiles/normalize.h"
# include "../../LibraryFiles/AtomUtils.h"
# include "../../LibraryFiles/minimize.h"
# include "../../LibraryFiles/MinMax.h"
# include "../../LibraryFiles/MonteCarlo.h"
# include "../../LibraryFiles/MonteCarloTemplate.h"
# include "../../LibraryFiles/VectorManip.h"
# include "../../LibraryFiles/time.h"
# include "../../LibraryFiles/MemoryUsage.h"
# include "../../LibraryFiles/PrintPdb.h"
# include "../../LibraryFiles/IOUtils.h"
# include "../../LibraryFiles/StringUtils.h"
# include "../../LibraryFiles/rmsd.h"
# include "../../LibraryFiles/ReadPdb.h"
# include "../../LibraryFiles/TMScore.h"
# include "SpaceGroup.h"
# include "Params.h"
# include "LatticeStruct.h"
# include "XRay.h"
# include "RMSD.h"
# include "ReadExperimentalDataFile.h"
# include "Initialize.h"
# include "LatticeUtils.h"
# include "molecularReplacement.h"
# include "MaximumLikelihood.h"
# include "CalcClashScore.h"

void iterativeElectronDensityModification(XRayStruct &expXRay, XRayStruct &xray, LatticeStruct &lattice, XRayParamStruct &params)
{
	Real r;
	for (int i=0;i<params.MRIterations;i++)
	{
		calcScattering(lattice, xray);
		calcDensity(xray, expXRay, lattice);
		r=RFactor(xray.i, expXRay.i);
	}
}

Real quantifyMatch(XRayStruct &xray, XRayStruct &expXRay, string xRayScoreType)
{
	int npoint=xray.i.size();
	int nexp=expXRay.i.size();
	Real score=0;
	if (npoint!=nexp || npoint==0)
	{
		string errorStr="The number of calculated points ";
		errorStr+=IntToStr(npoint)+" does not match the number of ";
		errorStr+="experimental points "+IntToStr(nexp);
		error(errorStr, __LINE__, __FILE__);
	}
	if (xRayScoreType=="RFactor") score=RFactorFast(xray.i, expXRay.F);
	else if (xRayScoreType=="DiffOfSquares") score=calcBestSqrDiff(xray.i, expXRay.i);
	else if (xRayScoreType=="Pearson") score=1.0-calcPearson(xray.i, expXRay.i);
	else if (xRayScoreType=="Product") score=calcDotProduct(xray.i, expXRay.i);
	else if (xRayScoreType=="MLRF0") score=calcMLRF0(xray, expXRay);
	else if (xRayScoreType=="MLTF") score=calcMLTF(xray, expXRay, 0);
	else
	{
		string errorStr="Invalid option for XRayScoreType,"+xRayScoreType+". ";
		errorStr+="The valid values are RFactor, DiffOfSquares, ";
		errorStr+="Pearson, and MLRF0.";
		error(errorStr, __LINE__, __FILE__);
	}
	return score;
}

Real calcLogSqrDiff(Vector &calci, Vector &expi)
{
	int Size1=calci.size();
	int Size2=expi.size();
	Real logCalc, logObs, scale;
	Vector logCalcs, logObss;

	if (Size1!=Size2 || Size1==0)
	{
		string errorStr="Size1= "+toStr(Size1)+" Size2= "+toStr(Size2);
		error(errorStr, __LINE__, __FILE__);
	}
	scale=calcScaleForSqrDiff(calci, expi);
	for (int i=0;i<Size1;i++)
	{
		if (calci[i]>0 && expi[i]>0)
		{
			logCalc=log(calci[i]*scale);
			logObs=log(expi[i]);
			SafePushBack(logCalcs, logCalc, "logCalcs");
			SafePushBack(logObss, logObs, "logObs");
		}
	}
	return calcBestSqrDiff(logCalcs, logObss);
}

Real quantifyWeightedMatch(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	int scoreIndex;
	int npoint=xray.i.size();
	int nexp=expXRay.i.size();
	Real score=0;
	Real r, sqrDiff, pearson, p;
 
	//timeval start, end;

	if (npoint!=nexp || npoint==0)
	{
		string errorStr="The number of calculated points ";
		errorStr+=IntToStr(npoint)+" does not match the number of ";
		errorStr+="experimental points "+IntToStr(nexp);
		error(errorStr, __LINE__, __FILE__);
	}
	if (params.RFactorWeight!=0) 
	{
		//gettimeofday(&start, NULL);
		scoreIndex=xray.scoreIndex[RFACTOR];
		r=RFactorFast(xray.i, expXRay.F);
		score+=r*params.RFactorWeight/xray.stdDev[scoreIndex];
		xray.score[scoreIndex]=r;
		//gettimeofday(&end, NULL);
		if (isnan(score))
		{
			error("score= nan", __LINE__, __FILE__);
		}
	}
	if (params.DiffOfSquaresWeight!=0) 
	{
		//Real sqrDiff=calcBestSqrDiff(xray.i, expXRay.i)*params.DiffOfSquaresWeight*3000.0;
		//gettimeofday(&start, NULL);
		scoreIndex=xray.scoreIndex[SQR_DIFF];
		sqrDiff=calcBestSqrDiff(xray.i, expXRay.i);
		score+=sqrDiff*params.DiffOfSquaresWeight/xray.stdDev[scoreIndex];
		xray.score[scoreIndex]=sqrDiff;
		//gettimeofday(&end, NULL);
		//score+=sqrDiff/Real(npoint);
		if (isnan(score))
		{
			error("score= nan", __LINE__, __FILE__);
		}
	}
	if (params.PearsonWeight!=0) 
	{
		//gettimeofday(&start, NULL);
		scoreIndex=xray.scoreIndex[PEARSON];
		pearson=1.0-calcPearson(xray.i, expXRay.i);
		score+=pearson*params.PearsonWeight/xray.stdDev[scoreIndex];
		xray.score[scoreIndex]=pearson;
		//gettimeofday(&end, NULL);
		if (isnan(score))
		{
			error("score= nan", __LINE__, __FILE__);
		}
	}
	if (params.ProductWeight!=0) 
	{
		//Real p=calcDotProduct(xray.i, expXRay.i)*params.ProductWeight*3000.0;
		//score+=p/Real(npoint);
		//gettimeofday(&start, NULL);
		scoreIndex=xray.scoreIndex[PRODUCT];
		p=calcDotProduct(xray.i, expXRay.i);
		score+=p*params.ProductWeight/xray.stdDev[scoreIndex];
		xray.score[scoreIndex]=p;
		//gettimeofday(&end, NULL);
		if (isnan(score))
		{
			error("score= nan", __LINE__, __FILE__);
		}
	}
	if (params.MLRF0Weight!=0) 
	{
		//gettimeofday(&start, NULL);
		score+=calcMLRF0(xray, expXRay)*params.MLRF0Weight;
		//gettimeofday(&end, NULL);
		if (isnan(score))
		{
			error("score= nan", __LINE__, __FILE__);
		}
	}
	if (params.RotationWeight!=0)
	{
		//gettimeofday(&start, NULL);
		score+=calcRotationPossibleScore(xray, expXRay)*params.RotationWeight;
		//gettimeofday(&end, NULL);
		if (isnan(score))
		{
			error("score= nan", __LINE__, __FILE__);
		}
	}
	if (params.LogSquareDiffWeight!=0)
	{
		//gettimeofday(&start, NULL);
		score+=calcLogSqrDiff(xray.i, expXRay.i)*params.LogSquareDiffWeight;
		//gettimeofday(&end, NULL);
		if (isnan(score))
		{
			error("score= nan", __LINE__, __FILE__);
		}
	}
	return score;
}

Real quantifyMatch(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	return quantifyMatch(xray, expXRay, params.XRayScoreType);
}

Real phase(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	Real initRFactor, finalRFactor;
	calcScatteringComplex(Proteins, xray);
	initRFactor=quantifyMatch(xray, expXRay, params);
	calcDensity(xray, expXRay, lattice);
	iterativeElectronDensityModification(expXRay, xray, lattice, params);
	calcScattering(lattice, xray);
	finalRFactor=quantifyMatch(xray, expXRay, params);
	return finalRFactor;
}

Real phase(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	vector<ProteinStruct> Proteins;

	splitIntoChains(Protein, Proteins);
	return phase(Proteins, xray, expXRay, lattice, params);
}

Real phase(Vector degOfFreedom, MRStruct args)
{
	placeProteins(args.Proteins, degOfFreedom);
	return phase(args.Proteins, args.xray, args.expXRay, args.lattice, args.params);
}

Real calcFitness(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	Real score=0;
	if (params.TranslationRotationMethod=="Patterson") 
	{
		score=-calcPattersonOverlap(Proteins, xray, expXRay);
	}
	else if (params.TranslationRotationMethod=="Phase")
	{
		score=phase(Proteins, xray, expXRay, lattice, params);
	}
	else if (params.TranslationRotationMethod=="RFactor")
	{
		score=calcScattering(Proteins, xray, expXRay, params);
	}
	else
	{
		string errorStr="Unknown TranslationRotationMethod ";
		errorStr+=params.TranslationRotationMethod;
		errorStr+=" Acceptable values are Patterson Phase RFactor.";
		error(errorStr, __LINE__, __FILE__);
	}
	return score;
}

Real calcFitness(Vector degOfFreedom, MRStruct args)
{
	Real score=0;
	placeProteins(args.Proteins, degOfFreedom);
	degOfFreedomToRotMat(args.Proteins, degOfFreedom);
	score=calcFitness(args.Proteins, args.xray, args.expXRay, args.lattice, args.params);
	return score;
}

Real calcFitness(Vector &degOfFreedom, vector<ProteinStruct> Proteins, XRayStruct &XRay, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	Real score=0;

	placeProteins(Proteins, degOfFreedom);
	degOfFreedomToRotMat(Proteins, degOfFreedom);
	setSegmentOccupancies(Proteins[0], degOfFreedom, 1);
	score=calcFitness(Proteins, XRay, expXRay, lattice, params);
	return score;
}


void calcScattering(Vector &degOfFreedom, ProteinStruct Protein, XRayStruct &xray, LatticeStruct &lattice, XRayParamStruct &params)
{
	int nprot, nsym, npoint;
	Vector dx, translation;
	Matrix rotationMatrix;
	Get3DVectorSize(xray.complexSymReal, nprot, nsym, npoint);

	if (nprot!=1)
	{
		Safe3DAlloc(xray.complexSymReal, 1, nsym, npoint, "real");
		Safe3DAlloc(xray.complexSymImag, 1, nsym, npoint, "imag");
		Safe3DAlloc(xray.complexSymTempReal, 1, nsym, npoint, "tempReal");
		Safe3DAlloc(xray.complexSymTempImag, 1, nsym, npoint, "tempImag");
	}

	if (!params.UseFastTranslationRotation)
	{
		rotateAtoms(Protein.Atoms, 0, 0, degOfFreedom[Z]);
		rotateAtoms(Protein.Atoms, 0, degOfFreedom[Y], 0);
		rotateAtoms(Protein.Atoms, degOfFreedom[X], 0, 0);
	}

	SafeAlloc(dx, 3, "dx");
	if (params.UseFastTranslationRotation)
	{
		calcRotationMatrix(rotationMatrix, degOfFreedom[X], degOfFreedom[Y], degOfFreedom[Z]);
		SafeAlloc(translation, 3, "translation");
		calcFastTranslationRotation(Protein, xray, rotationMatrix, 0, params);
		xray.complexSymTempReal=xray.complexSymReal;
		xray.complexSymTempImag=xray.complexSymImag;
	}
	else
	{
		calcSymmetryScattering(xray, Protein, 0);
	}
	dx[X]=degOfFreedom[3+X];
	dx[Y]=degOfFreedom[3+Y];
	dx[Z]=degOfFreedom[3+Z];
	translateXRay2(xray, Protein, dx, 0);
	calcTotalFormFactor(xray);

}

void initializeXRayProtein(ProteinStruct &Protein, XRayStruct &xray)
{
	int nsym=Protein.symmetryOperations.size();

	ProteinStruct tempProtein;
	tempProtein=Protein;
	SafeAlloc(tempProtein.center, 3, "center");
	removeNonCA(tempProtein.Atoms);
	SafeAlloc(xray.Proteins, tempProtein, nsym*27, "Proteins");
}

Real calcFitnessChange(Vector &degOfFreedom, ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params, int pick)
{
	int prot=pick/6;
	int scoreIndex;
	Real score=0;
	Real clashScore;
	Real x, y, z;
	Real clashCC, clashCS, clashSS;
	//timeval start, end;

	int nprot, nsym, npoint;
	Get3DVectorSize(xray.complexSymReal, nprot, nsym, npoint);

	if (!params.UseFastTranslationRotation)
	{
		rotateAtoms(Protein.Atoms, 0, 0, degOfFreedom[prot*6+Z]);
		rotateAtoms(Protein.Atoms, 0, degOfFreedom[prot*6+Y], 0);
		rotateAtoms(Protein.Atoms, degOfFreedom[prot*6+X], 0, 0);
	}
	//moveAtoms(Protein.Atoms, degOfFreedom[prot*6+3], degOfFreedom[prot*6+4], degOfFreedom[prot*6+5]);
	//degOfFreedomToRotMat(Protein, degOfFreedom);
	//gettimeofday(&start, NULL);
	calcXRayScatteringChange(Protein, xray, pick, degOfFreedom, params);
	//gettimeofday(&end, NULL);
	if (params.CorrectionFile!="") correctIntensity(xray);
	score=calcRFactor(xray, expXRay, params);
	x=Protein.XBoxLength;
	y=Protein.YBoxLength;
	z=Protein.ZBoxLength;
	if (params.ClashScoreUse=="DuringMR")
	{
		int nXrayProt=xray.Proteins.size();
		if (nXrayProt==0)
		{
			initializeXRayProtein(Protein, xray);
		}
		//Real clashScoreSlow;
		clashScore=calcClashScoreFast(degOfFreedom, xray.Proteins, params, clashCC, clashCS, clashSS);
/*
		if (pick<6)
		{
			clashScore=calcClashScoreFast(degOfFreedom, xray.Proteins, params);
			xray.clashScore=clashScore;
		}
		else
		{
			clashScore=xray.clashScore;
		}
*/
		scoreIndex=xray.scoreIndex[CLASH_CC];
		score+=clashCC*params.ClashWeightCoreCore/xray.stdDev[scoreIndex];
		xray.score[scoreIndex]=clashCC;

		scoreIndex=xray.scoreIndex[CLASH_CS];
		score+=clashCS*params.ClashWeightCoreSurface/xray.stdDev[scoreIndex];
		xray.score[scoreIndex]=clashCS;

		scoreIndex=xray.scoreIndex[CLASH_SS];
		score+=clashSS*params.ClashWeightSurfaceSurface/xray.stdDev[scoreIndex];
		xray.score[scoreIndex]=clashSS;

		if (isnan(score))
		{
			error("ClashScore is nan ", __LINE__, __FILE__);
		}
	}
	//score+=penalizeOutOfBounds(degOfFreedom, x, y, z);
	return score;

}

Real calcFitnessTranslation(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &dPos)
{
	Real score;
	vector<ProteinStruct> tempProteins;
	if (params.UseFastTranslation)
	{
		translateXRay2(xray, Proteins[0], dPos, 0);
		calcTotalFormFactor(xray);
		calcIntensity(xray);
		score=calcRFactor(xray, expXRay, params);
	}
	else
	{
		SafePushBack(tempProteins, Proteins[0], "tempProteins");
		moveAtoms(Proteins[0].Atoms, tempProteins[0].Atoms, dPos[X], dPos[Y], dPos[Z]);
		score=calcScattering(tempProteins, xray, expXRay, params);
	}
	return score;
}

Real calcFitnessTranslation(Vector dPos, MRStruct args)
{
	Real score;
	score=calcFitnessTranslation(args.Proteins, args.xray, args.expXRay, args.params, dPos);
	return score;
}

void degOfFreedomToRotMat(Matrix &rotationMatrix, Real rx, Real ry, Real rz)
{
	Matrix rotationMatrixXY;
	Matrix rotationMatrixX, rotationMatrixY, rotationMatrixZ;

	calcRotationMatrixAroundX(rx, rotationMatrixX);
	calcRotationMatrixAroundY(ry, rotationMatrixY);	
	calcRotationMatrixAroundZ(rz, rotationMatrixZ);	
	matrixMultiply(rotationMatrixY, rotationMatrixX, rotationMatrixXY);
	//matrixMultiply(rotationMatrixX, rotationMatrixY, rotationMatrixXY);
	matrixMultiply(rotationMatrixZ, rotationMatrixXY, rotationMatrix);
	//Protein.rotationMatrix=rotationMatrix;
	//SafeAlloc(Protein.translation, 3, "Protein.translation");
	//Protein.translation[X]=0;
	//Protein.translation[Y]=0;
	//Protein.translation[Z]=0;
}

void degOfFreedomToRotMat(vector<ProteinStruct> &Proteins, Vector &degOfFreedom)
{
	int nprot=Proteins.size();
	Vector subV;

	for (int i=0;i<nprot;i++)
	{
		subV=getSubVector(degOfFreedom, i*6, i*6+5);
		degOfFreedomToRotMat(Proteins[i].rotationMatrix, degOfFreedom[i*6+X], degOfFreedom[i*6+Y], degOfFreedom[i*6+Z]);
	}
}

void degOfFreedomToRotMatAndTranslation(ProteinStruct &Protein, Vector &degOfFreedom)
{
	Matrix rotationMatrix, rotationMatrixXY;
	Matrix rotationMatrixX, rotationMatrixY, rotationMatrixZ;

	calcRotationMatrixAroundX(degOfFreedom[0], rotationMatrixX);	
	calcRotationMatrixAroundY(degOfFreedom[1], rotationMatrixY);	
	calcRotationMatrixAroundZ(degOfFreedom[2], rotationMatrixZ);	
	matrixMultiply(rotationMatrixX, rotationMatrixY, rotationMatrixXY);
	matrixMultiply(rotationMatrixZ, rotationMatrixXY, rotationMatrix);
	rotationMatrix[X][X]=1.0;rotationMatrix[X][Y]=0;rotationMatrix[X][Z]=0;
	rotationMatrix[Y][X]=0;rotationMatrix[Y][Y]=1.0;rotationMatrix[Y][Z]=0;
	rotationMatrix[Z][X]=0;rotationMatrix[Z][Y]=0;rotationMatrix[Z][Z]=1.0;
	Protein.rotationMatrix=rotationMatrix;
	SafeAlloc(Protein.translation, 3, "Protein.translation");
	//Protein.translation[X]=degOfFreedom[3];
	//Protein.translation[Y]=degOfFreedom[4];
	//Protein.translation[Z]=degOfFreedom[5];
}

void rotateProteins(vector<ProteinStruct> Proteins1, Vector &degOfFreedom, vector<ProteinStruct> &Proteins2)
{
	int nprot=Proteins1.size();
	for (int i=0;i<nprot;i++)
	{
		rotateAtoms(Proteins1[i].Atoms, 0, 0, degOfFreedom[6*i+2]);
		rotateAtoms(Proteins1[i].Atoms, 0, degOfFreedom[6*i+1], 0);
		rotateAtoms(Proteins1[i].Atoms, degOfFreedom[6*i], 0, 0);
		moveIntoUnitCell(Proteins1[i]);
	}
	Proteins2=Proteins1;
}

void placeProtein(ProteinStruct &Protein1, Vector &degOfFreedom, ProteinStruct &Protein2)
{
	Protein2=Protein1;
	rotateAtoms(Protein2.Atoms, 0, 0, degOfFreedom[Z]);
	rotateAtoms(Protein2.Atoms, 0, degOfFreedom[Y], 0);
	rotateAtoms(Protein2.Atoms, degOfFreedom[X], 0, 0);
	moveAtoms(Protein2.Atoms, degOfFreedom[3+X], degOfFreedom[3+Y], degOfFreedom[3+Z]);
	moveIntoUnitCell(Protein2);
}

void placeProtein(ProteinStruct &Protein, Vector &degOfFreedom)
{
	rotateAtoms(Protein.Atoms, 0, 0, degOfFreedom[Z]);
	rotateAtoms(Protein.Atoms, 0, degOfFreedom[Y], 0);
	rotateAtoms(Protein.Atoms, degOfFreedom[X], 0, 0);
	moveAtoms(Protein.Atoms, degOfFreedom[3+X], degOfFreedom[3+Y], degOfFreedom[3+Z]);
	moveIntoUnitCell(Protein);
}

void placeProteins(vector<ProteinStruct> &Proteins, Vector &degOfFreedom)
{
	int nprot=Proteins.size();
	int natom;
	for (int i=0;i<nprot;i++)
	{
		rotateAtoms(Proteins[i].Atoms, 0, 0, degOfFreedom[6*i+Z]);
		rotateAtoms(Proteins[i].Atoms, 0, degOfFreedom[6*i+Y], 0);
		rotateAtoms(Proteins[i].Atoms, degOfFreedom[6*i+X], 0, 0);
		moveAtoms(Proteins[i].Atoms, degOfFreedom[6*i+3], degOfFreedom[6*i+4], degOfFreedom[6*i+5]);
		natom=Proteins[i].Atoms.size();
		for (int j=0;j<natom;j++)
		{
			Proteins[i].Atoms[j].ChainName=ALPHABET[i];
		}
		moveIntoUnitCell(Proteins[i]);
	}
}

void placeProteins(Vector &degOfFreedom, ProteinStruct &Protein, vector<ProteinStruct> &Proteins)
{
	int nprot=degOfFreedom.size()/6;

	SafeAlloc(Proteins, Protein, nprot, "Proteins");
	placeProteins(Proteins, degOfFreedom);
}

void placeProteins(vector<ProteinStruct> &Proteins1, Vector &degOfFreedom, vector<ProteinStruct> &Proteins2)
{
	Proteins2=Proteins1;
	placeProteins(Proteins2, degOfFreedom);
}

void placeProteins(vector<ProteinStruct> Proteins, Vector &degOfFreedom, ProteinStruct &Complex, XRayStruct xray)
{
	Complex=Proteins[0];
	Complex.Atoms.clear();
	placeProteins(Proteins, degOfFreedom);
	mergeProteins(Proteins, Complex.Atoms);
}

void placeProteins(vector<ProteinStruct> Proteins, Vector &degOfFreedom, ProteinStruct &Complex)
{
	Complex=Proteins[0];
	Complex.Atoms.clear();
	placeProteins(Proteins, degOfFreedom);
	mergeProteins(Proteins, Complex.Atoms);
}

void makeUnitCell(ProteinStruct &Protein, vector<ProteinStruct> &allProteins)
{
	int noperations=Protein.symmetryOperations.size();
	ProteinStruct tempProtein;
	tempProtein=Protein;
	if (noperations==0) {setSpaceGroupP1(Protein); noperations=1;}
	SafeAlloc(allProteins, tempProtein, noperations, "allProteins");
	for (int j=0;j<noperations;j++)
	{
		allProteins[j].Atoms=Protein.Atoms;
		applySymmetryOperation(allProteins[j].Atoms, Protein.symmetryOperations[j]);
		//tempProtein.Atoms=Protein.Atoms;
		//applySymmetryOperation(tempProtein.Atoms, Protein.symmetryOperations[j]);
		//SafePushBack(allProteins, tempProtein, "allProteins");
	}
}

void makeImages(vector<ProteinStruct> &allProteins)
{
	int nprot=allProteins.size();
	int count=nprot;
	int totalProt;
	Real dx, dy, dz;
	VectorStruct a, b, c;
	ProteinStruct tempProtein;

	a=allProteins[0].a;
	b=allProteins[0].b;
	c=allProteins[0].c;

	totalProt=27*nprot;
	allProteins.resize(totalProt);

	for (Real h=-1.0;h<1.1;h+=1.0)
	{
		for (Real k=-1.0;k<1.1;k+=1.0)
		{
			for (Real l=-1.0;l<1.1;l+=1.0)
			{
				if (!(h==0 && k==0 && l==0))
				{
					for (int i=0;i<nprot;i++)
					{
						dx=h*a.x+k*b.x+l*c.x;
						dy=h*a.y+k*b.y+l*c.y;
						dz=h*a.z+k*b.z+l*c.z;
						allProteins[count]=allProteins[i];
						moveAtoms(allProteins[count].Atoms, dx, dy, dz);
						count++;
						//SafePushBack(allProteins, tempProtein, "allProteins");
					}
				}
			}
		}
	}
}

void makeUnitCellAndImages(ProteinStruct &Protein, vector<ProteinStruct> &allProteins)
{
	makeUnitCell(Protein, allProteins);
	makeImages(allProteins);
}

void makeUnitCellAndImages(vector<ProteinStruct> &Proteins, vector<ProteinStruct> &allProteins)
{
	int nprot=Proteins.size();
	for (int i=0;i<nprot;i++) makeUnitCell(Proteins[i], allProteins);
	makeImages(allProteins);
}

void makeUnitCellAndImagesFast(vector<ProteinStruct> &Proteins, vector<ProteinStruct> &allProteins)
{
	//Slower not faster.  Why?
	int nprot=Proteins.size();
	int noperations=Proteins[0].symmetryOperations.size();
	int index=0;
	Real dx, dy, dz;
	VectorStruct a, b, c;

	a=Proteins[0].a;
	b=Proteins[0].b;
	c=Proteins[0].c;
	if (noperations==0) {setSpaceGroupP1(Proteins[0]); noperations=1;}
	//SafeAlloc(allProteins, Proteins[0], noperations*nprot*27, "allProteins");
	allProteins.resize(noperations*nprot*27);
	for (int i=0;i<nprot;i++)
	{
		SafeAlloc(Proteins[i].center, 3, "center");
		calcCenter(Proteins[i].Atoms, Proteins[i].center[X], Proteins[i].center[Y], Proteins[i].center[Z]);
		for (int j=0;j<noperations;j++)
		{
			allProteins[index].Atoms=Proteins[i].Atoms;
			applySymmetryOperation(allProteins[index].Atoms, Proteins[i].symmetryOperations[j]);
			SafeAlloc(allProteins[index].center, 3, "center");
			calcCenter(allProteins[index].Atoms, allProteins[index].center[X], allProteins[index].center[Y], allProteins[index].center[Z]);
			//printVector(allProteins[index].center);
			index++;
		}
	}
	for (Real h=-1.0;h<1.1;h+=1.0)
	{
		for (Real k=-1.0;k<1.1;k+=1.0)
		{
			for (Real l=-1.0;l<1.1;l+=1.0)
			{
				if (!(h==0 && k==0 && l==0))
				{
					dx=h*a.x+k*b.x+l*c.x;
					dy=h*a.y+k*b.y+l*c.y;
					dz=h*a.z+k*b.z+l*c.z;
					for (int i=0;i<nprot*noperations;i++)
					{
						allProteins[index]=allProteins[i];
						moveAtoms(allProteins[index].Atoms, dx, dy, dz);
						//printVector(allProteins[index].center);
						allProteins[index].center[X]+=dx;
						allProteins[index].center[Y]+=dy;
						allProteins[index].center[Z]+=dz;
						//printVector(allProteins[index].center);
						index++;
					}
				}
			}
		}
	}
/*
	for (int i=0;i<nprot*noperations*27;i++)
	{
		printVector(allProteins[i].center);
	}
*/
}

void makeUnitCellAndImagesFast2(vector<ProteinStruct> &Proteins, vector<ProteinStruct> &allProteins)
{
	int nprot=Proteins.size();
	int noperations=Proteins[0].symmetryOperations.size();
	int index=0;
	Real r, r2, cutoff;
	Real x, y, z;
	Real dx, dy, dz;
	Real centerx, centery, centerz;
	VectorStruct a, b, c;
	timeval start, end;

	a=Proteins[0].a;
	b=Proteins[0].b;
	c=Proteins[0].c;
	if (noperations==0) {setSpaceGroupP1(Proteins[0]); noperations=1;}
	gettimeofday(&start, NULL);
	allProteins.reserve(noperations*nprot*27);
	SafeAlloc(allProteins, Proteins[0], noperations*nprot, "allProteins");
	gettimeofday(&end, NULL);
	gettimeofday(&start, NULL);
	for (int i=0;i<nprot;i++)
	{
		SafeAlloc(Proteins[i].center, 3, "center");
		calcCenter(Proteins[i].Atoms, Proteins[i].center[X], Proteins[i].center[Y], Proteins[i].center[Z]);
		for (int j=0;j<noperations;j++)
		{
			allProteins[index].Atoms=Proteins[i].Atoms;
			applySymmetryOperation(allProteins[index].Atoms, Proteins[i].symmetryOperations[j]);
			SafeAlloc(allProteins[index].center, 3, "center");
			calcCenter(allProteins[index].Atoms, allProteins[index].center[X], allProteins[index].center[Y], allProteins[index].center[Z]);
			index++;
		}
	}
	gettimeofday(&end, NULL);
	gettimeofday(&start, NULL);
	r=calcFurthestDistance(Proteins[0].Atoms);
	cutoff=4.0*r*r;
	gettimeofday(&start, NULL);
	for (Real h=-1.0;h<1.1;h+=1.0)
	{
		for (Real k=-1.0;k<1.1;k+=1.0)
		{
			for (Real l=-1.0;l<1.1;l+=1.0)
			{
				if (!(h==0 && k==0 && l==0))
				{
					dx=h*a.x+k*b.x+l*c.x;
					dy=h*a.y+k*b.y+l*c.y;
					dz=h*a.z+k*b.z+l*c.z;
					for (int i=0;i<nprot*noperations;i++)
					{
						centerx=allProteins[i].center[X]+dx;
						centery=allProteins[i].center[Y]+dy;
						centerz=allProteins[i].center[Z]+dz;
						for (int j=0;j<nprot;j++)
						{
							x=centerx-allProteins[j].center[X];
							y=centery-allProteins[j].center[Y];
							z=centerz-allProteins[j].center[Z];
							r2=x*x+y*y+z*z;
							if (r2<cutoff)
							{
								SafePushBack(allProteins, allProteins[j], "allProteins");
								moveAtoms(allProteins[index].Atoms, dx, dy, dz);
								allProteins[index].center[X]=centerx;
								allProteins[index].center[Y]=centery;
								allProteins[index].center[Z]=centerz;
								index++;
							}
						}
					}
				}
			}
		}
	}
	gettimeofday(&end, NULL);
}

void makeUnitCellAndImagesFast3(vector<ProteinStruct> &Proteins, vector<ProteinStruct> &allProteins)
{
	int nprot=Proteins.size();
	int noperations=Proteins[0].symmetryOperations.size();
	int index=0;
	Real r, r2, cutoff, buffer=3.0;
	Real x, y, z;
	Real dx, dy, dz;
	Real centerx, centery, centerz;
	VectorStruct a, b, c;

	a=Proteins[0].a;
	b=Proteins[0].b;
	c=Proteins[0].c;
	if (noperations==0) {setSpaceGroupP1(Proteins[0]); noperations=1;}
	//SafeAlloc(allProteins, Proteins[0], noperations*nprot*27, "allProteins");
	allProteins.resize(noperations*nprot*27);
	for (int i=0;i<noperations*nprot*27;i++)
	{
		SafeAlloc(allProteins[i].center, 3, "center");
	}
	for (int i=0;i<nprot;i++)
	{
		SafeAlloc(Proteins[i].center, 3, "center");
		calcCenter(Proteins[i].Atoms, Proteins[i].center[X], Proteins[i].center[Y], Proteins[i].center[Z]);
		for (int j=0;j<noperations;j++)
		{
			allProteins[index].Atoms=Proteins[i].Atoms;
			applySymmetryOperation(allProteins[index].Atoms, Proteins[i].symmetryOperations[j]);
			SafeAlloc(allProteins[index].center, 3, "center");
			calcCenter(allProteins[index].Atoms, allProteins[index].center[X], allProteins[index].center[Y], allProteins[index].center[Z]);
			index++;
		}
	}
	r=calcFurthestDistance(Proteins[0].Atoms);
	r+=buffer;
	cutoff=4.0*r*r;
	for (Real h=-1.0;h<1.1;h+=1.0)
	{
		for (Real k=-1.0;k<1.1;k+=1.0)
		{
			for (Real l=-1.0;l<1.1;l+=1.0)
			{
				if (!(h==0 && k==0 && l==0))
				{
					dx=h*a.x+k*b.x+l*c.x;
					dy=h*a.y+k*b.y+l*c.y;
					dz=h*a.z+k*b.z+l*c.z;
					for (int i=0;i<nprot*noperations;i++)
					{
						for (int j=0;j<nprot;j++)
						{
							centerx=allProteins[i].center[X]+dx;
							centery=allProteins[i].center[Y]+dy;
							centerz=allProteins[i].center[Z]+dz;
							x=centerx-allProteins[j].center[X];
							y=centery-allProteins[j].center[Y];
							z=centerz-allProteins[j].center[Z];
							//r2=sqrt(x*x+y*y+z*z)+2.0*buffer;
							//r2=r2*r2;
							r2=x*x+y*y+z*z;
							if (r2<cutoff)
							{
								//SafeAlloc(allProteins[index].center, 3, "center");
								allProteins[index]=allProteins[i];
								allProteins[index].center[X]=centerx;
								allProteins[index].center[Y]=centery;
								allProteins[index].center[Z]=centerz;
								moveAtoms(allProteins[index].Atoms, dx, dy, dz);
								index++;
							}
						}
					}
				}
			}
		}
	}
}

Real calcRotationTranslationFitness(Vector &degOfFreedom, vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	Real score=0;
	ProteinStruct Complex;

	placeProteins(Proteins, degOfFreedom);
	//degOfFreedomToRotMatAndTranslation(args.Protein, degOfFreedom);
	score=calcFitness(Proteins, xray, expXRay, lattice, params);

	return score;
}

Real calcRotationTranslationFitness(Vector degOfFreedom, MRStruct args)
{
	return calcRotationTranslationFitness(degOfFreedom, args.Proteins, args.xray, args.expXRay, args.lattice, args.params);
}

Real penalizeOutOfBounds(Vector &degOfFreedom, Real x, Real y, Real z)
{
	int nprot;
	Real penalty=0, penaltyWeight=1.0e15;

	nprot=degOfFreedom.size()/6;

	//penalty+=(degOfFreedom[0]-0.604154)*(degOfFreedom[0]-0.604154)*100.0;
	//penalty+=(degOfFreedom[1]-1.007720)*(degOfFreedom[1]-1.007720)*100.0;
	//penalty+=(degOfFreedom[2]-1.390380)*(degOfFreedom[2]-1.390380)*100.0;

	for (int j=0;j<nprot;j++)
	{
		for (int i=0;i<3;i++)
		{
			if (abs(degOfFreedom[j*6+i])>pi) 
			{
				penalty+=(abs(degOfFreedom[j*6+i])-pi)*(abs(degOfFreedom[j*6+i])-pi);
			}
		}
		if (abs(degOfFreedom[j*6+3])>x)
		{
			penalty+=(abs(degOfFreedom[j*6+3])-x)*(abs(degOfFreedom[j*6+3])-x);
		}
		if (abs(degOfFreedom[j*6+4])>y)
		{
			penalty+=(abs(degOfFreedom[j*6+4])-y)*(abs(degOfFreedom[j*6+4])-y);
		}
		if (abs(degOfFreedom[j*6+5])>z)
		{
			penalty+=(abs(degOfFreedom[j*6+5])-z)*(abs(degOfFreedom[j*6+5])-z);
		}
	}
	penalty*=penaltyWeight;

	return penalty;
}

Real calcRotationTranslationFitness2(Vector degOfFreedom, vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	Real score=0;
	Real x, y, z;
	ProteinStruct Complex;	

	placeProteins(Proteins, degOfFreedom);
	degOfFreedomToRotMat(Proteins, degOfFreedom); //Need to fix
	score=calcFitness(Proteins, xray, expXRay, lattice, params);
	x=Proteins[0].XBoxLength;
	y=Proteins[0].YBoxLength;
	z=Proteins[0].ZBoxLength;
	if (params.ClashScoreUse=="DuringMR")
	{
		score+=calcClashScore(Proteins, params);
	}
	//penalty=penalizeOutOfBounds(degOfFreedom, x, y, z);
	//score+=penalty;
	return score;
}

Real calcRotationTranslationFitness2(Vector degOfFreedom, MRStruct &args)
{
	Real score;

	score=calcRotationTranslationFitness2(degOfFreedom, args.Proteins, args.xray, args.expXRay, args.lattice, args.params);
	return score;
}


Real calcBFactorFitness(Vector bFactorScale, MRStruct args)
{
	int nprot=args.Proteins.size();
	Real score=0;

	for (int i=0;i<nprot;i++)
	{
		scaleBFactors(args.Proteins[i], bFactorScale[0]);	
	}
	score=calcScattering(args.Proteins, args.xray, args.expXRay, args.params);
	bulkSolventCorrection(args.xray, args.expXRay, args.params);
	return score;
}

Real calcFormFactorFitness(Vector formFactorScale, MRStruct args)
{
	Real score=0;

	args.params.scale=formFactorScale[0];
	solventCorrectedScatteringFactor(args.xray.f, args.xray.q, args.params);
	score=calcScattering(args.Proteins, args.xray, args.expXRay, args.params);
	return score;
}

Real MRRotateTranslate(ProteinStruct Protein, Vector degOfFreedom, XRayStruct xray, XRayStruct expXRay, LatticeStruct lattice, XRayParamStruct params)
{
	Real dScoreCriteria=0.0000001;
	Real maxHours=2.0;
	MRStruct args;
	Real (*func_to_minimize)(Vector, MRStruct)=calcFitness;

	args.Protein=Protein;
	args.xray=xray;
	args.expXRay=expXRay;
	args.lattice=lattice;
	args.params=params;

	return nonlinearConjugateGradient(degOfFreedom, args, func_to_minimize, dScoreCriteria, maxHours);
	//return OptimizeCoefficients(degOfFreedom, args, func_to_minimize, dScoreCriteria);
}

Real MRRotationalGridSearch(Vector &bestGrid, Real rinc, int npoint, ProteinStruct Protein, XRayStruct xray, XRayStruct &expXRay, LatticeStruct lattice, XRayParamStruct &params, vector<Array3D> &degOfFreedoms, Array3D &scores)
{
	Real r, rmin;
	Real xstart, ystart, zstart;
	Real xend, yend, zend;
	Vector degOfFreedom;
	vector<ProteinStruct> Proteins;
	timeval start, end;

	SafePushBack(Proteins, Protein, "Proteins");
	xstart=bestGrid[0]-rinc*(npoint-1)*0.5;
	ystart=bestGrid[1]-rinc*(npoint-1)*0.5;
	zstart=bestGrid[2]-rinc*(npoint-1)*0.5;
	xend=bestGrid[0]+rinc*(npoint-1)*0.5;
	yend=bestGrid[1]+rinc*(npoint-1)*0.5;
	zend=bestGrid[2]+rinc*(npoint-1)*0.5;
	SafeAlloc(degOfFreedom, 3, "defOfFreedom");
	degOfFreedom[0]=xstart;
	degOfFreedom[1]=ystart;
	degOfFreedom[2]=zstart;
	rmin=calcFitness(degOfFreedom, Proteins, xray, expXRay, lattice, params);
	Safe3DAlloc(scores, npoint, npoint, npoint, "scores");
	Safe3DAlloc(degOfFreedoms, degOfFreedom, npoint, npoint, npoint, "degOfFreedoms");
	bestGrid=degOfFreedom;
	for (int i=0;i<npoint;i++)
	{
		degOfFreedom[X]=xstart+rinc*Real(i);
		for (int j=0;j<npoint;j++)
		{
			degOfFreedom[Y]=ystart+rinc*Real(j);
			for (int k=0;k<npoint;k++)
			{
				degOfFreedom[Z]=zstart+rinc*Real(k);
				gettimeofday(&start, NULL);
				r=calcFitness(degOfFreedom, Proteins, xray, expXRay, lattice, params);
				gettimeofday(&end, NULL);
				scores[i][j][k]=r;
				degOfFreedoms[i][j][k]=degOfFreedom;
				if (r<rmin)
				{
					rmin=r;
					bestGrid=degOfFreedom;
				}
			}
		}
	}
	printVector(bestGrid, "bestGrid");
	return rmin;
}

Real MRRotationalGridSearch(Vector &bestGrid, ProteinStruct Protein, XRayStruct xray, XRayStruct &expXRay, LatticeStruct lattice, XRayParamStruct &params, vector<Array3D> &degOfFreedoms, Array3D &scores)
{
	int npoint=params.RotationGridSearchPoints;
	int iterations=1;
	Real rinc=2.0*pi/Real(npoint);
	Real rmin;
	Vector  oldGrid=bestGrid;
	for (int i=0;i<iterations;i++)
	{
		rmin=MRRotationalGridSearch(bestGrid, rinc, npoint, Protein, xray, expXRay, lattice, params, degOfFreedoms, scores);
		if (vectorEqual(oldGrid, bestGrid, rinc*0.5)) 
		{
			rinc*=0.25;
		}
		oldGrid=bestGrid;
	}
	return rmin;
}

Real minimizeRotation(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom, Real maxHours)
{
	Real dScoreCriteria=0.0000001;
	MRStruct args;
	Real (*func_to_minimize)(Vector, MRStruct)=calcFitness;
	args.Protein=Protein;
	args.xray=xray;
	args.expXRay=expXRay;
	args.params=params;

	//nonlinearConjugateGradient(degOfFreedom, args, func_to_minimize, dScoreCriteria);
	//conjugateGradient(degOfFreedom, args, func_to_minimize, dScoreCriteria, maxHours);
	optimizeUsingEstimatedGrad(degOfFreedom, args, func_to_minimize, dScoreCriteria, maxHours);
	//OptimizeCoefficients(degOfFreedom, args, func_to_minimize, dScoreCriteria);
	return calcFitness(degOfFreedom, args);
}

Real minimizeRotation(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Real maxHours)
{
	Real score;
	Vector degOfFreedom;

	SafeAlloc(degOfFreedom, 3, "degOfFreedom");
	calcUnitCellVectors(Protein, params);
	score=minimizeRotation(Protein, xray, expXRay, params, degOfFreedom, maxHours);
	rotateAtoms(Protein.Atoms, 0, 0, degOfFreedom[Z]);	
	rotateAtoms(Protein.Atoms, 0, degOfFreedom[Y], 0);	
	rotateAtoms(Protein.Atoms, degOfFreedom[X], 0, 0);	

	return score;
}

Real minimizeRotationOnly(ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom, Real maxHours)
{
	//minimizeRotationOnly takes in a vector, degOfFreedom, containing
	//both rotational and translational degrees of freedom, performs the
	//translation, and then optimizes the rotational degrees of freedom.
	//minimizeRotation takes in a vector containing only the 
	//rotational degrees of freedom and optimizes them.
	Real score;
	Vector dRot;
	timeval start, end;

	gettimeofday(&start, NULL);
	SafeAlloc(dRot, 3, "dPos");
	dRot[X]=degOfFreedom[X];
	dRot[Y]=degOfFreedom[Y];
	dRot[Z]=degOfFreedom[Z];
	moveAtoms(Protein.Atoms, degOfFreedom[3+X], degOfFreedom[3+Y], degOfFreedom[3+Z]);
	score=minimizeRotation(Protein, xray, expXRay, params, dRot, maxHours);
	degOfFreedom[X]=dRot[X];
	degOfFreedom[Y]=dRot[Y];
	degOfFreedom[Z]=dRot[Z];
	gettimeofday(&end, NULL);

	return score;
}

bool minimizeRotationUsingEstimatedScore(ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom)
{
	int maxSteps=20;
	Real length=1e-4;
	Real score, scoreOld;
	Matrix cartGrad;
	Vector rotGrad, backUpDOF;
	ProteinStruct backupProtein;
	vector<ProteinStruct> Proteins;
	EstimateScoreStruct args;
	Real (*func_to_minimize)(Vector, EstimateScoreStruct)=estimateXRayScore;

	backUpDOF=degOfFreedom;

	SafeAlloc(rotGrad, 3, "rotGrad");
	backupProtein=Protein;
	rotateAtoms(Protein.Atoms, 0, 0, degOfFreedom[Z]);	
	rotateAtoms(Protein.Atoms, 0, degOfFreedom[Y], 0);	
	rotateAtoms(Protein.Atoms, degOfFreedom[X], 0, 0);	
	SafePushBack(Proteins, Protein, "Proteins");
	scoreOld=calcScattering(Proteins, xray, expXRay, params);
	calcScatteringGradient(Protein, xray, expXRay, cartGrad);
	calcRotationalDerivative(backupProtein, degOfFreedom[X], degOfFreedom[Y], degOfFreedom[Z], cartGrad, rotGrad[X], rotGrad[Y], rotGrad[Z]);
	setLength(rotGrad, -length);	

	args.Protein=backupProtein;
	args.referenceProtein=Protein;
	args.cartGrad=cartGrad;

	minimizeAlongVector(degOfFreedom, rotGrad, args, func_to_minimize, maxSteps);
	Protein=backupProtein;
	rotateAtoms(Protein.Atoms, 0, 0, degOfFreedom[Z]);	
	rotateAtoms(Protein.Atoms, 0, degOfFreedom[Y], 0);	
	rotateAtoms(Protein.Atoms, degOfFreedom[X], 0, 0);	
	Proteins[0]=Protein;
	score=calcScattering(Proteins, xray, expXRay, params);

	if (score>scoreOld)
	{
		degOfFreedom=backUpDOF;
		return false;
	}
	else
	{
		return true;
	}
}

void minimizeRotationFast(ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom)
{
	bool worked;
	Vector dRot;
	timeval start, end;

	gettimeofday(&start, NULL);
	SafeAlloc(dRot, 3, "dRot");

	dRot[X]=degOfFreedom[X];
	dRot[Y]=degOfFreedom[Y];
	dRot[Z]=degOfFreedom[Z];
	moveAtoms(Protein.Atoms, degOfFreedom[3+X], degOfFreedom[3+Y], degOfFreedom[3+Z]);
	while (true)
	{
		worked=minimizeRotationUsingEstimatedScore(Protein, xray, expXRay, params, dRot);
		if (!worked) break;
	}
	degOfFreedom[X]=dRot[X];
	degOfFreedom[Y]=dRot[Y];
	degOfFreedom[Z]=dRot[Z];
	gettimeofday(&end, NULL);

}

Real calcAnalyticalRotationalGrad(Vector &degOfFreedom, ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &rotGrad)
{
	Matrix cartGrad;
	ProteinStruct tempProtein;

	SafeAlloc(rotGrad, 3, "rotGrad");
	tempProtein=Protein;

	rotateAtoms(Protein.Atoms, 0, 0, degOfFreedom[Z]);	
	rotateAtoms(Protein.Atoms, 0, degOfFreedom[Y], 0);	
	rotateAtoms(Protein.Atoms, degOfFreedom[X], 0, 0);	

	calcScatteringBfactor(Protein, xray);
	calcScatteringGradient(Protein, xray, expXRay, cartGrad);
	calcRotationalDerivative(tempProtein, degOfFreedom[X], degOfFreedom[Y], degOfFreedom[Z], cartGrad, rotGrad[X], rotGrad[Y], rotGrad[Z]);

	return 0.001;
}

Real minimizeRotationAnalytical(vector<Real> &coefficient, MRStruct args, Real dScoreCriteria, Real maxHours)
{
	//Algorithm from en.wikipedia.org/wiki/Energy_minimization
	int Size=coefficient.size();
	Real dScore, score, score_old, length=0.1;
	vector <Real> step, step_old, derivative, derivative_old;
	Real (*func_to_minimize)(vector<Real>, MRStruct)=calcFitness;
	score_old=func_to_minimize(coefficient, args);
	timeval start, end;

	gettimeofday(&start, NULL);
	SafeAlloc(step, Size, "step");
	SafeAlloc(derivative, Size, "derivative");
	SafeAlloc(step_old, Size, "step_old");
	calcAnalyticalRotationalGrad(coefficient, args.Protein, args.xray, args.expXRay, args.params, derivative);
	step=derivative;
	scaleVector(step, length);
	while (true)
	{
		score=minimizeAlongVector(coefficient, step, args, func_to_minimize);
		dScore=score-score_old;
		if (abs(dScore)<abs(score)*dScoreCriteria && score!=0)
		{
			ZeroVector(step_old);
			calcAnalyticalRotationalGrad(coefficient, args.Protein, args.xray, args.expXRay, args.params, derivative);
			step=derivative;
			score=minimizeAlongVector(coefficient, step, args, func_to_minimize);
		}
		dScore=score-score_old;
		score_old=score;
		gettimeofday(&end, NULL);
		score_old=score;
		gettimeofday(&end, NULL);
		if (abs(dScore)<abs(score)*dScoreCriteria && score!=0) break;
		else if (abs(dScore)<dScoreCriteria) break;
		if (calcTimeDiff(start, end)/3600.0>maxHours) break;
		derivative_old=derivative;
		length=calcAnalyticalRotationalGrad(coefficient, args.Protein, args.xray, args.expXRay, args.params, derivative);
		updateStep(derivative, derivative_old, step, step_old);
		step_old=step;
		scaleVector(step, length);
	}
	return score;
}

Real minimizeRotationAnalytical(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	Real score;
	Real dScoreCriteria=0.0000001;
	Real maxHours=0.5;
	Vector degOfFreedom;
	MRStruct args;

	args.Protein=Protein;
	args.xray=xray;
	args.expXRay=expXRay;
	args.params=params;

	SafeAlloc(degOfFreedom, 3, "degOfFreedom");
	calcUnitCellVectors(Protein, params);
	score=minimizeRotationAnalytical(degOfFreedom, args, dScoreCriteria, maxHours);
	rotateAtomsNoCenter(Protein.Atoms, 0, 0, degOfFreedom[Z]);	
	rotateAtomsNoCenter(Protein.Atoms, 0, degOfFreedom[Y], 0);	
	rotateAtomsNoCenter(Protein.Atoms, degOfFreedom[X], 0, 0);	

	return score;
}

Real minimizeRotationTranslationAnalytical(Vector &degOfFreedom, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Real maxHours)
{
	Real score;
	Real dScoreCriteria=10e-7;
	MRStruct args, gradArgs;
	Real (*func_to_minimize)(Vector, MRStruct)=calcRotationTranslationFitness;	
	Real (*grad_func)(Vector, MRStruct, Vector &)=calcAnalyticalRotationTranslationGrad;	

	args.Proteins=Proteins;
	args.xray=xray;
	args.expXRay=expXRay;
	args.params=params;

	gradArgs=args;

	score=optimizeCoefficients(degOfFreedom, args, gradArgs, func_to_minimize, grad_func, dScoreCriteria, maxHours);
	//score=optimizeCoefficients(degOfFreedom, args, func_to_minimize, dScoreCriteria, maxHours);
	//score=nonlinearConjugateGradient(degOfFreedom, args, gradArgs, func_to_minimize, grad_func, dScoreCriteria, maxHours);

	return score;
}

bool isLocalMinima(Array3D &scores, int i, int j, int k)
{
	int Size1, Size2, Size3;
	int startI, startJ, startK;
	int endI, endJ, endK;
	int realI, realK, realJ;

	Get3DVectorSize(scores, Size1, Size2, Size3, "scores");
	startI=i-1;
	startJ=j-1;
	startK=k-1;
	endI=i+1;
	endJ=j+1;
	endK=k+1;
	for (int a=startI;a<=endI;a++)
	{
		for (int b=startJ;b<=endJ;b++)
		{
			for (int c=startK;c<=endK;c++)
			{
				if (a<0) realI=Size1-1;
				else if (a>=Size1) realI=0;
				else realI=a;	

				if (b<0) realJ=Size2-1;
				else if (b>=Size2) realJ=0;
				else realJ=b;

				if (c<0) realK=Size3-1;
				else if (c>=Size3) realK=0;
				else realK=c;

				if (scores[realI][realJ][realK]<scores[i][j][k]) return false;				
			}
		}
	}
	return true;
}

void findLocalMinima(vector< Array3D > &degOfFreedoms, Array3D &scores, Matrix &localMinima, Vector &localMinimaScores)
{
	bool minima;
	int Size1, Size2, Size3;
	int Size1b, Size2b, Size3b;


	Get3DVectorSize(scores, Size1, Size2, Size3);
	Get3DVectorSize(degOfFreedoms, Size1b, Size2b, Size3b);

	if (Size1!=Size1b || Size2!=Size2b || Size3!=Size3b)
	{
		string errorStr;
		errorStr="degOfFreedoms array and scores array are not the ";
		errorStr+"size Size1= "+toStr(Size1)+" Size1b= "+toStr(Size1b);
		errorStr+" Size2= "+toStr(Size2)+" Size2b= "+toStr(Size2b);
		errorStr+" Size3= "+toStr(Size3)+" Size3b= "+toStr(Size3b);
		error(errorStr, __LINE__, __FILE__);
	}

	for (int i=0;i<Size1;i++)
	{
		for (int j=0;j<Size2;j++)
		{
			for (int k=0;k<Size3;k++)
			{
				minima=isLocalMinima(scores, i, j, k);
				if (minima)
				{
					SafePushBack(localMinima, degOfFreedoms[i][j][k], "localMinima");
					SafePushBack(localMinimaScores, scores[i][j][k], "localMinimaScores");
					//printVector(degOfFreedoms[i][j][k], "dof");
				}
			}
		}
	}

}

void getBestResults(Matrix &localMinima, Vector &localMinimaScores)
{
	int Size1=localMinimaScores.size();
	int Size2=localMinima.size();
	Real ave, best, w=0;
	Vector tempScores;
	Matrix tempLocalMinima;

	if (Size1!=Size2 || Size1==0)
	{
		error("Vectors don't match", __LINE__, __FILE__);
	}

	ave=calcAverage(localMinimaScores);
	findMin(localMinimaScores, best);
	for (int i=0;i<Size1;i++)
	{
		if (localMinimaScores[i]<ave*(1.0-w)+best*w)
		{
			SafePushBack(tempScores, localMinimaScores[i], "tempScores");
			SafePushBack(tempLocalMinima, localMinima[i], "tempLocalMinima");
		}
	}
	localMinima=tempLocalMinima;
	localMinimaScores=tempScores;
}

Real MRRotate(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params, Matrix &localMinima, Vector &degOfFreedom)
{
	Real rfactor;
	Vector localMinimaScores;
	vector<Array3D> degOfFreedoms;
	Array3D scores;
	vector<ProteinStruct> Proteins;

	SafeAlloc(degOfFreedom, 3, "degOfFreedom");
	MRRotationalGridSearch(degOfFreedom, Protein, xray, expXRay, lattice, params, degOfFreedoms, scores);
	findLocalMinima(degOfFreedoms, scores, localMinima, localMinimaScores);
	getBestResults(localMinima, localMinimaScores);
	//minimizeRotation(Protein, xray, expXRay, params, degOfFreedom);
	for (int i=0;i<3;i++)
	{
	}
	//return phase(Protein, xray, expXRay, lattice, params);
	//params.XRayScoreType="RFactor";
	SafePushBack(Proteins, Protein, "Proteins");
	rfactor=calcFitness(Proteins, xray, expXRay, lattice, params);
	return rfactor;
}

Real MRTranslationGridSearch(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &dPos)
{
	int steps=params.TranslationGridSearchPoints;
	Real score, bestScore, oldScore;
	Real xinc, yinc, zinc;
	Vector bestDPos;
	ProteinStruct tempProtein;
	vector<ProteinStruct> Proteins, tempProteins;
	timeval start, end;

	gettimeofday(&start, NULL);
	xinc=Protein.XBoxLength/steps;
	yinc=Protein.YBoxLength/steps;
	zinc=Protein.ZBoxLength/steps;

	SafePushBack(Proteins, Protein, "Proteins");
	SafePushBack(tempProteins, Protein, "tempProteins");
	ZeroVector(dPos);
	tempProtein=Protein;
	xray.symReal.clear();
	xray.symImag.clear();
	oldScore=calcScattering(Proteins, xray, expXRay, params);
	bestScore=calcScattering(tempProteins, xray, expXRay, params);
	bestDPos=dPos;
	for (int i=0;i<steps;i++)
	{
		dPos[X]=xinc*Real(i);
		for (int j=0;j<steps;j++)
		{
			dPos[Y]=yinc*Real(j);
			for (int k=0;k<steps;k++)
			{
				dPos[Z]=zinc*Real(k);
				score=calcFitnessTranslation(Proteins, xray, expXRay, params, dPos);
				if (score<bestScore)
				{
					bestScore=score;
					bestDPos=dPos;
				}
			}
		}
	}
	dPos=bestDPos;
	gettimeofday(&end, NULL);
	return bestScore;
}

Real minimizeTranslation(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom)
{
	//minimizeTranslationOnly takes in a vector, degOfFreedom, containing
	//both rotational and translational degrees of freedom, performs the
	//rotations, and then optimizes the translational degrees of freedom.
	//minimizeTranslation takes in a vector containing only the 
	//translational degrees of freedom and optimizes them.
	Real dScoreCriteria=1e-5, maxHours=0.01, score;
	MRStruct args;
	Real (*func_to_minimize)(Vector, MRStruct)=calcFitnessTranslation;
	calcCartToFrac(Proteins[0]);
	args.Proteins=Proteins;
	args.xray=xray;
	args.expXRay=expXRay;
	args.params=params;
	score=conjugateGradient(degOfFreedom, args, func_to_minimize, dScoreCriteria, maxHours);
	//score=optimizeUsingEstimatedGrad(degOfFreedom, args, func_to_minimize, dScoreCriteria, maxHours);
	return score;
}

Real minimizeTranslationOnly(vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom)
{
	//minimizeTranslationOnly takes in a vector, degOfFreedom, containing
	//both rotational and translational degrees of freedom, performs the
	//rotations, and then optimizes the translational degrees of freedom.
	//minimizeTranslation takes in a vector containing only the
	//translational degrees of freedom and optimizes them.
	Real score;
	Vector dPos;
	Matrix rotationMatrix;
	timeval start, end;
	gettimeofday(&start, NULL);
	SafeAlloc(dPos, 3, "dPos");
	dPos[X]=degOfFreedom[3+X];
	dPos[Y]=degOfFreedom[3+Y];
	dPos[Z]=degOfFreedom[3+Z];
	rotateAtoms(Proteins[0].Atoms, 0, 0, degOfFreedom[Z]);
	rotateAtoms(Proteins[0].Atoms, 0, degOfFreedom[Y], 0);
	rotateAtoms(Proteins[0].Atoms, degOfFreedom[X], 0, 0);
	//degOfFreedom[3+X]=0;
	//degOfFreedom[3+Y]=0;
	//degOfFreedom[3+Z]=0;
	if (params.UseFastTranslationRotation)
	{
		degOfFreedomToRotMat(Proteins, degOfFreedom);
		calcRotationMatrix(rotationMatrix, degOfFreedom[X], degOfFreedom[Y], degOfFreedom[Z]);
		calcFastTranslationRotation(Proteins[0], xray, rotationMatrix, 0, params);
		xray.complexSymTempReal[0]=xray.complexSymReal[0];
		xray.complexSymTempImag[0]=xray.complexSymImag[0];
	}
	else calcScatteringComplex(Proteins, xray, degOfFreedom);
	/*
	   calcTotalFormFactor(xray);
	   calcIntensity(xray);
	   calcAmplitude(xray);
	   int nmiller=xray.miller.size();
	   for (int i=0;i<nmiller;i++)
	   {
	   }
	 */
	score=minimizeTranslation(Proteins, xray, expXRay, params, dPos);
	degOfFreedom[3+X]=dPos[X];
	degOfFreedom[3+Y]=dPos[Y];
	degOfFreedom[3+Z]=dPos[Z];
	gettimeofday(&end, NULL);

	return score;
}

Real minimizeRotationTranslation(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom)
{
	Real dScoreCriteria=0.0000001, maxHours=200.0, score;
	MRStruct args;
	Real (*func_to_minimize)(Vector, MRStruct)=calcRotationTranslationFitness;
	args.Protein=Protein;
	args.xray=xray;
	args.expXRay=expXRay;
	args.params=params;
	SafePushBack(args.Proteins, Protein, "Protein");
	score=conjugateGradient(degOfFreedom, args, func_to_minimize, dScoreCriteria, maxHours);
	//RotateAtomsNoCenter(Protein.Atoms, 0, 0, degOfFreedom[2]);
	//RotateAtomsNoCenter(Protein.Atoms, 0, degOfFreedom[1], 0);
	//RotateAtomsNoCenter(Protein.Atoms, degOfFreedom[0], 0, 0);
	//moveAtoms(Protein.Atoms, degOfFreedom[3], degOfFreedom[4], degOfFreedom[5]);
	return score;
}

Real MRTranslate(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &dPos)
{
	Real score;

	SafeAlloc(dPos, 3, "dPos");

	score=MRTranslationGridSearch(Protein, xray, expXRay, params, dPos);

	return score;
}

void makeRandomRotationTranslation(ProteinStruct &Protein)
{
	Real phi, psi, theta;
	Real dx, dy, dz;

	phi=randDouble(-pi, pi);
	psi=randDouble(-pi, pi);
	theta=randDouble(-pi, pi);
	dx=randDouble(Protein.XBoxLength);
	dy=randDouble(Protein.YBoxLength);
	dz=randDouble(Protein.ZBoxLength);

	//RotateAtoms(Protein.Atoms, phi, psi, theta);
	moveAtoms(Protein.Atoms, dx, dy, dz);
}

Real MRTranslate(ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom, Vector &dPos)
{
	Real score;

	SafeAlloc(dPos, 3, "dPos");
	printVector(degOfFreedom);
	rotateAtoms(Protein.Atoms, 0, 0, degOfFreedom[Z]);	
	rotateAtoms(Protein.Atoms, 0, degOfFreedom[Y], 0);	
	rotateAtoms(Protein.Atoms, degOfFreedom[X], 0, 0);	
	//RotateAtomsNoCenter(Protein.Atoms, 0, 0, degOfFreedom[X]);	
	//RotateAtomsNoCenter(Protein.Atoms, 0, degOfFreedom[Y], 0);	
	//RotateAtomsNoCenter(Protein.Atoms, degOfFreedom[Z], 0, 0);	
	score=MRTranslate(Protein, xray, expXRay, params, dPos);
	return score;
}

Real MRTranslate(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Matrix &degOfFreedoms, Vector &bestDegOfFreedom, Vector &bestDPos)
{
	string outputFile;
	int ncandidates=degOfFreedoms.size();
	int count=0;
	Real bestScore, score;
	vector<ProteinStruct> Proteins;
	Vector dPos, degOfFreedom;

	SafePushBack(Proteins, Protein, "Proteins");
	SafeAlloc(degOfFreedom, 6, "degOfFreedom");
	bestScore=MRTranslate(Protein, xray, expXRay, params, degOfFreedoms[0], dPos);
	outputFile=params.ReplicaStructuresOutput+"_0.pdb";
	printReplica(degOfFreedom, Proteins, outputFile);
	bestDegOfFreedom=degOfFreedoms[0];
	bestDPos=dPos;
	for (int i=1;i<ncandidates;i++)
	{
		score=MRTranslate(Protein, xray, expXRay, params, degOfFreedoms[i], dPos);
		if (score<bestScore)
		{
			bestScore=score;
			bestDegOfFreedom=degOfFreedoms[i];
			bestDPos=dPos;
			Proteins[0]=Protein;
			outputFile=params.ReplicaStructuresOutput+"_"+toStr(i)+".pdb";
			for (int j=0;j<3;j++) degOfFreedom[j]=bestDegOfFreedom[j];
			for (int j=0;j<3;j++) degOfFreedom[j]=dPos[j];
			//printReplica(degOfFreedom, Proteins, outputFile);
			count++;
		}
		Proteins[0]=Protein;
		outputFile=params.ReplicaStructuresOutput+"_"+toStr(i)+".pdb";
		for (int j=0;j<3;j++) degOfFreedom[j]=bestDegOfFreedom[j];
		for (int j=0;j<3;j++) degOfFreedom[j]=dPos[j];
		printReplica(degOfFreedom, Proteins, outputFile);
	}
	return bestScore;
}

Real MRRotateTranslate(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	int noperations=Protein.symmetryOperations.size();
	Real rfactor, maxHours=0.5;
	Vector bestDegOfFreedom, bestDPos;
	Matrix degOfFreedoms; //Contains possible solutions to rotation problem.
	vector<ProteinStruct> Proteins;

	SafeAlloc(bestDPos, 3, "bestDPos");
	makeRandomRotationTranslation(Protein);
	MRRotate(Protein, xray, expXRay, lattice, params, degOfFreedoms, bestDegOfFreedom);
	if (noperations>1) MRTranslate(Protein, xray, expXRay, params, degOfFreedoms, bestDegOfFreedom, bestDPos);
	else minimizeRotation(Protein, xray, expXRay, params, bestDegOfFreedom, maxHours);
	rotateAtoms(Protein.Atoms, 0, 0, bestDegOfFreedom[Z]);	
	rotateAtoms(Protein.Atoms, 0, bestDegOfFreedom[Y], 0);	
	rotateAtoms(Protein.Atoms, bestDegOfFreedom[X], 0, 0);	
	moveAtoms(Protein.Atoms, bestDPos[X], bestDPos[Y], bestDPos[Z]);
	SafePushBack(Proteins, Protein, "Proteins");
	rfactor=calcScattering(Proteins, xray, expXRay, params);
	return rfactor;
}

void printDofs(Matrix &degOfFreedoms)
{
	int nrep, ndof;

	Get2DVectorSize(degOfFreedoms, nrep, ndof, "degOfFreedoms");
	for (int i=0;i<nrep;i++)
	{
		for (int j=0;j<ndof;j++)
		{
		}
	}
}

void initializeDegOfFreedoms(Matrix &degOfFreedoms, string degOfFreedomFile)
{
	vector<string> lines;
	int nlines, nrep, ndegOfFreedom, degOfFreedom, rep;
	Real value;

	Get2DVectorSize(degOfFreedoms, nrep, ndegOfFreedom);
	readLines(degOfFreedomFile, lines, "degOfFreedomFile");
	nlines=lines.size();

	for (int i=0;i<nlines;i++)
	{
		rep=toInt(GetWord2(lines[i], 1));
		degOfFreedom=toInt(GetWord2(lines[i], 2));
		value=toReal(GetWord2(lines[i], 3));
		if (rep<nrep && degOfFreedom<ndegOfFreedom)
		{
			degOfFreedoms[rep][degOfFreedom]=value;
		}
	}
	printDofs(degOfFreedoms);
}

void fixDegOfFreedoms(Vector &degOfFreedoms, Real x, Real y, Real z)
{
	int Size=degOfFreedoms.size();
	int nprot=Size/6;
	int n;

	for (int i=0;i<nprot;i++)
	{
		for (int j=0;j<3;j++)
		{
			n=int((abs(degOfFreedoms[i*6+j])+pi)/2.0*pi);
			if (degOfFreedoms[i*6+j]>pi) 
			{
				degOfFreedoms[i*6+j]-=2.0*pi*Real(n);
			}
			if (degOfFreedoms[i*6+j]<-pi) 
			{
				degOfFreedoms[i*6+j]+=2.0*pi*Real(n);
			}
		}
		if (degOfFreedoms[i*6+3]>x*0.5) degOfFreedoms[i*6+3]-=x;
		if (degOfFreedoms[i*6+3]<-x*0.5) degOfFreedoms[i*6+3]+=x;
		if (degOfFreedoms[i*6+4]>y*0.5) degOfFreedoms[i*6+4]-=y;
		if (degOfFreedoms[i*6+4]<-y*0.5) degOfFreedoms[i*6+4]+=y;
		if (degOfFreedoms[i*6+5]>z*0.5) degOfFreedoms[i*6+5]-=z;
		if (degOfFreedoms[i*6+5]<-z*0.5) degOfFreedoms[i*6+5]+=z;
	}
}

void fixDegOfFreedoms(Matrix &degOfFreedoms, Real x, Real y, Real z)
{
	int Size1=degOfFreedoms.size();

	for (int i=0;i<Size1;i++)
	{
		fixDegOfFreedoms(degOfFreedoms[i], x, y, z);
	}


}

int pickRandom(MoveStruct &move)
{
	int Size=move.probSum.size();
	int pick=UNK_INT;
	Real rand;

	rand=randDouble(0, 1.0);
	for (int i=0;i<Size-1;i++)
	{
		if (move.probSum[i]<=rand && move.probSum[i+1]>=rand)
		{
			pick=i;
		}
	}
	if (pick==UNK_INT) pick=Size-1;

	return pick;
}

void fixSpaceGroupPick(MoveStruct &move, string spaceGroup, int &pick)
{
	if (spaceGroup=="P3" || spaceGroup=="P31" || spaceGroup=="P32" || spaceGroup=="R3" || spaceGroup=="H3" || spaceGroup=="P4" || spaceGroup=="P41" || spaceGroup=="P42" || spaceGroup=="P43" || spaceGroup=="I4" || spaceGroup=="I41" || spaceGroup=="P6" || spaceGroup=="P61" || spaceGroup=="P65" || spaceGroup=="P62" || spaceGroup=="P64" || spaceGroup=="P63")
	{
		while (pick==Z+3)
		{
			pick=pickRandom(move);
		}
	}
	else if (spaceGroup=="P2" || spaceGroup=="P21" || spaceGroup=="C2" || spaceGroup=="A2" || spaceGroup=="I2" || spaceGroup=="P1211" || spaceGroup=="C121")
	{
		while (pick==Y+3)
		{
			pick=pickRandom(move);
		}
	}
	if (spaceGroup=="P1")
	{
		while (pick==X+3 || pick==Y+3 || pick==Z+3)
		{
			pick=pickRandom(move);
		}
	}
}

int pickRandomDOF(MoveStruct &move, string spaceGroup, int nprot)
{
	int pick=pickRandom(move);
	if (pick<nprot*6) fixSpaceGroupPick(move, spaceGroup, pick);
	return pick;
}

int pickRandomDOF(MoveStruct &move, string spaceGroup)
{
	int pick=pickRandom(move);
	fixSpaceGroupPick(move, spaceGroup, pick);
	return pick;
}

void fixMoveSize(int pick, Real xlength, Real ylength, Real zlength, Real &moveSize, Vector &coefficient)
{
	int n;

	if (pick%6<3)
	{
		if (moveSize>2.0*pi)
		{
			n=int(moveSize/(2.0*pi));
			moveSize-=Real(n)*2.0*pi;
		}
		if (moveSize<-2.0*pi)
		{
			n=int(abs(moveSize)/(2.0*pi));
			moveSize+=Real(n)*2.0*pi;
		}
	}
	if (pick%6==3)
	{
		if (moveSize>xlength)
		{
			n=int(moveSize/xlength);
			moveSize-=Real(n)*xlength;
		}
		if (moveSize<xlength)
		{
			n=int(abs(moveSize)/xlength);
			moveSize+=Real(n)*xlength;
		}
	}
	if (pick%6==4)
	{
		if (moveSize>ylength)
		{
			n=int(moveSize/ylength);
			moveSize-=Real(n)*ylength;
		}
		if (moveSize<ylength)
		{
			n=int(abs(moveSize)/ylength);
			moveSize+=Real(n)*ylength;
		}
	}
	if (pick%6==5)
	{
		if (moveSize>zlength)
		{
			n=int(moveSize/zlength);
			moveSize-=Real(n)*zlength;
		}
		if (moveSize<zlength)
		{
			n=int(abs(moveSize)/zlength);
			moveSize+=Real(n)*zlength;
		}
	}

}

void fixDOF(Real xlength, Real ylength, Real zlength, Vector &coefficient)
{
	bool change=false;
	int mod;
	int ndof=coefficient.size();
	for (int i=0;i<ndof;i++)
	{
		if (i%6<3)
		{
			if (coefficient[i]<-pi) 
			{
				mod=int(abs(coefficient[i])/(2.0*pi))+1;
				coefficient[i]+=2.0*pi*Real(mod);
				change=true;
			}
			if (coefficient[i]>pi) 
			{
				mod=int(abs(coefficient[i])/(2.0*pi))+1;
				coefficient[i]-=2.0*pi*Real(mod);
				change=true;
			}
		}
		if (i%6==3)
		{
			if (coefficient[i]<-xlength*0.5) 
			{
				mod=int(abs(coefficient[i])/(xlength))+1;
				coefficient[i]+=xlength*Real(mod);
				if (coefficient[i]>xlength*0.5)
				{
					coefficient[i]-=xlength;
				}
				change=true;
			}
			if (coefficient[i]>xlength*0.5) 
			{
				mod=int(abs(coefficient[i])/(xlength))+1;
				coefficient[i]-=xlength*Real(mod);
				change=true;
			}
		}
		if (i%6==4)
		{
			if (coefficient[i]<-ylength*0.5) 
			{
				mod=int(abs(coefficient[i])/(ylength))+1;
				coefficient[i]+=ylength*Real(mod);
				if (coefficient[i]>ylength*0.5)
				{
					coefficient[i]-=ylength;
				}
				change=true;
			}
			if (coefficient[i]>ylength*0.5) 
			{
				mod=int(abs(coefficient[i])/(ylength))+1;
				coefficient[i]-=ylength*Real(mod);
				change=true;
			}
		}
		if (i%6==5)
		{
			if (coefficient[i]<-zlength*0.5) 
			{
				mod=int(abs(coefficient[i])/(zlength))+1;
				coefficient[i]+=zlength*Real(mod);
				if (coefficient[i]>zlength*0.5)
				{
					coefficient[i]-=zlength;
				}
				change=true;
			}
			if (coefficient[i]>zlength*0.5) 
			{
				mod=int(abs(coefficient[i])/(zlength))+1;
				coefficient[i]-=zlength*Real(mod);
				change=true;
			}
		}
	}		
	if (change) fixDOF(xlength, ylength, zlength, coefficient);
}

bool checkIfAllZero(Vector &coefficient, int nprot)
{
	int ndof=coefficient.size();
	for (int i=nprot*6;i<ndof;i++)
	{
		if (coefficient[i]!=0) return false;
	}
	return true;
}

void fixDOF(Vector &coefficient, int nprot)
{
	bool allZero;
	int ndof=coefficient.size();

	for (int i=nprot*6;i<ndof;i++)
	{
		if (coefficient[i]<0) coefficient[i]=0;
		if (coefficient[i]>1.0) coefficient[i]=1.0; 
	}
	allZero=checkIfAllZero(coefficient, nprot);
	if (allZero) coefficient[nprot*6]=1.0;
}

void fixMoveSizeOccupancy(Real &moveSize, Vector &coefficient)
{
	if (moveSize>1.0) moveSize=0.5;
	else if (moveSize<-1.0) moveSize=-0.5;
}

void fixMoveSize(int pick, Real xlength, Real ylength, Real zlength, Real &moveSize, Vector &coefficient, int nprot)
{
	if (pick<nprot*6) fixMoveSize(pick, xlength, ylength, zlength, moveSize, coefficient);
	else fixMoveSizeOccupancy(moveSize, coefficient);
}

void makeRandomChange(Vector &coefficient, MoveStruct &move, int &pick, Real xlength, Real ylength, Real zlength, string spaceGroup, int nprot)
{
	Real moveSize, moveSizeProb;

	pick=pickRandomDOF(move, spaceGroup, nprot);
	moveSize=randGaussian(move.width[pick]);
	fixMoveSize(pick, xlength, ylength, zlength, moveSize, coefficient, nprot);
	//coefficient[pick]+=moveSize;
	if (pick<6) coefficient[pick]+=moveSize;
	else
	{
		if (coefficient[pick]==0) coefficient[pick]=1.0;
		else coefficient[pick]=0;
	}
	fixDOF(xlength, ylength, zlength, coefficient);
	fixDOF(coefficient, nprot);
	moveSizeProb=calcGaussian(moveSize, move.width[pick]);
	SafePushBack(move.index, pick, "index");
	SafePushBack(move.size, moveSize, "size");
	SafePushBack(move.moveSizeProb, moveSizeProb, "moveSizeProb");
}

Real getRx(Matrix &rotMat, Real ry)
{
	Real rx, mzz, d=1.0e-6;


	rx=asin(-rotMat[Y][Z]/cos(ry));

	mzz=cos(ry)*cos(rx);

	if (abs(mzz+rotMat[Z][Z])<d) rx=pi-rx;
	else if (abs(mzz-rotMat[Z][Z])>d)
	{
		string errorStr="ry= "+toStr(ry)+" rx= "+toStr(rx);
		errorStr+=" mzz= "+toStr(mzz)+".  Cannot calculate angles";
		errorStr+=" for rotation matrix";
		printMatrix(rotMat);
		error(errorStr, __LINE__, __FILE__);
	}
	return rx;
}

Real getRz(Matrix &rotMat, Real ry)
{
	Real rz, mxy, d=1.0e-6;

	rz=acos(rotMat[X][X]/cos(ry));
	mxy=-sin(rz)*cos(ry);
	if (abs(mxy+rotMat[X][Y])<d) rz=-rz;
	else if (abs(mxy-rotMat[X][Y])>d)
	{
		string errorStr="ry= "+toStr(ry)+" rz= "+toStr(rz);
		errorStr+=" mxy= "+toStr(mxy)+".  Cannot calculate angles";
		errorStr+=" for rotation matrix";
		printMatrix(rotMat);
		error(errorStr, __LINE__, __FILE__);
	}
	return rz;
}

void rotationMatrixToCoefficients(Matrix &rotMat, Real &rx, Real &ry, Real &rz)
{
	//Rz(rx)Ry(ry)Rx(rz)
	Matrix rotMat2, rotMat3, rotMat4;
	Real myx, d=1.0e-6;

	ry=asin(rotMat[X][Z]);
	rx=getRx(rotMat, ry);
	rz=getRz(rotMat, ry);

	myx=cos(rx)*sin(rz)+sin(rx)*sin(ry)*cos(rz);

	if (abs(myx-rotMat[Y][X])>d) ry=pi-ry;

	rx=getRx(rotMat, ry);
	rz=getRz(rotMat, ry);

	myx=cos(rx)*sin(rz)+sin(rx)*sin(ry)*cos(rz);

	if (abs(myx-rotMat[Y][X])>d)
	{
		Safe2DAlloc(rotMat2, 3, 3, "rotMat2");
		rotMat2[X][X]=cos(rz)*cos(ry);
		rotMat2[X][Y]=-cos(ry)*sin(rz);
		rotMat2[X][Z]=sin(ry);
	
		rotMat2[Y][X]=cos(rx)*sin(rz)+sin(rx)*sin(ry)*cos(rz);
		rotMat2[Y][Y]=cos(rx)*cos(rz)-sin(rx)*sin(ry)*sin(rz);
		rotMat2[Y][Z]=-sin(rx)*cos(ry);

		rotMat2[Z][X]=sin(rx)*sin(rz)-cos(rx)*sin(ry)*cos(rz);
		rotMat2[Z][Y]=sin(rx)*cos(rz)+cos(rx)*sin(ry)*sin(rz);
		rotMat2[Z][Z]=cos(ry)*cos(rx);


		printMatrix(rotMat);
		printMatrix(rotMat2);
		endProgram(__LINE__, __FILE__);
	}


}

void calcNCSTranslation(Vector &degOfFreedomNCS, VectorStruct &v, Matrix &rotMat, Real &dx, Real &dy, Real &dz)
{
	Real a;
	Real xc, yc, zc;
	Vector dIn, dOut;

	SafeAlloc(dIn, 3, "dIn");
	SafeAlloc(dOut, 3, "dOut");

	a=v.x*(degOfFreedomNCS[3+X]-degOfFreedomNCS[8]);
	a+=v.y*(degOfFreedomNCS[3+Y]-degOfFreedomNCS[9]);
	a+=v.z*degOfFreedomNCS[3+Z];

	xc=degOfFreedomNCS[8]+a*v.x;
	yc=degOfFreedomNCS[9]+a*v.y;
	zc=a*v.z;

	xc=degOfFreedomNCS[8];
	yc=degOfFreedomNCS[9];
	zc=0;

	dIn[X]=degOfFreedomNCS[3+X]-xc;
	dIn[Y]=degOfFreedomNCS[3+Y]-yc;
	dIn[Z]=degOfFreedomNCS[3+Z]-zc;

	matrixMultiply(dIn, rotMat, dOut); 

	dx=xc+dOut[X];
	dy=yc+dOut[Y];
	dz=zc+dOut[Z];
}

void calcDOFFromNCS(Vector &degOfFreedomNCS, Vector &degOfFreedom)
{
	Matrix rotMat, rotMat1, totalRotMat;
	VectorStruct v;

	SafeAlloc(degOfFreedom, 12, "degOfFreedom");
	for (int i=0;i<6;i++)
	{
		degOfFreedom[i]=degOfFreedomNCS[i];
	}
	v.x=sin(degOfFreedomNCS[6])*cos(degOfFreedomNCS[7]);
	v.y=sin(degOfFreedomNCS[6])*sin(degOfFreedomNCS[7]);
	v.z=cos(degOfFreedomNCS[6]);

	calcRotationMatrix(v, pi, rotMat);
	calcRotationMatrix(rotMat1, degOfFreedom[X], degOfFreedom[Y], degOfFreedom[Z]);
	matrixMultiply(rotMat, rotMat1, totalRotMat);
	rotationMatrixToCoefficients(totalRotMat, degOfFreedom[6+X], degOfFreedom[6+Y], degOfFreedom[6+Z]);	
	calcNCSTranslation(degOfFreedomNCS, v, rotMat, degOfFreedom[9+X], degOfFreedom[9+Y], degOfFreedom[9+Z]);
}

void calcDOFFromNCS(Matrix &degOfFreedomNCS, Matrix &degOfFreedom)
{
	int nrep=degOfFreedomNCS.size();
	Vector tempDegOfFreedom;

	degOfFreedom.clear();
	for (int i=0;i<nrep;i++)
	{
		calcDOFFromNCS(degOfFreedomNCS[i], tempDegOfFreedom);
		SafePushBack(degOfFreedom, tempDegOfFreedom, "degOfFreedom");
	}
}

void calcDOFFromNCS(Matrix &degOfFreedomNCS)
{
	Matrix degOfFreedom;

	calcDOFFromNCS(degOfFreedomNCS, degOfFreedom);
	degOfFreedomNCS=degOfFreedom;
}

void fixMoveSizeNCS(int pick, Real xlength, Real ylength, Real zlength, Real &moveSize)
{
	int n;

	if (pick<3 || pick==8 || pick==9)
	{
		if (moveSize>2.0*pi)
		{
			n=int(moveSize/(2.0*pi));
			moveSize-=Real(n)*2.0*pi;
		}
		if (moveSize<-2.0*pi)
		{
			n=int(abs(moveSize)/(2.0*pi));
			moveSize+=Real(n)*2.0*pi;
		}
	}
	if (pick==3 || pick==6 || pick==7)
	{
		if (moveSize>xlength)
		{
			n=int(moveSize/xlength);
			moveSize-=Real(n)*xlength;
		}
		if (moveSize<xlength)
		{
			n=int(abs(moveSize)/xlength);
			moveSize+=Real(n)*xlength;
		}
	}
	if (pick==4)
	{
		if (moveSize>ylength)
		{
			n=int(moveSize/ylength);
			moveSize-=Real(n)*ylength;
		}
		if (moveSize<ylength)
		{
			n=int(abs(moveSize)/ylength);
			moveSize+=Real(n)*ylength;
		}
	}
	if (pick==5)
	{
		if (moveSize>zlength)
		{
			n=int(moveSize/zlength);
			moveSize-=Real(n)*zlength;
		}
		if (moveSize<zlength)
		{
			n=int(abs(moveSize)/zlength);
			moveSize+=Real(n)*zlength;
		}
	}
}

void makeRandomChangeNCS(Vector &degOfFreedom, Vector &degOfFreedomNCS, MoveStruct &move, int &pick, Real xlength, Real ylength, Real zlength, string spaceGroup)
{
	Real moveSize, moveSizeProb;

	pick=pickRandomDOF(move, spaceGroup);
	moveSize=randGaussian(move.width[pick]);
	fixMoveSizeNCS(pick, xlength, ylength, zlength, moveSize);
	degOfFreedomNCS[pick]+=moveSize;
	calcDOFFromNCS(degOfFreedomNCS, degOfFreedom);
	fixDOF(xlength, ylength, zlength, degOfFreedom);
	moveSizeProb=calcGaussian(moveSize, move.width[pick]);
	SafePushBack(move.index, pick, "index");
	SafePushBack(move.size, moveSize, "size");
	SafePushBack(move.moveSizeProb, moveSizeProb, "moveSizeProb");
}

Real calcScore(Vector &degOfFreedoms, vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	//int nprot=Proteins.size();
	int nprot, nsym=0, npoint;
	int scoreIndex;
	Real score=0;
	Real clashCC, clashCS, clashSS;
	Real clashScore;
	Vector degOfFreedomsNCS, dx;
	Matrix rotationMatrix;


	SafeAlloc(dx, 3, "dx");
	nprot=Proteins.size();
	if (nprot!=0)
	{
		nsym=Proteins[0].symmetryOperations.size();
	}
	else
	{
		error("No proteins ", __LINE__, __FILE__);
	}
	npoint=xray.miller.size();
	Safe3DAlloc(xray.complexSymReal, nprot, nsym, npoint, "xray.complexSymReal");
	Safe3DAlloc(xray.complexSymImag, nprot, nsym, npoint, "xray.complexSymImag");
	Safe3DAlloc(xray.complexSymTempReal, nprot, nsym, npoint, "xray.complexSymReal");
	Safe3DAlloc(xray.complexSymTempImag, nprot, nsym, npoint, "xray.complexSymImag");
	Get3DVectorSize(xray.complexSymReal, nprot, nsym, npoint);
/*
	if (params.UseNCS)
	{
		calcDOFFromNCS(degOfFreedoms, degOfFreedomsNCS);
		degOfFreedomToRotMat(Proteins, degOfFreedomsNCS);
		score=calcFitness(degOfFreedomsNCS, Proteins, xray, expXRay, lattice, params);
	}
	else 
	{
		if (params.UseFastTranslationRotation)
		{
			//sumScatteringContinuousSegments(xray, degOfFreedoms, params.NumCopies*6);
			if (params.OptimizeOccupancies) sumScatteringContinuousSegments(xray, Proteins[0], degOfFreedoms, params.NumCopies*6);
			degOfFreedomToRotMat(Proteins, degOfFreedoms);
			calcRotationMatrix(rotationMatrix, degOfFreedoms[X], degOfFreedoms[Y], degOfFreedoms[Z]);
			//calcFastTranslationRotation(Proteins, xray, params);
			calcFastTranslationRotation(Proteins[0], xray, rotationMatrix, 0, params);
			dx[X]=degOfFreedoms[X];
			dx[Y]=degOfFreedoms[Y];
			dx[Z]=degOfFreedoms[Z];
			translateXRay2(xray, Proteins[0], dx, 0);
			calcTotalFormFactor(xray);
			//score=calcFitness(degOfFreedoms, Proteins, xray, expXRay, lattice, params);
		}
		else
		{
			score=calcFitness(degOfFreedoms, Proteins, xray, expXRay, lattice, params);
		}
	}
	score=calcFitness(degOfFreedoms, Proteins, xray, expXRay, lattice, params);
*/
	if (params.OptimizeOccupancies)
	{
		//sumScatteringContinuousSegments(xray, degOfFreedom, nprot*6);
		sumScatteringContinuousSegments(xray, Proteins[0], degOfFreedoms, nprot*6);
	}
	for (int prot=0;prot<nprot;prot++)
	{
		score=calcFitnessChange(degOfFreedoms, Proteins[prot], xray, expXRay, lattice, params, prot*6);
	}
	if (params.ClashScoreUse=="DuringMR")
	{
		initializeXRayProtein(Proteins[0], xray);
		clashScore=calcClashScoreFast(degOfFreedoms, xray.Proteins, params, clashCC, clashCS, clashSS);
		scoreIndex=xray.scoreIndex[CLASH_CC];
		score+=clashCC*params.ClashWeightCoreCore/xray.stdDev[scoreIndex];
		xray.score[scoreIndex]=clashCC;

		scoreIndex=xray.scoreIndex[CLASH_CS];
		score+=clashCS*params.ClashWeightCoreSurface/xray.stdDev[scoreIndex];
		xray.score[scoreIndex]=clashCS;


		scoreIndex=xray.scoreIndex[CLASH_SS];
		score+=clashSS*params.ClashWeightSurfaceSurface/xray.stdDev[scoreIndex];
		xray.score[scoreIndex]=clashSS;
		//score+=clashScore;
	}

	//score=calcRotationTranslationFitness2(degOfFreedoms, Proteins, xray, expXRay, lattice, params);
	//score=calcFitness(degOfFreedoms, Proteins, xray, expXRay, lattice, params);
	//for (int prot=0;prot<nprot;prot++)
	//{
	//	score=calcFitnessChange(degOfFreedoms, Proteins[prot], xray, expXRay, lattice, params, prot*6);
	//}
	//calcTotalFormFactor(xray);
	//calcIntensity(xray);
	//score=calcRFactor(xray, expXRay, params);
	return score;
}

void calcScores(Matrix &degOfFreedoms, Vector &score, Matrix &scores, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	int nrep=degOfFreedoms.size();
	for (int i=0;i<nrep;i++)
	{
		score[i]=calcScore(degOfFreedoms[i], Proteins, xray, expXRay, lattice, params);
		scores[i]=xray.score;
	}
}

void calcStandardDeviations(Matrix &scores, Vector &stdDev)
{
	int Size1, Size2;
	Vector ave, aveSqr;

	Get2DVectorSize(scores, Size1, Size2, "scores");
	
	SafeAlloc(stdDev, Size2, "stdDev");
	SafeAlloc(ave, Size2, "ave");
	SafeAlloc(aveSqr, Size2, "aveSqr");
	for (int i=0;i<Size1;i++)
	{
		for (int j=0;j<Size2;j++)
		{
			ave[j]+=scores[i][j];		
			aveSqr[j]+=scores[i][j]*scores[i][j];
		}
	}

	for (int i=0;i<Size2;i++)
	{
		ave[i]/=Real(Size1);
		aveSqr[i]/=Real(Size1);
		stdDev[i]=sqrt(aveSqr[i]-ave[i]*ave[i]);
		if (stdDev[i]==0) stdDev[i]=1.0;
	}
}

Real monteCarloNCS(Vector &degOfFreedomNCS, vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params, Real temp, int nsteps, MoveStruct &move)
{
	int pick, nmove;
	int updateMoveSizeInterval=50;
	int nMoveTypes=degOfFreedomNCS.size();
	int nprot=Proteins.size();
	Real score, score_old, score_initial=0, dScore, randNum;
	Real acceptProb, acceptRate, time;
	Vector degOfFreedom;
	//Real scoreSlow;
	Vector old_degOfFreedom, subV;
	XRayStruct oldXRay;
	//XRayStruct slowXRay;
	timeval start, end;


	calcDOFFromNCS(degOfFreedomNCS, degOfFreedom);
	fixDOF(Proteins[0].XBoxLength, Proteins[0].YBoxLength, Proteins[0].ZBoxLength, degOfFreedom);

	if (params.UseFastTranslationRotation)
	{
		degOfFreedomToRotMat(Proteins, degOfFreedom);
		calcFastTranslationRotation(Proteins, xray, params);
	}
	else
	{
		calcScatteringComplex(Proteins, xray, degOfFreedom);
	}
	//slowXRay=xray;
	score_initial=calcFitness(degOfFreedom, Proteins, xray, expXRay, lattice, params);
	for (int prot=0;prot<nprot;prot++)
	{
		//score_initial=calcFitnessChange(degOfFreedom, ProteinsOriginal[prot], xray, expXRay, lattice, params, prot*6);
		score_initial=calcFitnessChange(degOfFreedom, Proteins[prot], xray, expXRay, lattice, params, prot*6);
		score_initial=calcFitnessChange(degOfFreedom, Proteins[prot], xray, expXRay, lattice, params, prot*6+3);
	}
	//params.UseFastTranslationRotation=false;
	//scoreSlow=calcFitness(degOfFreedom, Proteins, slowXRay, expXRay, lattice, params);
	//params.UseFastTranslationRotation=true;
/*
	int npoint=xray.complexSymReal[0][0].size();
	int nsym=xray.complexSymReal[0].size();

	for (int i=0;i<nprot;i++)
	{
		for (int j=0;j<nsym;j++)
		{
			for (int k=0;k<npoint;k++)
			{
			}
		}
	}
	endProgram(__LINE__, __FILE__);
*/
	score=score_initial;
	oldXRay.complexSymReal=xray.complexSymReal;
	oldXRay.complexSymImag=xray.complexSymImag;
	oldXRay.complexSymTempReal=xray.complexSymTempReal;
	oldXRay.complexSymTempImag=xray.complexSymTempImag;
	score_old=score_initial;
	old_degOfFreedom=degOfFreedomNCS;
	for (int i=0;i<nsteps;i++)
	{
		gettimeofday(&start, NULL);
		makeRandomChangeNCS(degOfFreedom, degOfFreedomNCS, move, pick, Proteins[0].XBoxLength, Proteins[0].YBoxLength, Proteins[0].ZBoxLength, Proteins[0].spaceGroup);
		for (int j=0;j<nprot;j++)
		{
			score=calcFitnessChange(degOfFreedom, Proteins[j], xray, expXRay, lattice, params, j*6);
			score=calcFitnessChange(degOfFreedom, Proteins[j], xray, expXRay, lattice, params, j*6+3);
		}
		//params.UseFastTranslationRotation=false;
		//scoreSlow=calcFitness(degOfFreedom, Proteins, slowXRay, expXRay, lattice, params);
		//params.UseFastTranslationRotation=true;
		//printVector(degOfFreedom);
		gettimeofday(&end, NULL);
		time=calcTimeDiff(start, end);
		dScore=score-score_old;
		acceptProb=calcAcceptProb(dScore, temp);
		SafePushBack(move.score, score, "score");
		SafePushBack(move.dScore, dScore, "dScore");
		SafePushBack(move.temp, temp, "temp");
		SafePushBack(move.acceptProb, acceptProb, "acceptProb");
		SafePushBack(move.time, time, "time");
		randNum=randDouble(0,1.0);
		if (acceptProb>randNum)
		{
			score_old=score;
			old_degOfFreedom=degOfFreedomNCS;
			//oldXRay.complexSymReal[prot]=xray.complexSymReal[prot];
			//oldXRay.complexSymImag[prot]=xray.complexSymImag[prot];
			//oldXRay.complexSymTempReal[prot]=xray.complexSymTempReal[prot];
			//oldXRay.complexSymTempImag[prot]=xray.complexSymTempImag[prot];
			oldXRay.complexSymReal=xray.complexSymReal;
			oldXRay.complexSymImag=xray.complexSymImag;
			oldXRay.complexSymTempReal=xray.complexSymTempReal;
			oldXRay.complexSymTempImag=xray.complexSymTempImag;
		}
		else
		{
			degOfFreedomNCS=old_degOfFreedom;
			//xray.complexSymReal[prot]=oldXRay.complexSymReal[prot];
			//xray.complexSymImag[prot]=oldXRay.complexSymImag[prot];
			//xray.complexSymTempReal[prot]=oldXRay.complexSymTempReal[prot];
			//xray.complexSymTempImag[prot]=oldXRay.complexSymTempImag[prot];
			xray.complexSymReal=oldXRay.complexSymReal;
			xray.complexSymImag=oldXRay.complexSymImag;
			xray.complexSymTempReal=oldXRay.complexSymTempReal;
			xray.complexSymTempImag=oldXRay.complexSymTempImag;
		}
		nmove=move.score.size();
		if (nmove%updateMoveSizeInterval==0 && i!=nmove)
		{
			updateMoveSizes(move, nMoveTypes);
			updateMoveProbabilities(move, nMoveTypes);
		}
	}
	acceptRate=calcAcceptRate(move);
		<<" temp= "<<temp<<" acceptRate= "<<acceptRate<<endl;
	return score_old;
}

Real monteCarlo(Vector &degOfFreedom, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params, Real temp, int nsteps, MoveStruct &move)
{
	int prot=0, pick, nmove;
	int updateMoveSizeInterval=50;
	int nMoveTypes=degOfFreedom.size();
	int nprot=Proteins.size();
	Real score, score_old, score_initial=0, dScore, randNum;
	Real acceptProb, acceptRate, time;
	//Real scoreSlow;
	Vector old_degOfFreedom, subV, scores;
	XRayStruct oldXRay;
	//XRayStruct slowXRay;
	timeval start, end;

	//slowXRay=xray;
	if (params.OptimizeOccupancies)
	{
		//sumScatteringContinuousSegments(xray, degOfFreedom, nprot*6);
		sumScatteringContinuousSegments(xray, Proteins[0], degOfFreedom, nprot*6);
	}
	for (int prot=0;prot<nprot;prot++)
	{
		score_initial=calcFitnessChange(degOfFreedom, Proteins[prot], xray, expXRay, lattice, params, prot*6);
	}
	//params.UseFastTranslationRotation=false;
	//scoreSlow=calcFitness(degOfFreedom, Proteins, slowXRay, expXRay, lattice, params);
	//params.UseFastTranslationRotation=true;
/*
	int nsym, npoint;
	Get3DVectorSize(xray.complexSymTempReal, nprot, nsym, npoint);
	for (int n=0;n<nprot;n++)
	{
		for (int j=0;j<nsym;j++)
		{
			for (int k=0;k<npoint;k++)
			{
				//	<<" imag= "<<xray.complexSymTempImag[n][j][k]<<" imagslow= "<<slowXRay.complexSymTempImag[n][j][k]
					<<" imag= "<<xray.complexSymTempImag[n][j][k]
					<<" n= "<<n<<" j= "<<j<<" k= "<<k
					<<" miller= "<<xray.miller[k].Pos[X]
					<<" "<<xray.miller[k].Pos[Y]
					<<" "<<xray.miller[k].Pos[Z]<<endl;
			}
		}
	}
	for (int n=0;n<npoint;n++)
	{
	}
*/
	/*
	   params.UseFastTranslationRotation=false;
	   scoreSlow=calcFitness(degOfFreedom, ProteinsOriginal, xraySlow, expXRay, lattice, params);
	   int npoint=xray.miller.size();
	   for (int i=0;i<npoint;i++)
	   {
	   }
	 */
	score=score_initial;
	oldXRay.complexSymReal=xray.complexSymReal;
	oldXRay.complexSymImag=xray.complexSymImag;
	oldXRay.complexSymTempReal=xray.complexSymTempReal;
	oldXRay.complexSymTempImag=xray.complexSymTempImag;
	oldXRay.clashScore=xray.clashScore;
	score_old=score_initial;
	old_degOfFreedom=degOfFreedom;
	for (int i=0;i<nsteps;i++)
	{
		gettimeofday(&start, NULL);
		makeRandomChange(degOfFreedom, move, pick, Proteins[0].XBoxLength, Proteins[0].YBoxLength, Proteins[0].ZBoxLength, Proteins[0].spaceGroup, nprot);
		//printVector(degOfFreedom);
		prot=pick/6;
		if (prot>=nprot) prot=nprot-1;
		score=calcFitnessChange(degOfFreedom, Proteins[prot], xray, expXRay, lattice, params, pick);
		//params.UseFastTranslationRotation=false;
		//scoreSlow=calcFitness(degOfFreedom, Proteins, slowXRay, expXRay, lattice, params);
		//params.UseFastTranslationRotation=true;
		//printVector(degOfFreedom);
		//if (abs(score-scoreSlow)>1.0) endProgram(__LINE__, __FILE__);
		gettimeofday(&end, NULL);
		time=calcTimeDiff(start, end);
		dScore=score-score_old;
		acceptProb=calcAcceptProb(dScore, temp);
		SafePushBack(move.score, score, "score");
		SafePushBack(move.dScore, dScore, "dScore");
		SafePushBack(move.temp, temp, "temp");
		SafePushBack(move.acceptProb, acceptProb, "acceptProb");
		SafePushBack(move.time, time, "time");
		randNum=randDouble(0,1.0);
		if (acceptProb>randNum)
		{
			score_old=score;
			old_degOfFreedom=degOfFreedom;
			oldXRay.complexSymReal[prot]=xray.complexSymReal[prot];
			oldXRay.complexSymImag[prot]=xray.complexSymImag[prot];
			oldXRay.clashScore=xray.clashScore;
			if (params.NumCopies>1)
			{
				oldXRay.complexSymTempReal[prot]=xray.complexSymTempReal[prot];
				oldXRay.complexSymTempImag[prot]=xray.complexSymTempImag[prot];
			}
		}
		else
		{
			degOfFreedom=old_degOfFreedom;
			xray.complexSymReal[prot]=oldXRay.complexSymReal[prot];
			xray.complexSymImag[prot]=oldXRay.complexSymImag[prot];
			//xray.clashScore=oldXRay.clashScore;
			if (params.NumCopies>1)
			{
				xray.complexSymTempReal[prot]=oldXRay.complexSymTempReal[prot];
				xray.complexSymTempImag[prot]=oldXRay.complexSymTempImag[prot];
			}
			if (params.OptimizeOccupancies)
			{
				//sumScatteringContinuousSegments(xray, degOfFreedom, nprot*6);
				sumScatteringContinuousSegments(xray, Proteins[0], degOfFreedom, nprot*6);
			}
		}
		nmove=move.score.size();
		if (nmove%updateMoveSizeInterval==0 && i!=nmove)
		{
			updateMoveSizes(move, nMoveTypes);
			updateMoveProbabilities(move, nMoveTypes);
		}
	}
	acceptRate=calcAcceptRate(move);
	if (params.Verbose) 
	{
		<<" score_final= "<<score_old
		<<" temp= "<<temp<<" acceptRate= "<<acceptRate<<endl;
	}
	return score_old;
}



void calcScoresWithClash(Matrix &degOfFreedoms, Vector &score, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	int nrep=degOfFreedoms.size();
	for (int i=0;i<nrep;i++)
	{
		score[i]=calcRotationTranslationFitness2(degOfFreedoms[i], Proteins, xray, expXRay, lattice, params);
		score[i]+=calcClashScore(degOfFreedoms[i], Proteins, params);
	}
}

void printReplica(Vector &degOfFreedom, vector<ProteinStruct> Proteins, string outputFile)
{
	ProteinStruct Complex;
	Vector subVector;
	int ndof=degOfFreedom.size();
	for (int i=0;i<ndof;i++)
	{
	}
	setSegmentOccupancies(Proteins[0], degOfFreedom, 1);
	subVector=getSubVector(degOfFreedom, 0, 5);
	placeProteins(Proteins, subVector, Complex);
	printPdb(outputFile, Complex);
}

void printReplicas(Matrix &degOfFreedoms, vector<ProteinStruct> &Proteins, XRayParamStruct &params, string replicaStructuresOutput)
{
	//Outputs the coordinates of the replicas.
	string outputFile;
	string suffix;
	int Size=degOfFreedoms.size();

	if (params.ClashScoreUse=="DuringMR") suffix="Clash_";
	if (params.OptimizeOccupancies) suffix+="OptimizeOccupancies_";	
	for (int i=0;i<Size;i++)
	{
		fixDOF(Proteins[0].XBoxLength, Proteins[0].YBoxLength, Proteins[0].ZBoxLength, degOfFreedoms[i]);
		outputFile=replicaStructuresOutput+"_"+suffix+toStr(i)+".pdb";
		printReplica(degOfFreedoms[i], Proteins, outputFile);
	}
}

void roundCoordinates(vector<AtomStruct> &Atoms)
{
	//Rounds coordinates to nearest 0.001 A, to match pdb format, so better
	//comparisons can be made from structures read from pdb files.
	int natom=Atoms.size();

	for (int i=0;i<natom;i++)
	{
		Atoms[i].x=roundReal(Atoms[i].x, 0.001);
		Atoms[i].y=roundReal(Atoms[i].y, 0.001);
		Atoms[i].z=roundReal(Atoms[i].z, 0.001);
	}
}

void calcScores(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &scores, vector<string> &scoreNames)
{
	Real r, diffOfSqrs, pearson, product, logSqrDiff, rotationPossibleScore;
	Real densityScore, phaseScore, mltf, nativeDensityCorrelation, emll;
	Real estimatedRmsd;
	scores.clear();
	scoreNames.clear();

	//r=RFactor(xray.i, expXRay.i);
	r=RFactorFast(xray.i, expXRay.F);
	SafePushBack(scores, r, "scores");
	SafePushBack(scoreNames, string("RFactor"), "scoreNames");

	diffOfSqrs=calcRFactor(xray, expXRay, "DiffOfSquares");
	SafePushBack(scores, diffOfSqrs, "scores");
	SafePushBack(scoreNames, string("DiffOfSquares"), "scoreNames");

	pearson=calcRFactor(xray, expXRay, "Pearson");
	SafePushBack(scores, pearson, "scores");
	SafePushBack(scoreNames, string("Pearson"), "scoreNames");

	product=-calcRFactor(xray, expXRay, "Product");
	SafePushBack(scores, product, "scores");
	SafePushBack(scoreNames, string("Product"), "scoreNames");

	logSqrDiff=calcLogSqrDiff(xray.i, expXRay.i);
	SafePushBack(scores, logSqrDiff, "scores");
	SafePushBack(scoreNames, string("LogSqrDiff"), "scoreNames");

	if (params.CalcRotationPossibleScore)
	{
		rotationPossibleScore=calcRotationPossibleScore(xray, expXRay);
	}
	else rotationPossibleScore=0;
	SafePushBack(scores, rotationPossibleScore, "scores");
	SafePushBack(scoreNames, string("RotationPossibleScore"), "scoreNames");
	
	densityScore=0;
	SafePushBack(scores, densityScore, "scores");
	SafePushBack(scoreNames, string("DensityScore"), "scoreNames");

	phaseScore=0;
	SafePushBack(scores, phaseScore, "scores");
	SafePushBack(scoreNames, string("PhaseScore"), "scoreNames");

	mltf=0;
	SafePushBack(scores, mltf, "scores");
	SafePushBack(scoreNames, string("mltf"), "scoreNames");

	if (params.NativePdbFile!="")
	{
		//nativeDensityCorrelation=calcNativeDensityCorrelation(xray, lattice, Protein, params);
		nativeDensityCorrelation=0;
	}
	else nativeDensityCorrelation=0;
	SafePushBack(scores, nativeDensityCorrelation, "scores");
	SafePushBack(scoreNames, string("nativeDensityCorrelation"), "scoreNames");

	//emll=calcEmpiricalMaximumLikelihood(xray, expXRay, Protein, params);
	emll=0;
	SafePushBack(scores, emll, "scores");
	SafePushBack(scoreNames, string("emll"), "scoreNames");

	estimatedRmsd=0;
	SafePushBack(scores, estimatedRmsd, "scores");
	SafePushBack(scoreNames, string("estimatedRmsd"), "scoreNames");
}

void printScores(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	vector<string> scoreNames;
	int nscore;
	Vector scores;

	calcScores(xray, expXRay, params, scores, scoreNames);
	nscore=scores.size();

	for (int i=0;i<nscore;i++)
	{
	}	
}

void printScores(Vector &degOfFreedoms, vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	//Prints all of the components of XRayScore for a replica.
	int nprot=Proteins.size();
	Real total;
	Real clashScore, clashCC, clashCS, clashSS;
	Real r, diffOfSqrs, pearson, product;
	Real logSqrDiff, rotationPossibleScore;
	ProteinStruct Complex;
	vector<ProteinStruct> OriginalProteins;
	LatticeStruct lattice;

	OriginalProteins=Proteins;
	printClashScore(degOfFreedoms, Proteins[0], params);
	clashScore=calcClashScore(degOfFreedoms, Proteins[0], params, clashCC, clashCS, clashSS);
	//Safe3DAlloc(xray.complexSymReal, 0, 0, 0, "xray.complexSymReal");
	//Safe3DAlloc(xray.complexSymImag, 0, 0, 0, "xray.complexSymImag");
	for (int i=0;i<nprot;i++)
	{
		rotateAtoms(Proteins[i].Atoms, 0, 0, degOfFreedoms[i*6+Z]);
		rotateAtoms(Proteins[i].Atoms, 0, degOfFreedoms[i*6+Y], 0);
		rotateAtoms(Proteins[i].Atoms, degOfFreedoms[i*6+X], 0, 0);
	}
	//calcScatteringComplex(Proteins, xray, degOfFreedoms);
	placeProteins(OriginalProteins, degOfFreedoms, Complex, xray);
	roundCoordinates(Complex.Atoms);
	calcSymmetryScattering(xray, Complex);
	calcIntensity(xray);
	//params.UseFastTranslationRotation=false;
	//calcFitness(degOfFreedoms, OriginalProteins, xray, expXRay, lattice, params);
	//calcFitnessChange(degOfFreedoms, OriginalProteins[0], xray, expXRay, lattice, params, 0);
	if (params.CorrectionFile!="") correctIntensity(xray);
	calcAmplitude(xray);

	//r=calcRFactor(xray, expXRay, "RFactor");
	r=RFactor(xray.i, expXRay.i);
	diffOfSqrs=calcRFactor(xray, expXRay, "DiffOfSquares");
	pearson=calcRFactor(xray, expXRay, "Pearson");
	product=calcRFactor(xray, expXRay, "Product");
	logSqrDiff=calcLogSqrDiff(xray.i, expXRay.i);
	rotationPossibleScore=calcRotationPossibleScore(xray, expXRay);

	total=diffOfSqrs*params.DiffOfSquaresWeight;
	total+=r*params.RFactorWeight;
	total+=pearson*params.PearsonWeight;
	total+=product*params.ProductWeight;
	total+=rotationPossibleScore*params.RotationWeight;
	total+=clashScore;


}

void printScores(Matrix &degOfFreedoms, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	//Prints all of the components of XRayScore for all of the replicas.
	int nrep=degOfFreedoms.size();

	for (int i=0;i<nrep;i++)
	{
		//if (i==0) xray.boolVerbose=true;
		//else xray.boolVerbose=false;
		printScores(degOfFreedoms[i], Proteins, xray, expXRay, params);
	}
}

void printXRay(Vector &degOfFreedoms, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayParamStruct &params, string xrayOutFile)
{
	//Rotates and translates a protein, calculates the scattering, and 
	//outputs it to a file.
	//int nprot=Proteins.size();
	ProteinStruct Complex;

	//for (int i=0;i<nprot;i++)
	//{
	//	rotateAtoms(Proteins[i].Atoms, 0, 0, degOfFreedoms[i*6+Z]);
	//	rotateAtoms(Proteins[i].Atoms, 0, degOfFreedoms[i*6+Y], 0);
	//	rotateAtoms(Proteins[i].Atoms, degOfFreedoms[i*6+X], 0, 0);
	//}
	//calcScatteringComplex(Proteins, xray, degOfFreedoms);

	placeProteins(Proteins, degOfFreedoms, Complex, xray);
	roundCoordinates(Complex.Atoms);
	calcSymmetryScattering(xray, Complex);
	calcIntensity(xray);
	if (params.CorrectionFile!="") correctIntensity(xray);
	printXRay(xrayOutFile, xray, Proteins[0]);
}

void printXRay(Matrix degOfFreedoms, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayParamStruct &params)
{
	//Prints the calculated scattering of all of the replicas.
	int nrep=degOfFreedoms.size();
	string xrayOutFile;

	if (params.ReplicaXRayOutput!="")
	{
		for (int i=0;i<nrep;i++)
		{
			xrayOutFile=params.ReplicaXRayOutput+"_"+toStr(i)+".txt";
			printXRay(degOfFreedoms[i], Proteins, xray, params, xrayOutFile);
		}
	}
}

void optimizeAll(Matrix &degOfFreedoms, Vector &score, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	//Selects which method to use to perform rigid body optimization.
	int nmin, nrep;
	Real maxHours=0.02;
	Real dScoreCriteria=1.0e-15;
	//Real (*func_to_minimize)(Vector, MRStruct)=calcRotationTranslationFitness;
	Real (*func_to_minimize)(Vector, MRStruct)=calcFitness;
	MRStruct args;

	nrep=degOfFreedoms.size();
	nmin=nrep;
	//SafeAlloc(degOfFreedoms[0], 3, "degOfFreedoms");
	//setSpaceGroupP1(Proteins[0]);
	args.Proteins=Proteins;
	args.xray=xray;
	args.expXRay=expXRay;
	args.params=params;
	args.params.UseFastTranslationRotation=false;
	for (int i=0;i<nmin;i++)
	{
		//for (int j=0;j<3;j++)
		//{
		args.params.RotationWeight=0;
		minimizeTranslationOnly(Proteins, xray, expXRay, params, degOfFreedoms[i]);
		//args.params.RotationWeight=params.RotationWeight;
		//minimizeRotationFast(Proteins[0], xray, expXRay, params, degOfFreedoms[i]);
		//}
		if (params.MinimizeRotationTranslation=="Analytical")
		{
			minimizeRotationTranslationAnalytical(degOfFreedoms[i], Proteins, xray, expXRay, params, maxHours);	
		}
		else if (params.MinimizeRotationTranslation=="RotationOnly")
		{
			minimizeRotationOnly(Proteins[0], xray, expXRay, params, degOfFreedoms[i], maxHours);
		}
		else if (params.MinimizeRotationTranslation=="EstimatedGradient")
		{
			//optimizeUsingEstimatedGrad(degOfFreedoms[i], args, func_to_minimize, dScoreCriteria, maxHours);
		}
		else if (params.MinimizeRotationTranslation=="ConjugateGradient")
		{
			conjugateGradient(degOfFreedoms[i], args, func_to_minimize, dScoreCriteria, maxHours);
		}
		else if (params.MinimizeRotationTranslation=="OneAtATime")
		{
			minimizeOneAtATime(degOfFreedoms[i], args, func_to_minimize, maxHours);
		}
		else
		{
			string errorStr="Unrecognized option for MinimizeRotationTranslation "+params.MinimizeRotationTranslation;
			errorStr+=". Acceptable values are Analytical, RotationOnly, EstimatedGradient, ConjugateGradient, and OneAtATime";
		}
		//ProteinStruct Complex;
		//placeProteins(Proteins, degOfFreedoms[0], Complex);
		//calcScatteringBfactor(Proteins[0], xray);
		//endProgram(__LINE__, __FILE__);
	}
}

void initializeMove(int ndof, MoveStruct &move, ProteinStruct &Protein, int index, int nrep, XRayParamStruct &params)
{
	//Initializes the move sizes for the replica-exchange Monte Carlo 
	//simulation when non-crystallographic symmetry is not used.
	Real dAngle;
	Real minAngle=params.MinAngleChange;
	Real maxAngle=params.MaxAngleChange;
	Real dFrac;
	Real minFrac=params.MinTranslationFraction;
	Real maxFrac=params.MaxTranslationFraction;
	int nprot=ndof/6;

	dAngle=maxAngle-minAngle;
	dFrac=maxFrac-minFrac;
	SafeAlloc(move.width, ndof, "width");
	initializeMoveProbabilities(move, ndof);
	for (int i=0;i<nprot;i++)
	{
		for (int j=0;j<3;j++)
		{
			move.width[i*6+j]=minAngle+dAngle*Real(index)/Real(nrep);
		}
		move.width[i*6+3]=Protein.XBoxLength*(minFrac+dFrac*Real(index)/Real(nrep));
		move.width[i*6+4]=Protein.YBoxLength*(minFrac+dFrac*Real(index)/Real(nrep));
		move.width[i*6+5]=Protein.ZBoxLength*(minFrac+dFrac*Real(index)/Real(nrep));
	}

	for (int i=nprot*6;i<ndof;i++)
	{
		move.width[i]=0.1;
	}
}

void initializeMoveNCS(int ndof, MoveStruct &move, ProteinStruct &Protein, int index, int nrep, XRayParamStruct &params)
{
	//Sets the initial move sizes for when there are multiple proteins in
	//the asymmetric unit cell and there is non-crystallographic symmetry
	Real dAngle;
	Real minAngle=params.MinAngleChange;
	Real maxAngle=params.MaxAngleChange;
	Real dFrac;
	Real minFrac=params.MinTranslationFraction;
	Real maxFrac=params.MaxTranslationFraction;

	dAngle=maxAngle-minAngle;
	dFrac=maxFrac-minFrac;
	SafeAlloc(move.width, ndof, "width");
	initializeMoveProbabilities(move, ndof);
	for (int j=0;j<3;j++)
	{
		move.width[j]=minAngle+dAngle*Real(index)/Real(nrep);
	}
	move.width[3]=Protein.XBoxLength*(minFrac+dFrac*Real(index)/Real(nrep));
	move.width[4]=Protein.YBoxLength*(minFrac+dFrac*Real(index)/Real(nrep));
	move.width[5]=Protein.ZBoxLength*(minFrac+dFrac*Real(index)/Real(nrep));
	move.width[6]=minAngle+dAngle*Real(index)/Real(nrep);
	move.width[7]=minAngle+dAngle*Real(index)/Real(nrep);
	move.width[8]=Protein.XBoxLength*(minFrac+dFrac*Real(index)/Real(nrep));
	move.width[9]=Protein.YBoxLength*(minFrac+dFrac*Real(index)/Real(nrep));
}

void initializeMoves(int nrep, int ndof, vector<MoveStruct> &moves, ProteinStruct &Protein, XRayParamStruct &params)
{
	//Sets the initial move sizes of the replica-exchange Monte Carlo search
	MoveStruct tempMove;
	SafeAlloc(moves, tempMove, nrep, "moves");
	for (int i=0;i<nrep;i++)
	{
		if (params.UseNCS)
		{
			initializeMoveNCS(ndof, moves[i], Protein, i, nrep, params);	
		}
		else
		{
			initializeMove(ndof, moves[i], Protein, i, nrep, params);	
		}
	}
}

bool checkForGap(Vector &scores, Real stdDevBest)
{
	//Checks if there is a gap in the XRayScores of the replicas.
	int nrep=scores.size();
	int start, end;

	start=nrep/10;
	end=nrep/5;
	for (int i=start;i<end;i++)
	{
		if (scores[i]-scores[i-1]>stdDevBest*2.0) return true;
	}
	return false;
}

Real calcRMSD(vector<ProteinStruct> &Proteins1, Vector &degOfFreedom, vector<ProteinStruct> Proteins2)
{
	//Rotates and translates a protein according to the degrees of freedom
	//and then calculates the crystallographic RMSD.
	placeProteins(Proteins2, degOfFreedom);
	return calcRMSD(Proteins1, Proteins2, "NONE");
}

bool checkSimilarity(Matrix &degOfFreedom, vector<ProteinStruct> Proteins)
{
	//Checks if low temperature replicas are near to the lowest temperature
	//replica.  If they are the seach is terminated early.
	int nrep=degOfFreedom.size();
	int nbest=5;
	int nprot=Proteins.size();
	Real rmsd, cutoff=0.5;
	vector<ProteinStruct> copyProteins;

	for (int i=0;i<nprot;i++) removeNonCA(Proteins[i].Atoms);

	if (nbest>nrep) nbest=nrep;
	copyProteins=Proteins;

	placeProteins(copyProteins, degOfFreedom[0]);
	for (int i=1;i<nbest;i++)
	{
		rmsd=calcRMSD(copyProteins, degOfFreedom[i], Proteins);
		if (rmsd>cutoff) return false;
	}
	return true;
}

bool checkForConvergence(Vector &scores, Matrix &degOfFreedom, vector<ProteinStruct> &Proteins)
{
	//Checks if low temperature replicas are all the same.  If they are
	//the search is ended early.
	bool converged=false;
	int nrep=scores.size();
	int nbest=int(nrep*0.1);
	Real stdDev, stdDevBest;;
	Vector best;

	stdDev=calcStandardDeviation(scores);
	best=getSubVector(scores, 0, nbest);
	stdDevBest=calcStandardDeviation(best);
	//if (true)
	if (stdDevBest<stdDev*0.1)
	{
		converged=checkSimilarity(degOfFreedom, Proteins);
	}
	return converged;
}

void calcRMSDMatrix(Matrix &degOfFreedoms, ProteinStruct &Protein, Matrix &rmsd)
{
	//Calculates a crystallographic RMSD matrix using all replicas.
	int ndof, nrep, nprot, ntotal, rep, prot;
	int natom=Protein.Atoms.size();
	Vector subV;
	ProteinStruct tempProtein1, tempProtein2;
	vector<ProteinStruct> tempProteins;

	tempProtein1=Protein;
	tempProtein1.Atoms.clear();
	for (int i=0;i<natom;i++)
	{
		if (Protein.Atoms[i].AtomName=="CA")
		{
			SafePushBack(tempProtein1.Atoms, Protein.Atoms[i], "tempProtein");
			break;
		}		
	}
	Get2DVectorSize(degOfFreedoms, nrep, ndof, "degOfFreedoms");

	nprot=ndof/6;
	ntotal=nprot*nrep;
	Safe2DAlloc(rmsd, ntotal, ntotal, "rmsd");

	for (int i=0;i<ntotal;i++)
	{
		rep=i/ndof;
		prot=i%nprot;
		subV=getSubVector(degOfFreedoms[rep], prot*6, prot*6+5);
		placeProtein(tempProtein1, subV, tempProtein2);
		SafePushBack(tempProteins, tempProtein2, "tempProteins");
	}	

	for (int i=0;i<ntotal;i++)
	{
		for (int j=i+1;j<ntotal;j++)
		{
			rmsd[i][j]=calcRMSD_NoAlign(tempProteins[i], tempProteins[j]);
			rmsd[j][i]=rmsd[i][j];
		}
	}

}

void calcScatteringForAllModels(Vector &degOfFreedom, XRayStruct &xray, LatticeStruct &lattice, vector<ProteinStruct> &Proteins, XRayParamStruct &params, Matrix &real, Matrix &imag)
{
	int ndof=degOfFreedom.size();
	int nprot=ndof/6;
	Vector subV;

	for (int prot=0;prot<nprot;prot++)
	{
		subV=getSubVector(degOfFreedom, prot*6, prot*6+5);
		calcScattering(subV, Proteins[prot], xray, lattice, params);
		SafePushBack(real, xray.real, "real");
		SafePushBack(imag, xray.imag, "imag");
	}
}

void calcScatteringForAllModels(Matrix &degOfFreedoms, XRayStruct &xray, LatticeStruct &lattice, vector<ProteinStruct> &Proteins, XRayParamStruct &params, Matrix &real, Matrix &imag)
{
	int nrep=degOfFreedoms.size();

	for (int i=0;i<nrep;i++)
	{
		calcScatteringForAllModels(degOfFreedoms[i], xray, lattice, Proteins, params, real, imag);
	}
}

Real getScoreForPermutation(Matrix &real, Matrix &imag, XRayStruct &xray, XRayStruct &expXRay, vector<int> &permutation, Matrix &rmsds, XRayParamStruct &params)
{
	//This function is used when there are multiple proteins in the 
	//asymmetric unit.  Some replicas might find the correct solution for
	//one protein and another replica might find the correct solution for 
	//the other.  Thus, all combinations of proteins are tested.  This finds
	//the scores for a permutation efficiently.
	int nreal, npoint;
	int nmodel=permutation.size();
	Real rotationWeight=params.RotationWeight;
	Real score=0, overlapPenalty=1000.0, cutoff=4.0;

	ZeroVector(xray.real);
	ZeroVector(xray.imag);

	Get2DVectorSize(real, nreal, npoint, "real");
	for (int i=0;i<nmodel;i++)
	{
		if (permutation[i]==1)
		{
			for (int j=0;j<npoint;j++)
			{
				xray.real[j]+=real[i][j];
				xray.imag[j]+=imag[i][j];
			}
		}
	}
	for (int i=0;i<nmodel;i++)
	{
		if (permutation[i]==1)
		{
			for (int j=i+1;j<nmodel;j++)
			{
				if (permutation[j]==1)
				{
					if (rmsds[i][j]<cutoff)
					{
						score+=overlapPenalty;
					}
				}
			}
		}
	}
	calcIntensity(xray);
	//score+=RFactorFast(intensity, expXRay.F);
	params.RotationWeight=0;
	score+=calcRFactor(xray, expXRay, params);
	params.RotationWeight=rotationWeight;
	return score;
}

void initializePermutation(vector<int> &permutation, int nmodel, int nprot)
{
	SafeAlloc(permutation, nmodel, "permutation");
	for (int i=0;i<nprot;i++)
	{
		permutation[i]=1;
	}
}

void calcScoresForPermutations(Matrix &real, Matrix &imag, XRayStruct &xray, XRayStruct expXRay, int nprot, Vector &scores, vector< vector<int> > &permutations, Matrix &rmsds, int nuse, XRayParamStruct &params)
{
	//This function is used when there are multiple proteins in the 
	//asymmetric unit.  Some replicas might find the correct solution for
	//one protein and another replica might find the correct solution for 
	//the other.  Thus, all combinations of proteins are tested.
	int count=0;
	int max=100000000;
	int nreal, npoint;
	int reserveSize;
	vector<int> permutation;
	Real score, maxTime=600;
	Vector realTotal, imagTotal, intensity;
	timeval start, end;

	gettimeofday(&start, NULL);
	reserveSize=int(nChooseK(nuse, nprot));
	if (reserveSize>max) reserveSize=max;
	scores.reserve(reserveSize);
	permutations.reserve(reserveSize);	
	Get2DVectorSize(real, nreal, npoint, "real");
	SafeAlloc(realTotal, npoint, "realTotal");
	SafeAlloc(imagTotal, npoint, "imagTotal");
	initializePermutation(permutation, nuse, nprot);
	while (true)
	{
		if (isLastPermutation(permutation)) break;
		score=getScoreForPermutation(real, imag, xray, expXRay, permutation, rmsds, params);
		SafePushBack(scores, score, "scores");
		SafePushBack(permutations, permutation, "permutations");
		getNextPermutation(permutation);
		count++;
		if (count%1000000==0)
		{
			for (int i=0;i<nuse;i++)
			{
			}
			gettimeofday(&end, NULL);
			if (calcTimeDiff(start, end)>maxTime) break;
		}
		if (count>max) break;
	}
}

void getNewDegOfFreedoms(Vector &scores, vector< vector<int> > &permutations, Matrix &degOfFreedoms)
{
	//When there are multiple proteins in the asymmetric unit the proteins
	//are mixed and matched to produce lower scoring placements.  This 
	//creates the new degrees of freedom resulting from the mixing.
	int npermutation, nmodel;
	int nrep, nprot, prot;
	Real cutoff;
	Vector tempScores, temp;
	Matrix newDegOfFreedoms;

	Get2DVectorSize(permutations, npermutation, nmodel, "permutations");
	Get2DVectorSize(degOfFreedoms, nrep, nprot, "degOfFreedoms");
	nprot/=6;

	tempScores=scores;
	sort(tempScores);
	cutoff=tempScores[nrep-1];

	for (int i=0;i<npermutation;i++)
	{
		if (scores[i]<=cutoff)
		{
			temp.clear();
			for (int j=0;j<nmodel;j++)
			{
			}
			for (int j=0;j<nmodel;j++)
			{
				if (permutations[i][j]==1)
				{
					nrep=j/(nprot);
					prot=j%nprot;
					for (int k=0;k<6;k++)
					{
						SafePushBack(temp, degOfFreedoms[nrep][prot*6+k], "temp");
					}	
				}
			}
			SafePushBack(newDegOfFreedoms, temp, "newDegOfFreedoms");
		}
	}
	degOfFreedoms=newDegOfFreedoms;
}

void calcPermutationsOfProteins(Matrix &degOfFreedoms, XRayStruct xray, XRayStruct &expXRay, LatticeStruct &lattice, vector<ProteinStruct> &Proteins, XRayParamStruct &params)
{
	//When there are multiple proteins in the assymetric unit cell
	//this function calculates all possible combinations of proteins
	//from the replicas.
	int nprot=Proteins.size();
	int nrep=degOfFreedoms.size();
	int nuse;
	vector< vector<int> > permutations;
	Vector scores;
	Matrix real, imag, rmsds;

	if (params.NumCopies>3) nuse=nrep*nprot/2;
	else nuse=nrep*nprot;

	calcRMSDMatrix(degOfFreedoms, Proteins[0], rmsds);
	calcScatteringForAllModels(degOfFreedoms, xray, lattice, Proteins, params, real, imag);
	calcScoresForPermutations(real, imag, xray, expXRay, nprot, scores, permutations, rmsds, nuse, params);
	getNewDegOfFreedoms(scores, permutations, degOfFreedoms);
}

void shuffleReplicas(Matrix &degOfFreedoms)
{
	//Overwrites proteins in the high temperature replicas with proteins
	//in the low temperature replicas.  This can be used when there are
	//multiple proteins in the assymetric unit.  The idea is that if
	//there are two replicas which have part of the solution they can be
	//combined into one correct solution.
	vector<bool> used;
	int nrep, ndof, nprot;
	int pick, rep, prot;

	Get2DVectorSize(degOfFreedoms, nrep, ndof, "degOfFreedoms");
	nprot=ndof/6;

	SafeAlloc(used, false, nrep*nprot/2, "used");
	for (int i=nrep/2;i<nrep;i++)
	{
		for (int j=0;j<nprot;j++)
		{
			while (true)
			{
				pick=randInt(0, nrep*nprot/2);
				if (!used[pick])
				{
					used[pick]=true;
					rep=pick/nprot;
					prot=pick%nprot;
					for (int k=0;k<6;k++)
					{
						degOfFreedoms[i][j*6+k]=degOfFreedoms[rep][prot*6+k];
					}
					break;
				}
			}
		}
	}
}

void printContinuationFile(string continuationFile, Matrix &degOfFreedoms)
{
	//Prints out the final degrees of freedom so that the run can be 
	//continued later.
	ofstream file;
	int Size1, Size2;

	Get2DVectorSize(degOfFreedoms, Size1, Size2, "degOfFreedoms");
	OpenFileForWriting(continuationFile, file);

	for (int i=0;i<Size1;i++)
	{
		for (int j=0;j<Size2;j++)
		{
			file <<i<<" "<<j<<" "<<degOfFreedoms[i][j]<<endl;
		}
	}
	file.close();
}

void outputDegreesOfFreedom(Matrix &degOfFreedoms, int cycle, XRayParamStruct &params)
{
	//Prints the degrees of freedom to a trajectory file, because printing
	//actual structures takes up too much memory.
	ofstream file;
	int nreplica, ndof;

	Get2DVectorSize(degOfFreedoms, nreplica, ndof, "degOfFreedoms");
	if (params.PrintTrajectoryMax>nreplica)
	{
		params.PrintTrajectoryMax=nreplica;
	}
	OpenFileForWriting(params.DegOfFreedomTrajectoryFile, file);
	for (int i=0;i<params.PrintTrajectoryMax;i++)
	{
		for (int j=0;j<ndof;j++)
		{
			file <<degOfFreedoms[i][j]<<"\t";
		}
		file <<cycle<<"\t"<<i<<endl;
	}
	file.close();
}

void outputDegreesOfFreedomHeader(XRayParamStruct &params)
{
	//Prints a header for the degree of freedom file.
	ofstream file;
	OpenFileForWriting(params.DegOfFreedomTrajectoryFile, file);
	for (int i=0;i<params.NumCopies;i++)
	{
		file <<"RX\tRY\tRZ\tDX\tDY\tDZ\t";
	}
	if (params.UseNCS) file <<"NCS\tNCS\tNCS\tNCS\t";
	for (int i=0;i<params.NumOccupancySegments;i++)
	{
		file <<"O\t";
	}
	file <<"Cycle\tReplica\t"<<"NumCopies= "<<params.NumCopies<<"\tNumOccupancySegments= "<<params.NumOccupancySegments<<endl;
	file.close();
}

Real calcOrientation(ProteinStruct &Protein1, ProteinStruct &Protein2)
{
	//Calculates by how much the model must be rotate to match the target.
	vector<bool> use;
	int natom=Protein1.Atoms.size();
	int nsym=Protein1.symmetryOperations.size();
	Real theta, phi, psi, bestPsi;
	ProteinStruct tempProtein;

	
	SafeAlloc(use, true, natom, "use");
	tempProtein=Protein2;
	RMSD(Protein1.Atoms, tempProtein.Atoms, theta, phi, psi, use);
	psi=calcDeltaAngle(0, psi);
	bestPsi=psi;

	for (int i=1;i<nsym;i++)
	{
		tempProtein=Protein2;
		applyRotationMatrixToAtoms(tempProtein.Atoms, Protein1.symmetryOperations[i].rotationMat);
		RMSD(Protein1.Atoms, tempProtein.Atoms, theta, phi, psi, use);
		psi=calcDeltaAngle(0, psi);
		if (psi<bestPsi) bestPsi=psi;
		
	}
	return bestPsi;
}

void roundOccupancies(vector<AtomStruct> &Atoms)
{
	//Rounds occupancies for the purpose of making it equivalent to what
	//is read in from the pdb so a more accurate comparison can be made
	//between when occupancies are read in from a pdb file and when they
	//are generated by the program.
	int natom=Atoms.size();

	for (int i=0;i<natom;i++)
	{
		roundReal(Atoms[i].Occupancy, 0.01);
	}
}

void outputTrajectorySnapShot(Vector &degOfFreedom, vector<ProteinStruct> &Proteins, ProteinStruct &Native, XRayStruct &xray, XRayStruct &expXRay, string name, ofstream &file, XRayParamStruct &params)
{
	int nmiller=xray.miller.size();
	Real total;
	Real clashScore, clashCC, clashCS, clashSS;
	Real r, diffOfSqrs, pearson, product;
	Real xtalRmsd, pointXtalRmsd, orientation;
	Vector subVector;
	ProteinStruct Complex;
	vector<ProteinStruct> newProteins;
	LatticeStruct lattice;

	printClashScore(degOfFreedom, Proteins[0], params);
	clashScore=calcClashScore(degOfFreedom, Proteins[0], params, clashCC, clashCS, clashSS);
	printVector(degOfFreedom, "degOfFreedom");
	subVector=getSubVector(degOfFreedom, 0, 6*params.NumCopies-1);
	setSegmentOccupancies(Proteins[0], degOfFreedom, 1);
	string outputFile;
	outputFile="/nfs/amino-home/jouko/project/mr/run/Decoys_CB_2_hours_3Drobot_clash_2_hours_CalphaRadius_4.2/1DVO/pdb_out/"+name+".pdb";
	printReplica(degOfFreedom, Proteins, outputFile);
	placeProteins(subVector, Proteins[0], newProteins);
	calcSymmetryScattering(xray, newProteins[0]);
	calcIntensity(xray);
	if (params.CorrectionFile!="") correctIntensity(xray);
	calcAmplitude(xray);

	r=RFactor(xray.i, expXRay.i);
	diffOfSqrs=calcRFactor(xray, expXRay, "DiffOfSquares")*3000.0/Real(nmiller);
	pearson=calcRFactor(xray, expXRay, "Pearson");
	product=calcRFactor(xray, expXRay, "Product")*3000.0/Real(nmiller);

	total=diffOfSqrs*params.DiffOfSquaresWeight;
	total+=r*params.RFactorWeight;
	total+=pearson*params.PearsonWeight;
	total+=product*params.ProductWeight;
	total+=clashScore;

	file <<diffOfSqrs<<"\t"<<clashCC<<"\t"<<clashCS<<"\t"<<clashSS<<"\t"
	<<r<<"\t"<<pearson<<"\t"<<product<<"\t";
	if (params.NativePdbFile!="")
	{
		xtalRmsd=calcRMSD_NoAlign(Native, newProteins[0]);
		pointXtalRmsd=calcPointXtalRmsd(Native, newProteins[0]);
		orientation=calcOrientation(Native, newProteins[0]);
		file <<xtalRmsd<<"\t"<<pointXtalRmsd<<"\t"<<orientation<<"\t";
	}
	file <<name<<endl;
}

void outputTrajectorySnapShot(Matrix &degOfFreedoms, int cycle, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	ofstream file;
	string name;
	ProteinStruct Native, tempProtein;

	OpenFileForWriting(params.TrajectoryOutputFile, file);
	if (params.NativePdbFile!="") readPdb(params.NativePdbFile, Native);
	tempProtein=Proteins[0];
	calcTMScore(Native.Atoms, tempProtein.Atoms);
	if (params.PrintTrajectoryMax>params.NumReplicas)
	{
		params.PrintTrajectoryMax=params.NumReplicas;
	}
	for (int i=0;i<params.PrintTrajectoryMax;i++)
	{
		name="replica_"+toStr(i)+"_"+toStr(cycle);
		outputTrajectorySnapShot(degOfFreedoms[i], Proteins, tempProtein, xray, expXRay, name, file, params);
	}
}

void printCorrectConformation(vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	//If the native structure is given, calculates the XRayScore and other
	//things and outputs them to the score file, for the purpose of 
	//comparison.
	ofstream file;
	int nmiller=xray.miller.size();
	Real clashScore, clashCC, clashCS, clashSS;
	Real r, diffOfSqrs, pearson, product;
	Real xtalRmsd, pointXtalRmsd, orientation;
	ProteinStruct Native;

	OpenFileForWriting(params.TrajectoryOutputFile, file);
	readPdb(params.NativePdbFile, Native);
	Native.spaceGroup=Proteins[0].spaceGroup;
	Native.a=Proteins[0].a;
	Native.b=Proteins[0].b;
	Native.c=Proteins[0].c;
	Native.cartToFrac=Proteins[0].cartToFrac;
	Native.fracToCart=Proteins[0].fracToCart;
	Native.XBoxLength=Proteins[0].XBoxLength;
	Native.YBoxLength=Proteins[0].YBoxLength;
	Native.ZBoxLength=Proteins[0].ZBoxLength;
	Native.symmetryOperations=Proteins[0].symmetryOperations;
	calcTMScore(Native.Atoms, Proteins[0].Atoms);
	clashScore=calcClashScore(Proteins[0], params, clashCC, clashCS, clashSS);
	calcSymmetryScattering(xray, Proteins[0]);
	calcIntensity(xray);
	if (params.CorrectionFile!="") correctIntensity(xray);
	calcAmplitude(xray);

	r=RFactor(xray.i, expXRay.i);
	diffOfSqrs=calcRFactor(xray, expXRay, "DiffOfSquares")*3000.0/Real(nmiller);
	pearson=calcRFactor(xray, expXRay, "Pearson");
	product=calcRFactor(xray, expXRay, "Product")*3000.0/Real(nmiller);
	xtalRmsd=calcRMSD_NoAlign(Native, Proteins[0]);
	pointXtalRmsd=calcPointXtalRmsd(Native, Proteins[0]);
	int natom2=Proteins[0].Atoms.size();
	for (int i=0;i<natom2;i++)
	{
		printAtomInfo(Proteins[0].Atoms[i]);
	}
	orientation=calcOrientation(Native, Proteins[0]);


	file <<diffOfSqrs<<"\t"<<clashCC<<"\t"<<clashCS<<"\t"<<clashSS<<"\t"
	<<r<<"\t"<<pearson<<"\t"<<product<<"\t"<<xtalRmsd<<"\t"<<pointXtalRmsd
	<<"\t"<<orientation<<"\tcorrect"<<endl;
}

void printTrajectoryHeader(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	//Outputs a header for the trajectory file.
	ofstream file;

	OpenFileForWriting(params.TrajectoryOutputFile, file);
	file <<"XRayScore\tClashCC\tClashCS\tClashSS\tRFactor\tpearson\tproduct\t"
	<<"Xtal_rmsd\tDistance\torientation\tname"<<endl;
	file.close();
	if (params.NativePdbFile!="")
	{
		printCorrectConformation(Proteins, xray, expXRay, params);
	}
}

Real replicaExchange(Matrix &degOfFreedoms, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params, Real minTemp, Real maxTemp, int nrep, int nsteps, int ncycles)
{
	//Performs the replica-exchange Monte Carlo search. There is the 
	//possibility of optimizing the temperatures during REMC.
	bool converged;
	int best, Size, ndof;
	int minCyclesForTemperatureUpdate=50;
	int checkConvergenceInterval=50;
	vector<int> hotOrCold;	//Did the replica last visit the lowest or highest temperature replica?  See J. Chem. Phys. 134, 244111 (2011)
	Real maxHours=params.MaxReplicaExchangeTime;
	Real bestScore=0, finalScore;
	Real minScore, maxScore;
	//Real targetSwapRate=0.965;
	Vector degOfFreedom, score, temps, totalSwaps, swapAttempts;
	Vector hotHistogram, coldHistogram, fractionHot, fractionCold; // See J. Chem. Phys. 134, 244111 (2011)
	Matrix scores;
	vector<MoveStruct> moves;
	timeval start, end;

	Safe2DAlloc(scores, nrep, NUM_SCORE_TYPES, "scores");
	gettimeofday(&start, NULL);
	SafeAlloc(score, nrep, "score");
	SafeAlloc(totalSwaps, nrep, "totalSwaps");
	SafeAlloc(swapAttempts, nrep, "swapAttempts");
	setTemps(temps, minTemp, maxTemp, nrep);
	initializeHotAndCold(hotOrCold, hotHistogram, coldHistogram, nrep);
	degOfFreedom=degOfFreedoms[0];
	Size=degOfFreedom.size();
	//printScores(degOfFreedoms, Proteins, xray, expXRay, params);
	if (params.UseNCS) ndof=10;
	else if (params.OptimizeOccupancies) ndof=Proteins.size()*6+params.NumOccupancySegments;
	else ndof=Proteins.size()*6;
	//initializeMoves(degOfFreedoms, moves);
	initializeMoves(nrep, ndof, moves, Proteins[0], params);
	calcScores(degOfFreedoms, score, scores, Proteins, xray, expXRay, lattice, params);
	calcStandardDeviations(scores, xray.stdDev);
	calcScores(degOfFreedoms, score, scores, Proteins, xray, expXRay, lattice, params);
	MinMax(score, minScore, maxScore);
	if (params.SetTempsBasedOnScores)
	{
		minTemp=abs(minScore*params.MinTempFrac);
		maxTemp=abs(maxScore*params.MaxTempFrac);
	}
	setTemps(temps, minTemp, maxTemp, nrep);
	multipleSwaps(degOfFreedoms, score, totalSwaps, swapAttempts, temps, hotOrCold, hotHistogram, coldHistogram);
	ZeroVector(totalSwaps); ZeroVector(swapAttempts);
	bestScore=minScore;

	for (int i=0;i<ncycles;i++)
	{
		calcStandardDeviations(scores, xray.stdDev);
		for (int j=0;j<nrep;j++)
		{
			if (params.UseNCS)
			{
				score[j]=monteCarloNCS(degOfFreedoms[j], Proteins, xray, expXRay, lattice, params, temps[j], nsteps, moves[j]);
			}
			else
			{
				score[j]=monteCarlo(degOfFreedoms[j], Proteins, xray, expXRay, lattice, params, temps[j], nsteps, moves[j]);
			}
			if (score[j]<bestScore)
			{
				bestScore=score[j];
				degOfFreedom=degOfFreedoms[j];
			}
			scores[j]=xray.score;
		}
		if (params.PrintTrajectory && i%params.PrintTrajectoryInterval==0)
		{
			//outputTrajectorySnapShot(degOfFreedoms, i, Proteins, xray, expXRay, params);
			outputDegreesOfFreedom(degOfFreedoms, i, params);
		}
		if (params.EndIfConverged && i%checkConvergenceInterval==0)
		{
			converged=checkForConvergence(score, degOfFreedoms, Proteins);
			if (converged) 
			{
				break;
			}
		}
		if (params.NumCopies>1 && params.ShuffleDegreesOfFreedom && i%params.ShuffleInterval==0 && !params.UseNCS)
		{
			//shuffleReplicas(degOfFreedoms);
			calcPermutationsOfProteins(degOfFreedoms, xray, expXRay, lattice, Proteins, params);
			calcScores(degOfFreedoms, score, scores, Proteins, xray, expXRay, lattice, params);
		}
		multipleSwaps(degOfFreedoms, score, totalSwaps, swapAttempts, temps, hotOrCold, hotHistogram, coldHistogram);
		best=findMin(score, finalScore);
		for (int k=0;k<Size;k++)
		{
		}
		gettimeofday(&end, NULL);
			<<" maxHours= "<<maxHours<<endl;
		if (calcTimeDiff(start, end)/3600.0>maxHours) break;
		if (params.UseDynamicTemperature && i%minCyclesForTemperatureUpdate==0 && i>0)
		{
			//updateTemperatures(moves, temps, targetSwapRate, maxTemp);
			calcFractionHot(hotHistogram, coldHistogram, fractionHot);
			fixFractionHot(fractionHot);
			updateTemperatureUsingHotHistogram(fractionHot, temps);
		}
	}
	if (params.UseNCS)
	{
		calcDOFFromNCS(degOfFreedoms);
	}

	if (params.ClashScoreUse=="AtEnd")
	{
		calcScoresWithClash(degOfFreedoms, score, Proteins, xray, expXRay, lattice, params);
	}

	if (!params.UseNCS && params.RigidBodyOptimization) optimizeAll(degOfFreedoms, score, Proteins, xray, expXRay, lattice, params);
	degOfFreedoms[nrep-1]=degOfFreedom;
	multipleSwaps(degOfFreedoms, score, totalSwaps, swapAttempts, temps, hotOrCold, hotHistogram, coldHistogram);
	printSwapRates(totalSwaps, swapAttempts);
	printMoveSizes(moves);
	printMoveProbabilities(moves);
	if (!params.UseNCS) bestScore=monteCarlo(degOfFreedoms[0], Proteins, xray, expXRay, lattice, params, temps[0], nsteps*nrep, moves[0]);
	else bestScore=score[0];
	//SafePushBack(degOfFreedoms, degOfFreedom, "degOfFreedoms");
	printReplicas(degOfFreedoms, Proteins, params, params.ReplicaStructuresOutput);
	printScores(degOfFreedoms, Proteins, xray, expXRay, params);
	printXRay(degOfFreedoms, Proteins, xray, params);
	for (int i=0;i<Size;i++)
	{
	}
	//printClashScore(degOfFreedoms[0], Proteins[0], params);		
	if (params.ContinuationFile!="")
	{
		printContinuationFile(params.ContinuationFile, degOfFreedoms);
	}
	if (params.PrintTrajectory)
	{
		//outputTrajectorySnapShot(degOfFreedoms, i, Proteins, xray, expXRay, params);
		outputDegreesOfFreedom(degOfFreedoms, ncycles, params);
	}
	return bestScore;
}

void initializeDegOfFreedomNoNCS(Matrix &degOfFreedoms, ProteinStruct &Protein, XRayParamStruct &params)
{
	//Initializes the degrees of freedom for when non-crystallographic 
	//symmetry is not used.  First 3 degrees of freedom are rotational 
	//last three degrees of freedom are translational.  For some space
	//groups translation along 1 or more axis is arbitrary and the dof
	//are set to 0.
	string spaceGroup=Protein.spaceGroup;
	int nrep=params.NumReplicas;

	if (params.NumOccupancySegments!=0 && !params.OptimizeOccupancies)
	{
		params.NumOccupancySegments=0;
	}
	Safe2DAlloc(degOfFreedoms, nrep, params.NumCopies*6+params.NumOccupancySegments, "degOfFreedoms");
	for (int i=0;i<params.NumCopies;i++)
	{
		for (int j=0;j<nrep;j++)
		{
			degOfFreedoms[j][6*i]=randDouble(0, pi*0.5);
			degOfFreedoms[j][6*i+1]=randDouble(0, pi*0.5);
			degOfFreedoms[j][6*i+2]=randDouble(0, pi*0.5);
			degOfFreedoms[j][6*i+3]=randDouble(0, Protein.XBoxLength);
			degOfFreedoms[j][6*i+4]=randDouble(0, Protein.YBoxLength);
			degOfFreedoms[j][6*i+5]=randDouble(0, Protein.ZBoxLength);
			if (spaceGroup=="P3" || spaceGroup=="P31" || spaceGroup=="P32" || spaceGroup=="R3" || spaceGroup=="H3" || spaceGroup=="P4" || spaceGroup=="P41" || spaceGroup=="P42" || spaceGroup=="P43" || spaceGroup=="I4" || spaceGroup=="I41" || spaceGroup=="P6" || spaceGroup=="P61" || spaceGroup=="P65" || spaceGroup=="P62" || spaceGroup=="P64" || spaceGroup=="P63")
			{
				degOfFreedoms[j][Z+3]=0;
			}
			else if (spaceGroup=="P2" || spaceGroup=="P21" || spaceGroup=="C2" || spaceGroup=="A2" || spaceGroup=="I2" || spaceGroup=="P1211" || spaceGroup=="C121")
			{
				degOfFreedoms[j][Y+3]=0;
			}
			else if (spaceGroup=="P1")
			{
				degOfFreedoms[j][X+3]=0;
				degOfFreedoms[j][Y+3]=0;
				degOfFreedoms[j][Z+3]=0;
			}
		}
	}
	for (int i=0;i<nrep;i++)
	{
		for (int j=0;j<params.NumOccupancySegments;j++)
		{
			degOfFreedoms[i][params.NumCopies*6+j]=1.0;
		}
	}
	if (params.DegOfFreedomFile!="" && exist(params.DegOfFreedomFile))
	{
		initializeDegOfFreedoms(degOfFreedoms, params.DegOfFreedomFile);
		fixDegOfFreedoms(degOfFreedoms, Protein.XBoxLength, Protein.YBoxLength, Protein.ZBoxLength);
	}

}

void initializeDegOfFreedomNCS(Matrix &degOfFreedoms, ProteinStruct &Protein, XRayParamStruct &params)
{
	//Inititializes the degrees of freedom from REMC when 
	//non-crystallographic symmetry is used.  First 3 degrees of freedom
	//are rotational degrees of freedom.  Second 3 degrees of freedom are
	//translational degrees of freedom.  The last 4 degrees of freedom
	//define the axis of rotation relating the two proteins.  The first 2
	//of these define the orientation and the last two define the location.
	//For some space groups translation along one or more axis is arbitrary
	//so the degree of freedom is set to 0.
	string spaceGroup=Protein.spaceGroup;
	int nrep=params.NumReplicas;

	Safe2DAlloc(degOfFreedoms, nrep, 10, "degOfFreedoms");
	for (int j=0;j<nrep;j++)
	{
		degOfFreedoms[j][0]=randDouble(0, pi*0.5);
		degOfFreedoms[j][1]=randDouble(0, pi*0.5);
		degOfFreedoms[j][2]=randDouble(0, pi*0.5);
		degOfFreedoms[j][3]=randDouble(0, Protein.XBoxLength);
		degOfFreedoms[j][4]=randDouble(0, Protein.YBoxLength);
		degOfFreedoms[j][5]=randDouble(0, Protein.ZBoxLength);
		degOfFreedoms[j][6]=randDouble(0, pi);
		degOfFreedoms[j][7]=randDouble(0, pi);
		degOfFreedoms[j][8]=randDouble(0, Protein.XBoxLength);
		degOfFreedoms[j][9]=randDouble(0, Protein.YBoxLength);

		if (spaceGroup=="P3" || spaceGroup=="P31" || spaceGroup=="P32" || spaceGroup=="R3" || spaceGroup=="H3" || spaceGroup=="P4" || spaceGroup=="P41" || spaceGroup=="P42" || spaceGroup=="P43" || spaceGroup=="I4" || spaceGroup=="I41" || spaceGroup=="P6" || spaceGroup=="P61" || spaceGroup=="P65" || spaceGroup=="P62" || spaceGroup=="P64" || spaceGroup=="P63")
		{
			degOfFreedoms[j][Z+3]=0;
		}
		else if (spaceGroup=="P2" || spaceGroup=="P21" || spaceGroup=="C2" || spaceGroup=="A2" || spaceGroup=="I2" || spaceGroup=="P1211" || spaceGroup=="C121")
		{
			degOfFreedoms[j][Y+3]=0;
		}
		else if (spaceGroup=="P1")
		{
			degOfFreedoms[j][X+3]=0;
			degOfFreedoms[j][Y+3]=0;
			degOfFreedoms[j][Z+3]=0;
		}
	}
}

void initializeDegOfFreedom(Matrix &degOfFreedoms, ProteinStruct &Protein, XRayParamStruct &params)
{
	//Decides which type of degrees of freedom you have. With no ncs
	//there are 6n degrees of freedom. 3 rotation, 3 translation.
	//With ncs there are 2 proteins. For first protein there are 3
	//rotation and 3 translation degrees of freedom.  For the second there
	//are 4, dealing with the axis of rotation.
	//NCS means non-crystallographic symmetry.
	if (params.UseNCS)
	{
		initializeDegOfFreedomNCS(degOfFreedoms, Protein, params);
	}
	else
	{
		initializeDegOfFreedomNoNCS(degOfFreedoms, Protein, params);
	}
}

void addOccupanciesToDegOfFreedom(Matrix &degOfFreedoms, int NumOccupancySegments)
{
	//When occupanicies are not optimized there are 6n degrees of freedom,
	//but when occupancies are optimized there are 6n+nsegments dof.
	int nrep, ndof;

	Get2DVectorSize(degOfFreedoms, nrep, ndof, "degOfFreedoms");

	for (int i=0;i<nrep;i++)
	{
		for (int j=0;j<NumOccupancySegments;j++)
		{
			SafePushBack(degOfFreedoms[i], 1.0, "degOfFreedoms");
		}
	}
}

void turnOnOptimizeOccupancies(Matrix &degOfFreedoms, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayParamStruct &params, XRayParamStruct &copyParams)
{
	//Used when occupancies are not optimized in the first phase, but not
	//in the second.
	copyParams.NumOccupancySegments=10;
	copyParams.OptimizeOccupancies=true;
	calcScatteringContinuousSegments(Proteins[0], xray, copyParams);	
	copyParams.MaxReplicaExchangeTime=params.MaxReplicaExchangeOptimizeOccupancyTime;
	addOccupanciesToDegOfFreedom(degOfFreedoms, copyParams.NumOccupancySegments);
}

Real MRReplicaExchange(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	//Performs replica-exchange without clash score, followed by 
	//replica-exchange with clash score, followed by optional gradient
	//based minimization.
	string spaceGroup=Protein.spaceGroup;
	int nrep, nstep, ncycles;
	Real minTemp, maxTemp;
	Vector degOfFreedom;
	Matrix degOfFreedoms;
	MRStruct args;
	Real dScoreCriteria=1.0e-5;
	Real maxHours=1.0;
	//Real (*func_to_minimize)(Vector, MRStruct&)=calcRotationTranslationFitness2;
	Real (*func_to_minimize2)(Vector, MRStruct)=calcRotationTranslationFitness;
	XRayParamStruct copyParams;
	ProteinStruct Complex;
	vector<ProteinStruct> Proteins, Proteins2;

	SafeAlloc(Proteins, Protein, params.NumCopies, "Proteins");
	//center(Protein.Atoms);
	args.Proteins=Proteins;
	args.xray=xray;
	args.expXRay=expXRay;
	args.params=params;
	nrep=params.NumReplicas;
	nstep=params.StepsPerCycle;
	ncycles=params.NumCycles;
	minTemp=params.MinTemp;
	maxTemp=params.MaxTemp;
	initializeDegOfFreedom(degOfFreedoms, Protein, params);
	//if (params.PrintTrajectory) printTrajectoryHeader(Proteins, xray, expXRay, params);
	if (params.PrintTrajectory) outputDegreesOfFreedomHeader(params);
	replicaExchange(degOfFreedoms, Proteins, xray, expXRay, lattice, params, minTemp, maxTemp, nrep, nstep, ncycles);
	copyParams=params;
	copyParams.ClashScoreUse="DuringMR";
	copyParams.MaxReplicaExchangeTime=params.MaxReplicaExchangeClashTime;
	replicaExchange(degOfFreedoms, Proteins, xray, expXRay, lattice, copyParams, minTemp, maxTemp, nrep, nstep, ncycles);
	//turnOnOptimizeOccupancies(degOfFreedoms, Proteins, xray, params, copyParams);
	//replicaExchange(degOfFreedoms, Proteins, xray, expXRay, lattice, copyParams, minTemp, maxTemp, nrep, nstep, ncycles);
	degOfFreedom=degOfFreedoms[0];
	//nonlinearConjugateGradient(degOfFreedom, args, func_to_minimize2, dScoreCriteria);
	if (params.RigidBodyOptimization)
	{
		conjugateGradient(degOfFreedom, args, func_to_minimize2, dScoreCriteria, maxHours);
	}
	//Proteins=args.Proteins;
	//placeProteins(Proteins, degOfFreedom, Complex, xray);
	//return calcScattering(Proteins, xray, expXRay, params);
	return 0;
}

Real optimizeBFactorScale(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	//Attempts to optimize a scale factor applied to the B-factors of all of 
	//the atoms.
	Real dScoreCriteria=0.0000001;
	Real maxHours=1.0;
	Vector bFactorScale;
	MRStruct args;
	Real (*func_to_minimize)(Vector, MRStruct)=calcBFactorFitness;
	vector<ProteinStruct> Proteins;

	args.Protein=Protein;
	args.xray=xray;
	args.expXRay=expXRay;
	args.params=params;

	SafeAlloc(bFactorScale, Real(1.0), 1, "bFactorScale");
	optimizeCoefficients(bFactorScale, args, func_to_minimize, dScoreCriteria, maxHours);
	scaleBFactors(Protein, bFactorScale[0]);
	SafePushBack(Proteins, Protein, "Proteins");
	return calcScattering(Proteins, xray, expXRay, params);

}

Real optimizeFormFactorScale(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	//Attempts to optimize the structure factors of atom types, with a scale
	//factor for the width of the gaussians used to represent them.
	Real dScoreCriteria=0.0000001;
	Real maxHours=1.0;
	Vector formFactorScale;
	MRStruct args;
	Real (*func_to_minimize)(Vector, MRStruct)=calcFormFactorFitness;
	vector<ProteinStruct> Proteins;

	args.Protein=Protein;
	args.xray=xray;
	args.expXRay=expXRay;
	args.params=params;

	SafeAlloc(formFactorScale, Real(1.9142), 1, "formFactorScale");
	optimizeCoefficients(formFactorScale, args, func_to_minimize, dScoreCriteria, maxHours);
	solventCorrectedScatteringFactor(args.xray.f, args.xray.q, args.params);
	SafePushBack(Proteins, Protein, "Proteins");
	return calcScattering(Proteins, xray, expXRay, params);

}

void molecularReplacement(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params)
{
	//There are a few options here most of which are not really mr.
	//Some of them have not been tested in a long time and may be broken
	//The main option is MRReplicaExchange which performs mr using replica 
	//exchange monte carlo simulations.  MRRotateTranslate is the more 
	//traditional approach which performs a systematic rotational search
	//followed by systematic translational searches. optimizeBFactorScale 
	//is not mr.  It attempts to find a scale factor that when applied to
	//B-factors optimizes the fit to experimental data.  
	//optimizeFormFactorScale is similar except that it optimizes the 
	//scattering factors of atom types. phase(...) assumes that the model is
	//correctly placed, calculates the phases from the model, calculates the
	//electron density from the phases and experimental data, and iteratively
	//attempts to improve the electron density.  
	vector<ProteinStruct> Proteins;

	SafePushBack(Proteins, Protein, "Proteins");
	if (params.RotateTranslate) MRRotateTranslate(Protein, xray, expXRay, lattice, params);
	if (params.ReplicaExchange) MRReplicaExchange(Protein, xray, expXRay, lattice, params);
	if (params.OptimizeBFactorScale) optimizeBFactorScale(Protein, xray, expXRay, params);
	if (params.OptimizeFormFactorScale) optimizeFormFactorScale(Protein, xray, expXRay, params);
	if (params.Phase) phase(Proteins, xray, expXRay, lattice, params);

}

void molecularReplacement(string inputPdbFile, XRayParamStruct &params)
{
	LatticeStruct lattice;
	ProteinStruct Protein, UnitCell;
	XRayStruct xray, expXRay;
	parsePdb(inputPdbFile, Protein, params.RemoveWaters, params.RemoveUnknownAtoms);
	calcUnitCellVectors(Protein, params);
	findCoreAtoms(Protein.Atoms, params.CalphaRadius, params.AsaCutoff);
	if (params.SpaceGroup!="") setSpaceGroup(Protein, params.SpaceGroup);
	readExperimentalDataFile(params.XRayInput, expXRay);
	//removeNegativeMiller(expXRay);
	removeNegativeIntensity(expXRay);
	//combineMillerIndexes(expXRay);
	initialize(lattice, xray, expXRay, params);
	removeSystematicAbsences(xray, Protein);
	removeNonoverlapping(xray, expXRay);
	initializeScatteringFactors(xray, Protein, params);
	normalizeXRay(expXRay.i);
	if (params.UseFastTranslationRotation)
	{
		if (params.OptimizeOccupancies)
		{
			calcScatteringContinuousSegments(Protein, xray, params);	
		}
		else calcPlaceProteinAndScatteringContinuous(Protein, xray, xray.realContinuous, xray.imagContinuous, params);
	}
	millerToQ(expXRay, lattice);
	if (params.TranslationRotationMethod=="Patterson") 
	{
		expXRay.patterson.real=xray.patterson.real;
		expXRay.patterson.imag=xray.patterson.imag;
		expXRay.patterson.u=xray.patterson.u;
		calcPatterson(expXRay);
	}
	Protein.intSpaceGroup=spaceGroupStrToSpaceGroupInt(Protein.spaceGroup);
	molecularReplacement(Protein, xray, expXRay, lattice, params);
	//makeUnitCell(Protein, UnitCell);
	printPdb(params.DensityOutput, Protein);
	//printDensity(params.DensityOutput, lattice);
	if (params.PrintOutput) printXRay(params.XRayOutput, xray, Protein);
}

