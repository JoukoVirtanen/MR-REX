# include <vector>
# include <iostream>
# include <algorithm>
# include <sys/time.h>

using namespace std;

# include "../../LibraryFiles/Structures.h"
# include "../../LibraryFiles/MathUtils.h"
# include "../../LibraryFiles/ReadPdb.h"
# include "../../LibraryFiles/AtomUtils.h"
# include "../../LibraryFiles/time.h"
# include "../../LibraryFiles/ConvertStructures.h"
# include "../../LibraryFiles/Cube.h"
# include "../../LibraryFiles/PrintPdb.h"
# include "../../LibraryFiles/normalize.h"
# include "../../LibraryFiles/rmsd.h"
# include "../../LibraryFiles/TMScore.h"

# include "XRayStruct.h"
# include "Params.h"
# include "XRay.h"
# include "LatticeStruct.h"
# include "molecularReplacement.h"
# include "CalcXRay.h"
# include "ReadExperimentalDataFile.h"
# include "ReadParameters.h"
# include "Initialize.h"
# include "SpaceGroup.h"
# include "MaximumLikelihood.h"
# include "CalcClashScore.h"
# include "FlexibleEditing.h"
# include "EmpiricalMaximumLikelihood.h"


void printXRayTerminal(XRayStruct &xray)
{
	int nmiller=xray.miller.size();

	for (int i=0;i<nmiller;i++)
	{
		cout <<xray.miller[i].Pos[X]<<"\t"<<xray.miller[i].Pos[Y]<<"\t"<<xray.miller[i].Pos[Z]<<"\t"<<xray.i[i]<<endl;
	}
}

void printXRayTerminal(XRayStruct &xray, XRayStruct &expXRay)
{
	int nmiller=xray.miller.size();

	for (int i=0;i<nmiller;i++)
	{
		cout <<xray.miller[i].Pos[X]<<"\t"<<xray.miller[i].Pos[Y]<<"\t"<<xray.miller[i].Pos[Z]<<"\t"<<xray.i[i]<<"\t"<<expXRay.miller[i].Pos[X]<<"\t"<<expXRay.miller[i].Pos[Y]<<"\t"<<expXRay.miller[i].Pos[Z]<<endl;
	}
}

void setVdwByAtomType(vector<AtomStruct> &Atoms)
{
	int natom=Atoms.size();

	for (int i=0;i<natom;i++)
	{
		if (Atoms[i].atomid==0) Atoms[i].vdw=0.5;
		else Atoms[i].vdw=1.5;
	}
}

void calcScattering(string pdbFile, ProteinStruct &Protein, vector<CubeStruct> &Cubes, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &scores, Real estimatedRmsd)
{
	string name, outputPdbFile;
	int npoint=xray.miller.size();
	Real phaseScore, densityScore;
	Real r, diffOfSqrs, pearson, product, mltf;
	Real rotationPossibleScore=UNK_REAL, logSqrDiff;
	Real nativeDensityCorrelation=0;
	Real clashCC, clashCS, clashSS;
	Real estimatedRmsdProtein=0;
	Real emll;
	Real aveOccupancy=1.0;
	Vector translation;
	Matrix rotationMatrix;
	XRayStruct xrayDensity, xrayPhase;
	LatticeStruct lattice, expLattice;
	timeval start, end;

	cout <<"In calcScattering"<<endl;

	SafeAlloc(scores, MLTF+1, "scores");
	SafeAlloc(translation, 3, "translation");
	Safe2DAlloc(rotationMatrix, 3, 3, "rotationMatrix");
	rotationMatrix[0][0]=1.0;
	rotationMatrix[1][1]=1.0;
	rotationMatrix[2][2]=1.0;

	initializeLattice(lattice, params);
	expLattice=lattice;

	scaleBFactors(Protein, params.BFactorScale);
	setVdwByAtomType(Protein.Atoms);
	gettimeofday(&start, NULL);
	findCoreAtoms(Protein.Atoms, params.CalphaRadius, params.AsaCutoff);
	gettimeofday(&end, NULL);
	//cout <<"findCoreAtoms took "<<calcTimeDiff(start, end)<<endl;

	//cout <<"Before calcClashScore"<<endl;
	//setSurfaceAtomsClashingWithCoreAtomsToZero(Protein.Atoms);
	//flexibleScattering(Protein, xray, expXRay, params);
	//setIncorrectAtomsToZero(Protein, params);
	if (params.ProbabilitiesFile!="") 
	{
		normalizeXRay(expXRay.i);
		cout <<"estimatedMinRmsd= "<<estimatedRmsd<<endl;
		estimatedRmsdProtein=estimateRmsd(Protein, xray, expXRay, params);
		if (estimatedRmsd<0 || estimatedRmsd>8) estimatedRmsd=7.5;
		setOccupancies(Protein, xray, expXRay, params, estimatedRmsd);
		name=getBase(pdbFile);
		outputPdbFile=params.OccupanciesPdbDir+"/"+name;
		//printPdb(outputPdbFile, Protein);
	}
	if (params.CubeFile!="")
	{
		readPdb(params.CubeFile, Cubes);
		subtractDensity(Cubes, params.BulkDensity);
		calcScatteringCubeProtein(Protein, Cubes, xray);
	}
	else if (params.DcdFilePaths!="") 
	{
		aveOverTrajectory(params.DcdFilePaths, Protein, xray, params);
	}
	else if (params.UseFastTranslationRotation)
	{
		if (params.FastTranslationRotationInterpolation=="Linear")
		{
			cout <<"Before calcScatteringContinuous"<<endl;
			calcScatteringContinuous(Protein, xray, xray.realContinuous, xray.imagContinuous, params);
		}
/*
		else if (params.FastTranslationRotationInterpolation=="Hessian")
		{
			cout <<"Before calcScatteringContinuousGrad"<<endl;
			calcScatteringContinuousGrad(Protein, xray, params);
		}
*/
		else
		{
			string errorStr="Unrecognized ";
			errorStr+="FastTranslationRotationInterpolation ";
			errorStr+=params.FastTranslationRotationInterpolation;
			errorStr+="\nacceptable values are Linear and Hessian";
			error(errorStr, __LINE__, __FILE__);
		}
		calcFastTranslationRotation(Protein, xray, Protein.rotationMatrix, 0, params);
		xray.complexSymTempReal=xray.complexSymReal;
		xray.complexSymTempImag=xray.complexSymImag;
		translateXRay2(xray, Protein, Protein.translation, 0);
		calcTotalFormFactor(xray);
		calcIntensity(xray);
	}
	else 
	{
		cout <<"Before calcSymmetryScattering"<<endl;
		cout <<"ExcludedVolumeMethod= "<<params.ExcludedVolumeMethod<<endl;
		if (params.ExcludedVolumeMethod=="Cube") 
		{
			cout <<"Before calcExcludedScattering"<<endl;
			calcExcludedScattering(Protein, xray, params); 
		}
		calcSymmetryScattering(xray, Protein);
		cout <<"After calcSymmetryScattering"<<endl;
		calcIntensity(xray);
	}
	aveOccupancy=calcAverageOccupancy(Protein.Atoms);
	if (params.CorrectionFile!="") correctIntensity(xray);
	if (params.OccupancyCorrectionFile!="" && params.OptimizeOccupancies) correctOccupancyIntensity(xray, aveOccupancy);
	fixBFactor(xray, expXRay);
	//printXRayTerminal(xray);
	//printXRayTerminal(xray, expXRay);
	calcAmplitude(xray);
	if (params.XRayInput!="")
	{
		normalizeXRay(expXRay.i);
		if (params.OptimizeSegmentOccupanciesAfterMR)
		{
			if (params.SegmentsFromFile)
			{
				optimizeSegmentsFromFileAfterMR(xray, expXRay, Protein, params);
			}
			else
			{
				optimizeSegmentOccupanciesAfterMR(xray, expXRay, Protein, params);
			}
			aveOccupancy=calcAverageOccupancy(Protein.Atoms);
		}
		//readExperimentalDataFile(params.XRayInput, expXRay);
		//cout <<"After readCifFile"<<endl;
		//millerToQ(expXRay, Protein.XBoxLength, Protein.YBoxLength, Protein.ZBoxLength);
		millerToQ(expXRay, params);
		//cout <<"expXRay.q.size()= "<<expXRay.q.size()<<endl;
		//compareAverageIntensities(xray, expXRay);
		//cout <<"After compareAverageIntensities"<<endl;
		if (params.BulkSolventCorrection)
		{
			cout <<"Before bulkSolventCorrection"<<endl;
			bulkSolventCorrection(xray, expXRay, params);
			compareAverageIntensities(xray, expXRay);
		}
		//cout <<"Before RFactor"<<endl;
		gettimeofday(&start, NULL);
		//r=calcRFactor(xray, expXRay, "RFactor");
		r=RFactor(xray.i, expXRay.i);
		//r=RFactorFast(xray.i, expXRay.F);
		gettimeofday(&end, NULL);
		//cout <<"calcRFactor took "<<calcTimeDiff(start, end)<<endl;

		gettimeofday(&start, NULL);
		//diffOfSqrs=calcRFactor(xray, expXRay, "DiffOfSquares")*3000.0/Real(npoint);
		diffOfSqrs=calcRFactor(xray, expXRay, "DiffOfSquares");
		gettimeofday(&end, NULL);
		//cout <<"diffOfSqrs took "<<calcTimeDiff(start, end)<<endl;

		gettimeofday(&start, NULL);
		pearson=calcRFactor(xray, expXRay, "Pearson");
		gettimeofday(&end, NULL);
		//cout <<"pearson took "<<calcTimeDiff(start, end)<<endl;

		gettimeofday(&start, NULL);
		//product=-calcRFactor(xray, expXRay, "Product")*3000.0/Real(npoint);
		product=-calcRFactor(xray, expXRay, "Product");
		gettimeofday(&end, NULL);
		//cout <<"product took "<<calcTimeDiff(start, end)<<endl;


		gettimeofday(&start, NULL);
		//logSqrDiff=calcLogSqrDiff(xray.i, expXRay.i)/Real(npoint);
		logSqrDiff=calcLogSqrDiff(xray.i, expXRay.i);
		gettimeofday(&end, NULL);
		//cout <<"calcLogSqrDiff took "<<calcTimeDiff(start, end)<<endl;

		gettimeofday(&start, NULL);
		if (params.CalcRotationPossibleScore)
		{
			rotationPossibleScore=calcRotationPossibleScore(xray, expXRay);
		}
		gettimeofday(&end, NULL);
		//cout <<"rotationPossibleScore took "<<calcTimeDiff(start, end)<<endl;

		gettimeofday(&start, NULL);
		//xrayDensity=xray;
		//calcExperimentalDensity(xrayDensity, expXRay, expLattice);
		//densityScore=calcProteinDensityMatch(Protein, xrayDensity, lattice, expLattice);
		densityScore=0;

		if (params.NativePdbFile!="")
		{
			cout <<"Before calcNativeDensityCorrelation"<<endl;
			//nativeDensityCorrelation=calcNativeDensityCorrelation(xray, lattice, Protein, params);
			nativeDensityCorrelation=0;
		}
		cout <<"nativeDensityCorrelation= "<<nativeDensityCorrelation<<endl;

		//densityScore=0;
		//string densityExpOutput="/home/joukov/project/mr/pdb_out/Decoys_CB_8_hours_3Drobot_nrep_300_4_1DVO/1DVO/1DVO_decoy0_1_mr/ExpDensity.pdb";
		//string densityOutput="/home/joukov/project/mr/pdb_out/Decoys_CB_8_hours_3Drobot_nrep_300_4_1DVO/1DVO/1DVO_decoy0_1_mr/Density.pdb";
		//printDensity(densityExpOutput, expLattice);
		//printDensity(densityOutput, lattice);
		//endProgram(__LINE__, __FILE__);
		cout <<"densityScore= "<<densityScore<<endl;
		//densityScore=0;
		gettimeofday(&end, NULL);
		//cout <<"calcLogSqrDiff took "<<calcTimeDiff(start, end)<<endl;

		gettimeofday(&start, NULL);
		cout <<"Before phase"<<endl;
		//xrayPhase=xray;
		//phaseScore=phase(Protein, xrayPhase, expXRay, lattice, params);
		phaseScore=0;
		gettimeofday(&end, NULL);
		//cout <<"phase took "<<calcTimeDiff(start, end)<<endl;


		gettimeofday(&start, NULL);
		cout <<"expXRay.Ferror.size()= "<<expXRay.Ferror.size()<<endl;
		//mltf=calcMLTF(xray, expXRay, 0);
		mltf=0;
		gettimeofday(&end, NULL);

		//emll=calcEmpiricalMaximumLikelihood(xray, expXRay, Protein, params);
		emll=0;

		cout <<"RFactor= "<<r<<endl;
		cout <<"DiffOfSquares= "<<diffOfSqrs<<endl;
		cout <<"Pearson= "<<pearson<<endl;
		cout <<"Product= "<<product<<endl;
		cout <<"RotationPossibleScore= "<<rotationPossibleScore<<endl;
		cout <<"LogSqrDiff= "<<logSqrDiff<<endl;
		cout <<"DensityScore= "<<densityScore<<endl;
		cout <<"PhaseScore= "<<phaseScore<<endl;
		cout <<"mltf= "<<mltf<<endl;
		cout <<"estimatedRmsd= "<<estimatedRmsdProtein<<endl;
		cout <<"emll= "<<emll<<endl;
		cout <<"random= "<<npoint*log(0.51)<<endl;	
		cout <<"aveOccupancy= "<<aveOccupancy<<endl;

		scores[RFACTOR]=r;
		scores[SQR_DIFF]=diffOfSqrs;
		scores[PEARSON]=pearson;
		scores[PRODUCT]=product;
		scores[ROTATION]=rotationPossibleScore;
		scores[LOG_SQR_DIFF]=logSqrDiff;
		scores[PHASE]=phaseScore;
		scores[DENSITY_CC]=densityScore;
		scores[MLTF]=mltf;
	}
	gettimeofday(&start, NULL);
	calcClashScoreRecord(Protein, params, clashCC, clashCS, clashSS);
	calcClashScore(Protein, params, clashCC, clashCS, clashSS);
	gettimeofday(&end, NULL);
	cout <<"printClashScore took "<<calcTimeDiff(start, end)<<endl;

	if (params.PrintOutput) printXRay(params.XRayOutput, xray, Protein);
	cout <<"After printXRay"<<endl;
	scores[CLASH_CC]=clashCC;
	scores[CLASH_CS]=clashCS;
	scores[CLASH_SS]=clashSS;
	cout <<"clashCC= "<<clashCC<<endl;
	cout <<"clashCS= "<<clashCS<<endl;
	cout <<"clashSS= "<<clashSS<<endl<<endl;
}

void calcMultipleScattering(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	ofstream file;
	int nline;
	vector<string> lines, logs;
	Real estimatedRmsd=UNK_REAL;
	Vector scores;
	ProteinStruct Protein;
	vector<CubeStruct> Cubes;

	cout <<"In calcMultipleScattering"<<endl;

	readLines(params.PdbPathsFile, lines, "PdbPathsFile");
	readLines(params.LogFiles, logs, "LogFiles");
	nline=lines.size();

	if (params.ProbabilitiesFile!="") 
	{
		normalizeXRay(expXRay.i);
		estimatedRmsd=findMinEstimatedRmsd(params.PdbPathsFile, xray, expXRay, params);
	}

	for (int i=0;i<nline;i++)
	{
		//if (i==1) xray.boolVerbose=true;
		//else xray.boolVerbose=false;
		parsePdb(lines[i], Protein, params.RemoveWaters, params.RemoveUnknownAtoms);
		//setOccupancies(Protein.Atoms, 1.0);
		cout <<"natom= "<<Protein.Atoms.size()<<endl;
		//cout <<"Before calcUnitCellVectors"<<endl;
		calcUnitCellVectors(Protein, params);
		//cout <<"Before scaleBFactors"<<endl;
		//scaleBFactors(Protein, params.BFactorScale);	
		//cout <<"Before setSpaceGroup"<<endl;
		if (params.SpaceGroup!="") setSpaceGroup(Protein, params.SpaceGroup);
		//cout <<"Before calcScattering"<<endl;
		//if (i!=0)
		//{
			calcScattering(lines[i], Protein, Cubes, xray, expXRay, params, scores, estimatedRmsd);
		//}
		//else
		//{
		//	optimizeSegmentOccupanciesAfterMR(xray, expXRay, Protein, params);
		//	break;
		//}
		//cout <<"After calcScattering"<<endl;
/*
 
		OpenFileForWriting(logs[i], file);
		file <<"clashCC= "<<scores[CLASH_CC]<<endl;
		file <<"clashCS= "<<scores[CLASH_CS]<<endl;
		file <<"clashSS= "<<scores[CLASH_SS]<<endl;
		file <<"RFactor= "<<scores[RFACTOR]<<endl;
		file <<"DiffOfSquares= "<<scores[SQR_DIFF]<<endl;
		file <<"Pearson= "<<scores[PEARSON]<<endl;
		file <<"Product= "<<scores[PRODUCT]<<endl;
		file <<"RotationPossibleScore= "<<scores[ROTATION]<<endl;
		file.close();
*/
	}
} 

void calcMultipleScatteringFromDOF(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	int nprot;
	Real estimatedRmsd=UNK_REAL;
	Vector scores;
	vector<CubeStruct> Cubes;
	vector< vector<ProteinStruct> > Proteins;

	dofToProteins(params.CenteredAndRotatedOutputPdb, params.DegOfFreedomTrajectoryFile, Proteins);
	nprot=Proteins.size();
	for (int i=0;i<nprot;i++)
	{
		assignAtomIDs(Proteins[i][0].Atoms);
		calcUnitCellVectors(Proteins[i][0], params);
		//scaleBFactors(Proteins[i][0], params.BFactorScale);	
		if (params.SpaceGroup!="") setSpaceGroup(Proteins[i][0], params.SpaceGroup);
		calcScattering("none", Proteins[i][0], Cubes, xray, expXRay, params, scores, estimatedRmsd);
	}
}

void calcScattering(string inputPdbFile, XRayParamStruct &params)
{
	Real estimatedRmsd=UNK_REAL;
	Vector scores;
	LatticeStruct lattice;
	ProteinStruct Protein;
	vector<CubeStruct> Cubes;
	XRayStruct xray, expXRay;
	timeval start, end;


	parsePdb(inputPdbFile, Protein, params.RemoveWaters, params.RemoveUnknownAtoms);

	cout <<"Before calcUnitCellVectors"<<endl;
	calcUnitCellVectors(Protein, params);
	cout <<"Before setSpaceGroup"<<endl;
	if (params.SpaceGroup!="") setSpaceGroup(Protein, params.SpaceGroup);
	gettimeofday(&start, NULL);
	if (params.XRayInput!="") 
	{
		readExperimentalDataFile(params.XRayInput, expXRay);
		//removeNegativeMiller(expXRay);
		removeNegativeIntensity(expXRay);
	}
	gettimeofday(&end, NULL);
	cout <<"readExperimentalDataFile took "<<calcTimeDiff(start, end)<<endl;

	params.PdbFile=inputPdbFile;
	
	gettimeofday(&start, NULL);
	initialize(lattice, xray, expXRay, params);
	gettimeofday(&end, NULL);
	cout <<"initialize took "<<calcTimeDiff(start, end)<<endl;

	//scaleBFactors(Protein, params.BFactorScale);	
	gettimeofday(&start, NULL);
	initializeScatteringFactors(xray, Protein, params);
	gettimeofday(&end, NULL);
	cout <<"initializeScatteringFactors took "<<calcTimeDiff(start, end)<<endl;
	//removeSystematicAbsences(xray, Protein);
	removeNonoverlapping(xray, expXRay);
	initializeCentric(xray, Protein);
	initializeCentric(expXRay, Protein);

	//xray.boolVerbose=false;
	cout <<"DegOfFreedomTrajectoryFile= "<<params.DegOfFreedomTrajectoryFile<<endl;
	if (params.PdbPathsFile!="") 
	{
		calcMultipleScattering(xray, expXRay, params);
	}
	else if (params.DegOfFreedomTrajectoryFile!="")
	{
		cout <<"Before calcMultipleScatteringFromDOF"<<endl;
		calcMultipleScatteringFromDOF(xray, expXRay, params);	
	}
	else
	{
		calcScattering("none", Protein, Cubes, xray, expXRay, params, scores, estimatedRmsd);
	}
}

int main(int argc, char *argv[])
{
	string inputPdbFile, paramsFile;
	XRayParamStruct params;
	timeval seconds1, seconds2;

	cout <<"STARTING"<<endl;
	if (argc!=3)
	{
		cout <<"Usage: "<<argv[0]<<" inputPdbFile paramsFile"<<endl;
		exit(EXIT_FAILURE);
	}

	inputPdbFile=argv[1];
	paramsFile=argv[2];

	cout <<"Before SetDefaultParams"<<endl;
	setDefaultParams(params);
	cout <<"params.UseExpLookUp= "<<params.UseExpLookUp<<endl;
	cout <<"Before ReadParameterFile"<<endl;
	ReadParameterFile(paramsFile, params);
	cout <<"params.UseExpLookUp= "<<params.UseExpLookUp<<endl;

	gettimeofday(&seconds1, NULL);
	cout <<"About to enter calcScattering"<<endl;
	calcScattering(inputPdbFile, params);
	gettimeofday(&seconds2, NULL);
	cout <<"Run time: "<<calcTimeDiff(seconds1, seconds2)<<endl;
	cout <<"DONE"<<endl;

	return EXIT_SUCCESS;
}
