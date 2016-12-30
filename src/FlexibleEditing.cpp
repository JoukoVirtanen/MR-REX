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
# include "../../LibraryFiles/PrintPdb.h"
# include "../../LibraryFiles/normalize.h"
# include "../../LibraryFiles/rmsd.h"
# include "../../LibraryFiles/TMScore.h"
# include "../../LibraryFiles/MinMax.h"

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

const int CC=0;
const int CS=1;
const int SS=2;

void readSegmentsFile(string segmentsFile, vector<int> &segmentStart, vector<int> &segmentEnd)
{
	vector<string> lines, str;
	int start, end, nstr;
	int nline, maxSegments=13;

	readLines(segmentsFile, lines, "SegmentsFile");
	nline=lines.size();

	if (nline>maxSegments) nline=maxSegments;

	for (int i=0;i<nline;i++)
	{
		Tokenize2(lines[i], "-", str);
		nstr=str.size();
		if (nstr==2)
		{
			start=toInt(str[0]);
			end=toInt(str[1]);
			SafePushBack(segmentStart, start, "segmentStart");
			SafePushBack(segmentEnd, end, "segmentEnd");
		}
	}
}

void intToBinary(int integer, Vector &binary)
{
	int Size=binary.size();
	int power;

	for (int i=Size-1;i>=0;i--)
	{
		power=int(calcPower(2.0, i));
		if (integer>=power)
		{
			binary[i]=1.0;
			integer-=power;
		}
		else binary[i]=0;
	}
}

void residueSegmentToAtomSegment(vector<AtomStruct> &Atoms, int &segmentStart, int &segmentEnd)
{
	bool first=false;
	int natom=Atoms.size();
	int tempSegmentEnd=0;

	for (int i=0;i<natom;i++)
	{
		if (Atoms[i].ResidueNum==segmentStart && !first)
		{
			segmentStart=i;
			first=true;
		}
		if (Atoms[i].ResidueNum==segmentEnd)
		{
			tempSegmentEnd=i;
		}
	}
	segmentEnd=tempSegmentEnd;
}

void residueSegmentsToAtomSegments(vector<AtomStruct> &Atoms, vector<int> &segmentStart, vector<int> &segmentEnd)
{
	int nsegment=segmentStart.size();

	for (int i=0;i<nsegment;i++)
	{
		residueSegmentToAtomSegment(Atoms, segmentStart[i], segmentEnd[i]);
	}
}

void setSegmentOccupancies(ProteinStruct &Protein, Vector &segmentOccupancies, vector<int> &segmentStart, vector<int> &segmentEnd)
{
	int nsegment=segmentStart.size();

	setOccupancies(Protein.Atoms, 0.0);
	for (int i=0;i<nsegment;i++)
	{
		for (int j=segmentStart[i];j<=segmentEnd[i];j++)
		{
			Protein.Atoms[j].Occupancy=segmentOccupancies[i];
		}
	}
}

void setSegmentOccupancies(ProteinStruct &Protein, Vector &segmentOccupancies, vector<int> &segmentStart, vector<int> &segmentEnd, Real defaultOccupancy)
{
	int nsegment=segmentStart.size();

	setOccupancies(Protein.Atoms, defaultOccupancy);
	for (int i=0;i<nsegment;i++)
	{
		for (int j=segmentStart[i];j<=segmentEnd[i];j++)
		{
			Protein.Atoms[j].Occupancy=segmentOccupancies[i];
		}
	}
/*
	int natom=Protein.Atoms.size();
	for (int i=0;i<natom;i++)
	{
		cout <<"i= "<<i<<" occupancy= "<<Protein.Atoms[i].Occupancy<<endl;
	}
*/
}

void setSegmentOccupancies(ProteinStruct &Protein, Vector &degOfFreedom, XRayParamStruct &params)
{
	int nsegment;
	int ndof=degOfFreedom.size();
	vector<int> segmentStart, segmentEnd;
	Vector segmentOccupancies;

	readSegmentsFile(params.SegmentsFile, segmentStart, segmentEnd);	
	residueSegmentsToAtomSegments(Protein.Atoms, segmentStart, segmentEnd);
	nsegment=segmentStart.size();
	segmentOccupancies=getSubVector(degOfFreedom, ndof-nsegment-1, ndof-2);
	setSegmentOccupancies(Protein, segmentOccupancies, segmentStart, segmentEnd, degOfFreedom[ndof-1]);
}

void setSegmentOccupancies(ProteinStruct &Protein, Vector &degOfFreedom)
{
	int nsegment;
	int ndof=degOfFreedom.size();
	Vector segmentOccupancies;

	nsegment=Protein.segmentStart.size();
	segmentOccupancies=getSubVector(degOfFreedom, ndof-nsegment, ndof-1);
	setSegmentOccupancies(Protein, segmentOccupancies, Protein.segmentStart, Protein.segmentEnd, 1.0);
}

void setSegmentOccupancyToOne(vector<AtomStruct> &Atoms, int segmentStart, int segmentEnd)
{
	int natom=Atoms.size();

	for (int i=0;i<natom;i++)
	{
		if (i>=segmentStart && i<=segmentEnd)
		{
			Atoms[i].Occupancy=1.0;
		}
		else Atoms[i].Occupancy=0;
	}
}	

void setFixedAndSegmentOccupancyToOne(vector<AtomStruct> &Atoms, vector<int> &segmentStart, vector<int> &segmentEnd, int segment)
{
	int nsegment=segmentStart.size();

	setOccupancies(Atoms, 1.0);
	for (int i=0;i<nsegment;i++)
	{
		for (int j=segmentStart[i];j<=segmentEnd[i];j++)
		{
			if (j!=segment) Atoms[j].Occupancy=0;
		}
	}
}

void setNonFixedOccupancyToZero(ProteinStruct &Protein, vector<int> &segmentStart, vector<int> &segmentEnd)
{
	int nsegment=segmentStart.size();

	setOccupancies(Protein.Atoms, 1.0);
	for (int i=0;i<nsegment;i++)
	{
		for (int j=segmentStart[i];j<=segmentEnd[i];j++)
		{
			Protein.Atoms[j].Occupancy=0;
		}
	}
}

void calcSegmentScattering(XRayStruct &xray, ProteinStruct &Protein, XRayParamStruct &params, int segment, Vector &real, Vector &imag)
{
	Vector segmentOccupancies;

	SafeAlloc(segmentOccupancies, params.NumOccupancySegments, "segmentOccupancies");
	segmentOccupancies[segment]=1.0;
	setSegmentOccupancies(Protein, segmentOccupancies, 0);
	calcSymmetryScattering(xray, Protein);
	real=xray.real;
	imag=xray.imag;
}

void calcSegmentScattering(XRayStruct &xray, ProteinStruct &Protein, XRayParamStruct &params, int segmentStart, int segmentEnd, Vector &real, Vector &imag)
{
	setSegmentOccupancyToOne(Protein.Atoms, segmentStart, segmentEnd);
	calcSymmetryScattering(xray, Protein);
	real=xray.real;
	imag=xray.imag;
}

void addSegmentScattering(Vector &segmentOccupancies, Matrix &real, Matrix &imag, XRayStruct &xray)
{
	int nsegment=segmentOccupancies.size();
	int npoint=xray.real.size();

	zeroVector(xray.real);
	zeroVector(xray.imag);
	for (int i=0;i<nsegment;i++)
	{
		for (int j=0;j<npoint;j++)
		{
			xray.real[j]+=real[i][j]*segmentOccupancies[i];
			xray.imag[j]+=imag[i][j]*segmentOccupancies[i];
		}
	}
}

void addSegmentScattering(Vector &segmentOccupancies, Matrix &real, Matrix &imag, Vector &fixedReal, Vector &fixedImag, XRayStruct &xray)
{
	int nsegment=segmentOccupancies.size();
	int npoint=xray.real.size();

	xray.real=fixedReal;
	xray.imag=fixedImag;
	for (int i=0;i<nsegment;i++)
	{
		for (int j=0;j<npoint;j++)
		{
			xray.real[j]+=real[i][j]*segmentOccupancies[i];
			xray.imag[j]+=imag[i][j]*segmentOccupancies[i];
		}
	}
}

void calcSegmentClashScores(ProteinStruct &Protein, XRayParamStruct &params, Matrix &intraSegmentClashes, Array3D &interSegmentClashes)
{
	int nsegment=params.NumOccupancySegments;
	Real clashCC, clashCS, clashSS;
	Vector segmentOccupancies;
	vector<ProteinStruct> Proteins;

	SafeAlloc(segmentOccupancies, nsegment, "segmentOccupancies");
	Safe2DAlloc(intraSegmentClashes, nsegment, 3, "intraSegmentClashes");
	Safe3DAlloc(interSegmentClashes, nsegment, nsegment, 3, "interSegmentClashes");
	SafePushBack(Proteins, Protein, "Proteins");

	for (int i=0;i<nsegment;i++)
	{
		zeroVector(segmentOccupancies);
		segmentOccupancies[i]=1.0;
		setSegmentOccupancies(Proteins[0], segmentOccupancies, 0);
		calcClashScore(Proteins[0], params, clashCC, clashCS, clashSS);
		intraSegmentClashes[i][CC]=clashCC;
		intraSegmentClashes[i][CS]=clashCS;
		intraSegmentClashes[i][SS]=clashSS;
	}

	for (int i=0;i<nsegment;i++)
	{
		for (int j=i+1;j<nsegment;j++)
		{
			zeroVector(segmentOccupancies);
			segmentOccupancies[i]=1.0;
			segmentOccupancies[j]=1.0;
			setSegmentOccupancies(Proteins[0], segmentOccupancies, 0);
			calcClashScore(Proteins[0], params, clashCC, clashCS, clashSS);
			interSegmentClashes[i][j][CC]=clashCC-intraSegmentClashes[i][CC]-intraSegmentClashes[j][CC];
			interSegmentClashes[i][j][CS]=clashCS-intraSegmentClashes[i][CS]-intraSegmentClashes[j][CS];
			interSegmentClashes[i][j][SS]=clashSS-intraSegmentClashes[i][SS]-intraSegmentClashes[j][SS];
		}
	}
}

void calcSegmentClashScores(ProteinStruct &Protein, XRayParamStruct &params, vector<int> &segmentStart, vector<int> &segmentEnd, Real &fixedClash, Vector &fixedNonFixedClashes, Vector &intraSegmentClashes, Matrix &interSegmentClashes)
{
	int nsegment=segmentStart.size();
	Real clashCC, clashCS, clashSS;
	Vector segmentOccupancies;

	SafeAlloc(fixedNonFixedClashes, nsegment, "fixedNonFixedClashes");
	SafeAlloc(intraSegmentClashes, nsegment, "intraSegmentClashes");
	SafeAlloc(segmentOccupancies, nsegment, "segmentOccupancies");
	Safe2DAlloc(interSegmentClashes, nsegment, nsegment, "interSegmentClashes");

	for (int i=0;i<nsegment;i++)
	{
		setSegmentOccupancyToOne(Protein.Atoms, segmentStart[i], segmentEnd[i]);
		intraSegmentClashes[i]=calcClashScore(Protein, params, clashCC, clashCS, clashSS);
		//cout <<"i= "<<i<<" intraSegmentClashes= "<<intraSegmentClashes[i]<<endl;
	}

	for (int i=0;i<nsegment;i++)
	{
		setFixedAndSegmentOccupancyToOne(Protein.Atoms, segmentStart, segmentEnd, i);
		fixedNonFixedClashes[i]=calcClashScore(Protein, params, clashCC, clashCS, clashSS);
		fixedNonFixedClashes[i]-=fixedClash;
		//cout <<"i= "<<i<<" fixedNonFixedClashes= "<<fixedNonFixedClashes[i]<<endl;
	}

	for (int i=0;i<nsegment;i++)
	{
		for (int j=i+1;j<nsegment;j++)
		{
			zeroVector(segmentOccupancies);
			segmentOccupancies[i]=1.0;
			segmentOccupancies[j]=1.0;
			setSegmentOccupancies(Protein, segmentOccupancies, segmentStart, segmentEnd);
			interSegmentClashes[i][j]=calcClashScore(Protein, params, clashCC, clashCS, clashSS);
			interSegmentClashes[i][j]-=intraSegmentClashes[i];
			interSegmentClashes[i][j]-=intraSegmentClashes[j];
			//cout <<"i= "<<i<<" j= "<<j<<" interSegmentClashes= "<<interSegmentClashes[i][j]<<endl;
		}
	}
}

Real calcClashScore(Vector &segmentOccupancies, Matrix &intraSegmentClashes, Array3D &interSegmentClashes, XRayParamStruct &params, Real &clashCC, Real &clashCS, Real &clashSS)
{
	int nsegment=segmentOccupancies.size();
	Real clashScore=0;

	clashCC=clashCS=clashSS=0;
	for (int i=0;i<nsegment;i++)
	{
		clashCC+=intraSegmentClashes[i][CC]*segmentOccupancies[i];
		clashCS+=intraSegmentClashes[i][CS]*segmentOccupancies[i];
		clashSS+=intraSegmentClashes[i][SS]*segmentOccupancies[i];
		for (int j=i+1;j<nsegment;j++)
		{
			clashCC+=interSegmentClashes[i][j][CC]*segmentOccupancies[i]*segmentOccupancies[j];
			clashCS+=interSegmentClashes[i][j][CS]*segmentOccupancies[i]*segmentOccupancies[j];
			clashSS+=interSegmentClashes[i][j][SS]*segmentOccupancies[i]*segmentOccupancies[j];
		}
	}
	clashScore=clashCC*params.ClashWeightCoreCore;
	clashScore+=clashCS*params.ClashWeightCoreSurface;
	clashScore+=clashSS*params.ClashWeightSurfaceSurface;
	return clashScore;
}

Real calcClashScore(Vector &segmentOccupancies, Real fixedClash, Vector &fixedNonFixedClashes, Vector &intraSegmentClashes, Matrix &interSegmentClashes)
{
	int nsegment=segmentOccupancies.size();
	int Size1, Size2;
	Real clashScore=0;

	Get2DVectorSize(interSegmentClashes, Size1, Size2, "interSegmentClashes");
	//cout <<"Size1= "<<Size1<<" Size2= "<<Size2<<endl;
	//cout <<"fixedNonFixedClashes.size()= "<<fixedNonFixedClashes.size()<<endl;
	//cout <<"intraSegmentClashes.size()= "<<intraSegmentClashes.size()<<endl;
	clashScore=fixedClash;
	for (int i=0;i<nsegment;i++)
	{
		clashScore+=fixedNonFixedClashes[i]*segmentOccupancies[i];
		clashScore+=intraSegmentClashes[i]*segmentOccupancies[i];
		for (int j=i+1;j<nsegment;j++)
		{
			clashScore+=interSegmentClashes[i][j]*segmentOccupancies[i]*segmentOccupancies[j];
		}
	}
	return clashScore;
}

Real calcAveOccupancy(Vector &segmentOccupancies)
{
	return calcAverage(segmentOccupancies);
}

Real optimizeSegmentOccupanciesAfterMR(XRayStruct &xray, XRayStruct &expXRay, ProteinStruct &Protein, XRayParamStruct &params)
{
	//string XRayOutput;
	int ncombinations;
	Real clashScore;
	Real clashCC, clashCS, clashSS;
	Real score=UNK_REAL, bestScore=UNK_REAL;
	//Real aveOccupancy;
	Vector segmentOccupancies;
	Vector bestOccupancies;
	Vector tempReal, tempImag, zero;
	Matrix intraSegmentClashes;
	Matrix real, imag;
	Array3D interSegmentClashes;
	vector<ProteinStruct> Proteins;
	LatticeStruct lattice;

	SafePushBack(Proteins, Protein, "Proteins");
	SafeAlloc(segmentOccupancies, params.NumOccupancySegments, "segmentOccupancies");
	ncombinations=int(calcPower(2.0, params.NumOccupancySegments));
	SafeAlloc(zero, 6*params.NumCopies, "zero");

	for (int i=0;i<params.NumOccupancySegments;i++)
	{
		calcSegmentScattering(xray, Protein, params, i, tempReal, tempImag);
		SafePushBack(real, tempReal, "real");
		SafePushBack(imag, tempImag, "imag");
	}

	calcSegmentClashScores(Protein, params, intraSegmentClashes, interSegmentClashes);
	for (int i=1;i<ncombinations;i++)
	{
		intToBinary(i, segmentOccupancies);
		//cout <<endl;
		//cout <<"i= "<<i<<endl;
		//printVector(segmentOccupancies);
		addSegmentScattering(segmentOccupancies, real, imag, xray);
		calcIntensity(xray);
		if (params.CorrectionFile!="") correctIntensity(xray);
		//score=calcRFactor(xray, expXRay, params);
		//cout <<"score= "<<score<<endl;
		score=0;
		clashScore=calcClashScore(segmentOccupancies, intraSegmentClashes, interSegmentClashes, params, clashCC, clashCS, clashSS);
		//cout <<"clashScore= "<<clashScore<<endl;
		score+=clashScore;
		//aveOccupancy=calcAveOccupancy(segmentOccupancies);
		//cout <<"aveOccupancy= "<<aveOccupancy<<endl;
		//printScores(xray, expXRay, params);
		//cout <<"clashCC= "<<clashCC<<endl;		
		//cout <<"clashCS= "<<clashCS<<endl;		
		//cout <<"clashSS= "<<clashSS<<endl;
		//XRayOutput=params.XRayOutput+"/"+toStr(i)+"_"+toStr(aveOccupancy)+".txt";	
		//printXRay(XRayOutput, xray, Protein);
		if (score<bestScore || i==1)
		{
			bestScore=score;
			bestOccupancies=segmentOccupancies;
		}
		//cout <<"score= "<<score<<" bestScore= "<<bestScore<<endl;
		//cout <<endl;
	}
	setSegmentOccupancies(Protein, bestOccupancies, 0);
	addSegmentScattering(bestOccupancies, real, imag, xray);
	calcIntensity(xray);
	if (params.CorrectionFile!="") correctIntensity(xray);
	//cout <<"bestOccupancies:"<<endl;
	//printVector(bestOccupancies);
	//gettimeofday(&end, NULL);
	//cout <<"Took "<<calcTimeDiff(start, end)<<endl;
	return bestScore;
}


Real optimizeSegmentOccupanciesAfterMR(XRayStruct &xray, XRayStruct &expXRay, ProteinStruct &Protein, XRayParamStruct &params, vector<int> &segmentStart, vector<int> &segmentEnd)
{
	int nsegment=segmentStart.size();
	int ncombinations;
	Real clashCC, clashCS, clashSS;
	Real fixedClash;
	Real score, clashScore;
	Real bestScore=UNK_REAL;
	Vector segmentOccupancies;
	Vector bestOccupancies;
	Vector tempReal, tempImag;
	Vector fixedReal, fixedImag;
	Vector fixedNonFixedClashes;
	Vector intraSegmentClashes;
	Matrix interSegmentClashes;
	Matrix real, imag;

	//cout <<"Before setNonFixedOccupancyToZero"<<endl;	
	SafeAlloc(segmentOccupancies, nsegment, "segmentOccupancies");
	setNonFixedOccupancyToZero(Protein, segmentStart, segmentEnd);
	//printLocation(__LINE__, __FILE__);
	calcSymmetryScattering(xray, Protein);
	//printLocation(__LINE__, __FILE__);
	fixedReal=xray.real;
	fixedImag=xray.imag;
	fixedClash=calcClashScore(Protein, params, clashCC, clashCS, clashSS);
	//cout <<"fixedClash= "<<fixedClash<<endl;
	//printLocation(__LINE__, __FILE__);
	ncombinations=int(calcPower(2.0, nsegment));

	for (int i=0;i<nsegment;i++)
	{
		calcSegmentScattering(xray, Protein, params, segmentStart[i], segmentEnd[i], tempReal, tempImag);
		SafePushBack(real, tempReal, "real");
		SafePushBack(imag, tempImag, "imag");
	}
	//printLocation(__LINE__, __FILE__);

	calcSegmentClashScores(Protein, params, segmentStart, segmentEnd, fixedClash, fixedNonFixedClashes, intraSegmentClashes, interSegmentClashes);
	//printLocation(__LINE__, __FILE__);
	//cout <<"ncombinations= "<<ncombinations<<endl;
	for (int i=0;i<ncombinations;i++)
	{
		intToBinary(i, segmentOccupancies);
		//printLocation(__LINE__, __FILE__);
		addSegmentScattering(segmentOccupancies, real, imag, fixedReal, fixedImag, xray);
		//printLocation(__LINE__, __FILE__);
		calcIntensity(xray);
		//printLocation(__LINE__, __FILE__);
		if (params.CorrectionFile!="") correctIntensity(xray);
		//printLocation(__LINE__, __FILE__);
		score=calcRFactor(xray, expXRay, params);
		//printLocation(__LINE__, __FILE__);
		clashScore=calcClashScore(segmentOccupancies, fixedClash, fixedNonFixedClashes, intraSegmentClashes, interSegmentClashes);
		//printLocation(__LINE__, __FILE__);
		score+=clashScore;
		//cout <<"score= "<<score<<" bestScore= "<<bestScore<<endl;

		if (score<bestScore || i==0)
		{
			bestScore=score;
			bestOccupancies=segmentOccupancies;
		}
		//printLocation(__LINE__, __FILE__);
	}
	//printLocation(__LINE__, __FILE__);
	setSegmentOccupancies(Protein, bestOccupancies, segmentStart, segmentEnd, 1.0);
	//printLocation(__LINE__, __FILE__);
	addSegmentScattering(bestOccupancies, real, imag, fixedReal, fixedImag, xray);
	calcIntensity(xray);
	if (params.CorrectionFile!="") correctIntensity(xray);

	return bestScore;
}

Real optimizeSegmentsFromFileAfterMR(XRayStruct &xray, XRayStruct &expXRay, ProteinStruct &Protein, XRayParamStruct &params)
{
	vector<int> segmentStart, segmentEnd;
	Real score;

	readSegmentsFile(params.SegmentsFile, segmentStart, segmentEnd);	
	residueSegmentsToAtomSegments(Protein.Atoms, segmentStart, segmentEnd);
	score=optimizeSegmentOccupanciesAfterMR(xray, expXRay, Protein, params, segmentStart, segmentEnd);

	return score;
}

void setIncorrectAtomsToZero(ProteinStruct &Protein, XRayParamStruct &params)
{
	int natom=Protein.Atoms.size();
	Real r, cutoff=2.0;
	ProteinStruct Native;

	readPdb(params.NativePdbFile, Native);

	calcTMScore(Protein.Atoms, Native.Atoms);

	for (int i=0;i<natom;i++)
	{
		r=atomDistance(Protein.Atoms[i], Native.Atoms[i]);
		if (r>cutoff)
		{
			Protein.Atoms[i].Occupancy=0;
		}
	}
}

void readProbabilitiesFile(string probabilitesFile, Real &binSize, Array3D &probabilities, Array3D &trimerProbability)
{
	bool probability=true;
	bool trimer=false;
	vector<string> lines, str;
	int nline, sign=0, bin, nstr;
	int feature1, feature2, feature3;
	Real clashScore;

	readLines(probabilitesFile, lines, "probabilitesFile");
	nline=lines.size();
	Safe3DAlloc(trimerProbability, 4, 4, 4, "trimerProbability");
	Safe3DAlloc(probabilities, 501, 2, 3, "probabilities");

	binSize=0.1;
	//binSize=getValueOfReal("binSize=", lines[0]);
	for (int i=1;i<nline;i++)
	{
		if (lines[i]=="PositiveXRayScore") sign=0;
		if (lines[i]=="NegativeXRayScore") sign=1;
		Tokenize2(lines[i], "\t", str);
		nstr=str.size();
		if (nstr==4 && isReal(str[0]) && probability)
		{
			clashScore=toReal(str[0]);
			bin=int(clashScore/binSize+0.5);
			for (int i=0;i<3;i++)
			{
				probabilities[bin][sign][i]=toReal(str[i+1]);
			}
		}
		if (nstr==4 && isReal(str[0]) && trimer)
		{
			feature1=toInt(str[0]);
			feature2=toInt(str[1]);
			feature3=toInt(str[2]);
			trimerProbability[feature1][feature2][feature3]=toReal(str[3]);
			//cout <<"feature= "<<feature1<<" "<<feature2<<" "<<feature3<<" "<<trimerProbability[feature1][feature2][feature3]<<endl;
		}
		if (nstr!=4 && trimer) trimer=false;
		if (lines[i]=="TrimerProbabilities") trimer=true;
	}
}

void readProbabilitiesFile(string probabilitesFile, Real &binSize, Array3D &probabilities, Array3D &trimerProbability, Real estimatedRmsd)
{
	bool probability=false;
	bool trimer=false;
	vector<string> lines, str;
	int nline, sign=0, bin, nstr;
	int feature1, feature2, feature3;
	Real clashScore;
	Real minRmsd, maxRmsd;

	readLines(probabilitesFile, lines, "probabilitesFile");
	nline=lines.size();
	Safe3DAlloc(trimerProbability, 4, 4, 4, "trimerProbability");
	Safe3DAlloc(probabilities, 501, 2, 3, "probabilities");

	binSize=0.1;
	for (int i=1;i<nline;i++)
	{
		if (lines[i]=="PositiveXRayScore") 
		{
			sign=0;
			probability=true;
			//binSize=getValueOfReal("binSize=", lines[i-1]);
		}
		if (lines[i]=="NegativeXRayScore") 
		{
			sign=1;
			probability=true;
			//binSize=getValueOfReal("binSize=", lines[i-1]);
		}
		Tokenize2(lines[i], "\t", str);
		nstr=str.size();
		if (nstr==4 && isReal(str[0]) && probability)
		{
			clashScore=toReal(str[0]);
			bin=int(clashScore/binSize+0.5);
			for (int i=0;i<3;i++)
			{
				probabilities[bin][sign][i]=toReal(str[i+1]);
			}
		}
		if (probability && !(nstr==1 || nstr==4)) probability=false;
		if (nstr!=4 && trimer && lines[i]!="TrimerProbabilities") trimer=false;
		//cout <<"nstr= "<<nstr<<" "<<lines[i]<<endl;
		if (nstr==2 && lines[i+1]=="TrimerProbabilities")
		{
			minRmsd=toReal(str[0]);
			maxRmsd=toReal(str[1]);
			if (estimatedRmsd>=minRmsd && estimatedRmsd<maxRmsd)
			{
				trimer=true;
			}
		}
		if (nstr==4 && isReal(str[0]) && trimer)
		{
			feature1=toInt(str[0]);
			feature2=toInt(str[1]);
			feature3=toInt(str[2]);
			trimerProbability[feature1][feature2][feature3]=toReal(str[3]);
			//cout <<"feature= "<<feature1<<" "<<feature2<<" "<<feature3<<" "<<trimerProbability[feature1][feature2][feature3]<<endl;
		}
		//if (lines[i]=="TrimerProbabilities") trimer=true;
	}
}

Real probabilityToOccupancy(Array3D &probabilities, Real binSize, Real clashScore, Real xrayScore)
{
	int maxBin=probabilities.size();
	int bin, sign;
	Real sum, occupancy;

	bin=int(clashScore/binSize+0.5);
	if (xrayScore>0) sign=0;
	else sign=1;
	if (bin<maxBin)
	{
		sum=probabilities[bin][sign][0]+probabilities[bin][sign][1];
		if (sum==0) occupancy=0;
		else occupancy=probabilities[bin][sign][0]/sum;
	}
	else occupancy=0;

	return occupancy;
}

void probabilitesToOccupancies(Array3D &probabilities, Real binSize, Vector &residueClashScore, Vector &residueXRayScore, vector<AtomStruct> &Atoms)
{
	int natom=Atoms.size();
	int residue;

	for (int i=0;i<natom;i++)
	{
		residue=Atoms[i].ResidueNum-1;
		Atoms[i].Occupancy=probabilityToOccupancy(probabilities, binSize, residueClashScore[residue], residueXRayScore[residue]);
	}
}

int getFeature(Real clashScore, Real xrayScore)
{
	int feature=0;

	if (clashScore>0) feature+=2;
	if (xrayScore<0) feature++;

	return feature;
}

void getFeatureSequence(Vector &clashScore, Vector &xrayScore, vector<int> &feature)
{
	int nres=clashScore.size();

	SafeAlloc(feature, nres, "feature");
	for (int i=0;i<nres;i++)
	{
		feature[i]=getFeature(clashScore[i], xrayScore[i]);
	}
}

void calcOccupancies(Array3D &probabilities, Array3D &trimerProbabilities, Real binSize, Vector &residueClashScore, Vector &residueXRayScore, Vector &occupancies)
{
	int nres=residueClashScore.size();
	vector<int> feature;

	SafeAlloc(occupancies, nres, "occupancies");

	occupancies[0]=probabilityToOccupancy(probabilities, binSize, residueClashScore[0], residueXRayScore[0]);
	occupancies[nres-1]=probabilityToOccupancy(probabilities, binSize, residueClashScore[nres-1], residueXRayScore[nres-1]);
	getFeatureSequence(residueClashScore, residueXRayScore, feature);

	for (int i=1;i<nres-1;i++)
	{
		occupancies[i]=trimerProbabilities[feature[i-1]][feature[i]][feature[i+1]];
		//cout <<"trimerProb["<<feature[i-1]<<"]["<<feature[i]<<"]["<<feature[i+1]<<"]= "<<trimerProbabilities[feature[i-1]][feature[i]][feature[i+1]]<<endl;
	}
}

void assignOccupancies(Vector &occupancies, vector<AtomStruct> &Atoms)
{
	int natom=Atoms.size();
	int nres=occupancies.size();
	int residue;

	for (int i=0;i<natom;i++)
	{
		residue=Atoms[i].ResidueNum-1;
		if (residue<nres)
		{
			Atoms[i].Occupancy=occupancies[residue];
		}
	}
}

void probabilitiesToOccupancies(Array3D &probabilities, Array3D &trimerProbabilities, Real binSize, Vector &residueClashScore, Vector &residueXRayScore, vector<AtomStruct> &Atoms)
{
	Vector occupancies;

	calcOccupancies(probabilities, trimerProbabilities, binSize, residueClashScore, residueXRayScore, occupancies);
	assignOccupancies(occupancies, Atoms);
}

void calcResidueClashScore(ProteinStruct Protein, XRayParamStruct &params, Vector &residueClashScore)
{
        calcClashScoreByAtom(Protein, params, residueClashScore);
}

void subtractResidueFromScattering(ProteinStruct &Protein, XRayStruct &xray, int residue)
{
        int natom=Protein.Atoms.size();

        for (int i=0;i<natom;i++)
        {
                if (Protein.Atoms[i].ResidueNum==residue)
                {
                        subtractAtomFromScattering(Protein.Atoms[i], Protein, xray);
                }
        }
        calcIntensity(xray);
}

Real calcResidueXRayScore(ProteinStruct &Protein, XRayStruct xray, XRayStruct &expXRay, int residue, XRayParamStruct &params)
{
        subtractResidueFromScattering(Protein, xray, residue);
        if (params.CorrectionFile!="") correctIntensity(xray);
        calcAmplitude(xray);
        return calcRFactor(xray, expXRay, params);
}

void calcResidueXRayScore(ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &residueXRayScore)
{
        int nres=getNumResidues(Protein.Atoms);
        Real xRayScore1, xRayScore2;

        SafeAlloc(residueXRayScore, nres, "residueXRayScore");
        calcCartToFrac(Protein);
        calcSymmetryScattering(xray, Protein);
        calcIntensity(xray);
        if (params.CorrectionFile!="") correctIntensity(xray);
        xRayScore1=calcRFactor(xray, expXRay, params);

        for (int i=1;i<=nres;i++)
        {
                xRayScore2=calcResidueXRayScore(Protein, xray, expXRay, i, params);
                residueXRayScore[i-1]=xRayScore2-xRayScore1;
        }
}

Real calcResidueDensity(ProteinStruct &Protein, int residue, LatticeStruct &lattice, LatticeStruct &expLattice)
{
	int natom=Protein.Atoms.size();
	int nsym=Protein.symmetryOperations.size();
	int xbin, ybin, zbin;
	int maxXBin, maxYBin, maxZBin;
	int pick;
	Real x, y, z;
	Real density=0, count=0;
	Real densityCorrelation;
	Vector densities;
	ProteinStruct tempProtein;

	Get3DVectorSize(lattice.density, maxXBin, maxYBin, maxZBin, "lattice.density");

	for (int i=0;i<nsym;i++)
	{
		tempProtein=Protein;
		applySymmetryOperation(tempProtein.Atoms, Protein.symmetryOperations[i]);
		convertToFractionalCoordinates(tempProtein);
		density=0;
		count=0;
		for (int i=0;i<natom;i++)
		{
			if (Protein.Atoms[i].ResidueNum==residue)
			{
				x=tempProtein.Atoms[i].frac.x;
				y=tempProtein.Atoms[i].frac.y;
				z=tempProtein.Atoms[i].frac.z;
				fractionalPositionToIndex(lattice, x, y, z, xbin, ybin, zbin);
				density+=(expLattice.density[xbin][ybin][zbin]-lattice.density[xbin][ybin][zbin]);
				count+=1.0;
			}
		}
		density/=count;
		SafePushBack(densities, density, "densities");
	}
	densityCorrelation=calcDensityCorrelation(lattice, expLattice);
	cout <<"densityCorrelation= "<<densityCorrelation<<endl;
	sort(densities);
	pick=nsym/4;
	if (pick>=nsym) pick=nsym-1;
	else if (pick<0) pick=0;
	return densities[pick];
}

Real calcResidueDensity(ProteinStruct &Protein, int residue, XRayStruct &xray, XRayStruct &combinedXRay, Real calcScale, Real expScale)
{
	int natom=Protein.Atoms.size();
	int nsym=Protein.symmetryOperations.size();
	int pick;
	Real x, y, z;
	Real density=0, count=0;
	Real expDensity, modelDensity;
	Vector densities;
	ProteinStruct tempProtein;

	for (int i=0;i<nsym;i++)
	{
		tempProtein=Protein;
		applySymmetryOperation(tempProtein.Atoms, Protein.symmetryOperations[i]);
		convertToFractionalCoordinates(tempProtein);
		density=0;
		count=0;
		for (int i=0;i<natom;i++)
		{
			if (Protein.Atoms[i].ResidueNum==residue)
			{
				x=tempProtein.Atoms[i].frac.x;
				y=tempProtein.Atoms[i].frac.y;
				z=tempProtein.Atoms[i].frac.z;
				expDensity=calcDensity(combinedXRay, x, y, z);
				modelDensity=calcDensity(xray, x, y, z);
				density+=(expDensity*expScale-modelDensity*calcScale);
				count+=1.0;
			}
		}
		density/=count;
		SafePushBack(densities, density, "densities");
	}
	sort(densities);
	pick=nsym/4;
	if (pick>=nsym) pick=nsym-1;
	else if (pick<0) pick=0;
	return densities[pick];
}

void calcResidueDensity(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, LatticeStruct &expLattice, Vector &residueDensity)
{
	int nres=getNumResidues(Protein.Atoms);
	Real calcScale, expScale;
	Real densityCorrelation;
	Real numElectrons, numLatticeElectrons, numExpLatticeElectrons;
	XRayStruct combinedXRay;

	SafeAlloc(residueDensity, nres, "residueDensity");
	combinePhaseAndIntensity(xray, expXRay, combinedXRay);
	calcDensity(xray, lattice);
	calcDensity(combinedXRay, expLattice);
	densityCorrelation=calcDensityCorrelation(lattice, expLattice);
	cout <<"densityCorrelation= "<<densityCorrelation<<endl;
	numElectrons=calcProteinElectrons(Protein.Atoms);
	numLatticeElectrons=calcLatticeElectrons(lattice);
	numExpLatticeElectrons=calcLatticeElectrons(expLattice);

	calcScale=numElectrons/numLatticeElectrons;
	expScale=numElectrons/numExpLatticeElectrons;

	for (int i=1;i<=nres;i++)
	{
		//residueDensity[i-1]=calcResidueDensity(Protein, i, lattice, expLattice);
		residueDensity[i-1]=calcResidueDensity(Protein, i, xray, combinedXRay, calcScale, expScale);
	}
}

void calcResidueDensityScore(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &residueDensity)
{
	LatticeStruct lattice, expLattice;

	initializeLattice(lattice, params);
	expLattice=lattice;
	calcDensity(Protein, xray, lattice);
	cout <<"Before calcExperimenalDensity"<<endl;
	calcExperimentalDensity(xray, expXRay, expLattice);
	cout <<"Before normalizeDensity"<<endl;
	normalizeDensity(Protein, lattice);
	cout <<"Before normalizeDensity"<<endl;
	normalizeDensity(Protein, expLattice);
	cout <<"Before calcResidueDensity"<<endl;
	calcResidueDensity(Protein, xray, expXRay, lattice, expLattice, residueDensity);
	cout <<"Leaving calcResidueDensityScore"<<endl;
}

void setOccupancies(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Real estimatedRmsd)
{
	Real binSize=0.1;
	Vector residueClashScore, residueXRayScore;
	Array3D probabilities, trimerProbabilities;

	readProbabilitiesFile(params.ProbabilitiesFile, binSize, probabilities, trimerProbabilities, estimatedRmsd);
	calcResidueClashScore(Protein, params, residueClashScore);
	calcResidueXRayScore(Protein, xray, expXRay, params, residueXRayScore);
	probabilitiesToOccupancies(probabilities, trimerProbabilities, binSize, residueClashScore, residueXRayScore, Protein.Atoms);
}

void setOccupancies(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	Real binSize;
	Vector residueClashScore, residueXRayScore;
	Array3D probabilities, trimerProbabilities;

	readProbabilitiesFile(params.ProbabilitiesFile, binSize, probabilities, trimerProbabilities);
	calcResidueClashScore(Protein, params, residueClashScore);
	calcResidueXRayScore(Protein, xray, expXRay, params, residueXRayScore);
	probabilitiesToOccupancies(probabilities, trimerProbabilities, binSize, residueClashScore, residueXRayScore, Protein.Atoms);
}

Real estimateRmsdFromAveOccupancy(Real aveOccupancy)
{
	return -27.3*aveOccupancy+26.4;
}

Real estimateRmsdFromAveOccupancy(Real aveOccupancy, Real residues)
{
	Real a, b;
	Real correction1, correction2, correction3;

	correction1=-27.3*aveOccupancy+26.4;
	a=-4.5254e-5*residues+0.018172;
	b=0.00097435*residues-0.0042486;
	correction2=a*exp(b*correction1);
	correction3=correction2*284.49-3.832;

	return correction3;
}

Real calcAverageOccupancy(vector<AtomStruct> &Atoms)
{
	int natom=Atoms.size();
	Real sum=0;

	for (int i=0;i<natom;i++)
	{
		sum+=Atoms[i].Occupancy;
	}
	return sum/Real(natom);
}

Real estimateRmsdFromAveOccupancy(vector<AtomStruct> &Atoms)
{
	Real aveOccupancy, rmsd;
	Real residues=Real(getNumResidues(Atoms));

	aveOccupancy=calcAverageOccupancy(Atoms);
	rmsd=estimateRmsdFromAveOccupancy(aveOccupancy, residues);

	return rmsd;
} 

Real estimateRmsd(ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	Real rmsd;

	setOccupancies(Protein, xray, expXRay, params);
	rmsd=estimateRmsdFromAveOccupancy(Protein.Atoms);

	return rmsd;
}

Real findMinEstimatedRmsd(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	int nprot=Proteins.size();
	Real rmsd, minRmsd;
	Vector rmsds;

	for (int i=1;i<nprot;i++)
	{
		assignAtomIDs(Proteins[i].Atoms);
		rmsd=estimateRmsd(Proteins[i], xray, expXRay, params);
		SafePushBack(rmsds, rmsd, "rmsds");
	}
	minRmsd=findMin(rmsds);
	//cout <<"rmsds: "<<endl;
	return minRmsd;
}

Real findMinEstimatedRmsd(string pdbPathsFile, XRayStruct xray, XRayStruct expXRay, XRayParamStruct &params)
{
	Real rmsd;
	vector<ProteinStruct> Proteins;

	readPdbs(pdbPathsFile, Proteins);
	rmsd=findMinEstimatedRmsd(Proteins, xray, expXRay, params);
	cout <<"minEstimatedRmsd= "<<rmsd<<endl;

	return rmsd;
}
