# include <vector>
# include <iostream>

using namespace std;

# include "../../LibraryFiles/MathUtils.h"
# include "../../LibraryFiles/normalize.h"
# include "../../LibraryFiles/MinMax.h"
# include "../../LibraryFiles/time.h"
# include "../../LibraryFiles/rmsd.h"
# include "../../LibraryFiles/ReadPdb.h"
# include "XRayStruct.h"
# include "XRay.h"
# include "EmpiricalMaximumLikelihood.h"

void readErrorHistogram(string errorHistogramFile, ErrorHistogramStruct &errorHistogram)
{
	bool allocated=false;
	vector<string> lines, str;
	int nline, nstr;
	int qBin, ratioBin;
	int maxQBin=0, maxRatioBin=0;

	readLines(errorHistogramFile, lines, "ErrorHistogramFile");
	nline=lines.size();

	for (int i=0;i<nline;i++)
	{
		Tokenize2(lines[i], " ", str);
		nstr=str.size();
		if (nstr==2)
		{
			if (str[0]=="maxQBin=") maxQBin=toInt(str[1]);
			if (str[0]=="maxRatioBin=") maxRatioBin=toInt(str[1]);
			if (str[0]=="qBinSize=") errorHistogram.qBinSize=toReal(str[1]);
			if (str[0]=="ratioBinSize=") errorHistogram.ratioBinSize=toReal(str[1]);
			if (str[0]=="minRatio=") errorHistogram.minRatio=toReal(str[1]);
			Safe2DAlloc(errorHistogram.histogram, maxQBin, maxRatioBin, "errorHistogram");
			allocated=true;
		}
		Tokenize2(lines[i], "\t", str);
		nstr=str.size();
		if (nstr==5 && isReal(str[0]) && allocated)
		{
			qBin=toInt(str[3]);
			ratioBin=toInt(str[4]);
			errorHistogram.histogram[qBin][ratioBin]=toReal(str[2]);
		}
	}
}

void addScatteringToDistribution(XRayStruct &xray, Matrix &distribution)
{
	int npoint=xray.i.size();

	for (int i=0;i<npoint;i++)
	{
		SafePushBack(distribution[i], xray.i[i], "distribution");
	}
}

void randomPlacement(ProteinStruct &Protein)
{
	Real randRX, randRY, randRZ;
	Real randX, randY, randZ;

	randRX=randDouble(0, 0.5*pi);
	randRY=randDouble(0, 0.5*pi);
	randRZ=randDouble(0, 0.5*pi);
	randX=randDouble(0, Protein.XBoxLength);
	randY=randDouble(0, Protein.YBoxLength);
	randZ=randDouble(0, Protein.ZBoxLength);

	rotateAroundX(Protein.Atoms, randRX);
	rotateAroundY(Protein.Atoms, randRY);
	rotateAroundZ(Protein.Atoms, randRZ);
	moveAtoms(Protein.Atoms, randX, 0, 0);
	moveAtoms(Protein.Atoms, 0, randY, 0);
	moveAtoms(Protein.Atoms, 0, 0, randZ);
}

void printCalculatedDistributionHistogram(Vector &histogram, Real binSize, Real min)
{
	int Size=histogram.size();

	cout <<"log(Icalc/Iobs)\tprobability"<<endl;
	for (int i=0;i<Size;i++)
	{
		cout <<Real(i)*binSize+min<<"\t"<<histogram[i]<<endl;
	}
	cout <<endl;
}

void printCalculatedDistribution(Vector &distribution, Real expI)
{
	int Size=distribution.size();
	Real binSize, min;
	Vector v=distribution;
	Vector histogram;

	for (int i=0;i<Size;i++)
	{
		v[i]=log(distribution[i]/expI);
	}
	calcHistogram(v, histogram, binSize, min);
	normalize(histogram);
	vectorMultiply(1.0/binSize, histogram);
	printCalculatedDistributionHistogram(histogram, binSize, min);
}

void printObservedDistribution(XRayStruct &expXRay, ErrorHistogramStruct &errorHistogram, int point)
{
	int Size;
	int qBin;	
        Real ratioBinSize=errorHistogram.ratioBinSize;
        Real minRatio=errorHistogram.minRatio;
        Real correctedObs;

	qBin=int(expXRay.qMag[point]/errorHistogram.qBinSize+0.5);	
	Size=errorHistogram.histogram[qBin].size();

	cout <<"correctedObs\tprobability"<<endl;
	for (int i=0;i<Size;i++)
	{
		correctedObs=Real(i)*ratioBinSize+minRatio;
		cout <<correctedObs<<"\t"<<errorHistogram.histogram[qBin][i]/errorHistogram.ratioBinSize<<endl;
	}
	cout <<endl;
}

void printCalculatedObservedDistributions(Matrix &distribution, XRayStruct &xray, XRayStruct &expXRay, ErrorHistogramStruct &errorHistogram, int point)
{
	Real logRatio;
	int qBin;	
        Real score;

	qBin=int(expXRay.qMag[point]/errorHistogram.qBinSize+0.5);	
	logRatio=log(xray.i[point]/expXRay.i[point]);
	cout <<"log(Icalc/Iobs)= "<<logRatio<<endl;
	printCalculatedDistribution(distribution[point], expXRay.i[point]);
	printObservedDistribution(expXRay, errorHistogram, point);
	score=calcEmpiricalMaximumLikelihood(xray.i[point], expXRay.i[point], errorHistogram, qBin, distribution[point]);

	cout <<"score= "<<score<<endl;
} 

void calcDistributionForEveryPoint(XRayStruct xray, XRayStruct &expXRay, ProteinStruct Protein, XRayParamStruct &params, Matrix &distribution)
{
	int nrandom=1000;
	int nmiller=xray.miller.size();
	Real scale;

	Safe2DAlloc(distribution, nmiller, 0, "distribution");
	for (int i=0;i<nrandom;i++)
	{
		randomPlacement(Protein);
		calcScatteringFast(Protein, xray);
		calcIntensity(xray);
		if (params.CorrectionFile!="") correctIntensity(xray);
		scale=calcScaleForSqrDiff(xray.i, expXRay.i);
		vectorMultiply(scale, xray.i);
		addScatteringToDistribution(xray, distribution);
	}
}

Real calcEmpiricalMaximumLikelihood(Real calc, Real correctedObs, Vector &distribution)
{
	int Size=distribution.size();
	Real prob=1.0;
	Real diff, randomDiff;

	diff=abs(calc-correctedObs);
	for (int i=0;i<Size;i++)
	{
		//cout <<"calc= "<<calc<<" correctedObs= "<<correctedObs<<" distribution= "<<distribution[i]<<endl;
		randomDiff=abs(correctedObs-distribution[i]);
		if (diff>randomDiff) prob+=1.0;
	}
	prob/=Real(Size);	
	return prob;
}

Real calcEmpiricalMaximumLikelihood(Real calc, Real obs, ErrorHistogramStruct &errorHistogram, int qBin, Vector &distribution)
{
	int Size=errorHistogram.histogram[qBin].size();
	Real score=0, pointScore;
	Real ratioBinSize=errorHistogram.ratioBinSize;
	Real minRatio=errorHistogram.minRatio;
	Real correctedObs;

	for (int i=0;i<Size;i++)
	{
		//cout <<"i= "<<i<<" Size= "<<Size<<endl;
		correctedObs=obs*exp(Real(i)*ratioBinSize+minRatio);
		//cout <<"correctedObs= "<<correctedObs<<endl;
		pointScore=calcEmpiricalMaximumLikelihood(calc, correctedObs, distribution);
		//cout <<"pointScore= "<<pointScore<<endl;
		score+=pointScore*errorHistogram.histogram[qBin][i];
	}
	//cout <<"score= "<<score<<endl;
	return log(score);
}

Real calcEmpiricalMaximumLikelihood(XRayStruct &xray, XRayStruct &expXRay, ErrorHistogramStruct &errorHistogram, Matrix &distribution)
{	
	int npoint=xray.i.size();
	int qBin;
	Real score=0;

	//for (int i=256;i<=256;i++)
	for (int i=0;i<npoint;i++)
	{
		qBin=int(xray.qMag[i]/errorHistogram.qBinSize+0.5);
		//cout <<"qBin= "<<qBin<<endl;
		score+=calcEmpiricalMaximumLikelihood(xray.i[i], expXRay.i[i], errorHistogram, qBin, distribution[i]);
		//cout <<"score= "<<score<<endl;
	}
	return score;
}

void readErrorPathsFile(string errorPathsFile, Real &rmsdBinSize, vector<string> &paths)
{
	vector<string> lines, str;
	int nline;

	readLines(errorPathsFile, lines, "errorPathsFile");
	nline=lines.size();
	Tokenize2(lines[0], "\t", str);
	rmsdBinSize=toReal(str[1]);
	for (int i=1;i<nline;i++)
	{
		SafePushBack(paths, lines[i], "paths");
	}
}

void getErrorHistogram(XRayParamStruct &params, ProteinStruct &Protein, ErrorHistogramStruct &errorHistogram)
{
	vector<string> paths;
	int bin;
	Real rmsd, rmsdBinSize;
	ProteinStruct Native;

	if (params.RmsdDependantELL)
	{
		readPdb(params.NativePdbFile, Native);
		rmsd=RMSD(Protein.Atoms, Native.Atoms);
		readErrorPathsFile(params.ErrorPathsFile, rmsdBinSize, paths);
		bin=int(rmsd/rmsdBinSize);
		readErrorHistogram(paths[bin], errorHistogram);
	}
	else
	{
		readErrorHistogram(params.ErrorHistogramFile, errorHistogram);
	}
}

Real calcEmpiricalMaximumLikelihood(XRayStruct &xray, XRayStruct &expXRay, ProteinStruct &Protein, XRayParamStruct &params)
{
	//int npoint=xray.i.size();
	int Size=xray.distribution.size();
	Real score, scale;
	ErrorHistogramStruct errorHistogram;

	//cout <<"In calcEmpiricalMaximumLikelihood"<<endl;

	calcQMagnitudes(xray);
	calcQMagnitudes(expXRay);
	scale=calcScaleForSqrDiff(xray.i, expXRay.i);
	vectorMultiply(scale, xray.i);
	//cout <<"Before readErrorHistogram"<<endl;
	//readErrorHistogram(params.ErrorHistogramFile, errorHistogram);
	getErrorHistogram(params, Protein, errorHistogram);
	//cout <<"Before calcDistributionForEveryPoint"<<endl;
	if (Size==0)
	{
		calcDistributionForEveryPoint(xray, expXRay, Protein, params, xray.distribution);
	}
	//int point=256;
	//cout <<"qMag= "<<xray.qMag[point]<<endl;
	//printCalculatedObservedDistributions(xray.distribution, xray, expXRay, errorHistogram, point);
	//for (int i=0;i<npoint;i++)
	//{
	//	cout <<"point= "<<i<<endl;
	//	printCalculatedObservedDistributions(distribution, xray, expXRay, errorHistogram, i);
	//	cout <<endl<<endl<<endl<<endl;
	//}
	//endProgram(__LINE__, __FILE__);
	//cout <<"Before calcEmpiricalMaximumLikelihood"<<endl;
	score=calcEmpiricalMaximumLikelihood(xray, expXRay, errorHistogram, xray.distribution);
	//cout <<"score= "<<score<<endl;

	return score;
}
