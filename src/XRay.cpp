# include <sys/time.h>

using namespace std;

# include "../../LibraryFiles/TypeDef.h"
# include "../../LibraryFiles/Structures.h"
# include "../../LibraryFiles/normalize.h"
# include "../../LibraryFiles/FormFactors.h"
# include "../../LibraryFiles/time.h"
# include "../../LibraryFiles/MemoryUsage.h"
# include "../../LibraryFiles/AtomUtils.h"
# include "../../LibraryFiles/MathStructures.h"
# include "../../LibraryFiles/MathUtils.h"
# include "../../LibraryFiles/VectorManip.h"
# include "../../LibraryFiles/IOUtils.h"
# include "../../LibraryFiles/ReadPdb.h"
# include "../../LibraryFiles/PrintPdb.h"
# include "../../LibraryFiles/GetNextFrame.h"
# include "../../LibraryFiles/MinMax.h"
# include "../../LibraryFiles/error.h"
# include "molecularReplacement.h"
# include "XRayStruct.h"
# include "XRay.h"
# include "Initialize.h"
# include "LatticeStruct.h"
# include "LatticeUtils.h"
# include "SpaceGroup.h"
# include "RMSD.h"
# include "FlexibleEditing.h"

XRayStruct&  XRayStruct::operator = (const XRayStruct &xray)
{
	centric=xray.centric;
	free=xray.free;
	i=xray.i;
	imag=xray.imag;
	real=xray.real;
	waterReal=xray.waterReal;
	waterImag=xray.waterImag;
	phase=xray.phase;
	error=xray.error;
	qMag=xray.qMag;
	F=xray.F;
	Ferror=xray.Ferror;
	Ffix=xray.Ffix;
	Fmove=xray.Fmove;
	realBlock=xray.realBlock;
	imagBlock=xray.imagBlock;
	correction=xray.correction;
	cosx=xray.cosx;
	cosy=xray.cosy;
	cosz=xray.cosz;
	sinx=xray.sinx;
	siny=xray.siny;
	sinz=xray.sinz;
	q=xray.q;
	miller=xray.miller;
	a=xray.a;
	b=xray.b;
	c=xray.c;
	alpha=xray.alpha;
	beta=xray.beta;
	gamma=xray.gamma;
	correctionBinSize=xray.correctionBinSize;
	hBinSize=xray.hBinSize;
	kBinSize=xray.kBinSize;
	lBinSize=xray.lBinSize;
	maxH=xray.maxH;
	maxK=xray.maxK;
	maxL=xray.maxL;
	minH=xray.minH;
	minK=xray.minK;
	minL=xray.minL;
	maxContinuousH=xray.maxContinuousH;
	maxContinuousK=xray.maxContinuousK;
	maxContinuousL=xray.maxContinuousL;
	realContinuous=xray.realContinuous;
	imagContinuous=xray.imagContinuous;
	//i3d=xray.i3d;
	//imag3d=xray.imag3d;
	//real3d=xray.real3d;
	//phase3d=xray.phase3d;
	q3d=xray.q3d;
	miller3d=xray.miller3d;
	dwf=xray.dwf;
	f=xray.f;
	fa=xray.fa;
	complexReal=xray.complexReal;
	complexImag=xray.complexImag;
	symReal=xray.symReal;
	symImag=xray.symImag;
	realBlockTranslated=xray.realBlockTranslated;
	imagBlockTranslated=xray.imagBlockTranslated;
	symTempReal=xray.symTempReal;
	symTempImag=xray.symTempImag;
	cosx1d=xray.cosx1d;
	cosy1d=xray.cosy1d;
	cosz1d=xray.cosz1d;
	sinx1d=xray.sinx1d;
	siny1d=xray.siny1d;
	sinz1d=xray.sinz1d;
	dwfx=xray.dwfx;
	dwfy=xray.dwfy;
	dwfz=xray.dwfz;
	matTrig1=xray.matTrig1;
	matTrig2=xray.matTrig2;
	complexSymReal=xray.complexSymReal;
	complexSymImag=xray.complexSymImag;
	complexSymTempReal=xray.complexSymTempReal;
	complexSymTempImag=xray.complexSymTempImag;
	patterson=xray.patterson;
	lookUp=xray.lookUp;
	clashScore=xray.clashScore;
	//boolVerbose=xray.boolVerbose;
	Proteins=xray.Proteins;
	OccupancyCorrection.occupancyBinSize=xray.OccupancyCorrection.occupancyBinSize;
	OccupancyCorrection.correctionBinSize=xray.OccupancyCorrection.correctionBinSize;
	OccupancyCorrection.correction=xray.OccupancyCorrection.correction;
	score=xray.score;
	stdDev=xray.stdDev;
	scoreIndex=xray.scoreIndex;
	
	return *this;
}

void millerToQ(Real h, Real k, Real l, Real a, Real b, Real c, PosStruct &q)
{
	q.Pos[X]=2.0*pi*h/a;
	q.Pos[Y]=2.0*pi*k/b;
	q.Pos[Z]=2.0*pi*l/c;
}

void millerToQ(XRayStruct &xray, Matrix &cartToFrac)
{
	int npoint=xray.q.size();
	int nmiller=xray.miller.size();
	int nqmag=xray.qMag.size();
	cout <<"In millerToQ npoint= "<<npoint<<" nmiller= "<<nmiller<<endl;
	if (nmiller==0) error("nmiller=0", __LINE__, __FILE__);
	if (npoint<nmiller) xray.q=xray.miller;
	if (nqmag!=nmiller) SafeAlloc(xray.qMag, nmiller, "qMag");
	for (int i=0;i<nmiller;i++)
	{
		matrixMultiply(xray.miller[i].Pos, cartToFrac, xray.q[i].Pos);
		for (int j=0;j<3;j++) xray.q[i].Pos[j]*=2.0*pi;
		xray.qMag[i]=calcMagnitude(xray.q[i]);
	}
	cout <<"xray.q.size()= "<<xray.q.size()<<endl;
	cout <<"Leaving millerToQ"<<endl;
}

void millerToQ(XRayStruct &xray, XRayParamStruct &params)
{
	Matrix fracToCart, cartToFrac;
	ProteinStruct tempProtein;

	calcUnitCellVectors(tempProtein, params);
	calcFracToCart(tempProtein.a, tempProtein.b, tempProtein.c, fracToCart);
	calcInverseMatrix(fracToCart, cartToFrac);

	millerToQ(xray, cartToFrac);
}

void millerToQ(XRayStruct &xray, Real a, Real b, Real c)
{
	int npoint=xray.q.size();
	int nmiller=xray.miller.size();
	int nqmag=xray.qMag.size();
	cout <<"In millerToQ npoint= "<<npoint<<" nmiller= "<<nmiller<<endl;
	cout <<"a= "<<a<<" b= "<<b<<" c= "<<c<<endl;
	if (nmiller==0) error("nmiller=0", __LINE__, __FILE__);
	if (npoint<nmiller) xray.q=xray.miller;
	if (nqmag!=nmiller) SafeAlloc(xray.qMag, nmiller, "qMag");
	for (int i=0;i<nmiller;i++)
	{
		xray.q[i].Pos[X]=2.0*pi*xray.miller[i].Pos[X]/a;
		xray.q[i].Pos[Y]=2.0*pi*xray.miller[i].Pos[Y]/b;
		xray.q[i].Pos[Z]=2.0*pi*xray.miller[i].Pos[Z]/c;
		xray.qMag[i]=calcMagnitude(xray.q[i]);
	}
	cout <<"xray.q.size()= "<<xray.q.size()<<endl;
	cout <<"Leaving millerToQ"<<endl;
}

void millerToQ(XRayStruct &xray, LatticeStruct &lattice)
{
	cout <<"In millerToQ(XRayStruct &xray, LatticeStruct &lattice)"<<endl;
	millerToQ(xray, lattice.a, lattice.b, lattice.c);
	cout <<"Leaving millerToQ2"<<endl;
}

void calcIntensity(Vector &real, Vector &imag, Vector &intensity)
{
	int nreal=real.size();
	int nimag=imag.size();
	int nintensity=intensity.size();

	if (nreal==0 || nreal!=nimag || (nintensity!=nreal && nintensity!=0))
	{
		string errorStr;
		errorStr="nintensity= "+IntToStr(nintensity);
		errorStr+=" nreal= "+IntToStr(nreal);
		errorStr+=" nimag= "+IntToStr(nimag);
		error(errorStr, __LINE__, __FILE__);
	}

	if (nintensity!=nreal) SafeAlloc(intensity, nreal, "intensity");

	for (int i=0;i<nreal;i++)
	{
		intensity[i]=real[i]*real[i]+imag[i]*imag[i];
		//cout <<"i= "<<i<<" intensity= "<<intensity[i]<<" real= "<<real[i]<<" imag= "<<imag[i]<<endl;
	}
}

void calcIntensity(XRayStruct &xray)
{
	calcIntensity(xray.real, xray.imag, xray.i);
}

void calcAmplitude(Vector &intensity, Vector &F)
{
	int nintensity=intensity.size();
	int namplitude=F.size();

	if (namplitude!=nintensity) SafeAlloc(F, nintensity, "F");

	for (int i=0;i<nintensity;i++)
	{
		F[i]=sqrt(intensity[i]);
	}
}

void calcAmplitude(Vector &real, Vector &imag, Vector &F)
{
	int nreal=real.size();
	int nimag=imag.size();
	int namplitude=F.size();	

	if (nimag!=nreal || nimag==0)
	{
		string errorStr;
		errorStr="nreal= "+IntToStr(nreal);
		errorStr+=" nimag= "+IntToStr(nimag);
		error(errorStr, __LINE__, __FILE__);
	}

	if (namplitude!=nreal) SafeAlloc(F, nreal, "F");

	for (int i=0;i<nreal;i++)
	{
		F[i]=sqrt(real[i]*real[i]+imag[i]*imag[i]);
	}

}

void calcAmplitude(Vector &intensity, Vector &real, Vector &imag, Vector &F)
{
	int nintensity=intensity.size();
	int nreal=real.size();
	int nimag=imag.size();

	if (nimag!=nreal || (nimag==0 && nintensity==0))
	{
		string errorStr;
		errorStr="nintensity= "+IntToStr(nintensity);
		errorStr+=" nreal= "+IntToStr(nreal);
		errorStr+=" nimag= "+IntToStr(nimag);
		error(errorStr, __LINE__, __FILE__);
	}

	if (nintensity!=0) calcAmplitude(intensity, F);
	else if (nreal!=0) calcAmplitude(real, imag, F);
	else error("This line should not be reached.", __LINE__, __FILE__);
}

void calcAmplitude(XRayStruct &xray)
{
	calcAmplitude(xray.i, xray.real, xray.imag, xray.F);
}

void checkPhase(XRayStruct &xray)
{
	int npoint=xray.q.size();
	int nphase=xray.phase.size();
	Real d=0.1;
	timeval start, end;
	gettimeofday(&start, NULL);

	if (npoint!=nphase)
	{
		cout <<"ERROR: nphase!=npoint"<<endl;
		exit(EXIT_FAILURE);
	}

	for (int i=0;i<npoint;i++)
	{
		if (abs(xray.real[i]-sqrt(xray.i[i])*cos(xray.phase[i]))<d || abs(xray.imag[i]-sqrt(xray.i[i])*sin(xray.phase[i]))<d)
		{
			cout <<"ERROR: Phases are not consistent"<<endl;
			cout <<"real= "<<xray.real[i]<<" imag= "<<xray.imag[i]<<" i= "<<xray.i[i]<<" phase= "<<xray.phase[i]<<endl;
			exit(EXIT_FAILURE);
		}
	}
	gettimeofday(&end, NULL);
	//cout <<"calcPhase took "<<calcTimeDiff(start, end)<<endl;
}

void calcPhase(XRayStruct &xray)
{
	int nreal=xray.real.size();
	int nimag=xray.imag.size();
	int nphase=xray.phase.size();
	timeval start, end;
	gettimeofday(&start, NULL);

	if (nreal!=nimag || nreal==0)
	{
		string errorStr;
		errorStr="nreal= "+IntToStr(nreal)+" nimag= "+IntToStr(nimag);
		error(errorStr, __LINE__, __FILE__);
	}

	if (nphase!=nreal)
	{
		SafeAlloc(xray.phase, nreal, "xray.phase");
	}

	for (int i=0;i<nreal;i++)
	{
		xray.phase[i]=atan2(xray.imag[i], xray.real[i]);
	}
	gettimeofday(&end, NULL);
	//cout <<"calcPhase took "<<calcTimeDiff(start, end)<<endl;
}

void calcTotalFormFactor(Matrix &realMat, Matrix &imagMat, Vector &real, Vector &imag)
{
	int nxray, npoint;

	Get2DVectorSize(realMat, nxray, npoint);

	zeroVector(real);
	zeroVector(imag);
	for (int i=0;i<nxray;i++)
	{
		for (int j=0;j<npoint;j++)
		{
			real[j]+=realMat[i][j];
			imag[j]+=imagMat[i][j];
		}
	}

}

void calcTotalFormFactor(Array3D &realMat, Array3D &imagMat, Vector &real, Vector &imag)
{
	int nprot, nxray, npoint;

	Get3DVectorSize(realMat, nprot, nxray, npoint);
	//cout <<"nprot= "<<nprot<<" nxray= "<<nxray<<" npoint= "<<npoint<<endl;
	zeroVector(real);
	zeroVector(imag);
	for (int i=0;i<nprot;i++)
	{
		for (int j=0;j<nxray;j++)
		{
			for (int k=0;k<npoint;k++)
			{
				real[k]+=realMat[i][j][k];
				imag[k]+=imagMat[i][j][k];
				//cout <<"realMat= "<<realMat[i][j][k]<<" real= "<<real[k]<<" i= "<<i<<" j= "<<j<<" k= "<<k<<endl;
			}
		}
	}

}

void addExcludedScattering(XRayStruct &xray)
{
	int nreal=xray.real.size();
	int nexcluded=xray.excludedReal.size();

	if (nreal==nexcluded)
	{
		for (int i=0;i<nreal;i++)
		{
			//cout <<"i= "<<i<<" real= "<<xray.real[i]<<" excludedReal= "<<xray.excludedReal[i]<<" imag= "<<xray.imag[i]<<" xray.excludedReal= "<<xray.excludedReal[i]<<endl;
			xray.real[i]+=xray.excludedReal[i];
			xray.imag[i]+=xray.excludedImag[i];
		}
	}
}

void calcTotalFormFactor(XRayStruct &xray)
{
	calcTotalFormFactor(xray.complexSymTempReal, xray.complexSymTempImag, xray.real, xray.imag);
	addExcludedScattering(xray);
}

void calcBlockFormFactor(Matrix &realIn, Matrix &imagIn, Vector &real, Vector &imag, int nunique)
{
	int nsym, npoint;
	int nreal=real.size();

	Get2DVectorSize(realIn, nsym, npoint, "realIn");
	if (nunique>nsym)
	{
		string errorStr="nsym= "+toStr(nsym);
		errorStr+=" nunique= "+toStr(nunique);
		error(errorStr, __LINE__, __FILE__);
	}
	if (nreal!=npoint)
	{
		SafeAlloc(real, npoint, "real");
		SafeAlloc(imag, npoint, "imag");
	}
	for (int i=0;i<nunique;i++)
	{
		for (int j=0;j<npoint;j++)
		{
			real[j]+=realIn[i][j];
			imag[j]+=imagIn[i][j];
		}
	}
}

void calcTotalTempSymFormFactor(XRayStruct &xray)
{
	int nxray=xray.complexReal.size();
	int npoint=xray.complexReal[0].size();

	zeroVector(xray.real);
	zeroVector(xray.imag);
	for (int i=0;i<nxray;i++)
	{
		for (int j=0;j<npoint;j++)
		{
			xray.real[j]+=xray.complexReal[i][j];
			xray.imag[j]+=xray.complexImag[i][j];
		}
	}
}

void printComponents(XRayStruct xray)
{
	string XRayOutput;
	int nxray=xray.complexReal.size();

	for (int i=0;i<nxray;i++)
	{
		xray.real=xray.complexReal[i];
		xray.imag=xray.complexImag[i];
		XRayOutput="/home/jouko/project/mr/placeProtein_"+IntToStr(i)+".txt";
		printXRayJV(XRayOutput,  xray);
	}
}

void convertToFractionalCoordinates(AtomStruct &Atom, VectorStruct &a, VectorStruct &b, VectorStruct &c)
{
	Vector variable;
	Matrix coefficient;

	SafeAlloc(variable, 3, "variable");
	Safe2DAlloc(coefficient, 4, 3, "coefficient");
	//cout <<"a= "<<vectorToStr(a)<<endl;
	//cout <<"b= "<<vectorToStr(b)<<endl;
	//cout <<"c= "<<vectorToStr(c)<<endl;
	coefficient[0][0]=a.x;
	coefficient[1][0]=b.x;
	coefficient[2][0]=c.x;
	coefficient[3][0]=Atom.x;

	coefficient[0][1]=a.y;
	coefficient[1][1]=b.y;
	coefficient[2][1]=c.y;
	coefficient[3][1]=Atom.y;

	coefficient[0][2]=a.z;
	coefficient[1][2]=b.z;
	coefficient[2][2]=c.z;
	coefficient[3][2]=Atom.z;

	equation(coefficient, variable);

	Atom.frac.x=variable[0];
	Atom.frac.y=variable[1];
	Atom.frac.z=variable[2];

	if (isinf(Atom.frac.x) || isnan(Atom.frac.x))
	{
		string errorStr;
		printMatrix(coefficient);
		errorStr="Unable to convert to fractional coordinates";
		error(errorStr, __LINE__, __FILE__);
	}
}

/*
   void convertToFractionalCoordinates(ProteinStruct &Protein)
   {
   int natom=Protein.Atoms.size();
   for (int i=0;i<natom;i++)
   {
   convertToFractionalCoordinates(Protein.Atoms[i], Protein.a, Protein.b, Protein.c);
   }
   }
 */

Real scatteringFactor(Real a[][NumAtomTypes], Real b[][NumAtomTypes], Real c[], int atom, Real q)
{
	Real f=0;
	for (int i=0;i<5;i++)
	{
		f+=a[i][atom]/exp(q*q*b[i][atom]);
	}
	f+=c[atom];

	return f;
}

Real scatteringFactorDer(Real h, Real k, Real l, Real a[][NumAtomTypes], Real b[][NumAtomTypes], int atom, Real unita, Real unitb, Real unitc, int index)
{
	Real f=0, q, exponent;

	q=(h/unita+k/unitb+l/unitc);
	for (int i=0;i<5;i++)
	{
		exponent=1.0/exp(4.0*pi*pi*b[i][atom]*q*q);
		f-=8.0*pi*a[i][atom]*b[i][atom]*q*exponent;
	}
	if (index==X) f/=unita;
	if (index==Y) f/=unitb;
	if (index==Z) f/=unitc;

	return f;
}

Real scatteringFactor(Real a[][NumAtomTypes], Real b[][NumAtomTypes], Real c[], int atom, PosStruct q)
{
	Real mag=0;

	for (int i=0;i<3;i++)
	{
		mag+=q.Pos[i]*q.Pos[i];
	}
	mag=sqrt(mag);

	return scatteringFactor(a, b, c, atom, mag);
}

Real scatteringFactor(int atom, Real q, Real scale)
{
	Real a[5][NumAtomTypes], b[5][NumAtomTypes], c[NumAtomTypes];

	FourGaussianParameters(a, b, c, scale);
	return scatteringFactor(a, b, c, atom, q);
}

void scatteringFactor(Matrix &f, vector<PosStruct> &q, Real scale)
{
	int npoint=q.size();
	Real a[5][NumAtomTypes], b[5][NumAtomTypes], c[NumAtomTypes];
	FourGaussianParameters(a, b, c, scale);
	Safe2DAlloc(f, NumAtomTypes, npoint, "freal");

	for (int i=0;i<NumAtomTypes;i++)
	{
		for (int j=0;j<npoint;j++)
		{
			f[i][j]=scatteringFactor(a, b, c, i, q[j]);
		}
	}
}

void setExcludedVolumeRadii(Vector &excludedVolumeRadii)
{
	SafeAlloc(excludedVolumeRadii, NumAtomTypes, "excludedVolumeRadii");
	excludedVolumeRadii[HYDROGEN]=0.902255;
	excludedVolumeRadii[CARBON]=1.43693;
	excludedVolumeRadii[NITROGEN]=1.24527;
	excludedVolumeRadii[OXYGEN]=1.22099;
	excludedVolumeRadii[SULFUR]=2.19596;	
}

Real gaussianScattering(Real r, Real density, Real q)
{
	return density*r*r*r*pi32/exp(q*q*r*r*0.25);
}

Real hardSphereScattering(Real r, Real density, Real q)
{
	//cout <<"In hardSphereScattering"<<endl;
	//endProgram(__LINE__, __FILE__);
	if (q==0) return density*4.0*pi/3.0;
	return density*4.0*pi*(sin(r*q)-r*q*cos(r*q))/(q*q*q);
}

Real gaussianScatteringDer(Real r, Real density, Real h, Real k, Real l, Real a, Real b, Real c, int index)
{
	Real g, q;

	q=h/a+k/b+l/c;
	g=density*r*r*r*r*r*pi32*pi*pi*q/exp(pi*pi*q*q*r*r);
	if (index==X) g/=a;
	if (index==Y) g/=b;
	if (index==Z) q/=c;

	return g;
}

Real calcMagnitude(PosStruct &q)
{
	Real mag=0;

	for (int i=0;i<3;i++)
	{
		mag+=q.Pos[i]*q.Pos[i];
	}
	mag=sqrt(mag);

	return mag;
}

void calcQMagnitudes(XRayStruct &xray)
{
	int nq=xray.q.size();

	cout <<"nq= "<<nq<<endl;
	SafeAlloc(xray.qMag, nq, "qMag");
	for (int i=0;i<nq;i++)
	{
		xray.qMag[i]=calcMagnitude(xray.q[i]);	
	}
}

Real gaussianScattering(Real r, Real density, PosStruct q)
{
	Real mag=calcMagnitude(q);

	return gaussianScattering(r, density, mag);
}

Real hardSphereScattering(Real r, Real density, PosStruct q)
{
	Real mag=calcMagnitude(q);

	return hardSphereScattering(r, density, mag);
}

void solventDummyAtomScatteringFactor(Matrix &g, vector<PosStruct> &q, Vector &excludedVolumeRadii, Real density)
{
	int npoint=q.size();
	Safe2DAlloc(g, NumAtomTypes, npoint, "freal");

	for (int i=0;i<NumAtomTypes;i++)
	{
		//cout <<endl;
		//cout <<"i= "<<i<<endl;
		for (int j=0;j<npoint;j++)
		{
			g[i][j]=gaussianScattering(excludedVolumeRadii[i], density, q[j]);
			//g[i][j]=hardSphereScattering(excludedVolumeRadii[i], density, q[j]);
			//cout <<calcMagnitude(q[j])<<"\t"<<g[i][j]<<endl;
		}
	}
	//endProgram(__LINE__, __FILE__);
}

void solventCorrectedScatteringFactor(Matrix &f, Matrix &g)
{
	int npoint, ntype;

	Get2DVectorSize(f, ntype, npoint, "f");

	for (int i=0;i<npoint;i++)
	{
		for (int j=0;j<ntype;j++)
		{
			f[j][i]-=g[j][i];
		}
	}
}

void solventCorrectedScatteringFactor(Matrix &f, vector<PosStruct> &q, Vector &excludedVolumeRadii, Real density, Real scale)
{
	Matrix g;
	scatteringFactor(f, q, scale);
	solventDummyAtomScatteringFactor(g, q, excludedVolumeRadii, density);
	solventCorrectedScatteringFactor(f, g);
}

void solventCorrectedScatteringFactor(Matrix &f, vector<PosStruct> &q, XRayParamStruct &params)
{
	int nq=q.size();
	Real density;
	Vector excludedVolumeRadii;

	if (nq==0) error("nq=0", __LINE__, __FILE__);

	if (params.CubeFile=="" && params.ExcludedVolumeMethod!="Cube") density=params.BulkDensity;
	else density=0;
	if (params.ExcludedVolumeRadii=="Default") setExcludedVolumeRadii(excludedVolumeRadii);
	else setExcludedVolumeRadii(excludedVolumeRadii, params);	
	solventCorrectedScatteringFactor(f, q, excludedVolumeRadii, density, params.scale);
}

Real solventCorrectedScatteringFactor(PosStruct &q, Real a[5][NumAtomTypes], Real b[5][NumAtomTypes], Real c[NumAtomTypes], int atom, Real excludedVolumeRadii, Real density)
{
	Real f, g;


	f=scatteringFactor(a, b, c, atom, q);
	g=gaussianScattering(excludedVolumeRadii, density, q);
	//g=hardSphereScattering(excludedVolumeRadii, density, q);

	return f-g;
}

Real solventCorrectedScatteringFactorDer(int h, int k, int l, Real a[5][NumAtomTypes], Real b[5][NumAtomTypes], int atom, Real excludedVolumeRadii, Real density, VectorStruct &unita, VectorStruct &unitb, VectorStruct &unitc, int index)
{
	Real f, g;
	Real alength, blength, clength;

	alength=vectorLength(unita);
	blength=vectorLength(unitb);
	clength=vectorLength(unitc);

	f=scatteringFactorDer(h, k, l, a, b, atom, alength, blength, clength, index);
	g=gaussianScatteringDer(h, k, l, excludedVolumeRadii, density, alength, blength, clength, index);

	return f-g;
}

Real calcDWF(PosStruct &q, Real bfactor)
{
	Real magnitude=calcMagnitude(q);
	bfactor/=(48.0*pi*pi);
	return 1.0/exp(magnitude*magnitude*bfactor);
}

Real calcDWF(PosStruct &q, Vector &anisou)
{
	Real sum=0;
	for (int i=0;i<3;i++) sum+=anisou[i]*q.Pos[i]*q.Pos[i];
	sum+=2.0*anisou[3]*q.Pos[X]*q.Pos[Y];
	sum+=2.0*anisou[4]*q.Pos[X]*q.Pos[Z];
	sum+=2.0*anisou[5]*q.Pos[Y]*q.Pos[Z];
	sum/=3.0;
	return 1.0/exp(sum);	
}

Real calcDWF(PosStruct &q, AtomStruct &Atom)
{
	int nanisou=Atom.anisou.size();

	if (nanisou==6) return calcDWF(q, Atom.anisou);
	else return calcDWF(q, Atom.BFactor);
}

Real calcDWF(PosStruct &q, Real bfactor, LookUpStruct &lookUp)
{
	Real magnitude=calcMagnitude(q);
	bfactor/=(48.0*pi*pi);
	//return 1.0/getExpLookUp(lookUp, magnitude*magnitude*bfactor);
	return getNegativeExpLookUp(lookUp, magnitude*magnitude*bfactor);
}

Real calcDWF(Real qMag, Real bfactor, LookUpStruct &lookUp)
{
	Real bfactorScale=0.002110858; // 1.0/(48.0*pi*pi)
	bfactor*=bfactorScale;
	return getNegativeExpLookUp(lookUp, qMag*qMag*bfactor);
}

Real calcDWF(PosStruct &q, Vector &anisou, LookUpStruct &lookUp)
{
	Real sum=0;
	for (int i=0;i<3;i++) sum+=anisou[i]*q.Pos[i]*q.Pos[i];
	sum+=2.0*anisou[3]*q.Pos[X]*q.Pos[Y];
	sum+=2.0*anisou[4]*q.Pos[X]*q.Pos[Z];
	sum+=2.0*anisou[5]*q.Pos[Y]*q.Pos[Z];
	sum/=3.0;
	return 1.0/getExpLookUp(lookUp, sum);	
}

Real calcDWF(PosStruct &q, AtomStruct &Atom, LookUpStruct &lookUp)
{
	int nanisou=Atom.anisou.size();

	if (nanisou==6) return calcDWF(q, Atom.anisou, lookUp);
	else return calcDWF(q, Atom.BFactor, lookUp);
}

Real calcDWF(XRayStruct &xray, int index, AtomStruct &Atom, LookUpStruct &lookUp)
{
	int nanisou=Atom.anisou.size();
	//cout <<"index= "<<index<<endl;
	//cout <<"xray.qMag= "<<xray.qMag[index]<<endl;
	//cout <<"Atom.BFactor= "<<Atom.BFactor<<endl;
	if (nanisou==6) return calcDWF(xray.q[index], Atom.anisou, lookUp);
	else return calcDWF(xray.qMag[index], Atom.BFactor, lookUp);
}

void calcScattering(ProteinStruct &Protein, XRayStruct &xray)
{
	int atomID;
	int npoint=xray.q.size();
	int nreal=xray.real.size();
	int natom=Protein.Atoms.size();
	Real dotproduct;
	timeval start, end;
	if (npoint!=nreal)
	{
		cout <<"ERROR: npoint!=nreal."<<endl;
		cout <<"npoint= "<<npoint<<" nreal= "<<nreal<<endl;
		exit(EXIT_FAILURE);
	}
	gettimeofday(&start, NULL);
	zeroVector(xray.i);
	zeroVector(xray.real);
	zeroVector(xray.imag);
	for (int i=0;i<npoint;i++)
	{
		for (int j=0;j<natom;j++)
		{
			atomID=Protein.Atoms[j].atomid;
			//freal*=Protein.Atoms[j].Occupancy;
			//freal=Protein.Atoms[j].BFactor;
			dotproduct=Protein.Atoms[j].x*xray.q[i].Pos[X];
			dotproduct+=Protein.Atoms[j].y*xray.q[i].Pos[Y];
			dotproduct+=Protein.Atoms[j].z*xray.q[i].Pos[Z];
			if (isnan(dotproduct))
			{
				cout <<"ERROR: Dotproduct is nan"<<endl;
				exit(EXIT_FAILURE);
			}
			//cout <<"f["<<i<<"]["<<atomID<<"]= "<<xray.f[i][atomID]<<endl;
			xray.real[i]+=xray.f[i][atomID]*cos(dotproduct);
			xray.imag[i]+=xray.f[i][atomID]*sin(dotproduct);
		}	
	}
	gettimeofday(&end, NULL);
	calcIntensity(xray);
	calcPhase(xray);
}

inline Real getCosLookUp(LookUpStruct &lookUp, Real &x)
{
	int bin;
	int maxBin=lookUp.cosLookUp.size();
	while (x>2.0*pi) x-=2.0*pi;
	while (x<0) x+=2.0*pi;	
	bin=int(x/lookUp.cosBin+0.5);
	if (bin>0 && bin<maxBin)
	{
		return lookUp.cosLookUp[bin];
	}
	else return cos(x);
}

inline Real getSinLookUp(LookUpStruct &lookUp, Real &x)
{
	int bin;
	int maxBin=lookUp.sinLookUp.size();
	while (x>2.0*pi) x-=2.0*pi;
	while (x<0) x+=2.0*pi;	
	bin=int(x/lookUp.sinBin+0.5);
	if (bin>0 && bin<maxBin)
	{
		return lookUp.sinLookUp[bin];
	}
	else return sin(x);
}

void calcScatteringBfactor(ProteinStruct &Protein, XRayStruct &xray)
{
	int atomID, natom;
	int npoint=xray.q.size();
	int nreal=xray.real.size();
	int nf1, nf2;
	int noperations;
	Real dotproduct, dwf;
	Real occupancy, product;
	Real cosdp, sindp;
	Real fracx, fracy, fracz;
	Vector real, imag;
	AtomStruct tempAtom;
	timeval fullStart, fullEnd;

	gettimeofday(&fullStart, NULL);
	Get2DVectorSize(xray.f, nf1, nf2);
	if (npoint!=nreal || npoint!=nf2 || npoint==0)
	{
		string errorStr="npoint!=nreal.";
		errorStr+=" npoint= "+toStr(npoint);
		errorStr+=" nreal= "+toStr(nreal);
		errorStr+=" nf= "+toStr(nf2);
		error(errorStr, __LINE__, __FILE__);
	}
	zeroVector(xray.i);
	zeroVector(xray.real);
	zeroVector(xray.imag);
	zeroVector(xray.waterReal);
	zeroVector(xray.waterImag);
	real=xray.real;
	imag=xray.imag;
	calcCartToFrac(Protein);
	natom=Protein.Atoms.size();
	noperations=Protein.symmetryOperations.size();
	if (noperations==0)
	{
		setSpaceGroupP1(Protein);
		noperations=1;
	}
	for (int j=0;j<natom;j++)
	{
		atomID=Protein.Atoms[j].atomid;
		occupancy=Protein.Atoms[j].Occupancy;
		for (int i=0;i<noperations;i++)
		{
			applySymmetryOperation(Protein.Atoms[j], Protein.symmetryOperations[i], tempAtom);
			convertToFractionalCoordinates(tempAtom.x, tempAtom.y, tempAtom.z, Protein.cartToFrac, fracx, fracy, fracz);
			for (int k=0;k<npoint;k++)
			{
				dotproduct=fracx*xray.miller[k].Pos[X];
				dotproduct+=fracy*xray.miller[k].Pos[Y];
				dotproduct+=fracz*xray.miller[k].Pos[Z];
				dotproduct*=2.0*pi;
				dwf=calcDWF(xray, k, tempAtom, xray.lookUp);
				getTrigLookUp(xray.lookUp, dotproduct, cosdp, sindp);
				product=xray.f[atomID][k]*dwf*occupancy;
				xray.real[k]+=product*cosdp;
				xray.imag[k]+=product*sindp;
				//printAtomInfo(Protein.Atoms[j]);
				//cout <<"fracx= "<<fracx<<" fracy= "<<fracy<<" fracz= "<<fracz<<endl;
				//cout <<"dwf= "<<dwf<<" f= "<<xray.f[atomID][k]<<endl;
				//cout <<"real= "<<xray.real[k]<<endl;
				/*			
							if (isnan(real) || isnan(imag))
							{
							cout <<"real= "<<xray.real[k]<<endl;
							cout <<"imag= "<<xray.imag[k]<<endl;
							cout <<"atomID= "<<atomID<<endl;
							cout <<"f= "<<xray.f[atomID][k]<<endl;
							cout <<"dotproduct= "<<dotproduct<<endl;
							printAtomInfo(Atoms[j]);
							cout <<"fracx= "<<fracx<<endl;
							cout <<"fracy= "<<fracy<<endl;
							cout <<"fracz= "<<fracz<<endl;
							cout <<"miller.x= "<<xray.miller[i].Pos[X]<<endl;
							cout <<"miller.y= "<<xray.miller[i].Pos[Y]<<endl;
							cout <<"miller.z= "<<xray.miller[i].Pos[Z]<<endl;
							cout <<"dwf= "<<dwf<<endl;
							cout <<"cosdp= "<<cosdp<<endl;
							cout <<"sindp= "<<sindp<<endl;
							cout <<"occupancy= "<<occupancy<<endl;
							error("Intensity is nan", __LINE__, __FILE__);
							}
				 */
			}
		}
	}
	gettimeofday(&fullEnd, NULL);
	cout <<"totalTime= "<<calcTimeDiff(fullStart, fullEnd)<<endl;
	calcIntensity(xray);
	calcAmplitude(xray);
	calcPhase(xray);
}

void calcRealImag(ProteinStruct &Protein, XRayStruct &xray, Matrix &cosx, Matrix &cosy, Matrix &cosz, Matrix &sinx, Matrix &siny, Matrix &sinz)
{
	int atomID;
	int h, k, l;
	int natom, npoint;
	Real occupancy, dwf;
	Real real, imag, product;
	//Real cxcy, cxsy, sxcy, sxsy;
	//timeval start, end;

	//gettimeofday(&start, NULL);
	natom=Protein.Atoms.size();
	npoint=xray.miller.size();
	for (int i=0;i<natom;i++)
	{
		atomID=Protein.Atoms[i].atomid;
		occupancy=Protein.Atoms[i].Occupancy;
		for (int j=0;j<npoint;j++)
		{
			real=imag=0;
			h=int(xray.miller[j].Pos[X]-xray.minH);
			k=int(xray.miller[j].Pos[Y]-xray.minK);
			l=int(xray.miller[j].Pos[Z]-xray.minL);

			real+=xray.cosx[i][h]*xray.cosy[i][k]*xray.cosz[i][l];
			real-=xray.cosx[i][h]*xray.siny[i][k]*xray.sinz[i][l];
			real-=xray.sinx[i][h]*xray.cosy[i][k]*xray.sinz[i][l];
			real-=xray.sinx[i][h]*xray.siny[i][k]*xray.cosz[i][l];

			imag-=xray.sinx[i][h]*xray.siny[i][k]*xray.sinz[i][l];
			imag+=xray.sinx[i][h]*xray.cosy[i][k]*xray.cosz[i][l];
			imag+=xray.cosx[i][h]*xray.siny[i][k]*xray.cosz[i][l];
			imag+=xray.cosx[i][h]*xray.cosy[i][k]*xray.sinz[i][l];

			//cout <<"alpha= "<<Protein.alpha<<" beta= "<<Protein.beta<<" gamma= "<<Protein.gamma<<endl;
			if (Protein.alpha==90.0 && Protein.beta==90.0 && Protein.gamma==90.0)
			{
				dwf=xray.dwfx[i][h]*xray.dwfy[i][k]*xray.dwfz[i][l];
				//cout <<"dwf= "<<dwf<<" dwfReal= "<<calcDWF(xray.q[j], Protein.Atoms[i].BFactor)<<endl;
			}
			else
			{
				dwf=calcDWF(xray.q[j], Protein.Atoms[i].BFactor, xray.lookUp);
			}
			product=xray.f[atomID][j]*dwf*occupancy;
			//if (j==0)
			//{
			//	cout <<"real= "<<real<<" imag= "<<imag<<" atomID= "<<atomID<<" f= "<<xray.f[atomID][j]<<" dwf= "<<dwf<<" occupancy= "<<occupancy<<endl;
			//	printAtomInfo(Protein.Atoms[i]);
			//}
			//cout <<"j= "<<j<<" real= "<<real<<" imag= "<<imag<<" atomID= "<<atomID<<" f= "<<xray.f[atomID][j]<<" dwf= "<<dwf<<" occupancy= "<<occupancy<<endl;
			xray.real[j]+=product*real;
			xray.imag[j]+=product*imag;
			/*
			   if (isnan(xray.real[j]) || isnan(xray.imag[j]))
			   {
			   cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<endl;
			   cout <<"f= "<<xray.f[atomID][j]<<endl;
			   cout <<"dwfx= "<<xray.dwfx[i][h]<<endl;
			   cout <<"dwfy= "<<xray.dwfy[i][k]<<endl;
			   cout <<"dwfz= "<<xray.dwfz[i][l]<<endl;
			   cout <<"occupancy= "<<occupancy<<endl;
			   cout <<"product= "<<product<<endl;
			   cout <<"real= "<<real<<endl;
			   cout <<"imag= "<<imag<<endl;
			   endProgram(__LINE__, __FILE__);
			   }
			 */
		}
	}
}

void calcTrigForScattering(ProteinStruct &Protein, XRayStruct &xray)
{
	int natom;
	int npoint=xray.q.size();
	int nreal=xray.real.size();
	int nf1, nf2;
	int Size1, Size2;
	Real fracx, fracy, fracz;
	Real fracxj, fracyj, fraczj;
	timeval start, end;


	//printAtomInfo(Protein.Atoms[0]);
	Get2DVectorSize(xray.f, nf1, nf2);
	if (npoint!=nreal || npoint!=nf2 || npoint==0)
	{
		string errorStr="npoint!=nreal.";
		errorStr+=" npoint= "+toStr(npoint);
		errorStr+=" nreal= "+toStr(nreal);
		errorStr+=" nf= "+toStr(nf2);
		error(errorStr, __LINE__, __FILE__);
	}
	setMinMiller(xray);
	natom=Protein.Atoms.size();
	gettimeofday(&start, NULL);
	Get2DVectorSize(xray.cosx, Size1, Size2);
	if (Size1!=natom || Size2!=xray.maxH-xray.minH+1)
	{
		Safe2DAlloc(xray.cosx, natom, xray.maxH-xray.minH+1, "cosx");
		Safe2DAlloc(xray.cosy, natom, xray.maxK-xray.minK+1, "cosy");
		Safe2DAlloc(xray.cosz, natom, xray.maxL-xray.minL+1, "cosz");
		Safe2DAlloc(xray.sinx, natom, xray.maxH-xray.minH+1, "sinx");
		Safe2DAlloc(xray.siny, natom, xray.maxK-xray.minK+1, "siny");
		Safe2DAlloc(xray.sinz, natom, xray.maxL-xray.minL+1, "sinz");
		Safe2DAlloc(xray.dwfx, natom, xray.maxH-xray.minH+1, "dwfx");
		Safe2DAlloc(xray.dwfy, natom, xray.maxK-xray.minK+1, "dwfy");
		Safe2DAlloc(xray.dwfz, natom, xray.maxL-xray.minL+1, "dwfz");

	}
	gettimeofday(&end, NULL);
	zeroVector(xray.i);
	zeroVector(xray.real);
	zeroVector(xray.imag);
	zeroVector(xray.waterReal);
	zeroVector(xray.waterImag);
	calcCartToFrac(Protein);
	gettimeofday(&start, NULL);
	for (int i=0;i<natom;i++)
	{
		convertToFractionalCoordinates(Protein.Atoms[i].x, Protein.Atoms[i].y, Protein.Atoms[i].z, Protein.cartToFrac, fracx, fracy, fracz);
		fracx*=2.0*pi;
		fracy*=2.0*pi;
		fracz*=2.0*pi;
		for (int j=xray.minH;j<=xray.maxH;j++)
		{
			fracxj=fracx*Real(j);
			getTrigLookUp(xray.lookUp, fracxj, xray.cosx[i][j-xray.minH], xray.sinx[i][j-xray.minH]);
			xray.dwfx[i][j-xray.minH]=calcDWF(2.0*pi*Real(j)/Protein.XBoxLength, Protein.Atoms[i].BFactor, xray.lookUp);
			//cout <<"fracxj= "<<fracxj<<" cosx= "<<xray.cosx[i][j-xray.minH]<<" sinx= "<<xray.sinx[i][j-xray.minH]<<endl;
			//cout <<"dwfx= "<<xray.dwfx[i][j-xray.minH]<<endl;
		}
		for (int j=xray.minK;j<=xray.maxK;j++)
		{
			fracyj=fracy*Real(j);
			getTrigLookUp(xray.lookUp, fracyj, xray.cosy[i][j-xray.minK], xray.siny[i][j-xray.minK]);
			xray.dwfy[i][j-xray.minK]=calcDWF(2.0*pi*Real(j)/Protein.YBoxLength, Protein.Atoms[i].BFactor, xray.lookUp);
		}
		for (int j=xray.minL;j<=xray.maxL;j++)
		{
			fraczj=fracz*Real(j);
			getTrigLookUp(xray.lookUp, fraczj, xray.cosz[i][j-xray.minL], xray.sinz[i][j-xray.minL]);
			xray.dwfz[i][j-xray.minL]=calcDWF(2.0*pi*Real(j)/Protein.ZBoxLength, Protein.Atoms[i].BFactor, xray.lookUp);
		}
	}
}

void calcScatteringFast(ProteinStruct &Protein, XRayStruct &xray)
{
	calcTrigForScattering(Protein, xray);
	calcRealImag(Protein, xray, xray.cosx, xray.cosy, xray.cosz, xray.sinx, xray.siny, xray.sinz);
}

Real estimateXRayScore(ProteinStruct &referenceProtein, ProteinStruct &Protein, Matrix &derivative, Real oldScore)
{
	//Estimates the X-ray score of Protein based on the X-ray score of 
	//referenceProtein and the gradient of the X-ray score with respect
	//to the cartesian coordinates.
	int natom1=referenceProtein.Atoms.size();
	int natom2=Protein.Atoms.size();
	Real dx, dy, dz;
	Real dScore=0;

	if (natom1!=natom2 || natom1==0)
	{
		string errorStr="natom1= "+toStr(natom1);
		errorStr+=" natom2= "+toStr(natom2);
		error(errorStr, __LINE__, __FILE__);
	}

	for (int i=0;i<natom1;i++)
	{
		dx=Protein.Atoms[i].x-referenceProtein.Atoms[i].x;
		dy=Protein.Atoms[i].y-referenceProtein.Atoms[i].y;
		dz=Protein.Atoms[i].z-referenceProtein.Atoms[i].z;
		dScore+=dx*derivative[i][X];
		dScore+=dy*derivative[i][Y];
		dScore+=dz*derivative[i][Z];
	}
	return oldScore+dScore;
}

Real estimateXRayScore(Vector dRot, EstimateScoreStruct args)
{
	rotateAtoms(args.Protein.Atoms, 0, 0, dRot[Z]);
	rotateAtoms(args.Protein.Atoms, 0, dRot[Y], 0);
	rotateAtoms(args.Protein.Atoms, dRot[X], 0, 0);
	return estimateXRayScore(args.referenceProtein, args.Protein, args.cartGrad, 0);
}

Real calcRotationalDerivative(vector<AtomStruct> &Atoms, int axis, Matrix &cartesian)
{
	int natom=Atoms.size();
	Real derivative=0;
	Real avex, avey, avez;
	Matrix rotationMatrix;
	VectorStruct v, diff;

	if (axis==X) {v.x=1.0; v.y=0; v.z=0;}
	if (axis==Y) {v.x=0; v.y=1.0; v.z=0;}
	if (axis==Z) {v.x=0; v.y=0; v.z=1.0;}
	Safe2DAlloc(rotationMatrix, 3, 3, "rotationMatrix");
	calcRotationMatrixDerivative(v, rotationMatrix);
	calcCenter(Atoms, avex, avey, avez);
	for (int i=0;i<natom;i++)
	{
		diff.x=Atoms[i].x-avex;
		diff.y=Atoms[i].y-avey;
		diff.z=Atoms[i].z-avez;
		applyRotationMatrix(diff.x, diff.y, diff.z, rotationMatrix);
		derivative+=diff.x*cartesian[i][X];
		derivative+=diff.y*cartesian[i][Y];
		derivative+=diff.z*cartesian[i][Z];
	}
	return derivative;
}

Real calcTranslationalDerivative(int axis, Matrix &cartesian)
{
	int Size1, Size2;
	Real derivative=0;

	Get2DVectorSize(cartesian, Size1, Size2);

	if (Size2!=3)
	{
		string errorStr="Size2= "+toStr(Size2)+", but should be 3";
		error(errorStr, __LINE__, __FILE__);
	}

	for (int i=0;i<Size1;i++)
	{
		derivative+=cartesian[i][axis];
	}
	return derivative;
}

Real calcDerivative(vector<AtomStruct> &Atoms, int angleType, int residue, Matrix &cartesian)
{
	int natom=Atoms.size();
	int atom1, atom2, atom3, atom4;
	Real derivative=0;
	Matrix rotationMatrix;
	VectorStruct v, diff;

	getDihedralAtoms(Atoms, atom1, atom2, atom3, atom4, angleType, residue);
	v=atomVector(Atoms[atom3], Atoms[atom2]);
	vectorNormalize(v);
	Safe2DAlloc(rotationMatrix, 3, 3, "rotationMatrix");
	calcRotationMatrixDerivative(v, rotationMatrix);
	for (int i=atom4;i<natom;i++)
	{
		diff=atomVector(Atoms[i], Atoms[atom3]);
		applyRotationMatrix(diff.x, diff.y, diff.z, rotationMatrix);
		derivative+=diff.x*cartesian[i][X];
		derivative+=diff.y*cartesian[i][Y];
		derivative+=diff.z*cartesian[i][Z];
	}
	return derivative;
}

void calcAnalyticalRotationalGrad(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &derivative)
{
	Real derx, dery, derz;
	Matrix gradient;

	calcUnitCellVectors(Protein, params);
	calcScatteringFast(Protein, xray);
	calcScatteringGradient(Protein, xray, expXRay, gradient);
	derx=calcRotationalDerivative(Protein.Atoms, X, gradient);
	dery=calcRotationalDerivative(Protein.Atoms, Y, gradient);
	derz=calcRotationalDerivative(Protein.Atoms, Z, gradient);
	SafePushBack(derivative, derx, "derivative");
	SafePushBack(derivative, dery, "derivative");
	SafePushBack(derivative, derz, "derivative");
}

void calcFiniteDifferenceRotationDer(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, int axis, Real &derLower, Real &derUpper)
{
	Real dr=1e-6;
	Real rfactor1, rfactor2, rfactor3;
	ProteinStruct tempProtein;
	vector<ProteinStruct> Proteins, tempProteins;

	SafePushBack(Proteins, Protein, "Proteins");
	SafePushBack(tempProteins, Protein, "tempProteins");
	rfactor1=calcScattering(Proteins, xray, expXRay, params);
	if (axis==X) rotateAtoms(tempProteins[0].Atoms, dr, 0, 0);
	if (axis==Y) rotateAtoms(tempProteins[0].Atoms, 0, dr, 0);
	if (axis==Z) rotateAtoms(tempProteins[0].Atoms, 0, 0, dr);
	rfactor2=calcScattering(tempProteins, xray, expXRay, params);
	derUpper=(rfactor2-rfactor1)/dr;
	tempProteins[0]=Protein;
	if (axis==X) rotateAtoms(tempProteins[0].Atoms, -dr, 0, 0);
	if (axis==Y) rotateAtoms(tempProteins[0].Atoms, 0, -dr, 0);
	if (axis==Z) rotateAtoms(tempProteins[0].Atoms, 0, 0, -dr);
	rfactor3=calcScattering(tempProteins, xray, expXRay, params);
	derLower=(rfactor1-rfactor3)/dr;
}

/*
   void calcFiniteDifferenceTranslationDer(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, int axis, Real &derLower, Real &derUpper)
   {
   Real dx=1e-6;
   Real rfactor1, rfactor2, rfactor3;
   ProteinStruct tempProtein;

   tempProtein=Protein;
   rfactor1=calcScattering(Protein, xray, expXRay, params);
   cout <<"rfactor1= "<<rfactor1<<endl;
   if (axis==X) moveAtoms(tempProtein.Atoms, dx, 0, 0);
   if (axis==Y) moveAtoms(tempProtein.Atoms, 0, dx, 0);
   if (axis==Z) moveAtoms(tempProtein.Atoms, 0, 0, dx);
   rfactor2=calcScattering(tempProtein, xray, expXRay, params);
   derUpper=(rfactor2-rfactor1)/dx;
   tempProtein=Protein;
   if (axis==X) moveAtoms(tempProtein.Atoms, -dx, 0, 0);
   if (axis==Y) moveAtoms(tempProtein.Atoms, 0, -dx, 0);
   if (axis==Z) moveAtoms(tempProtein.Atoms, 0, 0, -dx);
   rfactor3=calcScattering(tempProtein, xray, expXRay, params);
   derLower=(rfactor1-rfactor3)/dx;
   }

   void calcFiniteDifferenceRotationTranslationGrad(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &derivativeUpper, Vector &derivativeLower)
   {
   Real derrxUpper, derryUpper, derrzUpper;
   Real derrxLower, derryLower, derrzLower;
   Real derdxUpper, derdyUpper, derdzUpper;
   Real derdxLower, derdyLower, derdzLower;

   derivativeUpper.clear();
   derivativeLower.clear();

   calcFiniteDifferenceRotationDer(Protein, xray, expXRay, params, X, derrxUpper, derrxLower);
   calcFiniteDifferenceRotationDer(Protein, xray, expXRay, params, Y, derryUpper, derryLower);
   calcFiniteDifferenceRotationDer(Protein, xray, expXRay, params, Z, derrzUpper, derrzLower);
   calcFiniteDifferenceTranslationDer(Protein, xray, expXRay, params, X, derdxUpper, derdxLower);
   calcFiniteDifferenceTranslationDer(Protein, xray, expXRay, params, Y, derdyUpper, derdyLower);
   calcFiniteDifferenceTranslationDer(Protein, xray, expXRay, params, Z, derdzUpper, derdzLower);

   SafePushBack(derivativeUpper, derrxUpper, "derivativeUpper");
   SafePushBack(derivativeUpper, derryUpper, "derivativeUpper");
   SafePushBack(derivativeUpper, derrzUpper, "derivativeUpper");
   SafePushBack(derivativeUpper, derdxUpper, "derivativeUpper");
   SafePushBack(derivativeUpper, derdyUpper, "derivativeUpper");
   SafePushBack(derivativeUpper, derdzUpper, "derivativeUpper");

   SafePushBack(derivativeLower, derrxLower, "derivativeLower");
   SafePushBack(derivativeLower, derryLower, "derivativeLower");
   SafePushBack(derivativeLower, derrzLower, "derivativeLower");
   SafePushBack(derivativeLower, derdxLower, "derivativeLower");
   SafePushBack(derivativeLower, derdyLower, "derivativeLower");
   SafePushBack(derivativeLower, derdzLower, "derivativeLower");
   }
 */
void calcRotationalDerivativeX(Real phi, Real theta, Real psi, Matrix &rotMat)
{
	rotMat[X][X]=0;
	rotMat[X][Y]=0;
	rotMat[X][Z]=0;

	rotMat[Y][X]=cos(phi)*sin(theta)*cos(psi)-sin(phi)*sin(psi);
	rotMat[Y][Y]=-cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
	rotMat[Y][Z]=-cos(phi)*cos(theta);

	rotMat[Z][X]=sin(phi)*sin(theta)*cos(psi)+cos(phi)*sin(psi);
	rotMat[Z][Y]=-sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi);
	rotMat[Z][Z]=-sin(phi)*cos(theta);
}

void calcRotationalDerivativeY(Real phi, Real theta, Real psi, Matrix &rotMat)
{
	rotMat[X][X]=-sin(theta)*cos(psi);
	rotMat[X][Y]=sin(theta)*sin(psi);
	rotMat[X][Z]=cos(theta);

	rotMat[Y][X]=sin(phi)*cos(theta)*cos(psi);
	rotMat[Y][Y]=-sin(phi)*cos(theta)*sin(psi);
	rotMat[Y][Z]=sin(phi)*sin(theta);

	rotMat[Z][X]=-cos(phi)*cos(theta)*cos(psi);
	rotMat[Z][Y]=cos(phi)*cos(theta)*sin(psi);
	rotMat[Z][Z]=-cos(phi)*sin(theta);
}

void calcRotationalDerivativeZ(Real phi, Real theta, Real psi, Matrix &rotMat)
{
	rotMat[X][X]=-cos(theta)*sin(psi);
	rotMat[X][Y]=-cos(theta)*cos(psi);
	rotMat[X][Z]=0;

	rotMat[Y][X]=-sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi);
	rotMat[Y][Y]=-sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi);
	rotMat[Y][Z]=0;

	rotMat[Z][X]=cos(phi)*sin(theta)*sin(psi)+sin(phi)*cos(psi);
	rotMat[Z][Y]=cos(phi)*sin(theta)*cos(psi)-sin(phi)*sin(psi);
	rotMat[Z][Z]=0;
}

void calcRotationalDerivative(ProteinStruct &Protein, Real phi, Real theta, Real psi, Matrix &gradient, Real &derrx, Real &derry, Real &derrz)
{
	int natom=Protein.Atoms.size();
	Real avex, avey, avez;
	Matrix rotMatX, rotMatY, rotMatZ;
	VectorStruct diff;
	VectorStruct vx, vy, vz;	

	derrx=0;
	derry=0;
	derrz=0;

	Safe2DAlloc(rotMatX, 3, 3, "rotMatX");
	Safe2DAlloc(rotMatY, 3, 3, "rotMatY");
	Safe2DAlloc(rotMatZ, 3, 3, "rotMatZ");

	calcRotationalDerivativeX(phi, theta, psi, rotMatX);
	calcRotationalDerivativeY(phi, theta, psi, rotMatY);
	calcRotationalDerivativeZ(phi, theta, psi, rotMatZ);

	calcCenter(Protein.Atoms, avex, avey, avez);
	for (int i=0;i<natom;i++)
	{
		diff.x=Protein.Atoms[i].x-avex;
		diff.y=Protein.Atoms[i].y-avey;
		diff.z=Protein.Atoms[i].z-avez;
		vx=diff;
		vy=diff;
		vz=diff;
		applyRotationMatrix(vx.x, vx.y, vx.z, rotMatX);
		applyRotationMatrix(vy.x, vy.y, vy.z, rotMatY);
		applyRotationMatrix(vz.x, vz.y, vz.z, rotMatZ);

		derrx+=vx.x*gradient[i][X];
		derrx+=vx.y*gradient[i][Y];
		derrx+=vx.z*gradient[i][Z];

		derry+=vy.x*gradient[i][X];
		derry+=vy.y*gradient[i][Y];
		derry+=vy.z*gradient[i][Z];

		derrz+=vz.x*gradient[i][X];
		derrz+=vz.y*gradient[i][Y];
		derrz+=vz.z*gradient[i][Z];
	}
}

void checkRotationalGrad(ProteinStruct Protein, Real phi, Real psi, Real theta)
{
	Real d=1e-5;
	Real dx, dy, dz;
	Real avex, avey, avez;
	Matrix rotMat, rotMat2, rotMatXY;
	VectorStruct diff;
	ProteinStruct originalProtein, tempProtein;

	originalProtein=Protein;
	tempProtein=Protein;
	Safe2DAlloc(rotMat, 3, 3, "rotMat");
	Safe2DAlloc(rotMat2, 3, 3, "rotMat2");
	Safe2DAlloc(rotMatXY, 3, 3, "rotMatXY");
	calcCenter(Protein.Atoms, avex, avey, avez);

	rotMat2[X][X]=cos(theta)*cos(psi);
	rotMat2[X][Y]=-cos(theta)*sin(psi);
	rotMat2[X][Z]=sin(theta);

	rotMat2[Y][X]=sin(phi)*sin(theta)*cos(psi)+cos(phi)*sin(psi);
	rotMat2[Y][Y]=-sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi);
	rotMat2[Y][Z]=-sin(phi)*cos(theta);

	rotMat2[Z][X]=-cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
	rotMat2[Z][Y]=cos(phi)*sin(theta)*sin(psi)+sin(phi)*cos(psi);
	rotMat2[Z][Z]=cos(phi)*cos(theta);

	rotMatXY[X][X]=cos(theta);
	rotMatXY[X][Y]=0;
	rotMatXY[X][Z]=sin(theta);

	rotMatXY[Y][X]=sin(phi)*sin(theta);
	rotMatXY[Y][Y]=cos(phi);
	rotMatXY[Y][Z]=-sin(phi)*cos(theta);

	rotMatXY[Z][X]=-cos(phi)*sin(theta);
	rotMatXY[Z][Y]=sin(phi);
	rotMatXY[Z][Z]=cos(phi)*cos(theta);

	rotateAtomsNoCenter(Protein.Atoms, 0, 0, psi);
	rotateAtomsNoCenter(Protein.Atoms, 0, theta, 0);
	rotateAtomsNoCenter(Protein.Atoms, phi, 0, 0);

	printAtomInfo(Protein.Atoms[0]);

	diff.x=tempProtein.Atoms[0].x;
	diff.y=tempProtein.Atoms[0].y;
	diff.z=tempProtein.Atoms[0].z;
	applyRotationMatrix(diff.x, diff.y, diff.z, rotMat2);

	cout <<vectorToStr(diff)<<endl;
	endProgram(__LINE__, __FILE__);

	rotateAtoms(tempProtein.Atoms, 0, 0, theta);
	rotateAtoms(tempProtein.Atoms, 0, psi, 0);
	rotateAtoms(tempProtein.Atoms, phi+d, 0, 0);

	dx=(tempProtein.Atoms[0].x-Protein.Atoms[0].x)/d;
	dy=(tempProtein.Atoms[0].y-Protein.Atoms[0].y)/d;
	dz=(tempProtein.Atoms[0].z-Protein.Atoms[0].z)/d;

	diff.x=originalProtein.Atoms[0].x-avex;
	diff.y=originalProtein.Atoms[0].y-avex;
	diff.z=originalProtein.Atoms[0].z-avex;

	calcRotationalDerivativeX(phi, theta, psi, rotMat);
	applyRotationMatrix(diff.x, diff.y, diff.z, rotMat);

	cout <<"diff.x= "<<diff.x<<" dx= "<<dx<<endl;
	cout <<"diff.y= "<<diff.y<<" dy= "<<dy<<endl;
	cout <<"diff.z= "<<diff.z<<" dz= "<<dz<<endl;
	endProgram(__LINE__, __FILE__);
}

void calcAnalyticalRotationTranslationGrad(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom, Vector &derivative)
{
	Real derrx, derry, derrz;
	Real derdx, derdy, derdz;
	Vector derivativeUpper, derivativeLower;
	Matrix gradient;
	ProteinStruct tempProtein;
	timeval start, end;

	gettimeofday(&start, NULL);
	derivative.clear();
	tempProtein=Proteins[0];

	placeProteins(Proteins, degOfFreedom);
	calcUnitCellVectors(Proteins[0], params);
	calcScatteringBfactor(Proteins[0], xray);
	calcScatteringGradient(Proteins[0], xray, expXRay, gradient);
	calcRotationalDerivative(tempProtein, degOfFreedom[X], degOfFreedom[Y], degOfFreedom[Z], gradient, derrx, derry, derrz);
	derdx=calcTranslationalDerivative(X, gradient);
	derdy=calcTranslationalDerivative(Y, gradient);
	derdz=calcTranslationalDerivative(Z, gradient);

	SafePushBack(derivative, derrx, "derivative");
	SafePushBack(derivative, derry, "derivative");
	SafePushBack(derivative, derrz, "derivative");
	SafePushBack(derivative, derdx, "derivative");
	SafePushBack(derivative, derdy, "derivative");
	SafePushBack(derivative, derdz, "derivative");

	scaleVector(derivative, 0.01);
	gettimeofday(&end, NULL);

	cout <<"calcAnalyticalRotationTranslationGrad took "<<calcTimeDiff(start, end)<<endl;
}

Real calcAnalyticalRotationTranslationGrad(Vector degOfFreedom, MRStruct gradArgs, Vector &derivative)
{
	calcAnalyticalRotationTranslationGrad(gradArgs.Proteins, gradArgs.xray, gradArgs.expXRay, gradArgs.params, degOfFreedom, derivative);
	return 0.001;
}

void cartesianToDihedralGradient(ProteinStruct &Protein, Matrix &derivative, Vector &dihedralDer)
{
	int firstResidue, lastResidue;
	Real der;

	getFirstLastResidue(Protein.Atoms, firstResidue, lastResidue);
	dihedralDer.clear();
	for (int i=firstResidue;i<=lastResidue;i++)
	{
		if (i!=firstResidue)
		{
			der=calcDerivative(Protein.Atoms, PHI, i, derivative);
			SafePushBack(dihedralDer, der, "dihedralDer");
		}
		if (i!=lastResidue)
		{
			der=calcDerivative(Protein.Atoms, PSI, i, derivative);
			SafePushBack(dihedralDer, der, "dihedralDer");
		}
		if (i!=lastResidue)
		{
			der=calcDerivative(Protein.Atoms, OMEGA, i, derivative);
			SafePushBack(dihedralDer, der, "dihedralDer");
		}
	}
}

void calcScatteringGradientP1(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, Matrix &derivative)
{
	//Valid only for space group P1.
	int atomID, natom;
	int npoint=xray.q.size();
	int nreal=xray.real.size();
	Real dotproduct, dwf;
	Real occupancy;
	Real scale, product;
	Real imagDer, realDer;
	Real cosdp, sindp;
	Vector dRealdS, dImagdS; //The derivative of the xray score with respect to the real and imaginary components
	timeval fullStart, fullEnd;

	if (npoint!=nreal)
	{
		string errorStr="npoint!=nreal.";
		errorStr+="npoint= "+IntToStr(npoint);
		errorStr+="nreal= "+IntToStr(nreal);
		error(errorStr, __LINE__, __FILE__);
	}
	gettimeofday(&fullStart, NULL);
	calcIntensity(xray);
	convertToFractionalCoordinates(Protein);
	natom=Protein.Atoms.size();
	SafeAlloc(dRealdS, npoint, "dRealdS");
	SafeAlloc(dImagdS, npoint, "dImagdS");
	Safe2DAlloc(derivative, natom, 3, "derivative");
	scale=calcExperimentalScaleFactor(xray, expXRay);
	for (int i=0;i<npoint;i++)
	{
		dRealdS[i]=4.0*xray.real[i]*scale*(scale*xray.i[i]-expXRay.i[i]);
		dImagdS[i]=4.0*xray.imag[i]*scale*(scale*xray.i[i]-expXRay.i[i]);
	}
	for (int j=0;j<natom;j++)
	{
		atomID=Protein.Atoms[j].atomid;
		occupancy=Protein.Atoms[j].Occupancy;
		for (int i=0;i<npoint;i++)
		{
			dotproduct=Protein.Atoms[j].frac.x*xray.miller[i].Pos[X];
			dotproduct+=Protein.Atoms[j].frac.y*xray.miller[i].Pos[Y];
			dotproduct+=Protein.Atoms[j].frac.z*xray.miller[i].Pos[Z];
			dotproduct*=2.0*pi;
			dwf=calcDWF(xray, i, Protein.Atoms[j], xray.lookUp);
			getTrigLookUp(xray.lookUp, dotproduct, cosdp, sindp);
			for (int k=0;k<3;k++)
			{
				product=2.0*pi*xray.f[atomID][i]*dwf*xray.miller[i].Pos[k]*occupancy;
				realDer=-product*sindp;
				imagDer=product*cosdp;
				derivative[j][k]+=realDer*dRealdS[i];
				derivative[j][k]+=imagDer*dImagdS[i];
			}
			/*
			   if (isnan(real) || isnan(imag))
			   {
			   cout <<"real= "<<real<<endl;
			   cout <<"imag= "<<imag<<endl;
			   cout <<"atomID= "<<atomID<<endl;
			   cout <<"dotproduct= "<<dotproduct<<endl;
			   printAtomInfo(UnitCell.Atoms[j]);
			   cout <<"frac.x= "<<UnitCell.Atoms[j].frac.x<<endl;
			   cout <<"frac.y= "<<UnitCell.Atoms[j].frac.y<<endl;
			   cout <<"frac.z= "<<UnitCell.Atoms[j].frac.z<<endl;
			   cout <<"miller.x= "<<xray.miller[i].Pos[X]<<endl;
			   cout <<"miller.y= "<<xray.miller[i].Pos[Y]<<endl;
			   cout <<"miller.z= "<<xray.miller[i].Pos[Z]<<endl;
			   cout <<"dwf= "<<dwf<<endl;
			   cout <<"cosdp= "<<cosdp<<endl;
			   cout <<"sindp= "<<sindp<<endl;
			   error("Intensity is nan", __LINE__, __FILE__);
			   }
			 */
		}
	}
	for (int j=0;j<natom;j++)
	{
		derivative[j][X]/=Protein.XBoxLength;
		derivative[j][Y]/=Protein.YBoxLength;
		derivative[j][Z]/=Protein.ZBoxLength;
	}
	gettimeofday(&fullEnd, NULL);
	//cout <<"totalTime= "<<calcTimeDiff(fullStart, fullEnd)<<endl;
}

void calcScatteringGradient(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, Matrix &derivative)
{
	int atomID, natom;
	int npoint=xray.q.size();
	int nreal=xray.real.size();
	int noperations=Protein.symmetryOperations.size();
	int h, k, l;
	Real dwf;
	Real occupancy;
	Real scale, product;
	Real sum=0, sumSqr=0, crossTerm=0, term=0;
	Real dRealdScale, dImagdScale;
	Real real, imag;
	Vector dRealdS, dImagdS; //The derivative of the xray score with respect to the real and imaginary components
	Vector realDer, imagDer;
	Vector tempDerivative, tempDerivative2;
	Matrix fracToSym;
	AtomStruct tempAtom;
	ProteinStruct tempProtein;
	timeval fullStart, fullEnd;

	if (npoint!=nreal)
	{
		string errorStr="npoint!=nreal.";
		errorStr+="npoint= "+IntToStr(npoint);
		errorStr+="nreal= "+IntToStr(nreal);
		error(errorStr, __LINE__, __FILE__);
	}
	gettimeofday(&fullStart, NULL);
	calcIntensity(xray);
	calcCartToFrac(Protein);
	natom=Protein.Atoms.size();
	SafeAlloc(dRealdS, npoint, "dRealdS");
	SafeAlloc(dImagdS, npoint, "dImagdS");
	SafeAlloc(realDer, 3, "realDer");
	SafeAlloc(imagDer, 3, "imagDer");
	SafeAlloc(tempDerivative, 3, "tempDerivative");
	SafeAlloc(tempDerivative2, 3, "tempDerivative2");
	Safe2DAlloc(derivative, natom, 3, "derivative");
	Safe2DAlloc(fracToSym, 3, 3, "fracToSym");
	scale=calcExperimentalScaleFactor(xray, expXRay);
	cout <<"scale= "<<scale<<endl;
	for (int i=0;i<npoint;i++)
	{
		sum+=xray.i[i];
		sumSqr+=xray.i[i]*xray.i[i];
		crossTerm+=xray.i[i]*expXRay.i[i];
		term+=2.0*(scale*xray.i[i]-expXRay.i[i])*xray.i[i];
	}
	for (int i=0;i<npoint;i++)
	{

		dRealdScale=2.0*xray.real[i]*expXRay.i[i]*sumSqr;
		dRealdScale-=4.0*xray.real[i]*xray.i[i]*crossTerm;
		dRealdScale/=(sumSqr*sumSqr);
		dImagdScale=2.0*xray.imag[i]*expXRay.i[i]*sumSqr;
		dImagdScale-=4.0*xray.imag[i]*xray.i[i]*crossTerm;
		dImagdScale/=(sumSqr*sumSqr);
		dRealdS[i]=4.0*(scale*xray.i[i]-expXRay.i[i]);
		dRealdS[i]*=xray.real[i]*scale;
		dImagdS[i]=4.0*(scale*xray.i[i]-expXRay.i[i]);
		dImagdS[i]*=xray.imag[i]*scale;
		dRealdS[i]+=term*dRealdScale;
		dImagdS[i]+=term*dImagdScale;
	}
	for (int i=0;i<noperations;i++)
	{
		tempProtein=Protein;
		applySymmetryOperation(tempProtein.Atoms, Protein.symmetryOperations[i]);
		calcTrigForScattering(tempProtein, xray);
		for (int j=0;j<natom;j++)
		{
			atomID=Protein.Atoms[j].atomid;
			occupancy=Protein.Atoms[j].Occupancy;
			matrixMultiply(Protein.cartToFrac, Protein.symmetryOperations[i].rotationMat, fracToSym);
			zeroVector(tempDerivative);
			for (int n=0;n<npoint;n++)
			{
				real=imag=0;
				h=int(xray.miller[n].Pos[X]-xray.minH);
				k=int(xray.miller[n].Pos[Y]-xray.minK);
				l=int(xray.miller[n].Pos[Z]-xray.minL);

				real+=xray.cosx[j][h]*xray.cosy[j][k]*xray.cosz[j][l];
				real-=xray.cosx[j][h]*xray.siny[j][k]*xray.sinz[j][l];
				real-=xray.sinx[j][h]*xray.cosy[j][k]*xray.sinz[j][l];
				real-=xray.sinx[j][h]*xray.siny[j][k]*xray.cosz[j][l];

				imag-=xray.sinx[j][h]*xray.siny[j][k]*xray.sinz[j][l];
				imag+=xray.sinx[j][h]*xray.cosy[j][k]*xray.cosz[j][l];
				imag+=xray.cosx[j][h]*xray.siny[j][k]*xray.cosz[j][l];
				imag+=xray.cosx[j][h]*xray.cosy[j][k]*xray.sinz[j][l];

				if (Protein.alpha==90.0 && Protein.beta==90.0 && Protein.gamma==90.0)
				{
					dwf=xray.dwfx[j][h]*xray.dwfy[j][k]*xray.dwfz[j][l];
				}
				else
				{
					dwf=calcDWF(xray.q[n], Protein.Atoms[j].BFactor, xray.lookUp);
				}
				product=2.0*pi*xray.f[atomID][n]*dwf*occupancy;
				for (int m=0;m<3;m++)
				{
					realDer[m]=-product*xray.miller[k].Pos[m]*imag;
					imagDer[m]=product*xray.miller[k].Pos[m]*real;
				}
				for (int m=0;m<3;m++)
				{
					tempDerivative[m]+=realDer[m]*dRealdS[k];
					tempDerivative[m]+=imagDer[m]*dImagdS[k];	// dScore/dx=sigma(dScore/dF*dF/dx)
				}
			}
			matrixMultiply(tempDerivative, fracToSym, tempDerivative2);
			for (int m=0;m<3;m++)
			{
				derivative[j][m]+=tempDerivative2[m];
			} 
		}
	}
	gettimeofday(&fullEnd, NULL);
	//cout <<"totalTime= "<<calcTimeDiff(fullStart, fullEnd)<<endl;
}

void calcScatteringGradient_old(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, Matrix &derivative)
{
	int atomID, natom;
	int npoint=xray.q.size();
	int nreal=xray.real.size();
	int noperations=Protein.symmetryOperations.size();
	Real dotproduct, dwf;
	Real occupancy;
	Real scale, product;
	Real cosdp, sindp;
	Real fracx, fracy, fracz;
	Real sum=0, sumSqr=0, crossTerm=0, term=0;
	Real dRealdScale, dImagdScale;
	Vector dRealdS, dImagdS; //The derivative of the xray score with respect to the real and imaginary components
	Vector tempRealDer, tempImagDer;
	Vector tempCartRealDer, tempCartImagDer;
	Vector realDer, imagDer;
	Vector real, imag;
	Matrix fracToSym;
	AtomStruct tempAtom;
	timeval fullStart, fullEnd;

	if (npoint!=nreal)
	{
		string errorStr="npoint!=nreal.";
		errorStr+="npoint= "+IntToStr(npoint);
		errorStr+="nreal= "+IntToStr(nreal);
		error(errorStr, __LINE__, __FILE__);
	}
	gettimeofday(&fullStart, NULL);
	calcIntensity(xray);
	calcCartToFrac(Protein);
	natom=Protein.Atoms.size();
	SafeAlloc(dRealdS, npoint, "dRealdS");
	SafeAlloc(dImagdS, npoint, "dImagdS");
	SafeAlloc(tempRealDer, 3, "tempRealDer");
	SafeAlloc(tempImagDer, 3, "tempImagDer");
	SafeAlloc(tempCartRealDer, 3, "tempCartRealDer");
	SafeAlloc(tempCartImagDer, 3, "tempCartImagDer");
	SafeAlloc(realDer, 3, "realDer");
	SafeAlloc(imagDer, 3, "imagDer");
	Safe2DAlloc(derivative, natom, 3, "derivative");
	Safe2DAlloc(fracToSym, 3, 3, "fracToSym");
	scale=calcExperimentalScaleFactor(xray, expXRay);
	cout <<"scale= "<<scale<<endl;
	for (int i=0;i<npoint;i++)
	{
		sum+=xray.i[i];
		sumSqr+=xray.i[i]*xray.i[i];
		crossTerm+=xray.i[i]*expXRay.i[i];
		term+=2.0*(scale*xray.i[i]-expXRay.i[i])*xray.i[i];
	}
	for (int i=0;i<npoint;i++)
	{

		dRealdScale=2.0*xray.real[i]*expXRay.i[i]*sumSqr;
		dRealdScale-=4.0*xray.real[i]*xray.i[i]*crossTerm;
		dRealdScale/=(sumSqr*sumSqr);
		dImagdScale=2.0*xray.imag[i]*expXRay.i[i]*sumSqr;
		dImagdScale-=4.0*xray.imag[i]*xray.i[i]*crossTerm;
		dImagdScale/=(sumSqr*sumSqr);
		dRealdS[i]=4.0*(scale*xray.i[i]-expXRay.i[i]);
		dRealdS[i]*=xray.real[i]*scale;
		dImagdS[i]=4.0*(scale*xray.i[i]-expXRay.i[i]);
		dImagdS[i]*=xray.imag[i]*scale;
		dRealdS[i]+=term*dRealdScale;
		dImagdS[i]+=term*dImagdScale;
	}
	for (int j=0;j<natom;j++)
	{
		atomID=Protein.Atoms[j].atomid;
		occupancy=Protein.Atoms[j].Occupancy;
		for (int i=0;i<noperations;i++)
		{
			tempAtom=applySymmetryOperation(Protein.Atoms[j], Protein.symmetryOperations[i]);
			convertToFractionalCoordinates(tempAtom.x, tempAtom.y, tempAtom.z, Protein.cartToFrac, fracx, fracy, fracz);
			matrixMultiply(Protein.cartToFrac, Protein.symmetryOperations[i].rotationMat, fracToSym);
			fracx*=2.0*pi;
			fracy*=2.0*pi;
			fracz*=2.0*pi;
			for (int k=0;k<npoint;k++)
			{
				dotproduct=fracx*xray.miller[k].Pos[X];
				dotproduct+=fracy*xray.miller[k].Pos[Y];
				dotproduct+=fracz*xray.miller[k].Pos[Z];
				dwf=calcDWF(xray, k, Protein.Atoms[j], xray.lookUp);
				getTrigLookUp(xray.lookUp, dotproduct, cosdp, sindp);
				product=2.0*pi*xray.f[atomID][k]*dwf*occupancy;
				for (int m=0;m<3;m++)
				{
					tempRealDer[m]=-product*xray.miller[k].Pos[m]*sindp;
					tempImagDer[m]=product*xray.miller[k].Pos[m]*cosdp;
				}
				matrixMultiply(tempRealDer, fracToSym, realDer);	//Convert the derivative from fractional to cartesian coords
				matrixMultiply(tempImagDer, fracToSym, imagDer);	//and also take into account the symmetry operation
				for (int m=0;m<3;m++)
				{
					derivative[j][m]+=realDer[m]*dRealdS[k];
					derivative[j][m]+=imagDer[m]*dImagdS[k];	// dScore/dx=sigma(dScore/dF*dF/dx)
				}
			}
		}
	}
	gettimeofday(&fullEnd, NULL);
	//cout <<"totalTime= "<<calcTimeDiff(fullStart, fullEnd)<<endl;
}

void calcScatteringComplex(vector<ProteinStruct> &Proteins, XRayStruct &xray, Vector &degOfFreedom)
{
	int nprot=Proteins.size();
	int nsym=Proteins[0].symmetryOperations.size();
	int npoint=xray.real.size();
	Vector dx;
	/*
	   XRayStruct tempXRay;
	   tempXRay.q=xray.q;
	   tempXRay.qMag=xray.qMag;
	   tempXRay.miller=xray.miller;
	   tempXRay.real=xray.real;
	   tempXRay.imag=xray.imag;
	   tempXRay.f=xray.f;
	   tempXRay.lookUp=xray.lookUp;
	 */
	SafeAlloc(dx, 3, "dx");
	xray.complexSymReal.clear();
	xray.complexSymImag.clear();

	Safe3DAlloc(xray.complexSymReal, nprot, nsym, npoint, "complexSymReal");
	Safe3DAlloc(xray.complexSymImag, nprot, nsym, npoint, "complexSymImag");
	Safe3DAlloc(xray.complexSymTempReal, nprot, nsym, npoint, "complexSymTempReal");
	Safe3DAlloc(xray.complexSymTempImag, nprot, nsym, npoint, "complexSymTempImag");
	for (int i=0;i<nprot;i++)
	{
		//calcScatteringBfactor(Proteins[i], tempXRay);
		calcSymmetryScattering(xray, Proteins[i], i);
	}
	xray.complexSymTempReal=xray.complexSymReal;
	xray.complexSymTempImag=xray.complexSymImag;
	for (int i=0;i<nprot;i++)
	{
		dx[X]=degOfFreedom[i*6+3+X];
		dx[Y]=degOfFreedom[i*6+3+Y];
		dx[Z]=degOfFreedom[i*6+3+Z];
		translateXRay(xray, Proteins[i], dx, i);
		//translateXRay2(xray, Proteins[i], dx, i);
	} 
	calcTotalFormFactor(xray);
	calcIntensity(xray);
	calcPhase(xray);
}

void calcScatteringComplex(vector<ProteinStruct> &Proteins, XRayStruct &xray)
{
	int ndof=Proteins.size()*6;
	Vector degOfFreedom;

	SafeAlloc(degOfFreedom, ndof, "degOfFreedom");
	calcScatteringComplex(Proteins, xray, degOfFreedom);
}

void calcMove(XRayStruct &xray, int nmove)
{
	calcAmplitude(xray.complexReal[nmove], xray.complexImag[nmove], xray.Fmove);
}

void calcFixed(XRayStruct &xray, int nmove)
{
	int nxray=xray.complexReal.size();
	int npoint=xray.complexReal[0].size();
	Vector real, imag;	

	SafeAlloc(real, npoint, "real");
	SafeAlloc(imag, npoint, "imag");
	for (int i=0;i<npoint;i++)
	{
		for (int j=0;j<nxray;j++)
		{
			if (nxray!=nmove)
			{
				real[i]+=xray.complexReal[j][i];
				imag[i]+=xray.complexImag[j][i];
			}
		}
	}
	calcAmplitude(real, imag, xray.Ffix);
}

void calcRotationMatrix(Matrix &rotationMatrix, Real theta, Real phi, Real psi)
{
	Matrix xMat, yMat, zMat, yzMat;

	Safe2DAlloc(xMat, 3, 3, "xMat");
	Safe2DAlloc(yMat, 3, 3, "yMat");
	Safe2DAlloc(zMat, 3, 3, "zMat");
	Safe2DAlloc(yzMat, 3, 3, "yzMat");
	Safe2DAlloc(rotationMatrix, 3, 3, "rotationMatrix");

	calcRotationMatrixAroundX(theta, xMat);
	calcRotationMatrixAroundY(phi, yMat);
	calcRotationMatrixAroundZ(psi, zMat);

	matrixMultiply(yMat, zMat, yzMat);
	matrixMultiply(xMat, yzMat, rotationMatrix);
}

bool isRedundantSpaceGroup(string spaceGroup)
{
	//Returns true if space group has duplicate rotation matrixes.
	if (spaceGroup=="A1") return true;
	if (spaceGroup=="A2") return true;
	if (spaceGroup=="B112") return true;
	if (spaceGroup=="B2") return true;
	if (spaceGroup=="B2212") return true;
	if (spaceGroup=="C121") return true;
	if (spaceGroup=="C1211") return true;
	if (spaceGroup=="C21") return true;
	if (spaceGroup=="C222") return true;
	if (spaceGroup=="C2221") return true;
	if (spaceGroup=="F4132") return true;
	if (spaceGroup=="F222") return true;
	if (spaceGroup=="F23") return true;
	if (spaceGroup=="F432") return true;
	if (spaceGroup=="H3") return true;
	if (spaceGroup=="H32") return true;
	if (spaceGroup=="I121") return true;
	if (spaceGroup=="I1211") return true;
	if (spaceGroup=="I21") return true;
	if (spaceGroup=="I212121") return true;
	if (spaceGroup=="I213") return true;
	if (spaceGroup=="I222") return true;
	if (spaceGroup=="I23") return true;
	if (spaceGroup=="I-42D") return true;
	if (spaceGroup=="I-4C2") return true;
	if (spaceGroup=="I4") return true;
	if (spaceGroup=="I41") return true;
	if (spaceGroup=="I41/A") return true;
	if (spaceGroup=="I411") return true;
	if (spaceGroup=="I4122") return true;
	if (spaceGroup=="I4132") return true;
	if (spaceGroup=="I422") return true;
	if (spaceGroup=="I432") return true;

	return false;
}

bool isRedundantSpaceGroup(int spaceGroup)
{
	//Returns true if space group has duplicate rotation matrixes.

	switch (spaceGroup)
	{
		case A1: return true;
		case A2: return true;
		case B112: return true;
		case B2: return true;
		case B2212: return true;
		case C121: return true;
		case C1211: return true;
		case C21: return true;
		case C222: return true;
		case C2221: return true;
		case F4132: return true;
		case F222: return true;
		case F23: return true;
		case F432: return true;
		case H3: return true;
		case H32: return true;
		case I121: return true;
		case I1211: return true;
		case I21: return true;
		case I212121: return true;
		case I213: return true;
		case I222: return true;
		case I23: return true;
		case I_42D: return true;
		case I_4C2: return true;
		case I4: return true;
		case I41: return true;
		case I41_A: return true;
		case I411: return true;
		case I4122: return true;
		case I422: return true;
		case I432: return true;
		default: return false;
	}
	return false;
}

void calcXRayScatteringChange(ProteinStruct &Protein, XRayStruct &xray, int pick, Vector &degOfFreedom, XRayParamStruct &params)
{
	int prot=pick/6;
	int nprot=params.NumCopies;
	Vector dx, translation;
	Matrix rotationMatrix;
	//timeval start, end;
	SafeAlloc(dx, 3, "dx");
	//if (params.OptimizeOccupancies) sumScatteringContinuousSegments(xray, Protein, degOfFreedom, nprot*6);
	if (params.OptimizeOccupancies) 
	{
		if (pick>=nprot*6)
		{
			calcScatteringContinuousSegmentsChange(xray, Protein, degOfFreedom, nprot*6, pick);
		}
		if (pick<3)
		{
			sumScatteringContinuousSegments(xray, Protein, degOfFreedom, nprot*6);
		}
	}
	if (prot>=nprot)
	{
		if (params.UseFastTranslationRotation)
		{
			//sumScatteringContinuousSegments(xray, degOfFreedom, nprot*6);
			for (int i=0;i<nprot;i++)
			{
				calcRotationMatrix(rotationMatrix, degOfFreedom[i*6+X], degOfFreedom[i*6+Y], degOfFreedom[i*6+Z]);	
				//degOfFreedomToRotMat(rotationMatrix, degOfFreedom[i*6+X], degOfFreedom[i*6+Y], degOfFreedom[i*6+Z]);	
				calcFastTranslationRotation(Protein, xray, rotationMatrix, i, params);
				dx[X]=degOfFreedom[i*6+3+X];
				dx[Y]=degOfFreedom[i*6+3+Y];
				dx[Z]=degOfFreedom[i*6+3+Z];
				translateXRay2(xray, Protein, dx, i);
			}
			calcTotalFormFactor(xray);
		}
		else
		{
			setSegmentOccupancies(Protein, degOfFreedom, nprot);
			for (int i=0;i<nprot;i++)
			{
				calcSymmetryScattering(xray, Protein, i);
			}
		}
	}
	else
	{
		if (pick%6<3)	//Protein rotation
		{
			if (params.UseFastTranslationRotation)
			{
				calcRotationMatrix(rotationMatrix, degOfFreedom[prot*6+X], degOfFreedom[prot*6+Y], degOfFreedom[prot*6+Z]);	
				//gettimeofday(&start, NULL);
				calcFastTranslationRotation(Protein, xray, rotationMatrix, prot, params);
				//gettimeofday(&end, NULL);
				//cout <<"calcFastTranslationRotation took "<<calcTimeDiff(start, end)<<endl;
				//endProgram(__LINE__, __FILE__);
			}
			else
			{
				//gettimeofday(&start, NULL);
				calcSymmetryScattering(xray, Protein, prot);
				//gettimeofday(&end, NULL);
				//cout <<"calcSymmetryScattering took "<<calcTimeDiff(start, end)<<endl;
				//endProgram(__LINE__, __FILE__);
			}
		}
		dx[X]=degOfFreedom[prot*6+3+X];
		dx[Y]=degOfFreedom[prot*6+3+Y];
		dx[Z]=degOfFreedom[prot*6+3+Z];
		//if (false)
		//if (params.UseVeryFastTranslation && params.NumCopies==1 && isRedundantSpaceGroup(Protein.spaceGroup))
		if (params.UseVeryFastTranslation && params.NumCopies==1 && isRedundantSpaceGroup(Protein.intSpaceGroup))
		{
			translateXRayVeryFast(xray, Protein, dx, prot);
			if (params.RotationWeight!=0)
			{
				cout <<"Warning: Cannot use rotation score "
					<<"with VeryFastTranslation. Setting "
					<<"rotation weight to 0"<<endl;
				params.RotationWeight=0;
			}
		}
		else
		{
			translateXRay2(xray, Protein, dx, prot);
			calcTotalFormFactor(xray);
		}
	}
	//translateXRay2(xray, Protein, dx, prot);
	//calcTotalFormFactor(xray);
	calcIntensity(xray);
}

void rotateToFitUnitCell(ProteinStruct &Protein)
{
	int atom1, atom2;
	Real phi;
	Real length1, length2, length3;
	Vector translate;
	Matrix rotationMatrix;
	VectorStruct v, cp;
	//cout <<"In rotateToFitUnitCell"<<endl;
	SafeAlloc(translate, 3, "translate");
	Safe2DAlloc(rotationMatrix, 3, 3, "rotationMatrix");

	length1=vectorLength(Protein.a);
	length2=vectorLength(Protein.b);
	length3=vectorLength(Protein.c);

	v=findLongestVector(Protein.Atoms, atom1, atom2);
	translate[X]=-Protein.Atoms[atom1].x;
	translate[Y]=-Protein.Atoms[atom1].y;
	translate[Z]=-Protein.Atoms[atom1].z;
	moveAtoms(Protein.Atoms, translate[X], translate[Y], translate[Z]);
	if (length1>=length2 && length1>=length3) 
	{
		cp=calcCrossProduct(v, Protein.a);
		phi=calcAngle(Protein.a, v);
	}
	else if (length2>=length1 && length2>=length3) 
	{
		cp=calcCrossProduct(v, Protein.b);
		phi=calcAngle(Protein.b, v);
	}
	else if (length3>=length1 && length3>=length2) 
	{
		cp=calcCrossProduct(v, Protein.c);
		phi=calcAngle(Protein.c, v);
	}
	else
	{
		phi=0;
		error("Should not get here.", __LINE__, __FILE__);
	}
	vectorNormalize(cp);
	//cout <<"cp= "<<vectorToStr(cp)<<" phi= "<<phi<<endl;
	//printMatrix(rotationMatrix);
	calcRotationMatrix(cp, -phi, rotationMatrix);
	rotateAtomsAroundVector(Protein.Atoms, cp, 0.0, 0.0, 0.0, phi);
	cout <<"phi= "<<phi<<endl;
	cout <<"a= "<<vectorToStr(Protein.a)<<endl;
	cout <<"b= "<<vectorToStr(Protein.b)<<endl;
	cout <<"c= "<<vectorToStr(Protein.c)<<endl;
	cout <<"v= "<<vectorToStr(v)<<endl;
	v=atomVector(Protein.Atoms[atom1], Protein.Atoms[atom2]);
	cout <<"v= "<<vectorToStr(v)<<endl;
	//endProgram(__LINE__, __FILE__);
	moveAtoms(Protein.Atoms, -translate[X], -translate[Y], -translate[Z]);

	//matrixMultiply(translate, rotationMatrix, Protein.center);
	//Protein.principalAxes=rotationMatrix;
}

void rotateToFitUnitCell2(ProteinStruct &Protein)
{
	Real phi;
	Real length1, length2, length3;
	Vector eigenvector;
	Matrix momentOfInertia;
	VectorStruct cp, v;

	length1=vectorLength(Protein.a);
	length2=vectorLength(Protein.b);
	length3=vectorLength(Protein.c);

	calcMomentOfInertia(Protein.Atoms, momentOfInertia);
	calcEigenvector(momentOfInertia, eigenvector);
	v.x=eigenvector[X];
	v.y=eigenvector[Y];
	v.z=eigenvector[Z];
	if (length1>=length2 && length1>=length3) 
	{
		cp=calcCrossProduct(v, Protein.a);
		phi=calcAngle(Protein.a, v);
	}
	else if (length2>=length1 && length2>=length3) 
	{
		cp=calcCrossProduct(v, Protein.b);
		phi=calcAngle(Protein.b, v);
	}
	else if (length3>=length1 && length3>=length2) 
	{
		cp=calcCrossProduct(v, Protein.c);
		phi=calcAngle(Protein.c, v);
	}
	else
	{
		phi=0;
		error("Should not get here.", __LINE__, __FILE__);
	}
	cout <<"phi= "<<phi<<endl;
	vectorNormalize(cp);
	rotateAtomsAroundVector(Protein.Atoms, cp, 0.0, 0.0, 0.0, phi);
	center(Protein.Atoms);
}

void rotateToFitUnitCell3(ProteinStruct &Protein)
{
	Real phi, angle;
	Real length1, length2, length3;
	Vector eigenvector, majorAxis, minorAxis, dx;
	Matrix momentOfInertia, momentOfInertia2;
	Matrix rotMat1, rotMat2, rotMat3, rotMat4;
	Matrix rotMat43, rotMat432, rotMat4321, rotTemp, rotTemp2, rotMat21;
	VectorStruct cp, v, vz, minor, major;
	VectorStruct longestUCV, shortestUCV, newUCV; //UCV= Unit Cell Vector


	length1=vectorLength(Protein.a);
	length2=vectorLength(Protein.b);
	length3=vectorLength(Protein.c);

	calcMomentOfInertia(Protein.Atoms, momentOfInertia);
	calcEigenvector(momentOfInertia, eigenvector);
	v.x=eigenvector[X];
	v.y=eigenvector[Y];
	v.z=eigenvector[Z];

	vz.x=0;
	vz.y=0;
	vz.z=1.0;

	cp=calcCrossProduct(v, vz);
	phi=calcAngle(vz, v);

	cout <<"phi= "<<phi<<endl;
	vectorNormalize(cp);
	rotateAtomsAroundVector(Protein.Atoms, cp, 0.0, 0.0, 0.0, phi);
	calcMomentOfInertia(Protein.Atoms, momentOfInertia);
	calcRotationMatrix(cp, phi, rotMat1);

	cout <<"rotMat1:"<<endl;
	printMatrix(rotMat1);
	printMatrix(momentOfInertia);

	Safe2DAlloc(momentOfInertia2, 2, 2, "momentOfInertia");

	momentOfInertia2[X][X]=momentOfInertia[X][X];
	momentOfInertia2[X][Y]=momentOfInertia[X][Y];

	momentOfInertia2[Y][X]=momentOfInertia[Y][X];
	momentOfInertia2[Y][Y]=momentOfInertia[Y][Y];

	calcEigenvector(momentOfInertia2, eigenvector);

	printVector(eigenvector);
	Vector vx;
	SafePushBack(vx, Real(1.0), "vx");
	SafePushBack(vx, Real(0), "vx");
	angle=calcAngle(eigenvector, vx);

	rotateAtomsNoCenter(Protein.Atoms, 0, 0, angle);
	calcMomentOfInertia(Protein.Atoms, momentOfInertia);
	calcRotationMatrixAroundZ(angle, rotMat2);
	matrixMultiply(rotMat2, rotMat1, rotMat21);

	if (abs(momentOfInertia[X][Y])>100.0)
	{
		rotateAtomsNoCenter(Protein.Atoms, 0, 0, -2.0*angle);
		calcMomentOfInertia(Protein.Atoms, momentOfInertia);
		calcRotationMatrixAroundZ(-2.0*angle, rotTemp);
		matrixMultiply(rotTemp, rotMat21, rotTemp2);
		rotMat21=rotTemp2;
	}

	printMatrix(momentOfInertia);

	if (length1>=length2 && length1>=length3) 
	{
		longestUCV=Protein.a;
	}
	else if (length2>=length1 && length2>=length3) 
	{
		longestUCV=Protein.b;
	}
	else 
	{
		longestUCV=Protein.c;
	}

	if (length1<=length2 && length1<=length3) 
	{
		shortestUCV=Protein.a;
	}
	else if (length2<=length1 && length2<=length3) 
	{
		shortestUCV=Protein.b;
	}
	else  
	{
		shortestUCV=Protein.c;
	}

	if (momentOfInertia[X][X]>=momentOfInertia[Y][Y] && momentOfInertia[X][X]>=momentOfInertia[Z][Z])
	{
		major.x=1.0;
		major.y=0;
		major.z=0;
	}
	if (momentOfInertia[Y][Y]>=momentOfInertia[X][X] && momentOfInertia[Y][Y]>=momentOfInertia[Z][Z])
	{
		major.x=0;
		major.y=1.0;
		major.z=0;
	}
	if (momentOfInertia[Z][Z]>=momentOfInertia[X][X] && momentOfInertia[Z][Z]>=momentOfInertia[Y][Y])
	{
		major.x=0;
		major.y=0;
		major.z=1.0;
	}

	if (momentOfInertia[X][X]<=momentOfInertia[Y][Y] && momentOfInertia[X][X]<=momentOfInertia[Z][Z])
	{
		minor.x=1.0;
		minor.y=0;
		minor.z=0;
	}
	if (momentOfInertia[Y][Y]<=momentOfInertia[X][X] && momentOfInertia[Y][Y]<=momentOfInertia[Z][Z])
	{
		minor.x=0;
		minor.y=1.0;
		minor.z=0;
	}
	if (momentOfInertia[Z][Z]<=momentOfInertia[X][X] && momentOfInertia[Z][Z]<=momentOfInertia[Y][Y])
	{
		minor.x=0;
		minor.y=0;
		minor.z=1.0;
	}


	cp=calcCrossProduct(major, longestUCV);
	if (cp.x==0 && cp.y==0 && cp.z==0) cp.x=1.0;
	angle=calcAngle(major, longestUCV);

	cout <<"major= "<<vectorToStr(major)<<endl;
	cout <<"longestUCV= "<<vectorToStr(longestUCV)<<endl;
	cout <<"cp= "<<vectorToStr(cp)<<endl;
	cout <<"angle= "<<angle<<endl;

	rotateAtomsAroundVector(Protein.Atoms, cp, 0.0, 0.0, 0.0, -angle);
	calcRotationMatrix(cp, -angle, rotMat3);

	calcMomentOfInertia(Protein.Atoms, momentOfInertia);
	printMatrix(momentOfInertia);

	newUCV=rotateAroundVector(shortestUCV, -angle, cp);

	cp=calcCrossProduct(minor, newUCV);
	if (cp.x==0 && cp.y==0 && cp.z==0) cp.x=1.0;
	angle=calcAngle(minor, newUCV);

	cout <<"minor= "<<vectorToStr(minor)<<endl;
	cout <<"newUCV= "<<vectorToStr(newUCV)<<endl;
	cout <<"cp= "<<vectorToStr(cp)<<endl;
	cout <<"angle= "<<angle<<endl;

	rotateAtomsAroundVector(Protein.Atoms, cp, 0.0, 0.0, 0.0, -angle);
	calcRotationMatrix(cp, -angle, rotMat4);
	matrixMultiply(rotMat4, rotMat3, rotMat43);
	matrixMultiply(rotMat43, rotMat21, rotMat4321);
	calcInverseMatrix(rotMat4321, Protein.rotationMatrix);
	//calcInverseMatrix(rotMat1, Protein.rotationMatrix);

	calcMomentOfInertia(Protein.Atoms, momentOfInertia);
	printMatrix(momentOfInertia);

	SafeAlloc(dx, 3, "dx");
	calcCenter(Protein.Atoms, dx[X], dx[Y], dx[Z]);
	printVector(dx, "dx");
	//CenterAtomsInBox(Protein.Atoms, dx[X], dx[Y], dx[Z]);
	//dx[X]=-dx[X]; dx[Y]=-dx[Y]; dx[Z]=-dx[Z];
	matrixMultiply(Protein.rotationMatrix, dx, Protein.translation);
}

void getMinMaxFracCoordinates(ProteinStruct &Protein)
{
	int natom=Protein.Atoms.size();
	Real x,  y, z;
	Real minX, minY, minZ;
	Real maxX, maxY, maxZ;

	convertToFractionalCoordinates(Protein);

	minX=Protein.Atoms[0].frac.x;
	minY=Protein.Atoms[0].frac.y;
	minZ=Protein.Atoms[0].frac.z;

	maxX=Protein.Atoms[0].frac.x;
	maxY=Protein.Atoms[0].frac.y;
	maxZ=Protein.Atoms[0].frac.z;

	for (int i=0;i<natom;i++)
	{
		x=Protein.Atoms[i].frac.x;
		y=Protein.Atoms[i].frac.y;
		z=Protein.Atoms[i].frac.z;

		if (x<minX) minX=x;
		if (y<minY) minY=y;
		if (z<minZ) minZ=z;

		if (x>maxX) maxX=x;
		if (y>maxY) maxY=y;
		if (z>maxZ) maxZ=z;
	}

	cout <<"minX= "<<minX<<" maxX= "<<maxX<<endl;	
	cout <<"minY= "<<minY<<" maxY= "<<maxY<<endl;	
	cout <<"minZ= "<<minZ<<" maxZ= "<<maxZ<<endl;	
}

Real evaluateOrientation(ProteinStruct Protein, Real phi, Real psi, Real theta)
{
	int natom=Protein.Atoms.size();
	Real x, y, z;
	Real temp, sum=0;

	rotateAtoms(Protein.Atoms, phi, psi, theta);
	convertToFractionalCoordinates(Protein);

	for (int i=0;i<natom;i++)
	{
		x=Protein.Atoms[i].frac.x;
		y=Protein.Atoms[i].frac.y;
		z=Protein.Atoms[i].frac.z;
		temp=x*x*x+y*y*y+z*z*z;
		temp+=3.0*x*y*y+3.0*x*z*z;
		temp+=3.0*y*x*x+3.0*y*z*z;
		temp+=3.0*z*x*x+3.0*z*y*y;
		temp+=6.0*x*y*z;
		sum+=abs(temp);
	}
	return sum;
}

void findBestOrientation(ProteinStruct &Protein, Real &bestPhi, Real &bestPsi, Real &bestTheta)
{
	int increments=51;
	Real bestScore, score, inc;

	inc=2.0*pi/increments;
	bestScore=evaluateOrientation(Protein, 0, 0, 0);
	bestPhi=bestPsi=bestTheta=0;

	for (Real phi=-pi;phi<=pi;phi+=inc)
	{
		for (Real psi=-pi;psi<=pi;psi+=inc)
		{
			for (Real theta=-pi;theta<=pi;theta+=inc)
			{
				score=evaluateOrientation(Protein, phi, psi, theta);
				if (score<=bestScore)
				{
					bestPhi=phi;
					bestPsi=psi;
					bestTheta=theta;
				}
			}
		}
	}
	rotateAtoms(Protein.Atoms, bestPhi, bestPsi, bestTheta);
	getMinMaxFracCoordinates(Protein);
	printPdb("/home/jouko/pdbs/1UBQ/1UBQ_bestOrientation.pdb", Protein.Atoms);
	exit(EXIT_FAILURE);
}

void findBestOrientation(ProteinStruct &Protein)
{
	Real phi, psi, theta;

	findBestOrientation(Protein, phi, psi, theta);
}

void calcScatteringContinuous(ProteinStruct &Protein, XRayStruct &xray, Array3D &realContinuous, Array3D &imagContinuous, XRayParamStruct &params)
{
	//Used for FastRotationTranslation.  Calculate this once and then just
	//transform h, k, l to get f(h,k,l) instead of looping over atoms.
	int atomID, bestIndex;
	int npoint=xray.q.size();
	int nreal=xray.real.size();
	int natom=Protein.Atoms.size();
	int hBins, kBins, lBins;
	Real dist, max, dwf;
	Real product, occupancy, q;
	Real x, y, z;
	Real fracx, fracy, fracz;
	Real millerX, millerY, millerZ;
	Real real, imag, bulkSolventDensity, scale;
	Vector f, v, miller, miller2;
	Real a[5][NumAtomTypes], b[5][NumAtomTypes], c[NumAtomTypes];
	Vector excludedVolumeRadii;
	Matrix cosx, cosy, cosz;
	Matrix sinx, siny, sinz;
	Matrix dwfx, dwfy, dwfz;
	Matrix rotMat;
	PosStruct qContinuous;
	timeval start, end;

	if (npoint!=nreal)
	{
		cout <<"ERROR: npoint!=nreal."<<endl;
		cout <<"npoint= "<<npoint<<" nreal= "<<nreal<<endl;
		exit(EXIT_FAILURE);
	}
	gettimeofday(&start, NULL);
	zeroVector(xray.i);
	zeroVector(xray.real);
	zeroVector(xray.imag);
	zeroVector(xray.waterReal);
	zeroVector(xray.waterImag);
	SafeAlloc(qContinuous.Pos, 3, "qContinuous");
	SafeAlloc(miller, 3, "miller");
	SafeAlloc(miller2, 3, "miller2");
	scale=params.scale;
	bulkSolventDensity=params.BulkDensity;
	hBins=params.ContinuousHValues;
	kBins=params.ContinuousKValues*2;
	lBins=params.ContinuousLValues*2;
	cout <<"hBins= "<<hBins<<" kBins= "<<kBins<<" lBins= "<<lBins<<endl;
	Safe2DAlloc(cosx, natom, hBins, "cosx");
	Safe2DAlloc(cosy, natom, kBins, "cosx");
	Safe2DAlloc(cosz, natom, lBins, "cosx");
	Safe2DAlloc(sinx, natom, hBins, "sinx");
	Safe2DAlloc(siny, natom, kBins, "siny");
	Safe2DAlloc(sinz, natom, lBins, "sinz");
	Safe2DAlloc(dwfx, natom, hBins, "dwfx");
	Safe2DAlloc(dwfy, natom, kBins, "dwfy");
	Safe2DAlloc(dwfz, natom, lBins, "dwfz");
	SafeAlloc(f, NumAtomTypes, "f");

	int Size1, Size2;
	Get2DVectorSize(dwfx, Size1, Size2);
	cout <<"Size1= "<<Size1<<" Size2= "<<Size2<<endl;

	FourGaussianParameters(a, b, c, scale);
	setExcludedVolumeRadii(excludedVolumeRadii, params);
	//setExcludedVolumeRadii(excludedVolumeRadii);
	max=findMaxMillerForFastRotationAnalytical2(xray, params, bestIndex);
	for (int i=0;i<natom;i++)
	{
		fracx=Protein.Atoms[i].frac.x;
		fracy=Protein.Atoms[i].frac.y;
		fracz=Protein.Atoms[i].frac.z;
		for (int j=0;j<hBins;j++)
		{
			millerX=Real(j)*xray.hBinSize;
			q=2.0*pi*millerX/Protein.XBoxLength;
			x=2.0*pi*millerX*fracx;
			//cout <<"j= "<<j<<" millerX= "<<millerX<<" q= "<<q<<" x= "<<x<<endl;
			getTrigLookUp(xray.lookUp, x, cosx[i][j], sinx[i][j]);	
			dwfx[i][j]=calcDWF(q, Protein.Atoms[i].BFactor, xray.lookUp);
			cosx[i][j]=cos(x);
			sinx[i][j]=sin(x);
			//cout <<"cosx= "<<cosx[i][j]<<" sinx= "<<sinx[i][j]<<endl;
		}
		for (int j=0;j<kBins;j++)
		{
			millerY=Real(j)*xray.kBinSize;
			millerY-=Real(xray.maxContinuousK);
			q=2.0*pi*millerY/Protein.YBoxLength;
			y=2.0*pi*millerY*fracy;
			getTrigLookUp(xray.lookUp, y, cosy[i][j], siny[i][j]);	
			dwfy[i][j]=calcDWF(q, Protein.Atoms[i].BFactor, xray.lookUp);
			cosy[i][j]=cos(y);
			siny[i][j]=sin(y);
			//cout <<"j= "<<j<<" millerY= "<<millerY<<" q= "<<q<<" y= "<<y<<" cosy= "<<cosy[i][j]<<" siny= "<<siny[i][j]<<endl;
		}
		for (int j=0;j<lBins;j++)
		{
			millerZ=Real(j)*xray.lBinSize;
			millerZ-=Real(xray.maxContinuousL);
			q=2.0*pi*millerZ/Protein.ZBoxLength;
			//cout <<"j= "<<j<<" miller= "<<miller<<" q= "<<q<<endl;
			z=2.0*pi*millerZ*fracz;
			getTrigLookUp(xray.lookUp, z, cosz[i][j], sinz[i][j]);	
			dwfz[i][j]=calcDWF(q, Protein.Atoms[i].BFactor, xray.lookUp);
			cosz[i][j]=cos(z);
			sinz[i][j]=sin(z);
		}
	}
	//endProgram(__LINE__, __FILE__);
	cout <<"hBinSize= "<<xray.hBinSize<<endl;
	cout <<"kBinSize= "<<xray.kBinSize<<endl;
	cout <<"lBinSize= "<<xray.lBinSize<<endl;
	for (int h=0;h<hBins;h++)
	{
		cout <<"h= "<<h<<"\t";
		miller[X]=Real(h)*xray.hBinSize;
		//qContinuous.Pos[X]=2.0*pi*miller[X];
		//qContinuous.Pos[X]/=Protein.XBoxLength;
		for (int k=0;k<kBins;k++)
		{
			miller[Y]=Real(k)*xray.kBinSize-Real(xray.maxContinuousK);
			//qContinuous.Pos[Y]=2.0*pi*miller[Y];
			//qContinuous.Pos[Y]/=Protein.YBoxLength;
			for (int l=0;l<lBins;l++)
			{
				miller[Z]=Real(l)*xray.lBinSize-Real(xray.maxContinuousL);
				//qContinuous.Pos[Z]=2.0*pi*miller[Z];
				//qContinuous.Pos[Z]/=Protein.ZBoxLength;
				//Change the size of miller so that grid points just outside of boundary will be calculated
				//so that interpolation will work for points just within the boundary
				if (miller[X]>0) miller2[X]=miller[X]-xray.hBinSize;
				else miller2[X]=miller[X]+xray.hBinSize;
				if (miller[Y]>0) miller2[Y]=miller[Y]-xray.kBinSize;
				else miller2[Y]=miller[Y]+xray.kBinSize;
				if (miller[Z]>0) miller2[Z]=miller[Z]-xray.lBinSize;
				else miller2[Z]=miller[Z]+xray.lBinSize;
				matrixMultiply(miller2, Protein.cartToFrac, v);
				dist=v[X]*v[X]+v[Y]*v[Y]+v[Z]*v[Z];
				if (dist<max) //Only calculate intensity for miller indexes which can be reached by rotation
					//if (true) //Only calculate intensity for miller indexes which can be reached by rotation
				{
					matrixMultiply(miller, Protein.cartToFrac, qContinuous.Pos);
					for (int j=0;j<3;j++) qContinuous.Pos[j]*=2.0*pi;
					for (atomID=0;atomID<NumAtomTypes;atomID++)
					{
						f[atomID]=solventCorrectedScatteringFactor(qContinuous, a, b, c, atomID, excludedVolumeRadii[atomID], bulkSolventDensity);
					}
					//if (count%100==0)
					//{
					//	cout <<calcMagnitude(qContinuous)<<"\t"<<f[0]<<endl;
					//}
					for (int j=0;j<natom;j++)
					{

						atomID=Protein.Atoms[j].atomid;
						occupancy=Protein.Atoms[j].Occupancy;
						real=imag=0;

						real+=cosx[j][h]*cosy[j][k]*cosz[j][l];
						real-=cosx[j][h]*siny[j][k]*sinz[j][l];
						real-=sinx[j][h]*cosy[j][k]*sinz[j][l];
						real-=sinx[j][h]*siny[j][k]*cosz[j][l];

						imag-=sinx[j][h]*siny[j][k]*sinz[j][l];
						imag+=sinx[j][h]*cosy[j][k]*cosz[j][l];
						imag+=cosx[j][h]*siny[j][k]*cosz[j][l];
						imag+=cosx[j][h]*cosy[j][k]*sinz[j][l];
						if (Protein.alpha==90.0 && Protein.beta==90.0 && Protein.gamma==90.0)
						{
							dwf=dwfx[j][h]*dwfy[j][k]*dwfz[j][l];
						}
						else
						{
							dwf=calcDWF(qContinuous, Protein.Atoms[j].BFactor, xray.lookUp);
						}
						product=dwf*f[atomID]*occupancy;
						realContinuous[h][k][l]+=real*product;
						imagContinuous[h][k][l]+=imag*product;
						   //cout <<"cosx= "<<cosx[j][h]<<" cosy= "<<cosy[j][k]<<" cosz= "<<cosz[j][l]<<endl;
						   //cout <<"sinx= "<<sinx[j][h]<<" siny= "<<siny[j][k]<<" sinz= "<<sinz[j][l]<<endl;
						   //cout <<"product= "<<product<<" dwfx= "<<dwfx[j][h]<<" dwfy= "<<dwfy[j][k]<<" dwfz= "<<dwfz[j][l]<<endl;
						   //cout <<"f= "<<f[atomID]<<" occupancy= "<<occupancy<<endl;
						   //cout <<"real= "<<real<<" imag= "<<imag<<endl;
						/*
						   if (isnan(real) || isnan(imag))
						   {
						   cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<" j= "<<j<<endl;
						   cout <<"atomID= "<<atomID<<endl;
						   cout <<"dotproduct= "<<dp<<endl;
						   cout <<"h= "<<xray.continuousMiller[h][k][l].Pos[X]
						   <<" k= "<<xray.continuousMiller[h][k][l].Pos[Y]
						   <<" l= "<<xray.continuousMiller[h][k][l].Pos[Z]<<endl;
						   cout <<"x= "<<Protein.Atoms[j].frac.x
						   <<" y= "<<Protein.Atoms[j].frac.y
						   <<" z= "<<Protein.Atoms[j].frac.z<<endl;
						   cout <<"dwf= "<<dwf<<" real= "<<real<<" imag= "<<imag
						   <<" f= "<<f<<endl;
						   exit(EXIT_FAILURE);
						   }
						 */
					}
				}
			}
		}
	}
	//endProgram(__LINE__, __FILE__);
	gettimeofday(&end, NULL);
	cout <<"calcScatteringContinuous took "<<calcTimeDiff(start, end)<<endl;
}

void placeProteinForFastTranslationRotation(ProteinStruct &Protein)
{
	Real xave, yave, zave;

	SafeAlloc(Protein.center, 3, "Protein.center");
	rotateToFitUnitCell3(Protein);

	calcCenterAll(Protein.Atoms, xave, yave, zave);
	Protein.center[X]=xave;
	Protein.center[Y]=yave;
	Protein.center[Z]=zave;
	center(Protein.Atoms);
	calcCartToFrac(Protein);
	convertToFractionalCoordinates(Protein);
}

void calcPlaceProteinAndScatteringContinuous(ProteinStruct &Protein, XRayStruct &xray, Array3D &realContinious, Array3D &imagContinous, XRayParamStruct &params)
{
	placeProteinForFastTranslationRotation(Protein);
	printPdb(params.CenteredAndRotatedOutputPdb, Protein);
	calcScatteringContinuous(Protein, xray, xray.realContinuous, xray.imagContinuous, params);
}

void calcScatteringContinuousSegments(ProteinStruct &Protein, XRayStruct &xray, XRayParamStruct &params)
{
	int natom=Protein.Atoms.size();
	int nsegment=params.NumOccupancySegments;
	int start, end;
	ProteinStruct tempProtein;

	placeProteinForFastTranslationRotation(Protein);
	zeroContinuousScattering(xray);
	tempProtein=Protein;
	for (int i=0;i<nsegment;i++)
	{
		start=int(Real(natom)*Real(i)/Real(nsegment));
		end=int(Real(natom)*Real(i+1)/Real(nsegment))-1;
		cout <<"start= "<<start<<" end= "<<end<<endl;
		getSegment(Protein.Atoms, start, end, tempProtein.Atoms);
		SafePushBack(xray.realContinuousSegment, xray.realContinuous, "xray.realContinuousSegment");
		SafePushBack(xray.imagContinuousSegment, xray.imagContinuous, "xray.imagContinuousSegment");
		calcScatteringContinuous(tempProtein, xray, xray.realContinuousSegment[i], xray.imagContinuousSegment[i], params);
	}
}

void calcScatteringContinuousSegmentsFile(ProteinStruct &Protein, XRayStruct &xray, vector<int> &segmentStart, vector<int> &segmentEnd, XRayParamStruct &params)
{
	int nsegment=segmentStart.size();
	int natom=Protein.Atoms.size();
	int hBins, kBins, lBins;
	ProteinStruct tempProtein;

	cout <<"Protein.Atoms.size()= "<<Protein.Atoms.size()<<endl;
	placeProteinForFastTranslationRotation(Protein);
	cout <<"Protein.Atoms.size()= "<<Protein.Atoms.size()<<endl;
	zeroContinuousScattering(xray);
	cout <<"Protein.Atoms.size()= "<<Protein.Atoms.size()<<endl;
	tempProtein=Protein;
	cout <<"tempProtein.Atoms.size()= "<<tempProtein.Atoms.size()<<endl;
	for (int i=0;i<nsegment;i++)
	{
		getSegment(Protein.Atoms, segmentStart[i], segmentEnd[i], tempProtein.Atoms);
		SafePushBack(xray.realContinuousSegment, xray.realContinuous, "xray.realContinuousSegment");
		SafePushBack(xray.imagContinuousSegment, xray.imagContinuous, "xray.imagContinuousSegment");
		calcScatteringContinuous(tempProtein, xray, xray.realContinuousSegment[i], xray.imagContinuousSegment[i], params);
	}
	cout <<"tempProtein.Atoms.size()= "<<tempProtein.Atoms.size()<<endl;
	tempProtein=Protein;
	setNonFixedOccupancyToZero(tempProtein, segmentStart, segmentEnd);
	cout <<"tempProtein.Atoms.size()= "<<tempProtein.Atoms.size()<<endl;
	for (int i=0;i<natom;i++)
	{
		//cout <<"i= "<<i<<" occupancy= "<<tempProtein.Atoms[i].Occupancy<<endl;
	}
	SafePushBack(xray.realContinuousSegment, xray.realContinuous, "xray.realContinuousSegment");
	SafePushBack(xray.imagContinuousSegment, xray.imagContinuous, "xray.imagContinuousSegment");
	calcScatteringContinuous(tempProtein, xray, xray.realContinuousSegment[nsegment], xray.imagContinuousSegment[nsegment], params);
	hBins=params.ContinuousHValues;
	kBins=params.ContinuousKValues*2;
	lBins=params.ContinuousLValues*2;
	//Safe4DAlloc(xray.updatedContinuous, nsegment+1, hBins, kBins, lBins, "updatedContinuous");
	Safe3DAlloc(xray.updatedContinuous, hBins, kBins, lBins, "updatedContinuous");
	params.NumOccupancySegments=nsegment+1;
}

void calcScatteringContinuousSegmentsFile(ProteinStruct &Protein, XRayStruct &xray, XRayParamStruct &params)
{
	vector<int> segmentStart, segmentEnd;

	readSegmentsFile(params.SegmentsFile, segmentStart, segmentEnd);
	residueSegmentsToAtomSegments(Protein.Atoms, segmentStart, segmentEnd);
	calcScatteringContinuousSegmentsFile(Protein, xray, segmentStart, segmentEnd, params);
}

void calcScatteringContinuous(ProteinStruct &Protein, XRayStruct &xray, XRayParamStruct &params, Real h, Real k, Real l, Real &totalReal, Real &totalImag)
{
	int atomID;
	int npoint=xray.q.size();
	int nreal=xray.real.size();
	int natom=Protein.Atoms.size();
	//Real xave, yave, zave;
	Real dp, real, imag, dwf, f, bulkSolventDensity, scale;
	Real a[5][NumAtomTypes], b[5][NumAtomTypes], c[NumAtomTypes];
	Vector excludedVolumeRadii;
	Matrix rotMat;
	vector<AtomStruct> CopyAtoms;
	PosStruct q;
	timeval start, end;
	if (npoint!=nreal)
	{
		cout <<"ERROR: npoint!=nreal."<<endl;
		cout <<"npoint= "<<npoint<<" nreal= "<<nreal<<endl;
		exit(EXIT_FAILURE);
	}

	totalReal=0;
	totalImag=0;
	gettimeofday(&start, NULL);
	scale=params.scale;
	bulkSolventDensity=params.BulkDensity;
	CopyAtoms=Protein.Atoms;
	SafeAlloc(q.Pos, 3, "Pos");
	SafeAlloc(Protein.center, 3, "Protein.center");
	//rotateToFitUnitCell(Protein);
	//calcCenter(Protein.Atoms, xave, yave, zave);
	//Protein.center[X]=xave;
	//Protein.center[Y]=yave;
	//Protein.center[Z]=zave;
	//center(Protein.Atoms);
	convertToFractionalCoordinates(Protein);
	FourGaussianParameters(a, b, c, scale);
	setExcludedVolumeRadii(excludedVolumeRadii);
	millerToQ(h, k, l, Protein.XBoxLength, Protein.YBoxLength, Protein.ZBoxLength, q);
	//printAtomInfo(Protein.Atoms[0]);
	for (int j=0;j<natom;j++)
	{

		atomID=Protein.Atoms[j].atomid;
		dp=Protein.Atoms[j].frac.x*h; 
		dp+=Protein.Atoms[j].frac.y*k;
		dp+=Protein.Atoms[j].frac.z*l;
		dwf=calcDWF(q, Protein.Atoms[j], xray.lookUp);
		f=solventCorrectedScatteringFactor(q, a, b, c, atomID, excludedVolumeRadii[atomID], bulkSolventDensity);
		real=f*dwf*cos(2.0*pi*dp);
		imag=f*dwf*sin(2.0*pi*dp);
		real*=Protein.Atoms[j].Occupancy;
		imag*=Protein.Atoms[j].Occupancy;
		totalReal+=real;
		totalImag+=imag;
		if (isnan(real) || isnan(imag))
		{
			cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<" j= "<<j<<endl;
			cout <<"atomID= "<<atomID<<endl;
			cout <<"dotproduct= "<<dp<<endl;
			cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<endl;
			cout <<"x= "<<Protein.Atoms[j].frac.x
				<<" y= "<<Protein.Atoms[j].frac.y
				<<" z= "<<Protein.Atoms[j].frac.z<<endl;
			cout <<"dwf= "<<dwf<<" real= "<<real<<" imag= "<<imag
				<<" f= "<<f<<endl;
			exit(EXIT_FAILURE);
		}
	}
	//Protein.Atoms=CopyAtoms;
	gettimeofday(&end, NULL);
	//cout <<"calcScatteringBfactor took "<<calcTimeDiff(start, end)<<endl;
}

void calcScatteringContinousSegments(ProteinStruct &Protein, XRayStruct &xray, XRayParamStruct &params)
{

}

/*
   void calcScatteringContinuousGrad(ProteinStruct &Protein, XRayStruct &xray, XRayParamStruct &params)
   {
//Calculates the structure factor, gradient of the structure factor, 
//and hessian of the structure factor.
int atomID;
int npoint=xray.q.size();
int nreal=xray.real.size();
int natom=Protein.Atoms.size();
int hBins, kBins, lBins;
Real xn, xm;
Real xave, yave, zave;
Real dp, real, imag, dwf, f, df, bulkSolventDensity, scale;
Real realHes, imagHes;
Real a[5][NumAtomTypes], b[5][NumAtomTypes], c[NumAtomTypes];
Vector excludedVolumeRadii;
Matrix rotMat;
vector<AtomStruct> CopyAtoms;
Vector Pos;
timeval start, end;
if (npoint!=nreal)
{
cout <<"ERROR: npoint!=nreal."<<endl;
cout <<"npoint= "<<npoint<<" nreal= "<<nreal<<endl;
exit(EXIT_FAILURE);
}
gettimeofday(&start, NULL);
zeroVector(xray.i);
zeroVector(xray.real);
zeroVector(xray.imag);
zeroVector(xray.waterReal);
zeroVector(xray.waterImag);
SafeAlloc(Pos, 3, "Pos");
scale=params.scale;
bulkSolventDensity=params.BulkDensity;
CopyAtoms=Protein.Atoms;
SafeAlloc(Protein.center, 3, "Protein.center");
//rotateToFitUnitCell(Protein);
calcCenter(Protein.Atoms, xave, yave, zave);
Protein.center[X]=xave;
Protein.center[Y]=yave;
Protein.center[Z]=zave;
center(Protein.Atoms);
//findBestOrientation(Protein);
convertToFractionalCoordinates(Protein);
Get3DVectorSize(xray.continuousMiller, hBins, kBins, lBins, "continuousMiller");
FourGaussianParameters(a, b, c, scale);
setExcludedVolumeRadii(excludedVolumeRadii);
for (int h=0;h<hBins;h++)
{
for (int k=0;k<kBins;k++)
{
for (int l=0;l<lBins;l++)
{
for (int j=0;j<natom;j++)
{
Pos[X]=Protein.Atoms[j].frac.x;
Pos[Y]=Protein.Atoms[j].frac.y;
Pos[Z]=Protein.Atoms[j].frac.z;
atomID=Protein.Atoms[j].atomid;
dp=calcDotProduct(Protein.Atoms[j].frac, xray.continuousMiller[h][k][l].Pos);
dwf=calcDWF(xray.qContinuous[h][k][l], Protein.Atoms[j], xray.lookUp);
f=solventCorrectedScatteringFactor(xray.qContinuous[h][k][l], a, b, c, atomID, excludedVolumeRadii[atomID], bulkSolventDensity);
real=f*dwf*cos(2.0*pi*dp);
imag=f*dwf*sin(2.0*pi*dp);;
real*=Protein.Atoms[j].Occupancy;
imag*=Protein.Atoms[j].Occupancy;
xray.realContinuous[h][k][l]+=real;
xray.imagContinuous[h][k][l]+=imag;
for (int n=0;n<3;n++)
{
//This is only an approximation it assumes that d(f*dwf)/dh=0.  Fix later.
df=solventCorrectedScatteringFactorDer(h, k, l, a, b, atomID, excludedVolumeRadii[atomID], bulkSolventDensity, Protein.a, Protein.b, Protein.c, n);
xn=Pos[n];

xray.realContinuousDer[h][k][l].Pos[n]+=2.0*pi*xn*imag;
xray.realContinuousDer[h][k][l].Pos[n]+=df*dwf*cos(2.0*pi*dp);
xray.imagContinuousDer[h][k][l].Pos[n]-=2.0*pi*xn*real;
xray.realContinuousDer[h][k][l].Pos[n]+=df*dwf*sin(2.0*pi*dp);
for (int m=n;m<3;m++)
{
	xm=Pos[m];
	realHes=-4.0*pi*pi*xn*xm*real;
	imagHes=-4.0*pi*pi*xn*xm*imag;
	xray.realContinuousHes[h][k][l][n][m]+=realHes;
	xray.imagContinuousHes[h][k][l][n][m]+=imagHes;
	xray.realContinuousHes[h][k][l][m][n]=xray.realContinuousHes[h][k][l][n][m];
	xray.imagContinuousHes[h][k][l][m][n]=xray.imagContinuousHes[h][k][l][n][m];
}
}
if (isnan(real) || isnan(imag))
{
	cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<" j= "<<j<<endl;
	cout <<"atomID= "<<atomID<<endl;
	cout <<"dotproduct= "<<dp<<endl;
	cout <<"h= "<<xray.continuousMiller[h][k][l].Pos[X]
		<<" k= "<<xray.continuousMiller[h][k][l].Pos[Y]
		<<" l= "<<xray.continuousMiller[h][k][l].Pos[Z]<<endl;
	cout <<"x= "<<Protein.Atoms[j].frac.x
		<<" y= "<<Protein.Atoms[j].frac.y
		<<" z= "<<Protein.Atoms[j].frac.z<<endl;
	cout <<"dwf= "<<dwf<<" real= "<<real<<" imag= "<<imag
		<<" f= "<<f<<endl;
	exit(EXIT_FAILURE);
}
}
}
}
}
Protein.Atoms=CopyAtoms;
gettimeofday(&end, NULL);
cout <<"Leaving calcScatteringContinuous"<<endl;
//cout <<"calcScatteringBfactor took "<<calcTimeDiff(start, end)<<endl;
}
*/

Real interpolate3DFast(Real millerX, Real millerY, Real millerZ, Array3D &y, vector<Real> &pos, Real invDX, Real invDY, Real invDZ)
{
	//en.wikipedia.org/wiki/Trilinear_interpolation
	Real xd, yd, zd;
	Real c00, c10, c01, c11;
	Real c0, c1;

	xd=(pos[X]-millerX)*invDX;
	yd=(pos[Y]-millerY)*invDY;
	zd=(pos[Z]-millerZ)*invDZ;
	/*
	   for (int i=0;i<2;i++)
	   {
	   for (int j=0;j<2;j++)
	   {
	   for (int k=0;k<2;k++)
	   {
	   cout <<"y= "<<y[i][j][k]<<endl;
	   }
	   }
	   }
	 */
	c00=y[0][0][0]*(1.0-xd)+y[1][0][0]*xd;
	c10=y[0][1][0]*(1.0-xd)+y[1][1][0]*xd;
	c01=y[0][0][1]*(1.0-xd)+y[1][0][1]*xd;
	c11=y[0][1][1]*(1.0-xd)+y[1][1][1]*xd;

	c0=c00*(1.0-yd)+c10*yd;
	c1=c01*(1.0-yd)+c11*yd;

	return c0*(1.0-zd)+c1*zd;
}

Real interpolate3DFast(Real xd, Real yd, Real zd, Array3D &y)
{
	//en.wikipedia.org/wiki/Trilinear_interpolation
	Real c00, c10, c01, c11;
	Real c0, c1;

	/*
	   for (int i=0;i<2;i++)
	   {
	   for (int j=0;j<2;j++)
	   {
	   for (int k=0;k<2;k++)
	   {
	   cout <<"y= "<<y[i][j][k]<<endl;
	   }
	   }
	   }
	 */
	c00=y[0][0][0]*(1.0-xd)+y[1][0][0]*xd;
	c10=y[0][1][0]*(1.0-xd)+y[1][1][0]*xd;
	c01=y[0][0][1]*(1.0-xd)+y[1][0][1]*xd;
	c11=y[0][1][1]*(1.0-xd)+y[1][1][1]*xd;

	c0=c00*(1.0-yd)+c10*yd;
	c1=c01*(1.0-yd)+c11*yd;

	return c0*(1.0-zd)+c1*zd;
}

void interpolate3DFast(Real millerX, Real millerY, Real millerZ, Array3D &real, Array3D &imag, Vector &hNew, Real invDX, Real invDY, Real invDZ, Real &realOut, Real &imagOut)
{
	Real xd, yd, zd;
	Real c00, c10, c01, c11;
	Real c0, c1;

	xd=(hNew[X]-millerX)*invDX;
	yd=(hNew[Y]-millerY)*invDY;
	zd=(hNew[Z]-millerZ)*invDZ;

	//realOut=interpolate3DFast(xd, yd, zd, real);
	//imagOut=interpolate3DFast(xd, yd, zd, imag);
	c00=real[0][0][0]*(1.0-xd)+real[1][0][0]*xd;
	c10=real[0][1][0]*(1.0-xd)+real[1][1][0]*xd;
	c01=real[0][0][1]*(1.0-xd)+real[1][0][1]*xd;
	c11=real[0][1][1]*(1.0-xd)+real[1][1][1]*xd;

	c0=c00*(1.0-yd)+c10*yd;
	c1=c01*(1.0-yd)+c11*yd;

	realOut=c0*(1.0-zd)+c1*zd;

	c00=imag[0][0][0]*(1.0-xd)+imag[1][0][0]*xd;
	c10=imag[0][1][0]*(1.0-xd)+imag[1][1][0]*xd;
	c01=imag[0][0][1]*(1.0-xd)+imag[1][0][1]*xd;
	c11=imag[0][1][1]*(1.0-xd)+imag[1][1][1]*xd;

	c0=c00*(1.0-yd)+c10*yd;
	c1=c01*(1.0-yd)+c11*yd;

	imagOut=c0*(1.0-zd)+c1*zd;
}  

void interpolate3DVeryFast(Real millerX, Real millerY, Real millerZ, XRayStruct &xray, Vector &hNew, Real invDX, Real invDY, Real invDZ, Real &realOut, Real &imagOut, int hbin, int kbin, int lbin)
{
	Real xd, yd, zd;
	Real c00, c10, c01, c11;
	Real c0, c1;

	xd=(hNew[X]-millerX)*invDX;
	yd=(hNew[Y]-millerY)*invDY;
	zd=(hNew[Z]-millerZ)*invDZ;

	//realOut=interpolate3DFast(xd, yd, zd, real);
	//imagOut=interpolate3DFast(xd, yd, zd, imag);
	c00=xray.realContinuous[hbin][kbin][lbin]*(1.0-xd)+xray.realContinuous[hbin+1][kbin][lbin]*xd;
	c10=xray.realContinuous[hbin][kbin+1][lbin]*(1.0-xd)+xray.realContinuous[hbin+1][kbin+1][lbin]*xd;
	c01=xray.realContinuous[hbin][kbin][lbin+1]*(1.0-xd)+xray.realContinuous[hbin+1][kbin][lbin+1]*xd;
	c11=xray.realContinuous[hbin][kbin+1][lbin+1]*(1.0-xd)+xray.realContinuous[hbin+1][kbin+1][lbin+1]*xd;

	c0=c00*(1.0-yd)+c10*yd;
	c1=c01*(1.0-yd)+c11*yd;

	realOut=c0*(1.0-zd)+c1*zd;

	c00=xray.imagContinuous[hbin][kbin][lbin]*(1.0-xd)+xray.imagContinuous[hbin+1][kbin][lbin]*xd;
	c10=xray.imagContinuous[hbin][kbin+1][lbin]*(1.0-xd)+xray.imagContinuous[hbin+1][kbin+1][lbin]*xd;
	c01=xray.imagContinuous[hbin][kbin][lbin+1]*(1.0-xd)+xray.imagContinuous[hbin+1][kbin][lbin+1]*xd;
	c11=xray.imagContinuous[hbin][kbin+1][lbin+1]*(1.0-xd)+xray.imagContinuous[hbin+1][kbin+1][lbin+1]*xd;

	c0=c00*(1.0-yd)+c10*yd;
	c1=c01*(1.0-yd)+c11*yd;

	imagOut=c0*(1.0-zd)+c1*zd;
}  

void interpolate3DVeryFast2(Real millerX, Real millerY, Real millerZ, XRayStruct &xray, Vector &hNew, Real invDX, Real invDY, Real invDZ, Real &realOut, Real &imagOut, int hbin, int kbin, int lbin)
{
	Real xd, yd, zd;
	Real c00, c10, c01, c11;
	Real c0, c1;

	xd=(hNew[X]-millerX)*invDX;
	yd=(hNew[Y]-millerY)*invDY;
	zd=(hNew[Z]-millerZ)*invDZ;

	//cout <<"xd= "<<xd<<" yd= "<<yd<<" zd= "<<zd<<endl;
	//cout <<"hbin= "<<hbin<<" kbin= "<<kbin<<" lbin= "<<lbin<<endl;
	//cout <<"real= "<<xray.realContinuous[hbin][kbin][lbin]<<endl;

	c00=xray.realContinuous[hbin][kbin][lbin]*(1.0-zd)+xray.realContinuous[hbin][kbin][lbin+1]*zd;
	c10=xray.realContinuous[hbin][kbin+1][lbin]*(1.0-zd)+xray.realContinuous[hbin][kbin+1][lbin+1]*zd;
	c01=xray.realContinuous[hbin+1][kbin][lbin]*(1.0-zd)+xray.realContinuous[hbin+1][kbin][lbin+1]*zd;
	c11=xray.realContinuous[hbin+1][kbin+1][lbin]*(1.0-zd)+xray.realContinuous[hbin+1][kbin+1][lbin+1]*zd;

	//cout <<"c00= "<<c00<<endl;
	//cout <<"c10= "<<c10<<endl;
	//cout <<"c01= "<<c01<<endl;
	//cout <<"c11= "<<c11<<endl;

	c0=c00*(1.0-yd)+c10*yd;
	c1=c01*(1.0-yd)+c11*yd;

	//cout <<"c0= "<<c0<<endl;
	//cout <<"c1= "<<c1<<endl;

	realOut=c0*(1.0-xd)+c1*xd;
	//cout <<"realOut= "<<realOut<<endl<<endl;

	c00=xray.imagContinuous[hbin][kbin][lbin]*(1.0-zd)+xray.imagContinuous[hbin][kbin][lbin+1]*zd;
	c10=xray.imagContinuous[hbin][kbin+1][lbin]*(1.0-zd)+xray.imagContinuous[hbin][kbin+1][lbin+1]*zd;
	c01=xray.imagContinuous[hbin+1][kbin][lbin]*(1.0-zd)+xray.imagContinuous[hbin+1][kbin][lbin+1]*zd;
	c11=xray.imagContinuous[hbin+1][kbin+1][lbin]*(1.0-zd)+xray.imagContinuous[hbin+1][kbin+1][lbin+1]*zd;

	c0=c00*(1.0-yd)+c10*yd;
	c1=c01*(1.0-yd)+c11*yd;

	imagOut=c0*(1.0-xd)+c1*xd;
}  

void phaseAndInterpolate(XRayStruct &xray, Vector &hNew, Real invDX, Real invDY, Real invDZ, Real &realOut, Real &imagOut)
{
	//Corrects for rotation and translation and preforms interpolation
	//See calcFastTranslationRotation
	bool signFlip=false;
	int hbin, kbin, lbin;
	Real millerX, millerY, millerZ;

	if (hNew[X]<0)
	{
		signFlip=true;
		hNew[X]=-hNew[X]; hNew[Y]=-hNew[Y]; hNew[Z]=-hNew[Z];
	}
	hbin=int(hNew[X]*invDX);
	kbin=int((hNew[Y]+Real(xray.maxContinuousK))*invDY);
	lbin=int((hNew[Z]+Real(xray.maxContinuousL))*invDZ);
	//cout <<"hbin= "<<hbin<<" kbin= "<<kbin<<" lbin= "<<lbin<<endl;

	millerX=Real(hbin)*xray.hBinSize;
	millerY=Real(kbin)*xray.kBinSize-Real(xray.maxContinuousK);
	millerZ=Real(lbin)*xray.lBinSize-Real(xray.maxContinuousL);

	//cout <<endl;
	//cout <<"hNew:"<<endl;
	//printVector(hNew);
	//interpolate3DVeryFast(millerX, millerY, millerZ, xray, hNew, invDX, invDY, invDZ, realOut, imagOut, hbin, kbin, lbin);
	interpolate3DVeryFast2(millerX, millerY, millerZ, xray, hNew, invDX, invDY, invDZ, realOut, imagOut, hbin, kbin, lbin);
	if (signFlip) imagOut=-imagOut;
}



void calcFastTranslationRotation(SymmetryOperationStruct &symmetryOperation, Matrix &cartToFrac, Matrix &rotationMatrix, XRayStruct &xray, Vector &real, Vector &imag, XRayParamStruct &params, ProteinStruct &Protein)
{
	//From Acta Cryst. (1994). A50, 157-163 AMoRe: an Automated Package for
	//Molecular Replacement
	int npoint=xray.miller.size();
	//Real dp1, dp2, dp;
	Real realOut, imagOut;
	Real invDX, invDY, invDZ;
	Vector hMs, hMsD, hMsDR, hNew;	//hNew gives continuous Miller indexes in old frame.  x is translation in fractional ccordinates.
	Vector realContinuous, imagContinuous;
	Matrix fracToCart;
	Matrix RO, MsR;
	Matrix MsRO, DMsRO;
	PosStruct tempPos;
	//timeval start, end;
	//cout <<"In calcFastTranslationRotation"<<endl;
	zeroVector(xray.real);
	zeroVector(xray.imag);

	//gettimeofday(&start, NULL);
	calcInverseMatrix(cartToFrac, fracToCart);
	//SafeAlloc(ts, 3, "ts");
	SafeAlloc(hNew, 3, "hNew");

	//matrixMultiply(cartToFrac, symmetryOperation.translationVector, ts);
	matrixMultiply(rotationMatrix, fracToCart, RO);
	matrixMultiply(symmetryOperation.rotationMat, RO, MsRO);
	matrixMultiply(cartToFrac, MsRO, DMsRO);
	invDX=1.0/xray.hBinSize;
	invDY=1.0/xray.kBinSize;
	invDZ=1.0/xray.lBinSize;
	for (int i=0;i<npoint;i++)
	{
		matrixMultiplyFast(xray.miller[i].Pos, DMsRO, hNew);
		//matrixMultiply(xray.miller[i].Pos, DMsRO, hNew);
		if (params.FastTranslationRotationInterpolation=="Linear")
		{
			phaseAndInterpolate(xray, hNew, invDX, invDY, invDZ, realOut, imagOut);
			//calcScatteringContinuous(Protein, xray, params, hNew[X], hNew[Y], hNew[Z], realOut, imagOut); 
			//cout <<"realOut= "<<realOut<<" imagOut= "<<imagOut<<endl;
			//getttimeofday(&end, NULL);
			//cout <<"phaseAndInterpolate took "<<calcTimeDiff(start, end)<<endl;
		}
		//else if (params.FastTranslationRotationInterpolation=="Hessian")
		//{
		//phaseAndInterpolateHessian(xray, hNew, dp, realOut, imagOut, Protein, params);
		//}
		else
		{
			string errorStr="Unrecognized ";
			errorStr+="FastTranslationRotationInterpolation ";
			errorStr+=params.FastTranslationRotationInterpolation;
			errorStr+="\nacceptable values are Linear and Hessian";
			error(errorStr, __LINE__, __FILE__);
		}
		//xray.real[i]+=realOut;
		//xray.imag[i]+=imagOut;
		real[i]=realOut;
		imag[i]=imagOut;
	}
	//gettimeofday(&end, NULL);
	//cout <<"calcFastTranslationRotation took "<<calcTimeDiff(start, end)<<endl;
	//calcIntensity(xray);
	//calcPhase(xray);
	//cout <<"Leaving calcFastRotationTranslation"<<endl;
	//endProgram(__LINE__, __FILE__);
}

void calcTransformedMiller(PosStruct &miller, Matrix &rotationMatrix, Matrix &cartToFrac, Matrix &fracToCart, Vector &hNew)
{
	Matrix RO, DRO;

	matrixMultiply(rotationMatrix, fracToCart, RO);
	matrixMultiply(cartToFrac, RO, DRO);
	matrixMultiply(miller.Pos, DRO, hNew);
}

void findMaxMillerForFastRotation(XRayStruct &xray, Real theta, Real phi, Real psi, PosStruct &maxMiller, XRayParamStruct &params, int bestIndex)
{
	Vector hNew;
	Matrix fracToCart, cartToFrac, rotationMatrix;
	ProteinStruct tempProtein;
	//cout <<"In findMaxMillerForFastRotation"<<endl;
	calcUnitCellVectors(tempProtein, params);
	calcFracToCart(tempProtein.a, tempProtein.b, tempProtein.c, fracToCart);
	calcInverseMatrix(fracToCart, cartToFrac);
	calcRotationMatrix(rotationMatrix, theta, phi, psi);
	calcTransformedMiller(xray.miller[bestIndex], rotationMatrix, cartToFrac, fracToCart, hNew);
	if (abs(hNew[X])>maxMiller.Pos[X])
	{
		maxMiller.Pos[X]=abs(hNew[X]);
	}
	if (abs(hNew[Y])>maxMiller.Pos[Y])
	{
		maxMiller.Pos[Y]=abs(hNew[Y]);
	}
	if (abs(hNew[Z])>maxMiller.Pos[Z])
	{
		maxMiller.Pos[Z]=abs(hNew[Z]);
	}
}

void findMaxMillerForFastRotationSystematic(XRayStruct &xray, XRayParamStruct &params)
{
	int bestIndex;
	Real inc;
	PosStruct maxMiller;
	inc=pi/10.0;
	//inc=pi/100.0;
	SafeAlloc(maxMiller.Pos, 3, "maxMiller");
	findMaxMillerForFastRotationAnalytical2(xray, params, bestIndex);
	for (Real theta=0;theta<=pi;theta+=inc)
	{
		for (Real phi=0;phi<=pi;phi+=inc)
		{
			for (Real psi=0;psi<=pi;psi+=inc)
			{
				findMaxMillerForFastRotation(xray, theta, phi, psi, maxMiller, params, bestIndex);
			}
		}
	}
	cout <<"maxMiller.Pos[X]= "<<maxMiller.Pos[X]<<endl;
	cout <<"maxMiller.Pos[Y]= "<<maxMiller.Pos[Y]<<endl;
	cout <<"maxMiller.Pos[Z]= "<<maxMiller.Pos[Z]<<endl;
	xray.maxContinuousH=int(maxMiller.Pos[X]+2.0);
	xray.maxContinuousK=int(maxMiller.Pos[Y]+2.0);
	xray.maxContinuousL=int(maxMiller.Pos[Z]+2.0);
}

Real findMaxMillerForFastRotationAnalytical2(XRayStruct &xray, XRayParamStruct &params, int &bestIndex)
{
	int nmiller=xray.miller.size();
	Real max, dist;
	Vector v;
	Matrix fracToCart, cartToFrac;
	ProteinStruct tempProtein;

	calcUnitCellVectors(tempProtein, params);
	calcFracToCart(tempProtein.a, tempProtein.b, tempProtein.c, fracToCart);
	calcInverseMatrix(fracToCart, cartToFrac);

	matrixMultiply(xray.miller[0].Pos, cartToFrac, v);
	max=v[X]*v[X]+v[Y]*v[Y]+v[Z]*v[Z];
	bestIndex=0;
	for (int i=0;i<nmiller;i++)
	{
		matrixMultiply(xray.miller[i].Pos, cartToFrac, v);
		dist=v[X]*v[X]+v[Y]*v[Y]+v[Z]*v[Z];
		if (dist>max) 
		{
			max=dist;
			bestIndex=i;
		}
	}
	return max;
}

void findMaxMillerForFastRotationAnalytical(XRayStruct &xray, XRayParamStruct &params)
{
	int bestIndex;
	Real max;
	Real maxH, maxK, maxL;

	//max=Real(xray.maxH)*Real(xray.maxH)/(params.a*params.a);
	//max+=Real(xray.maxK)*Real(xray.maxK)/(params.b*params.b);
	//max+=Real(xray.maxL)*Real(xray.maxL)/(params.c*params.c);
	max=findMaxMillerForFastRotationAnalytical2(xray, params, bestIndex);
	maxH=sqrt(max*params.a*params.a);
	maxK=sqrt(max*params.b*params.b);
	maxL=sqrt(max*params.c*params.c);
	cout <<"maxH= "<<maxH<<endl;
	cout <<"maxK= "<<maxK<<endl;
	cout <<"maxL= "<<maxL<<endl;
	xray.maxContinuousH=int(maxH+2.0);
	xray.maxContinuousK=int(maxK+2.0);
	xray.maxContinuousL=int(maxL+2.0);
}

void findMaxMillerForFastRotation(XRayStruct &xray, XRayParamStruct &params)
{
	Real hBinSize, kBinSize, lBinSize;
	timeval start, end;

	cout <<"In findMaxMillerForFastRotation"<<endl;

	gettimeofday(&start, NULL);
	if (params.alpha!=90.0 || params.beta!=90.0 || params.gamma!=90.0)
	{
		cout <<"Before findMaxMillerForFastRotationSystematic"<<endl;
		findMaxMillerForFastRotationSystematic(xray, params);
		cout <<"After findMaxMillerForFastRotationSystematic"<<endl;
	}
	else
	{
		cout <<"Before findMaxMillerForFastRotationAnalytical"<<endl;
		findMaxMillerForFastRotationAnalytical(xray, params);
		cout <<"After findMaxMillerForFastRotationAnalytical"<<endl;
	}
	hBinSize=Real(xray.maxContinuousH)/Real(params.ContinuousHValues);
	kBinSize=Real(xray.maxContinuousK)/Real(params.ContinuousKValues);
	lBinSize=Real(xray.maxContinuousL)/Real(params.ContinuousLValues);
	xray.maxContinuousH=int(Real(xray.maxContinuousH)+hBinSize);
	xray.maxContinuousK=int(Real(xray.maxContinuousK)+kBinSize);
	xray.maxContinuousL=int(Real(xray.maxContinuousL)+lBinSize);
	gettimeofday(&end, NULL);
	cout <<"findMaxMillerForFastRotation took "<<calcTimeDiff(start, end)<<endl;
}

void findIdenticalRotations(vector<SymmetryOperationStruct> &symOps, vector< vector<int> > &identicalRotations)
{
	int nsym=symOps.size();

	Safe2DAlloc(identicalRotations, nsym, 0, "identicalRotations");
	for (int i=0;i<nsym;i++)
	{
		for (int j=i+1;j<nsym;j++)
		{
			if (matrixEqual(symOps[i].rotationMat, symOps[j].rotationMat))
			{
				SafePushBack(identicalRotations[j], i, "identicalRotations");
			}
		}
	}
}

void calcFastTranslationRotation(ProteinStruct &Protein, XRayStruct &xray, Matrix &rotationMatrix, int prot, XRayParamStruct &params)
{
	int noperations=Protein.symmetryOperations.size();
	//int nsym=xray.complexSymReal.size();
	int nidentical;
	int Size1, Size2, Size3, npoint, nprot;
	vector< vector<int> > identicalRotations;
	Matrix tempRotation;
	Vector tempReal, tempImag, ts;
	//timeval start, end;
	//cout <<"translation.size= "<<translation.size()<<" center.size= "<<Protein.center.size()<<endl;
	//cout <<"In calcFastTranslationRotation"<<endl;
	//cout <<"translation.size= "<<translation.size()<<" center.size= "<<Protein.center.size()<<endl;

	Safe2DAlloc(tempRotation, 3, 3, "tempRotation");
	//cout <<"translation.size= "<<translation.size()<<" center.size= "<<Protein.center.size()<<endl;
	//translation[X]+=Protein.center[X];
	//translation[Y]+=Protein.center[Y];
	//translation[Z]+=Protein.center[Z];
	//cout <<"Before Zero"<<endl;
	zeroVector(xray.real);
	zeroVector(xray.imag);
	//cout <<"Before calcCartToFrac"<<endl;
	calcCartToFrac(Protein);
	//matrixMultiply(translation, Protein.cartToFrac, x);
	if (noperations==0) 
	{
		setSpaceGroupP1(Protein);
		noperations=1;
	}
/*
	if (nsym<=prot)
	{
		string errorStr="nsym= "+toStr(nsym)+" prot= "+toStr(prot);
		error(errorStr, __LINE__, __FILE__);
	}
*/
	Get3DVectorSize(xray.complexSymReal, Size1, Size2, Size3);
	if (Size1==0 || Size2==0 || Size3==0)
	{
		nprot=params.NumCopies;
		npoint=xray.miller.size();
		Safe3DAlloc(xray.complexSymReal, nprot, noperations, npoint);
		Safe3DAlloc(xray.complexSymImag, nprot, noperations, npoint);
	}
	findIdenticalRotations(Protein.symmetryOperations, identicalRotations);
	for (int i=0;i<noperations;i++)
	{
		nidentical=identicalRotations[i].size();
		if (nidentical==0)
		{
			calcFastTranslationRotation(Protein.symmetryOperations[i], Protein.cartToFrac, rotationMatrix, xray, xray.complexSymReal[prot][i], xray.complexSymImag[prot][i], params, Protein);
			//xray.complexSymReal[prot][i]=xray.real;
			//xray.complexSymImag[prot][i]=xray.imag;
			//SafePushBack(xray.complexSymReal[prot], xray.real, "xray.complexSymReal");
			//SafePushBack(xray.complexSymImag[prot], xray.imag, "xray.complexSymImag");
			//cout <<"xray.real.size()= "<<xray.real.size()<<endl;
		}
		else
		{
			xray.complexSymReal[prot][i]=xray.complexSymReal[prot][identicalRotations[i][0]];
			xray.complexSymImag[prot][i]=xray.complexSymImag[prot][identicalRotations[i][0]];
			//SafePushBack(xray.complexSymReal[prot], xray.complexSymReal[prot][identicalRotations[i][0]], "xray.complexSymReal");
			//SafePushBack(xray.complexSymImag[prot], xray.complexSymImag[prot][identicalRotations[i][0]], "xray.complexSymImag");

		}
	}
	/*
	   for (int i=0;i<noperations;i++)
	   {
	   matrixMultiply(Protein.cartToFrac, Protein.symmetryOperations[i].translationVector, ts);
	   translateXRay(xray.lookUp, xray.complexSymReal[prot][i], xray.complexSymImag[prot][i], xray.miller, ts, tempReal, tempImag);
	   xray.complexSymReal[prot][i]=tempReal;
	   xray.complexSymImag[prot][i]=tempImag;
	   }
	 */
}

void zeroContinuousScattering(XRayStruct &xray)
{
	int maxH, maxK, maxL;

	Get3DVectorSize(xray.imagContinuous, maxH, maxK, maxL, "xray.imagContinuous");

	for (int h=0;h<maxH;h++)
	{
		for (int k=0;k<maxK;k++)
		{
			for (int l=0;l<maxL;l++)
			{
				xray.imagContinuous[h][k][l]=0;
				xray.realContinuous[h][k][l]=0;
			}
		}
	}
}

void sumScatteringContinuousSegments(XRayStruct &xray, Vector &degOfFreedom, int start)
{
	int nsegments=xray.realContinuousSegment.size();
	int maxH, maxK, maxL;
	timeval startt, endt;

	gettimeofday(&startt, NULL);
	zeroContinuousScattering(xray);
	Get3DVectorSize(xray.imagContinuous, maxH, maxK, maxL, "xray.imagContinuous");
	for (int i=0;i<nsegments;i++)
	{
		for (int h=0;h<maxH;h++)
		{
			for (int k=0;k<maxK;k++)
			{
				for (int l=0;l<maxL;l++)
				{
					xray.imagContinuous[h][k][l]+=xray.imagContinuousSegment[i][h][k][l]*degOfFreedom[i+start];
					xray.realContinuous[h][k][l]+=xray.realContinuousSegment[i][h][k][l]*degOfFreedom[i+start];
				}
			}
		}
	}
	gettimeofday(&endt, NULL);
	//cout <<"sumScatteringContinuousSegments took "<<calcTimeDiff(startt, endt)<<endl;
	//exit(EXIT_FAILURE);
}

void millerToBin(XRayStruct &xray, Vector &hNew, Real invDX, Real invDY, Real invDZ, int &hbin, int &kbin, int &lbin)
{
	if (hNew[X]<0)
	{
		hNew[X]=-hNew[X];
		hNew[Y]=-hNew[Y];
		hNew[Z]=-hNew[Z];
	}
	hbin=int(hNew[X]*invDX);
	kbin=int((hNew[Y]+Real(xray.maxContinuousK))*invDY);
	lbin=int((hNew[Z]+Real(xray.maxContinuousL))*invDZ);
}
/*
void sumScatteringContinuousSegments(XRayStruct &xray, Vector &degOfFreedom, int start, int hbin, int kbin, int lbin)
{
	int nsegments=xray.realContinuousSegment.size();
	for (int h=hbin;h<=hbin+1;h++)
	{
		for (int k=kbin;k<=kbin+1;k++)
		{
			for (int l=lbin;l<=lbin+1;l++)
			{
				xray.imagContinuous[h][k][l]=0;
				xray.realContinuous[h][k][l]=0;
			}
		}
	}
	for (int i=0;i<nsegments;i++)
	{
		for (int h=hbin;h<=hbin+1;h++)
		{
			for (int k=kbin;k<=kbin+1;k++)
			{
				for (int l=lbin;l<=lbin+1;l++)
				{
					xray.imagContinuous[h][k][l]+=xray.imagContinuousSegment[i][h][k][l]*degOfFreedom[i+start];
					xray.realContinuous[h][k][l]+=xray.realContinuousSegment[i][h][k][l]*degOfFreedom[i+start];
				}
			}
		}
	}
}
*/

void sumScatteringContinuousSegments(Array3D &continuous, const Array3D &continuousSegment, int hbin, int kbin, int lbin)
{
	for (int h=hbin;h<=hbin+1;h++)
	{
		for (int k=kbin;k<=kbin+1;k++)
		{
			for (int l=lbin;l<=lbin+1;l++)
			{
				continuous[h][k][l]+=continuousSegment[h][k][l];
			}
		}
	}
}

void sumScatteringContinuousSegments(XRayStruct &xray, Vector &degOfFreedom, int start, int hbin, int kbin, int lbin)
{
	int nsegments=xray.realContinuousSegment.size();
	for (int h=hbin;h<=hbin+1;h++)
	{
		for (int k=kbin;k<=kbin+1;k++)
		{
			for (int l=lbin;l<=lbin+1;l++)
			{
				xray.imagContinuous[h][k][l]=0;
				xray.realContinuous[h][k][l]=0;
			}
		}
	}
	for (int i=0;i<nsegments;i++)
	{
		if (degOfFreedom[i+start]==1.0)
		{
			for (int h=hbin;h<=hbin+1;h++)
			{
				for (int k=kbin;k<=kbin+1;k++)
				{
					for (int l=lbin;l<=lbin+1;l++)
					{
						xray.imagContinuous[h][k][l]+=xray.imagContinuousSegment[i][h][k][l];
						xray.realContinuous[h][k][l]+=xray.realContinuousSegment[i][h][k][l];
						//cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<" real= "<<xray.realContinuous[h][k][l]<<" realSegment= "<<xray.realContinuousSegment[i][h][k][l]<<endl;
					}
				}
			}
		}
	}
}

/*
void sumScatteringContinuousSegments(XRayStruct &xray, Vector &degOfFreedom, vector<bool> &binary, int start, int hbin, int kbin, int lbin)
{
	int nsegments=xray.realContinuousSegment.size();

	for (int i=0;i<nsegments;i++)
	{
		for (int h=hbin;h<=hbin+1;h++)
		{
			for (int k=kbin;k<=kbin+1;k++)
			{
				for (int l=lbin;l<=lbin+1;l++)
				{
					if (degOfFreedom[i+start]==1.0 && xray.updatedContinuous[i][h][k][l]==false)
					{
						xray.imagContinuous[h][k][l]+=xray.imagContinuousSegment[i][h][k][l];
						xray.realContinuous[h][k][l]+=xray.realContinuousSegment[i][h][k][l];
						xray.updatedContinuous[i][h][k][l]=true;
					}
					else if (degOfFreedom[i+start]==0.0 && xray.updatedContinuous[i][h][k][l]==true)
					{
                                                xray.imagContinuous[h][k][l]-=xray.imagContinuousSegment[i][h][k][l];
                                                xray.realContinuous[h][k][l]-=xray.realContinuousSegment[i][h][k][l];
                                                xray.updatedContinuous[i][h][k][l]=false;
					}
				}
			}
		}
	}
}
*/

void sumScatteringContinuousSegments(XRayStruct &xray, Vector &degOfFreedom, int occupancyNum, int start, int hbin, int kbin, int lbin)
{
	int nsegments=xray.realContinuousSegment.size();
	for (int h=hbin;h<=hbin+1;h++)
	{
		for (int k=kbin;k<=kbin+1;k++)
		{
			for (int l=lbin;l<=lbin+1;l++)
			{
				if (occupancyNum!=xray.updatedContinuous[h][k][l])
				{
					xray.imagContinuous[h][k][l]=0;
					xray.realContinuous[h][k][l]=0;
				}
			}
		}
	}
	for (int i=0;i<nsegments;i++)
	{
		if (degOfFreedom[i+start]==1.0)
		{
			for (int h=hbin;h<=hbin+1;h++)
			{
				for (int k=kbin;k<=kbin+1;k++)
				{
					for (int l=lbin;l<=lbin+1;l++)
					{
						if (occupancyNum!=xray.updatedContinuous[h][k][l])
						{
							xray.imagContinuous[h][k][l]+=xray.imagContinuousSegment[i][h][k][l];
							xray.realContinuous[h][k][l]+=xray.realContinuousSegment[i][h][k][l];
						}
					}
				}
			}
		}
	}
	for (int h=hbin;h<=hbin+1;h++)
	{
		for (int k=kbin;k<=kbin+1;k++)
		{
			for (int l=lbin;l<=lbin+1;l++)
			{
				xray.updatedContinuous[h][k][l]=occupancyNum;
			}
		}
	}
}

/*
void subtractScatteringContinuousSegment(XRayStruct &xray, int pick, int hbin, int kbin, int lbin)
{
	for (int h=hbin;h<=hbin+1;h++)
	{
		for (int k=kbin;k<=kbin+1;k++)
		{
			for (int l=lbin;l<=lbin+1;l++)
			{
				if (xray.updatedContinuous[pick][h][k][l]==true)
				{
					xray.imagContinuous[h][k][l]-=xray.imagContinuousSegment[pick][h][k][l];
					xray.realContinuous[h][k][l]-=xray.realContinuousSegment[pick][h][k][l];
					xray.updatedContinuous[pick][h][k][l]=false;
				}
			}
		}
	}
}

void addScatteringContinuousSegment(XRayStruct &xray, int pick, int hbin, int kbin, int lbin)
{
	//cout <<"pick= "<<pick<<" hbin= "<<hbin<<" kbin= "<<kbin<<" lbin= "<<lbin<<endl;
	for (int h=hbin;h<=hbin+1;h++)
	{
		for (int k=kbin;k<=kbin+1;k++)
		{
			for (int l=lbin;l<=lbin+1;l++)
			{
				if (xray.updatedContinuous[pick][h][k][l]==false)
				{
					xray.imagContinuous[h][k][l]+=xray.imagContinuousSegment[pick][h][k][l];
					xray.realContinuous[h][k][l]+=xray.realContinuousSegment[pick][h][k][l];
					xray.updatedContinuous[pick][h][k][l]=true;
				}
			}
		}
	}
}
*/

void subtractScatteringContinuousSegment(XRayStruct &xray, int pick, int occupancyNum, int hbin, int kbin, int lbin)
{
	for (int h=hbin;h<=hbin+1;h++)
	{
		for (int k=kbin;k<=kbin+1;k++)
		{
			for (int l=lbin;l<=lbin+1;l++)
			{
				//cout <<"occupancyNum= "<<occupancyNum<<" h= "<<h<<" k= "<<k<<" l= "<<l<<" update= "<<xray.updatedContinuous[h][k][l]<<endl;
				if (xray.updatedContinuous[h][k][l]!=occupancyNum)
				{
					xray.imagContinuous[h][k][l]-=xray.imagContinuousSegment[pick][h][k][l];
					xray.realContinuous[h][k][l]-=xray.realContinuousSegment[pick][h][k][l];
					xray.updatedContinuous[h][k][l]=occupancyNum;
				}
			}
		}
	}
}

void addScatteringContinuousSegment(XRayStruct &xray, int pick, int occupancyNum, int hbin, int kbin, int lbin)
{
	//cout <<"pick= "<<pick<<" hbin= "<<hbin<<" kbin= "<<kbin<<" lbin= "<<lbin<<endl;
	for (int h=hbin;h<=hbin+1;h++)
	{
		for (int k=kbin;k<=kbin+1;k++)
		{
			for (int l=lbin;l<=lbin+1;l++)
			{
				if (xray.updatedContinuous[h][k][l]!=occupancyNum)
				{
					xray.imagContinuous[h][k][l]+=xray.imagContinuousSegment[pick][h][k][l];
					xray.realContinuous[h][k][l]+=xray.realContinuousSegment[pick][h][k][l];
					xray.updatedContinuous[h][k][l]=occupancyNum;
				}
			}
		}
	}
}

void sumScatteringContinuousSegments(XRayStruct &xray, Vector &degOfFreedom, int start, int hbin, int kbin, int lbin, Matrix &real, Matrix &imag)
{
	int nsegments=xray.realContinuousSegment.size();

		for (int k=kbin;k<=kbin+1;k++)
		{
			for (int l=lbin;l<=lbin+1;l++)
			{
				imag[k][l]=0;
				real[k][l]=0;
			}
		}

	for (int i=0;i<nsegments;i++)
	{
			for (int k=kbin;k<=kbin+1;k++)
			{
				for (int l=lbin;l<=lbin+1;l++)
				{
					imag[k][l]+=xray.imagContinuousSegment[i][hbin][k][l]*degOfFreedom[i+start];
					real[k][l]+=xray.realContinuousSegment[i][hbin][k][l]*degOfFreedom[i+start];
				}
			}
	}
}

void sumScatteringContinuousSegmentsUnrolled(XRayStruct &xray, Vector &degOfFreedom, int start, int h, int k, int l, Real real[8], Real imag[8], Matrix &realSegment, Matrix &imagSegment)
{
	int nsegments=xray.realContinuousSegment.size();
	for (int i=0;i<nsegments;i++)
	{
		imag[0]=imagSegment[i][0]*degOfFreedom[i+start];
		imag[1]+=imagSegment[i][1]*degOfFreedom[i+start];
		imag[2]+=imagSegment[i][2]*degOfFreedom[i+start];
		imag[3]+=imagSegment[i][3]*degOfFreedom[i+start];
		imag[4]+=imagSegment[i][4]*degOfFreedom[i+start];
		imag[5]+=imagSegment[i][5]*degOfFreedom[i+start];
		imag[6]+=imagSegment[i][6]*degOfFreedom[i+start];
		imag[7]+=imagSegment[i][7]*degOfFreedom[i+start];
	}

	for (int i=0;i<nsegments;i++)
	{
		real[0]=realSegment[i][0]*degOfFreedom[i+start];
		real[1]+=realSegment[i][1]*degOfFreedom[i+start];
		real[2]+=realSegment[i][2]*degOfFreedom[i+start];
		real[3]+=realSegment[i][3]*degOfFreedom[i+start];
		real[4]+=realSegment[i][4]*degOfFreedom[i+start];
		real[5]+=realSegment[i][5]*degOfFreedom[i+start];
		real[6]+=realSegment[i][6]*degOfFreedom[i+start];
		real[7]+=realSegment[i][7]*degOfFreedom[i+start];
	}
}

void makeOccupanciesNonZero(Vector &degOfFreedom, int start)
{
	bool nonzero=false;
	int ndof=degOfFreedom.size();

	for (int i=start;i<ndof;i++)
	{
		if (degOfFreedom[i]==1.0) nonzero=true;
	}
	if (!nonzero) degOfFreedom[ndof-1]=1.0;
}

bool isDuplicate(vector<int> &v1, vector<int> &v2)
{
	int Size=v1.size();

	for (int i=0;i<Size;i++)
	{
		if (v1[i]!=v2[i]) return false;
	}
	return true;
}

bool isDuplicate(vector<int> &v1, vector< vector<int> > &v2)
{
	int Size=v2.size();

	for (int i=0;i<Size;i++)
	{
		if (isDuplicate(v1, v2[i])) return true;
	}
	return false;
}

int calcOccupancyNum(Vector &degOfFreedom, int start)
{
	int occupancyNum=0;
	int ndof=degOfFreedom.size();
	int pow=1;

	for (int i=start;i<ndof;i++)
	{
		if (degOfFreedom[i]==1.0) occupancyNum+=pow;
		pow*=2;
	}
	return occupancyNum;
}

void sumScatteringContinuousSegments(XRayStruct &xray, ProteinStruct &Protein, Vector &degOfFreedom, int start)
{
	//vector<bool> binary;
	int nmiller=xray.miller.size();
	int hbin, kbin, lbin;
	int prot=0;
	int nonredundant;
	//int nbinary;
	//int ndof=degOfFreedom.size();
	int occupancyNum=0;
	Real invDX, invDY, invDZ;
	Vector hNew;
	Matrix cartToFrac, fracToCart;
	Matrix rotationMatrix1, rotationMatrix2;
	Matrix RO, MsRO, DMsRO;

	//cout <<"In sumScatteringContinuousSegments"<<endl;
	//Safe3DArrayAlloc(imag, 2, 2, 2, "imag");
	//Safe3DArrayAlloc(real, 2, 2, 2, "real");
	invDX=1.0/xray.hBinSize;
	invDY=1.0/xray.kBinSize;
	invDZ=1.0/xray.lBinSize;
	cartToFrac=Protein.cartToFrac;
	fracToCart=Protein.fracToCart;

	//cout <<"Before calcRotationMatrix"<<endl;
	//degOfFreedomToRotMat(rotationMatrix1, degOfFreedom[prot*6+X], degOfFreedom[prot*6+Y], degOfFreedom[prot*6+Z]);	
	nonredundant=identifyRotationBlocks(Protein);
	makeOccupanciesNonZero(degOfFreedom, start);
	calcRotationMatrix(rotationMatrix1, degOfFreedom[prot*6+X], degOfFreedom[prot*6+Y], degOfFreedom[prot*6+Z]);	
	//nbinary=ndof-start+1;
	//SafeAlloc(binary, nbinary, "binary");
	occupancyNum=calcOccupancyNum(degOfFreedom, start);
	//for (int i=0;i<nbinary;i++)
	//{
	//	if (degOfFreedom[i+start]==1.0) binary[i]=true;
	//	else binary[i]=false;
	//}
	for (int i=0;i<nonredundant;i++)
	{	
		//cout <<"i= "<<i<<endl;
		matrixMultiply(rotationMatrix1, fracToCart, RO);
		matrixMultiply(Protein.symmetryOperations[i].rotationMat, RO, MsRO);
		matrixMultiply(cartToFrac, MsRO, DMsRO);
		for (int j=0;j<nmiller;j++)
		{
			//cout <<"j= "<<j<<endl;
			//calcTransformedMiller(xray.miller[j], rotationMatrix2, cartToFrac, fracToCart, hNew);
			matrixMultiply(xray.miller[j].Pos, DMsRO, hNew);
			//cout <<"Before millerToBin"<<endl;
			millerToBin(xray, hNew, invDX, invDY, invDZ, hbin, kbin, lbin);
			//cout <<"Before sumScatteringContinuousSegments(xray, degOfFreedom, start, hbin, kbin, lbin);"<<endl;
			//printVector(hNew, "hNew");
			//cout <<"hbin= "<<hbin<<" kbin= "<<kbin<<" lbin= "<<lbin<<endl;
			//**imag=&xray.imagContinuous[hbin][kbin][lbin];
			//**real=&xray.realContinuous[hbin][kbin][lbin];
			//sumScatteringContinuousSegments(xray, degOfFreedom, binary, start, hbin, kbin, lbin);
			sumScatteringContinuousSegments(xray, degOfFreedom, occupancyNum, start, hbin, kbin, lbin);
		}
	}
}

/*
void resetUpdateContinuous(XRayStruct &xray, ProteinStruct &Protein, Vector &degOfFreedom)
{
	int nmiller=xray.miller.size();
	int hbin, kbin, lbin;
	int prot=0;
	int nonredundant;
	Real invDX, invDY, invDZ;
	Vector hNew;
	Matrix cartToFrac, fracToCart;
	Matrix rotationMatrix1, rotationMatrix2;
	Matrix RO, MsRO, DMsRO;

	invDX=1.0/xray.hBinSize;
	invDY=1.0/xray.kBinSize;
	invDZ=1.0/xray.lBinSize;
	cartToFrac=Protein.cartToFrac;
	fracToCart=Protein.fracToCart;

	nonredundant=identifyRotationBlocks(Protein);
	calcRotationMatrix(rotationMatrix1, degOfFreedom[prot*6+X], degOfFreedom[prot*6+Y], degOfFreedom[prot*6+Z]);	
	for (int i=0;i<nonredundant;i++)
	{	
		matrixMultiply(rotationMatrix1, fracToCart, RO);
		matrixMultiply(Protein.symmetryOperations[i].rotationMat, RO, MsRO);
		matrixMultiply(cartToFrac, MsRO, DMsRO);
		for (int j=0;j<nmiller;j++)
		{
			matrixMultiply(xray.miller[j].Pos, DMsRO, hNew);
			millerToBin(xray, hNew, invDX, invDY, invDZ, hbin, kbin, lbin);
			for (int h=hbin;h<=hbin+1;h++)
			{
				for (int k=kbin;k<=kbin+1;k++)
				{
					for (int l=lbin;l<=lbin+1;l++)
					{
						xray.updatedContinuous[h][k][l]=false;
					}
				}
			}
		}
	}
}
*/

void calcScatteringContinuousSegmentsChange(XRayStruct &xray, ProteinStruct &Protein, Vector &degOfFreedom, int start, int pick)
{
	int nmiller=xray.miller.size();
	int hbin, kbin, lbin;
	int prot=0;
	int nonredundant;
	int ndof=degOfFreedom.size();
	int occupancyNum;
	Real invDX, invDY, invDZ;
	Vector hNew;
	Matrix cartToFrac, fracToCart;
	Matrix rotationMatrix1, rotationMatrix2;
	Matrix RO, MsRO, DMsRO;

	invDX=1.0/xray.hBinSize;
	invDY=1.0/xray.kBinSize;
	invDZ=1.0/xray.lBinSize;
	cartToFrac=Protein.cartToFrac;
	fracToCart=Protein.fracToCart;

	nonredundant=identifyRotationBlocks(Protein);
	if (pick==ndof-1)
	{
		if (degOfFreedom[pick]==0) degOfFreedom[pick]=1.0;
		else degOfFreedom[pick]=1.0;
		return;
	}
	//makeOccupanciesNonZero(degOfFreedom, start);
	occupancyNum=calcOccupancyNum(degOfFreedom, start);
	calcRotationMatrix(rotationMatrix1, degOfFreedom[prot*6+X], degOfFreedom[prot*6+Y], degOfFreedom[prot*6+Z]);	
	for (int i=0;i<nonredundant;i++)
	{	
		matrixMultiply(rotationMatrix1, fracToCart, RO);
		matrixMultiply(Protein.symmetryOperations[i].rotationMat, RO, MsRO);
		matrixMultiply(cartToFrac, MsRO, DMsRO);
		for (int j=0;j<nmiller;j++)
		{
			matrixMultiply(xray.miller[j].Pos, DMsRO, hNew);
			millerToBin(xray, hNew, invDX, invDY, invDZ, hbin, kbin, lbin);
			if (degOfFreedom[pick]==0)
			{
				//cout <<"Before subtractScatteringContinuousSegment"<<endl;
				subtractScatteringContinuousSegment(xray, pick-start, occupancyNum, hbin, kbin, lbin);
				//cout <<"After subtractScatteringContinuousSegment"<<endl;
			}
			else
			{
				//cout <<"Before addScatteringContinuousSegment"<<endl;
				addScatteringContinuousSegment(xray, pick-start, occupancyNum, hbin, kbin, lbin);
				//cout <<"After addScatteringContinuousSegment"<<endl;
			}
			//sumScatteringContinuousSegments(xray, degOfFreedom, start, hbin, kbin, lbin);
		}
	}
	//resetUpdateContinuous(xray, Protein, degOfFreedom);
}

void calcFastTranslationRotation(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayParamStruct &params)
{
	int nprot=Proteins.size();
	int nsym, npoint;
	int Size1, Size2, Size3;


	Get3DVectorSize(xray.complexSymReal, Size1, Size2, Size3);
	if (Size1==0 || Size2==0 || Size3==0)
	{
		nsym=Proteins[0].symmetryOperations.size();
		npoint=xray.miller.size();
		Safe3DAlloc(xray.complexSymReal, nprot, nsym, npoint, "complexSymReal");
		Safe3DAlloc(xray.complexSymImag, nprot, nsym, npoint, "complexSymImag");
		Safe3DAlloc(xray.complexSymTempReal, nprot, nsym, npoint, "complexSymTempReal");
		Safe3DAlloc(xray.complexSymTempImag, nprot, nsym, npoint, "complexSymTempImag");

	}

	for (int i=0;i<nprot;i++)
	{
		calcFastTranslationRotation(Proteins[i], xray, Proteins[i].rotationMatrix, i, params);
	}
}

void translateXRay(LookUpStruct &lookUp, Vector &symReal, Vector &symImag, vector<PosStruct> &miller, Vector &frac, Vector &symTempReal, Vector &symTempImag)
{
	int npoint=symReal.size();
	int ntemp=symTempReal.size();
	Real dp;
	Real cosdp, sindp;
	
	if (ntemp!=npoint)
	{
		SafeAlloc(symTempReal, npoint, "npoint");
		SafeAlloc(symTempImag, npoint, "npoint");
	}
	frac[X]*=2.0*pi;
	frac[Y]*=2.0*pi;
	frac[Z]*=2.0*pi;
	for (int i=0;i<npoint;i++)
	{
		dp=calcDotProduct(miller[i].Pos, frac);
		getTrigLookUp(lookUp, dp, cosdp, sindp);
		symTempReal[i]=cosdp*symReal[i]-sindp*symImag[i];
		symTempImag[i]=cosdp*symImag[i]+sindp*symReal[i];
	}
}

void translateXRayDecomposed(XRayStruct &xray, Vector &symReal, Vector &symImag, vector<PosStruct> &miller, Vector &frac, Vector &symTempReal, Vector &symTempImag)
{
	int Size=xray.cosx1d.size();
	int h, k, l;
	int npoint=symReal.size();
	int ntemp=symTempReal.size();
	Real x, y, z;
	Real cosdp, sindp;
	//Real dp, cosdp1, sindp1;
	
	if (ntemp!=npoint)
	{
		SafeAlloc(symTempReal, npoint, "npoint");
		SafeAlloc(symTempImag, npoint, "npoint");
	}
	frac[X]*=2.0*pi;
	frac[Y]*=2.0*pi;
	frac[Z]*=2.0*pi;

	if (Size==0)
	{
		setMinMiller(xray);
		setMaxMiller(xray);
		SafeAlloc(xray.cosx1d, xray.maxH-xray.minH+1, "cosx");
		SafeAlloc(xray.cosy1d, xray.maxK-xray.minK+1, "cosy");
		SafeAlloc(xray.cosz1d, xray.maxL-xray.minL+1, "cosz");
		SafeAlloc(xray.sinx1d, xray.maxH-xray.minH+1, "sinx");
		SafeAlloc(xray.siny1d, xray.maxK-xray.minK+1, "siny");
		SafeAlloc(xray.sinz1d, xray.maxL-xray.minL+1, "sinz");
	}
	//cout <<"minH= "<<xray.minH<<endl;
	//cout <<"maxH= "<<xray.maxH<<endl;

	for (int i=xray.minH;i<=xray.maxH;i++)
	{
		x=Real(i)*frac[X];
		getTrigLookUp(xray.lookUp, x, xray.cosx1d[i-xray.minH], xray.sinx1d[i-xray.minH]);
		//cout <<"x= "<<x<<" cosx1d= "<<xray.cosx1d[i-xray.minH]<<" sinx1d= "<<xray.sinx1d[i-xray.minH]<<endl;
	}
	for (int i=xray.minK;i<=xray.maxK;i++)
	{
		y=Real(i)*frac[Y];
		getTrigLookUp(xray.lookUp, y, xray.cosy1d[i-xray.minK], xray.siny1d[i-xray.minK]);
		//cout <<"y= "<<y<<" cosy1d= "<<xray.cosy1d[i-xray.minH]<<" siny1d= "<<xray.siny1d[i-xray.minH]<<endl;
	}
	for (int i=xray.minL;i<=xray.maxL;i++)
	{
		z=Real(i)*frac[Z];
		getTrigLookUp(xray.lookUp, z, xray.cosz1d[i-xray.minL], xray.sinz1d[i-xray.minL]);
		//cout <<"z= "<<z<<" cosz1d= "<<xray.cosz1d[i-xray.minH]<<" sinz1d= "<<xray.sinz1d[i-xray.minH]<<endl;
	}

	for (int i=0;i<npoint;i++)
	{
		//dp=calcDotProduct(miller[i].Pos, frac);
		//getTrigLookUp(xray.lookUp, dp, cosdp1, sindp1);
		h=int(xray.miller[i].Pos[X]-xray.minH);
		k=int(xray.miller[i].Pos[Y]-xray.minK);
		l=int(xray.miller[i].Pos[Z]-xray.minL);
		
		//cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<endl;

/*
		cout <<"cosx1d= "<<xray.cosx1d[h]<<endl;
		cout <<"cosy1d= "<<xray.cosy1d[k]<<endl;
		cout <<"cosz1d= "<<xray.cosz1d[l]<<endl;

		cout <<"sinx1d= "<<xray.sinx1d[h]<<endl;
		cout <<"siny1d= "<<xray.siny1d[k]<<endl;
		cout <<"sinz1d= "<<xray.sinz1d[l]<<endl;
*/
		
		//cosdp=xray.cosx1d[h]*xray.cosy1d[k]*xray.cosz1d[l];
		//cosdp-=xray.cosx1d[h]*xray.siny1d[k]*xray.sinz1d[l];
		//cosdp-=xray.sinx1d[h]*xray.cosy1d[k]*xray.sinz1d[l];
		//cosdp-=xray.sinx1d[h]*xray.siny1d[k]*xray.cosz1d[l];

		cosdp=xray.cosx1d[h]*(xray.cosy1d[k]*xray.cosz1d[l]-xray.siny1d[k]*xray.sinz1d[l]);
		cosdp-=xray.sinx1d[h]*(xray.cosy1d[k]*xray.sinz1d[l]+xray.siny1d[k]*xray.cosz1d[l]);

		//sindp=-xray.sinx1d[h]*xray.siny1d[k]*xray.sinz1d[l];
		//sindp+=xray.sinx1d[h]*xray.cosy1d[k]*xray.cosz1d[l];
		//sindp+=xray.cosx1d[h]*xray.siny1d[k]*xray.cosz1d[l];
		//sindp+=xray.cosx1d[h]*xray.cosy1d[k]*xray.sinz1d[l];

		sindp=xray.sinx1d[h]*(xray.cosy1d[k]*xray.cosz1d[l]-xray.siny1d[k]*xray.sinz1d[l]);
		sindp+=xray.cosx1d[h]*(xray.cosy1d[k]*xray.sinz1d[l]+xray.siny1d[k]*xray.cosz1d[l]);

		symTempReal[i]=cosdp*symReal[i]-sindp*symImag[i];
		symTempImag[i]=cosdp*symImag[i]+sindp*symReal[i];


		//cout <<"cosdp= "<<cosdp<<" sindp= "<<sindp<<" symReal= "
		//<<symReal[i]<<" symImag= "<<symImag[i]<<" symTempReal= "
		//<<symTempReal[i]<<" symTempImag= "<<symTempImag[i]<<endl;
	}
}

void translateXRayDecomposed2(XRayStruct &xray, Vector &symReal, Vector &symImag, vector<PosStruct> &miller, Vector &frac, Vector &symTempReal, Vector &symTempImag)
{
	int Size=xray.cosx1d.size();
	int h, k, l;
	int npoint=symReal.size();
	int ntemp=symTempReal.size();
	Real x, y, z;
	Real cosdp, sindp;
	//Real dp, cosdp1, sindp1;
	
	if (ntemp!=npoint)
	{
		SafeAlloc(symTempReal, npoint, "npoint");
		SafeAlloc(symTempImag, npoint, "npoint");
	}
	frac[X]*=2.0*pi;
	frac[Y]*=2.0*pi;
	frac[Z]*=2.0*pi;

	if (Size==0)
	{
		setMinMiller(xray);
		setMaxMiller(xray);
		SafeAlloc(xray.cosx1d, xray.maxH-xray.minH+1, "cosx");
		SafeAlloc(xray.cosy1d, xray.maxK-xray.minK+1, "cosy");
		SafeAlloc(xray.cosz1d, xray.maxL-xray.minL+1, "cosz");
		SafeAlloc(xray.sinx1d, xray.maxH-xray.minH+1, "sinx");
		SafeAlloc(xray.siny1d, xray.maxK-xray.minK+1, "siny");
		SafeAlloc(xray.sinz1d, xray.maxL-xray.minL+1, "sinz");
	}
	//cout <<"minH= "<<xray.minH<<endl;
	//cout <<"maxH= "<<xray.maxH<<endl;

	for (int i=xray.minH;i<=xray.maxH;i++)
	{
		x=Real(i)*frac[X];
		getTrigLookUp(xray.lookUp, x, xray.cosx1d[i-xray.minH], xray.sinx1d[i-xray.minH]);
		//cout <<"x= "<<x<<" cosx1d= "<<xray.cosx1d[i-xray.minH]<<" sinx1d= "<<xray.sinx1d[i-xray.minH]<<endl;
	}
	for (int i=xray.minK;i<=xray.maxK;i++)
	{
		y=Real(i)*frac[Y];
		getTrigLookUp(xray.lookUp, y, xray.cosy1d[i-xray.minK], xray.siny1d[i-xray.minK]);
		//cout <<"y= "<<y<<" cosy1d= "<<xray.cosy1d[i-xray.minH]<<" siny1d= "<<xray.siny1d[i-xray.minH]<<endl;
	}
	for (int i=xray.minL;i<=xray.maxL;i++)
	{
		z=Real(i)*frac[Z];
		getTrigLookUp(xray.lookUp, z, xray.cosz1d[i-xray.minL], xray.sinz1d[i-xray.minL]);
		//cout <<"z= "<<z<<" cosz1d= "<<xray.cosz1d[i-xray.minH]<<" sinz1d= "<<xray.sinz1d[i-xray.minH]<<endl;
	}

	for (int i=0;i<npoint;i++)
	{
		//dp=calcDotProduct(miller[i].Pos, frac);
		//getTrigLookUp(xray.lookUp, dp, cosdp1, sindp1);
		h=int(xray.miller[i].Pos[X]-xray.minH);
		k=int(xray.miller[i].Pos[Y]-xray.minK);
		l=int(xray.miller[i].Pos[Z]-xray.minL);
		
		//cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<endl;

/*
		cout <<"cosx1d= "<<xray.cosx1d[h]<<endl;
		cout <<"cosy1d= "<<xray.cosy1d[k]<<endl;
		cout <<"cosz1d= "<<xray.cosz1d[l]<<endl;

		cout <<"sinx1d= "<<xray.sinx1d[h]<<endl;
		cout <<"siny1d= "<<xray.siny1d[k]<<endl;
		cout <<"sinz1d= "<<xray.sinz1d[l]<<endl;
*/
		
		cosdp=xray.cosx1d[h]*xray.cosy1d[k]*xray.cosz1d[l];
		cosdp-=xray.cosx1d[h]*xray.siny1d[k]*xray.sinz1d[l];
		cosdp-=xray.sinx1d[h]*xray.cosy1d[k]*xray.sinz1d[l];
		cosdp-=xray.sinx1d[h]*xray.siny1d[k]*xray.cosz1d[l];

		sindp=-xray.sinx1d[h]*xray.siny1d[k]*xray.sinz1d[l];
		sindp+=xray.sinx1d[h]*xray.cosy1d[k]*xray.cosz1d[l];
		sindp+=xray.cosx1d[h]*xray.siny1d[k]*xray.cosz1d[l];
		sindp+=xray.cosx1d[h]*xray.cosy1d[k]*xray.sinz1d[l];

		symTempReal[i]=cosdp*symReal[i]-sindp*symImag[i];
		symTempImag[i]=cosdp*symImag[i]+sindp*symReal[i];


		//cout <<"cosdp= "<<cosdp<<" sindp= "<<sindp<<" symReal= "
		//<<symReal[i]<<" symImag= "<<symImag[i]<<" symTempReal= "
		//<<symTempReal[i]<<" symTempImag= "<<symTempImag[i]<<endl;
	}
}

void translateXRayDecomposed3(XRayStruct &xray, Vector &symReal, Vector &symImag, vector<PosStruct> &miller, Vector &frac, Vector &symTempReal, Vector &symTempImag)
{
	int Size=xray.cosx1d.size();
	int h, k, l;
	int npoint=symReal.size();
	int ntemp=symTempReal.size();
	Real a, b;
	Real x, y, z;
	Real cosdp, sindp;
	//Real dp, cosdp1, sindp1;
	
	if (ntemp!=npoint)
	{
		SafeAlloc(symTempReal, npoint, "npoint");
		SafeAlloc(symTempImag, npoint, "npoint");
	}
	frac[X]*=2.0*pi;
	frac[Y]*=2.0*pi;
	frac[Z]*=2.0*pi;

	if (Size==0)
	{
		setMinMiller(xray);
		setMaxMiller(xray);
		SafeAlloc(xray.cosx1d, xray.maxH-xray.minH+1, "cosx");
		SafeAlloc(xray.cosy1d, xray.maxK-xray.minK+1, "cosy");
		SafeAlloc(xray.cosz1d, xray.maxL-xray.minL+1, "cosz");
		SafeAlloc(xray.sinx1d, xray.maxH-xray.minH+1, "sinx");
		SafeAlloc(xray.siny1d, xray.maxK-xray.minK+1, "siny");
		SafeAlloc(xray.sinz1d, xray.maxL-xray.minL+1, "sinz");
	}
	//cout <<"minH= "<<xray.minH<<endl;
	//cout <<"maxH= "<<xray.maxH<<endl;

	for (int i=xray.minH;i<=xray.maxH;i++)
	{
		x=Real(i)*frac[X];
		getTrigLookUp(xray.lookUp, x, xray.cosx1d[i-xray.minH], xray.sinx1d[i-xray.minH]);
		//cout <<"x= "<<x<<" cosx1d= "<<xray.cosx1d[i-xray.minH]<<" sinx1d= "<<xray.sinx1d[i-xray.minH]<<endl;
	}
	for (int i=xray.minK;i<=xray.maxK;i++)
	{
		y=Real(i)*frac[Y];
		getTrigLookUp(xray.lookUp, y, xray.cosy1d[i-xray.minK], xray.siny1d[i-xray.minK]);
		//cout <<"y= "<<y<<" cosy1d= "<<xray.cosy1d[i-xray.minH]<<" siny1d= "<<xray.siny1d[i-xray.minH]<<endl;
	}
	for (int i=xray.minL;i<=xray.maxL;i++)
	{
		z=Real(i)*frac[Z];
		getTrigLookUp(xray.lookUp, z, xray.cosz1d[i-xray.minL], xray.sinz1d[i-xray.minL]);
		//cout <<"z= "<<z<<" cosz1d= "<<xray.cosz1d[i-xray.minH]<<" sinz1d= "<<xray.sinz1d[i-xray.minH]<<endl;
	}

	for (int i=0;i<npoint;i++)
	{
		//dp=calcDotProduct(miller[i].Pos, frac);
		//getTrigLookUp(xray.lookUp, dp, cosdp1, sindp1);
		h=int(xray.miller[i].Pos[X]-xray.minH);
		k=int(xray.miller[i].Pos[Y]-xray.minK);
		l=int(xray.miller[i].Pos[Z]-xray.minL);
		
		//cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<endl;

/*
		cout <<"cosx1d= "<<xray.cosx1d[h]<<endl;
		cout <<"cosy1d= "<<xray.cosy1d[k]<<endl;
		cout <<"cosz1d= "<<xray.cosz1d[l]<<endl;

		cout <<"sinx1d= "<<xray.sinx1d[h]<<endl;
		cout <<"siny1d= "<<xray.siny1d[k]<<endl;
		cout <<"sinz1d= "<<xray.sinz1d[l]<<endl;
*/
		
		a=xray.cosy1d[k]*xray.cosz1d[l]-xray.siny1d[k]*xray.sinz1d[l];
		b=xray.cosy1d[k]*xray.sinz1d[l]+xray.siny1d[k]*xray.cosz1d[l];

		cosdp=xray.cosx1d[h]*a;
		cosdp-=xray.sinx1d[h]*b;

		sindp=xray.sinx1d[h]*a;
		sindp+=xray.cosx1d[h]*b;

		symTempReal[i]=cosdp*symReal[i]-sindp*symImag[i];
		symTempImag[i]=cosdp*symImag[i]+sindp*symReal[i];


		//cout <<"cosdp= "<<cosdp<<" sindp= "<<sindp<<" symReal= "
		//<<symReal[i]<<" symImag= "<<symImag[i]<<" symTempReal= "
		//<<symTempReal[i]<<" symTempImag= "<<symTempImag[i]<<endl;
	}
}

void translateXRayDecomposed4(XRayStruct &xray, Vector &symReal, Vector &symImag, vector<PosStruct> &miller, Vector &frac, Vector &symTempReal, Vector &symTempImag)
{
	int Size=xray.cosx1d.size();
	int h, k, l;
	int npoint=symReal.size();
	int ntemp=symTempReal.size();
	Real x, y, z;
	Real cosdp, sindp;
	Matrix matTrig1, matTrig2;
	//Real dp, cosdp1, sindp1;
	
	if (ntemp!=npoint)
	{
		SafeAlloc(symTempReal, npoint, "npoint");
		SafeAlloc(symTempImag, npoint, "npoint");
	}
	frac[X]*=2.0*pi;
	frac[Y]*=2.0*pi;
	frac[Z]*=2.0*pi;

	if (Size==0)
	{
		setMinMiller(xray);
		setMaxMiller(xray);
		SafeAlloc(xray.cosx1d, xray.maxH-xray.minH+1, "cosx");
		SafeAlloc(xray.cosy1d, xray.maxK-xray.minK+1, "cosy");
		SafeAlloc(xray.cosz1d, xray.maxL-xray.minL+1, "cosz");
		SafeAlloc(xray.sinx1d, xray.maxH-xray.minH+1, "sinx");
		SafeAlloc(xray.siny1d, xray.maxK-xray.minK+1, "siny");
		SafeAlloc(xray.sinz1d, xray.maxL-xray.minL+1, "sinz");
	}
	
	Safe2DAlloc(matTrig1, xray.maxK-xray.minK+1, xray.maxL-xray.minL+1, "matTrig1");
	Safe2DAlloc(matTrig2, xray.maxK-xray.minK+1, xray.maxL-xray.minL+1, "matTrig2");

	//cout <<"minH= "<<xray.minH<<endl;
	//cout <<"maxH= "<<xray.maxH<<endl;

	for (int i=xray.minH;i<=xray.maxH;i++)
	{
		x=Real(i)*frac[X];
		getTrigLookUp(xray.lookUp, x, xray.cosx1d[i-xray.minH], xray.sinx1d[i-xray.minH]);
		//cout <<"x= "<<x<<" cosx1d= "<<xray.cosx1d[i-xray.minH]<<" sinx1d= "<<xray.sinx1d[i-xray.minH]<<endl;
	}
	for (int i=xray.minK;i<=xray.maxK;i++)
	{
		y=Real(i)*frac[Y];
		getTrigLookUp(xray.lookUp, y, xray.cosy1d[i-xray.minK], xray.siny1d[i-xray.minK]);
		//cout <<"y= "<<y<<" cosy1d= "<<xray.cosy1d[i-xray.minH]<<" siny1d= "<<xray.siny1d[i-xray.minH]<<endl;
	}
	for (int i=xray.minL;i<=xray.maxL;i++)
	{
		z=Real(i)*frac[Z];
		getTrigLookUp(xray.lookUp, z, xray.cosz1d[i-xray.minL], xray.sinz1d[i-xray.minL]);
		//cout <<"z= "<<z<<" cosz1d= "<<xray.cosz1d[i-xray.minH]<<" sinz1d= "<<xray.sinz1d[i-xray.minH]<<endl;
	}

	for (int k=0;k<=xray.maxK-xray.minK;k++)
	{
		for (int l=0;l<=xray.maxL-xray.minL;l++)
		{
			matTrig1[k][l]=xray.cosy1d[k]*xray.cosz1d[l]-xray.siny1d[k]*xray.sinz1d[l];
			matTrig2[k][l]=xray.cosy1d[k]*xray.sinz1d[l]+xray.siny1d[k]*xray.cosz1d[l];
		}
	}

	for (int i=0;i<npoint;i++)
	{
		//dp=calcDotProduct(miller[i].Pos, frac);
		//getTrigLookUp(xray.lookUp, dp, cosdp1, sindp1);
		h=int(xray.miller[i].Pos[X]-xray.minH);
		k=int(xray.miller[i].Pos[Y]-xray.minK);
		l=int(xray.miller[i].Pos[Z]-xray.minL);
		
		//cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<endl;

/*
		cout <<"cosx1d= "<<xray.cosx1d[h]<<endl;
		cout <<"cosy1d= "<<xray.cosy1d[k]<<endl;
		cout <<"cosz1d= "<<xray.cosz1d[l]<<endl;

		cout <<"sinx1d= "<<xray.sinx1d[h]<<endl;
		cout <<"siny1d= "<<xray.siny1d[k]<<endl;
		cout <<"sinz1d= "<<xray.sinz1d[l]<<endl;
*/
		

		cosdp=xray.cosx1d[h]*matTrig1[k][l];
		cosdp-=xray.sinx1d[h]*matTrig2[k][l];

		sindp=xray.sinx1d[h]*matTrig1[k][l];
		sindp+=xray.cosx1d[h]*matTrig2[k][l];

		symTempReal[i]=cosdp*symReal[i]-sindp*symImag[i];
		symTempImag[i]=cosdp*symImag[i]+sindp*symReal[i];


		//cout <<"cosdp= "<<cosdp<<" sindp= "<<sindp<<" symReal= "
		//<<symReal[i]<<" symImag= "<<symImag[i]<<" symTempReal= "
		//<<symTempReal[i]<<" symTempImag= "<<symTempImag[i]<<endl;
	}
}

void translateXRayDecomposed5(XRayStruct &xray, Vector &symReal, Vector &symImag, vector<PosStruct> &miller, Vector &frac, Vector &symTempReal, Vector &symTempImag)
{
	int Size=xray.cosx1d.size();
	int h, k, l;
	int npoint=symReal.size();
	int ntemp=symTempReal.size();
	Real x, y, z;
	Real cosdp, sindp;
	Matrix matTrig1, matTrig2;
	
	if (ntemp!=npoint)
	{
		SafeAlloc(symTempReal, npoint, "npoint");
		SafeAlloc(symTempImag, npoint, "npoint");
	}
	frac[X]*=2.0*pi;
	frac[Y]*=2.0*pi;
	frac[Z]*=2.0*pi;

	if (Size==0)
	{
		setMinMiller(xray);
		setMaxMiller(xray);
		SafeAlloc(xray.cosx1d, xray.maxH-xray.minH+1, "cosx");
		SafeAlloc(xray.cosy1d, xray.maxK-xray.minK+1, "cosy");
		SafeAlloc(xray.cosz1d, xray.maxL-xray.minL+1, "cosz");
		SafeAlloc(xray.sinx1d, xray.maxH-xray.minH+1, "sinx");
		SafeAlloc(xray.siny1d, xray.maxK-xray.minK+1, "siny");
		SafeAlloc(xray.sinz1d, xray.maxL-xray.minL+1, "sinz");

		Safe2DAlloc(xray.matTrig1, xray.maxK-xray.minK+1, xray.maxL-xray.minL+1, "matTrig1");
		Safe2DAlloc(xray.matTrig2, xray.maxK-xray.minK+1, xray.maxL-xray.minL+1, "matTrig2");
	}
	

	//cout <<"minH= "<<xray.minH<<endl;
	//cout <<"maxH= "<<xray.maxH<<endl;

	for (int i=xray.minH;i<=xray.maxH;i++)
	{
		x=Real(i)*frac[X];
		getTrigLookUp(xray.lookUp, x, xray.cosx1d[i-xray.minH], xray.sinx1d[i-xray.minH]);
		//cout <<"x= "<<x<<" cosx1d= "<<xray.cosx1d[i-xray.minH]<<" sinx1d= "<<xray.sinx1d[i-xray.minH]<<endl;
	}
	for (int i=xray.minK;i<=xray.maxK;i++)
	{
		y=Real(i)*frac[Y];
		getTrigLookUp(xray.lookUp, y, xray.cosy1d[i-xray.minK], xray.siny1d[i-xray.minK]);
		//cout <<"y= "<<y<<" cosy1d= "<<xray.cosy1d[i-xray.minH]<<" siny1d= "<<xray.siny1d[i-xray.minH]<<endl;
	}
	for (int i=xray.minL;i<=xray.maxL;i++)
	{
		z=Real(i)*frac[Z];
		getTrigLookUp(xray.lookUp, z, xray.cosz1d[i-xray.minL], xray.sinz1d[i-xray.minL]);
		//cout <<"z= "<<z<<" cosz1d= "<<xray.cosz1d[i-xray.minH]<<" sinz1d= "<<xray.sinz1d[i-xray.minH]<<endl;
	}
	for (int k=0;k<=xray.maxK-xray.minK;k++)
	{
		for (int l=0;l<=xray.maxL-xray.minL;l++)
		{
			xray.matTrig1[k][l]=xray.cosy1d[k]*xray.cosz1d[l]-xray.siny1d[k]*xray.sinz1d[l];
			xray.matTrig2[k][l]=xray.cosy1d[k]*xray.sinz1d[l]+xray.siny1d[k]*xray.cosz1d[l];
		}
	}
	for (int i=0;i<npoint;i++)
	{
		//dp=calcDotProduct(miller[i].Pos, frac);
		//getTrigLookUp(xray.lookUp, dp, cosdp1, sindp1);
		h=int(xray.miller[i].Pos[X]-xray.minH);
		k=int(xray.miller[i].Pos[Y]-xray.minK);
		l=int(xray.miller[i].Pos[Z]-xray.minL);
		
		//cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<endl;

		cosdp=xray.cosx1d[h]*xray.matTrig1[k][l];
		cosdp-=xray.sinx1d[h]*xray.matTrig2[k][l];

		sindp=xray.sinx1d[h]*xray.matTrig1[k][l];
		sindp+=xray.cosx1d[h]*xray.matTrig2[k][l];

		symTempReal[i]=cosdp*symReal[i]-sindp*symImag[i];
		symTempImag[i]=cosdp*symImag[i]+sindp*symReal[i];


		//cout <<"cosdp= "<<cosdp<<" sindp= "<<sindp<<" symReal= "
		//<<symReal[i]<<" symImag= "<<symImag[i]<<" symTempReal= "
		//<<symTempReal[i]<<" symTempImag= "<<symTempImag[i]<<endl;
	}
}


void calcSymmetryScattering(XRayStruct &xray, ProteinStruct &Protein, int prot)
{
	int noperations=Protein.symmetryOperations.size();
	int nsym=xray.complexSymReal.size();
	ProteinStruct tempProtein;

	/*
	   int npoint=xray.miller.size();
	   Real qMag, dwf;
	   if (npoint>1)
	   {
	   for (int i=3;i<=3;i++)
	   {
	   for (int j=0;j<npoint;j++)
	   {
	   qMag=calcMagnitude(xray.q[j]);
	   dwf=calcDWF(xray.q[j], Protein.Atoms[0].BFactor, xray.lookUp);
	   cout <<xray.qMag[j]<<"\t"<<xray.f[i][j]<<"\t"<<xray.f[i][j]*dwf<<"\t"<<dwf<<"\t"<<xray.miller[j].Pos[X]<<"\t"<<xray.miller[j].Pos[Y]<<"\t"<<xray.miller[j].Pos[Z]<<"\t"<<qMag<<endl;
	   }
	   }
	   }
	 */
	if (prot>=nsym)
	{
		string errorStr="prot= "+toStr(prot)+" nsym= "+toStr(nsym);
		error(errorStr, __LINE__, __FILE__);
	}
	tempProtein=Protein;
	setSpaceGroupP1(tempProtein);
	xray.complexSymReal[prot].clear();
	xray.complexSymImag[prot].clear();
	//int operation=1;
	for (int i=0;i<noperations;i++)
		//for (int i=operation;i<=operation;i++)
	{
		//SafeAlloc(Protein.symmetryOperations[i].translationVector, 3, "t");
		applySymmetryOperation(tempProtein.Atoms, Protein.symmetryOperations[i]);
		calcScatteringFast(tempProtein, xray);
		SafePushBack(xray.complexSymReal[prot], xray.real, "xray.symReal");
		SafePushBack(xray.complexSymImag[prot], xray.imag, "xray.symImag");
		tempProtein.Atoms=Protein.Atoms;
	}
/*
	int npoint=xray.complexSymReal[0][0].size();
	//if (xray.boolVerbose)
	//{
		cout <<"In boolVerbose= true"<<endl;
		for (int i=0;i<npoint;i++)
		{
			for (int j=0;j<noperations;j++)
			{
				cout <<"prot= "<<prot<<" "<<"i= "<<i<<" j= "<<j<<" "<<xray.complexSymReal[prot][j][i]<<"\t"<<xray.complexSymImag[prot][j][i]<<"\t"<<endl;
			}
			cout <<endl;
		}
	//}
*/
}

void calcSymmetryScattering(XRayStruct &xray, vector<ProteinStruct> &Proteins)
{
	int nprot=Proteins.size();
	int nsym;
	int npoint=xray.real.size();

	if (nprot==0)
	{
		error("No proteins.", __LINE__, __FILE__);
	}

	nsym=Proteins[0].symmetryOperations.size();

	xray.complexSymReal.clear();
	xray.complexSymImag.clear();

	Safe3DAlloc(xray.complexSymReal, nprot, nsym, npoint, "complexSymReal");
	Safe3DAlloc(xray.complexSymImag, nprot, nsym, npoint, "complexSymImag");
	Safe3DAlloc(xray.complexSymTempReal, nprot, nsym, npoint, "complexSymTempReal");
	Safe3DAlloc(xray.complexSymTempImag, nprot, nsym, npoint, "complexSymTempImag");

	//cout <<"After Safe3DAlloc"<<endl;

	for (int i=0;i<nprot;i++)
	{
		//cout <<"i= "<<i<<endl;
		calcSymmetryScattering(xray, Proteins[i], i);
	}
	//cout <<"After loop"<<endl;
	xray.complexSymTempReal=xray.complexSymReal;
	xray.complexSymTempImag=xray.complexSymImag;

	//cout <<"Before calcTotalFormFactor"<<endl;
	calcTotalFormFactor(xray);
	//cout <<"Before calcIntensity"<<endl;
	calcIntensity(xray);
	//cout <<"Before calcPhase"<<endl;
	calcPhase(xray);
	//cout <<"Leaving calcSymmetryScattering"<<endl;
}

void calcSymmetryScattering(XRayStruct &xray, ProteinStruct &Protein)
{
	vector<ProteinStruct> Proteins;

	splitIntoChains(Protein, Proteins);
	calcSymmetryScattering(xray, Proteins);
}

int identifyRotationBlocks(ProteinStruct &Protein)
{
	bool isEqual;
	int nunique;
	int nsym=Protein.symmetryOperations.size();

	nunique=nsym;
	for (int i=1;i<nsym;i++)
	{
		isEqual=matrixEqual(Protein.symmetryOperations[0].rotationMat, Protein.symmetryOperations[i].rotationMat);
		if (isEqual)
		{
			nunique=i;
			return nunique;
		}
	}
	return nunique;
}

void translateXRayBlock(XRayStruct &xray, ProteinStruct &Protein, Vector &dx, int prot, int nunique)
{
	int noperations=Protein.symmetryOperations.size();
	int nsym=xray.complexSymReal[0].size();
	Vector frac, dx2, ts;

	SafeAlloc(frac, 3, "frac");
	SafeAlloc(ts, 3, "ts");
	if (noperations!=nsym) calcSymmetryScattering(xray, Protein, prot);
	for (int i=0;i<nunique;i++)
	{
		matrixMultiply(Protein.symmetryOperations[i].rotationMat, dx, dx2);
		for (int j=0;j<3;j++)
		{
			dx2[j]+=Protein.symmetryOperations[i].translationVector[j];
		}
		convertToFractionalCoordinates(dx2[X], dx2[Y], dx2[Z], Protein, frac[X], frac[Y], frac[Z]);
		//translateXRay(xray.lookUp, xray.complexSymReal[prot][i], xray.complexSymImag[prot][i], xray.miller, frac, xray.complexSymTempReal[prot][i], xray.complexSymTempImag[prot][i]);
		translateXRayDecomposed5(xray, xray.complexSymReal[prot][i], xray.complexSymImag[prot][i], xray.miller, frac, xray.complexSymTempReal[prot][i], xray.complexSymTempImag[prot][i]);
	}
}

void translateXRayBlocks(ProteinStruct &Protein, XRayStruct &xray, int nunique, Vector &realBlock, Vector &imagBlock, Matrix &realTranslated, Matrix &imagTranslated)
{
	int nblock=realTranslated.size();
	Real t0, t1;
	Vector dx, frac;

	SafeAlloc(dx, 3, "dx");
	SafeAlloc(frac, 3, "frac");

	realTranslated[0]=realBlock;
	imagTranslated[0]=imagBlock;
	for (int i=1;i<nblock;i++)
	{
		for (int j=0;j<3;j++)
		{
			t0=Protein.symmetryOperations[0].translationVector[j];
			t1=Protein.symmetryOperations[i*nunique].translationVector[j];
			dx[j]=t1-t0;
		}
		convertToFractionalCoordinates(dx[X], dx[Y], dx[Z], Protein, frac[X], frac[Y], frac[Z]);
		//translateXRay(xray.lookUp, realBlock, imagBlock, xray.miller, frac, realTranslated[i], imagTranslated[i]);	
		translateXRayDecomposed5(xray, realBlock, imagBlock, xray.miller, frac, realTranslated[i], imagTranslated[i]);
	}
}

void translateXRayVeryFast(XRayStruct &xray, ProteinStruct &Protein, Vector &dx, int prot)
{
	//Accounts for translation very quickly by taking into account that 
	//in some space groups blocks of symmetry operations are related by
	//translation, and therefore not every symmetry mate has to be
	//translated separetely.  E.g I212121 has 8 symmetry operations but
	//only 4 unique rotation matrixes.  The two groups of 4 symmetry mates
	//are related by a translation.
	int nsym=xray.complexSymReal[0].size();
	int npoint=xray.miller.size();
	int nunique, nblock;
	int nrealBlock=xray.realBlock.size();
	Vector frac, dx2, ts;

	SafeAlloc(frac, 3, "frac");
	SafeAlloc(ts, 3, "ts");
	nunique=identifyRotationBlocks(Protein);
	translateXRayBlock(xray, Protein, dx, prot, nunique);
	nblock=nsym/nunique;
	if (nrealBlock!=npoint)
	{
		SafeAlloc(xray.realBlock, npoint, "realBlock");	
		SafeAlloc(xray.imagBlock, npoint, "imagBlock");
		Safe2DAlloc(xray.realBlockTranslated, nblock, npoint, "realTranslated");	
		Safe2DAlloc(xray.imagBlockTranslated, nblock, npoint, "imagTranslated");	
	}	
	zeroVector(xray.realBlock);
	zeroVector(xray.imagBlock);
	calcBlockFormFactor(xray.complexSymTempReal[prot], xray.complexSymTempImag[prot], xray.realBlock, xray.imagBlock, nunique);
	translateXRayBlocks(Protein, xray, nunique, xray.realBlock, xray.imagBlock, xray.realBlockTranslated, xray.imagBlockTranslated);
	calcTotalFormFactor(xray.realBlockTranslated, xray.imagBlockTranslated, xray.real, xray.imag);
}

void translateXRay2(XRayStruct &xray, ProteinStruct &Protein, Vector &dx, int prot)
{
	int noperations=Protein.symmetryOperations.size();
	int nsym=xray.complexSymReal[0].size();
	Vector frac, dx2, ts;
	SafeAlloc(frac, 3, "frac");
	SafeAlloc(ts, 3, "ts");
	
	if (noperations!=nsym) calcSymmetryScattering(xray, Protein, prot);
	for (int i=0;i<noperations;i++)
	{
		matrixMultiply(Protein.symmetryOperations[i].rotationMat, dx, dx2);
		for (int j=0;j<3;j++)
		{
			dx2[j]+=Protein.symmetryOperations[i].translationVector[j];
		}
		convertToFractionalCoordinates(dx2[X], dx2[Y], dx2[Z], Protein, frac[X], frac[Y], frac[Z]);
		//translateXRay(xray.lookUp, xray.complexSymReal[prot][i], xray.complexSymImag[prot][i], xray.miller, frac, xray.complexSymTempReal[prot][i], xray.complexSymTempImag[prot][i]);
		translateXRayDecomposed5(xray, xray.complexSymReal[prot][i], xray.complexSymImag[prot][i], xray.miller, frac, xray.complexSymTempReal[prot][i], xray.complexSymTempImag[prot][i]);
	}
	//calcTotalFormFactor(xray.symTempReal, xray.symTempImag, xray.real, xray.imag);
	//calcTotalFormFactor(xray.complexSymTempReal, xray.complexSymTempImag, xray.real, xray.imag);
	//calcIntensity(xray);
	//cout <<"Leaving translateXRay2"<<endl;
}

void translateXRay(XRayStruct &xray, ProteinStruct &Protein, Vector &dx, int prot)
{
	int noperations=Protein.symmetryOperations.size();
	int nsym=xray.complexSymReal[0].size();
	Vector frac, dx2, ts;
	//cout <<"In translateXRay"<<endl;
	SafeAlloc(frac, 3, "frac");
	SafeAlloc(ts, 3, "ts");
	if (noperations!=nsym) calcSymmetryScattering(xray, Protein, prot);
	for (int i=0;i<noperations;i++)
	{
		matrixMultiply(Protein.symmetryOperations[i].rotationMat, dx, dx2);
		convertToFractionalCoordinates(dx2[X], dx2[Y], dx2[Z], Protein, frac[X], frac[Y], frac[Z]);
		translateXRay(xray.lookUp, xray.complexSymReal[prot][i], xray.complexSymImag[prot][i], xray.miller, frac, xray.complexSymTempReal[prot][i], xray.complexSymTempImag[prot][i]);
	}
	//calcTotalFormFactor(xray.symTempReal, xray.symTempImag, xray.real, xray.imag);
	calcTotalFormFactor(xray.complexSymTempReal, xray.complexSymTempImag, xray.real, xray.imag);
	calcIntensity(xray);
	//cout <<"Leaving translateXRay"<<endl;
}

void trimProtein(ProteinStruct &Protein, XRayParamStruct &params)
{
	if (params.RemoveProtein) removeProteinAtoms(Protein.Atoms);
	if (params.RemoveWaters) removeWaters(Protein.Atoms);
}

void aveOverTrajectory(string dcdFilePaths, ProteinStruct &Protein, XRayStruct &xray, XRayParamStruct &params)
{
	bool isStructure=false;
	string DCDFile, PdbFilePaths;
	int nthstruct=1, NumStructures=0;
	int nInDcd, nthDcd=1;
	int maxStructures=params.MaxTrajectoryStructures;
	int Size=xray.real.size();
	Vector totalReal, totalImag;
	ProteinStruct tempProtein;
	cout <<"In aveOverTrajectory"<<endl;
	SafeAlloc(totalReal, Size);
	SafeAlloc(totalImag, Size);
	readPdb(params.PdbFile, Protein.Atoms);
	tempProtein=Protein;
	while (true)
	{
		isStructure=GetNextFrame(tempProtein, nthstruct, nInDcd, nthDcd, DCDFile, dcdFilePaths, PdbFilePaths);
		Protein=tempProtein;
		trimProtein(Protein, params);
		if (!isStructure) break;
		NumStructures++;
		if (NumStructures>=maxStructures) break;
		cout <<"Before calcSctteringBfactor"<<endl;
		calcScatteringBfactor(Protein, xray);
		cout <<"After calcScatteringBfactor"<<endl;
		totalReal=addVectors(xray.real, totalReal);
		totalImag=addVectors(xray.imag, totalImag);
	}
	xray.real=totalReal;
	xray.imag=totalImag;
	calcIntensity(xray);
	calcAmplitude(xray);
	calcPhase(xray);
}

void bulkSolventXRay(XRayStruct &xray, XRayStruct &xraySol, Real ksol, Real bsol)
{
	int npoint=xray.q.size();
	int nreal=xray.real.size();
	Real mag;

	if (npoint!=nreal)
	{
		cout <<"ERROR: nreal= "<<nreal<<" npoint= "<<npoint<<endl;
		exit(EXIT_FAILURE);
	}

	for (int i=0;i<npoint;i++)
	{
		mag=calcMagnitude(xray.q[i]);
		xraySol.real[i]=-ksol*xray.real[i]/exp(bsol*mag*mag);
		xraySol.imag[i]=-ksol*xray.imag[i]/exp(bsol*mag*mag);
	}	

}

void bulkSolventCorrection(XRayStruct &xray, XRayStruct &xraySol, XRayStruct &xrayCorrected, Real ksol, Real bsol)
{
	int npoint=xray.q.size();
	int nreal=xray.real.size();

	if (npoint!=nreal)
	{
		cout <<"ERROR: nreal= "<<nreal<<" npoint= "<<npoint<<endl;
		exit(EXIT_FAILURE);
	}
	bulkSolventXRay(xray, xraySol, ksol, bsol);
	for (int i=0;i<npoint;i++)
	{
		xrayCorrected.real[i]=xray.real[i]+xraySol.real[i]+xray.waterReal[i];
		xrayCorrected.imag[i]=xray.imag[i]+xraySol.imag[i]+xray.waterImag[i];
	}	
	calcIntensity(xrayCorrected);
	calcPhase(xrayCorrected);
}

Real calcBulkSolventCorrectionScore(XRayStruct &xray, XRayStruct &expXRay, XRayStruct &xraySol, XRayStruct &xrayCorrected, Real ksol, Real bsol, XRayParamStruct &params)
{
	bulkSolventCorrection(xray, xraySol, xrayCorrected, ksol, bsol);
	return calcRFactor(xrayCorrected, expXRay, params);
}

Real calcBulkSolventCorrectionScore(Vector solCoefficients, SolventStruct args)
{
	Real ksol=solCoefficients[0];
	Real bsol=solCoefficients[1];

	return calcBulkSolventCorrectionScore(args.xray, args.expXRay, args.xraySol, args.xrayCorrected, ksol, bsol, args.params);
}

Real bulkSolventCorrection(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	Real dScoreCriteria=0.00000001;
	Real RFactor, maxHours=200.0;
	Vector solCoefficients;
	Vector q, expQ, intensity, expIntensity;
	SolventStruct args;
	Real (*func_to_minimize)(Vector, SolventStruct)=calcBulkSolventCorrectionScore;


	args.xray=xray;
	args.expXRay=expXRay;
	args.xraySol=xray;
	args.xrayCorrected=xray;
	args.params=params;

	solCoefficients.resize(2);
	solCoefficients[0]=0.75;
	solCoefficients[1]=50.00;
	RFactor=conjugateGradient(solCoefficients, args, func_to_minimize, dScoreCriteria, maxHours);
	calcBulkSolventCorrectionScore(xray, expXRay, args.xraySol, args.xrayCorrected, solCoefficients[0], solCoefficients[1], params);
	xray=args.xrayCorrected;


	return RFactor;	
}

Real calcScattering(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	Real rfactor;
	int natom=Proteins[0].Atoms.size();

	for (int i=0;i<natom;i++)
	{
		printAtomInfo(Proteins[0].Atoms[i]);
	}
	if (params.UseFastTranslationRotation)
	{
		SafeAlloc(Proteins[0].translation, 3, "Protein.translation");
		calcFastTranslationRotation(Proteins, xray, params);
		xray.complexSymTempReal=xray.complexSymReal;
		xray.complexSymTempImag=xray.complexSymImag;
		//for (int i=0;i<nprot;i++)
		//{
		//	dx[X]=degOfFreedom[i*6+3+X];
		//	dx[Y]=degOfFreedom[i*6+3+Y];
		//	dx[Z]=degOfFreedom[i*6+3+Z];
		//	translateXRay2(xray, Protein, dx, i);
		//}
		calcTotalFormFactor(xray);
		calcIntensity(xray);
	}
	else 
	{
		//calcScatteringComplex(Proteins, xray);
		calcSymmetryScattering(xray, Proteins[0]);
		calcTotalFormFactor(xray);
		calcIntensity(xray);
	}
	if (params.BulkSolventCorrection)
	{
		bulkSolventCorrection(xray, expXRay, params);
	}
	if (params.CorrectionFile!="") correctIntensity(xray);
	if (params.OccupancyCorrectionFile!="" && params.OptimizeOccupancies) correctOccupancyIntensity(xray, Proteins[0], params);
	fixBFactor(xray, expXRay);
	//string xrayOutput="/home/joukov/project/mr/XRay_out/1D2N_1P_ModRefiner.txt";
	//printXRay(xrayOutput, xray, Proteins[0]);
	//endProgram(__LINE__, __FILE__);
	rfactor=calcRFactor(xray, expXRay, params);
	return rfactor;
}

Real calcScattering(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	vector<ProteinStruct> Proteins;

	SafePushBack(Proteins, Protein, "Proteins");
	return calcScattering(Proteins, xray, expXRay, params);
}

Real calcScattering(Vector degOfFreedom, MRStruct args)
{
	Real rfactor;
	moveAtoms(args.Protein.Atoms, degOfFreedom[3], degOfFreedom[4], degOfFreedom[5]);
	rotateAtoms(args.Protein.Atoms, degOfFreedom[0], degOfFreedom[1], degOfFreedom[2]);
	calcScatteringBfactor(args.Protein, args.xray);
	rfactor=RFactor(args.xray.i, args.expXRay.i);
	cout <<"RFactor= "<<rfactor<<endl;
	return rfactor;
}

void setResidueOccupancyToZero(vector<AtomStruct> &Atoms, int residue)
{
	int natom=Atoms.size();

	for (int i=0;i<natom;i++)
	{
		if (Atoms[i].ResidueNum==residue)
		{
			Atoms[i].Occupancy=0;
		}		
	}
}

void setSurfaceAtomsClashingWithCoreAtomsToZero(vector<AtomStruct> &Atoms)
{
	int natom=Atoms.size();
	int nclash;

	for (int i=0;i<natom;i++)
	{
		nclash=Atoms[i].clash.size();
		if (!Atoms[i].core && nclash>0)
		{
			for (int j=0;j<nclash;j++)
			{
				if (Atoms[Atoms[i].clash[j]].core)
				{
					setResidueOccupancyToZero(Atoms, Atoms[i].ResidueNum);
				}
			}
		}
	}
}

void addAtomToScattering(AtomStruct &Atom, ProteinStruct &Protein, XRayStruct &xray)
{
	int npoint=xray.q.size();
	int nsym=Protein.symmetryOperations.size();
	int atomID;
	Real fracx, fracy, fracz;
	Real dotproduct;
	Real cosdp, sindp;
	Real product, occupancy, dwf;
	AtomStruct tempAtom;

	atomID=Atom.atomid;
	occupancy=Atom.Occupancy;
	for (int i=0;i<nsym;i++)
	{
		applySymmetryOperation(Atom, Protein.symmetryOperations[i], tempAtom);
		convertToFractionalCoordinates(tempAtom.x, tempAtom.y, tempAtom.z, Protein.cartToFrac, fracx, fracy, fracz);
		fracx*=2.0*pi; fracy*=2.0*pi; fracz*=2.0*pi;	
		for (int j=0;j<npoint;j++)
		{
			dotproduct=fracx*xray.miller[j].Pos[X];
			dotproduct+=fracy*xray.miller[j].Pos[Y];
			dotproduct+=fracz*xray.miller[j].Pos[Z];
			dwf=calcDWF(xray, j, tempAtom, xray.lookUp);
			getTrigLookUp(xray.lookUp, dotproduct, cosdp, sindp);
			product=xray.f[atomID][j]*dwf*occupancy;
			xray.real[j]+=product*cosdp;
			xray.imag[j]+=product*sindp;
		}
	}
	calcIntensity(xray);
	calcAmplitude(xray);
	calcPhase(xray);
}

void subtractAtomFromScattering(AtomStruct &Atom, ProteinStruct &Protein, XRayStruct &xray)
{
	int npoint=xray.q.size();
	int nsym=Protein.symmetryOperations.size();
	int atomID;
	Real fracx, fracy, fracz;
	Real dotproduct;
	Real cosdp, sindp;
	Real product, occupancy, dwf;
	AtomStruct tempAtom;

	atomID=Atom.atomid;
	occupancy=Atom.Occupancy;
	for (int i=0;i<nsym;i++)
	{
		applySymmetryOperation(Atom, Protein.symmetryOperations[i], tempAtom);
		convertToFractionalCoordinates(tempAtom.x, tempAtom.y, tempAtom.z, Protein.cartToFrac, fracx, fracy, fracz);
		fracx*=2.0*pi; fracy*=2.0*pi; fracz*=2.0*pi;	
		for (int j=0;j<npoint;j++)
		{
			dotproduct=fracx*xray.miller[j].Pos[X];
			dotproduct+=fracy*xray.miller[j].Pos[Y];
			dotproduct+=fracz*xray.miller[j].Pos[Z];
			dwf=calcDWF(xray, j, tempAtom, xray.lookUp);
			getTrigLookUp(xray.lookUp, dotproduct, cosdp, sindp);
			product=xray.f[atomID][j]*dwf*occupancy;
			xray.real[j]-=product*cosdp;
			xray.imag[j]-=product*sindp;
		}
	}
	calcIntensity(xray);
	calcAmplitude(xray);
	calcPhase(xray);
}

void optimizeOccupanciesOfClashingAtoms(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, int atom1, int atom2, XRayParamStruct &params)
{
	//Assumes that all proteins are identical
	int best;
	int nprot=Proteins.size();
	Real score, bestScore;
	Real occupancy1, occupancy2;

	if (nprot==0)
	{
		error("No protein ", __LINE__, __FILE__);
	}

	occupancy1=Proteins[0].Atoms[atom1].Occupancy;
	occupancy2=Proteins[0].Atoms[atom2].Occupancy;

	//cout <<"atom1= "<<atom1<<" atom2= "<<atom2<<" occupancy1= "<<occupancy1<<" occupancy2= "<<occupancy2<<endl;

	for (int i=0;i<nprot;i++)
	{
		subtractAtomFromScattering(Proteins[i].Atoms[atom1], Proteins[i], xray);
		Proteins[i].Atoms[atom1].Occupancy=0;
	}
	if (params.CorrectionFile!="") correctIntensity(xray);
	bestScore=calcRFactor(xray, expXRay, params);
	best=0;

	for (int i=0;i<nprot;i++)
	{
		Proteins[i].Atoms[atom1].Occupancy=occupancy1;
		addAtomToScattering(Proteins[i].Atoms[atom1], Proteins[i], xray);
		subtractAtomFromScattering(Proteins[i].Atoms[atom2], Proteins[i], xray);
		Proteins[i].Atoms[atom2].Occupancy=0;
	}
	if (params.CorrectionFile!="") correctIntensity(xray);
	score=calcRFactor(xray, expXRay, params);
	if (score<bestScore)
	{
		bestScore=score;
		best=1;
	}
	for (int i=0;i<nprot;i++)
	{
		subtractAtomFromScattering(Proteins[i].Atoms[atom1], Proteins[i], xray);
		Proteins[i].Atoms[atom1].Occupancy=0;
		Proteins[i].Atoms[atom2].Occupancy=0;
	}
	if (params.CorrectionFile!="") correctIntensity(xray);
	score=calcRFactor(xray, expXRay, params);
	if (score<bestScore)
	{
		bestScore=score;
		best=2;
	}
	for (int i=0;i<nprot;i++)
	{
		if (best==0)
		{
			Proteins[i].Atoms[atom1].Occupancy=0;
			Proteins[i].Atoms[atom2].Occupancy=occupancy2;
			addAtomToScattering(Proteins[i].Atoms[atom2], Proteins[i], xray);
		}
		if (best==1)
		{
			Proteins[i].Atoms[atom1].Occupancy=occupancy1;
			addAtomToScattering(Proteins[i].Atoms[atom1], Proteins[i], xray);
			Proteins[i].Atoms[atom2].Occupancy=0;
		}
		if (best==2)
		{
			Proteins[i].Atoms[atom1].Occupancy=0;
			Proteins[i].Atoms[atom2].Occupancy=0;
		}
	}
}

void optimizeOccupanciesOfClashingAtoms(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	int natom=Proteins[0].Atoms.size();
	int nclash;

	calcScattering(Proteins, xray, expXRay, params);
	for (int i=0;i<natom;i++)
	{
		nclash=Proteins[0].Atoms[i].clash.size();
		//if (!Proteins[0].Atoms[i].core && nclash>0 && Proteins[0].Atoms[i].Occupancy!=0)
		//{
			for (int j=0;j<nclash;j++)
			{
				//if (!Proteins[0].Atoms[Proteins[0].Atoms[i].clash[j]].core)
				//{
				if (i<Proteins[0].Atoms[i].clash[j])
				{
					optimizeOccupanciesOfClashingAtoms(Proteins, xray, expXRay, i, Proteins[0].Atoms[i].clash[j], params);
					if (Proteins[0].Atoms[i].Occupancy==0) break;
				}
				//}
			}
		//}
	}
}

void flexibleScattering(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	int nprot=Proteins.size();
	for (int i=0;i<nprot;i++)
	{
		setSurfaceAtomsClashingWithCoreAtomsToZero(Proteins[i].Atoms);
	}
	optimizeOccupanciesOfClashingAtoms(Proteins, xray, expXRay, params);
}

void flexibleScattering(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	vector<ProteinStruct> Proteins;

	SafePushBack(Proteins, Protein, "Proteins");
	flexibleScattering(Proteins, xray, expXRay, params);
	Protein=Proteins[0];
}

void calcTrigForScattering(LatticeStruct &lattice, XRayStruct &xray)
{
	int Size1, Size2;
	int maxXBin, maxYBin, maxZBin;
	Real x, y, z;
	Real fracx, fracy, fracz;
	Real fracxj, fracyj, fraczj;
	timeval start, end;


	setMinMiller(xray);
	gettimeofday(&start, NULL);
	Get2DVectorSize(xray.cosx, Size1, Size2);
	Get3DVectorSize(lattice.density, maxXBin, maxYBin, maxZBin, "density");
	if (Size1!=maxXBin || Size2!=xray.maxH-xray.minH+1)
	{
		Safe2DAlloc(xray.cosx, maxXBin, xray.maxH-xray.minH+1, "cosx");
		Safe2DAlloc(xray.cosy, maxYBin, xray.maxK-xray.minK+1, "cosy");
		Safe2DAlloc(xray.cosz, maxZBin, xray.maxL-xray.minL+1, "cosz");
		Safe2DAlloc(xray.sinx, maxXBin, xray.maxH-xray.minH+1, "sinx");
		Safe2DAlloc(xray.siny, maxYBin, xray.maxK-xray.minK+1, "siny");
		Safe2DAlloc(xray.sinz, maxZBin, xray.maxL-xray.minL+1, "sinz");

	}
	gettimeofday(&end, NULL);
	gettimeofday(&start, NULL);
	for (int i=0;i<maxXBin;i++)
	{
		x=lattice.firstPos.Pos[X]+Real(i)*lattice.xCubeLength;
		fracx=x*2.0*pi;
		for (int j=xray.minH;j<=xray.maxH;j++)
		{
			fracxj=fracx*Real(j);
			getTrigLookUp(xray.lookUp, fracxj, xray.cosx[i][j-xray.minH], xray.sinx[i][j-xray.minH]);
		}	
	}
	for (int i=0;i<maxYBin;i++)
	{
		y=lattice.firstPos.Pos[Y]+Real(i)*lattice.yCubeLength;
		fracy=y*2.0*pi;
		for (int j=xray.minK;j<=xray.maxK;j++)
		{
			fracyj=fracy*Real(j);
			getTrigLookUp(xray.lookUp, fracyj, xray.cosy[i][j-xray.minK], xray.siny[i][j-xray.minK]);
		}	
	}
	for (int i=0;i<maxZBin;i++)
	{
		z=lattice.firstPos.Pos[Z]+Real(i)*lattice.zCubeLength;
		fracz=z*2.0*pi;
		for (int j=xray.minL;j<=xray.maxL;j++)
		{
			fraczj=fracz*Real(j);
			getTrigLookUp(xray.lookUp, fraczj, xray.cosz[i][j-xray.minL], xray.sinz[i][j-xray.minL]);
		}	
	}
}

void calcRealImag(LatticeStruct &lattice, XRayStruct &xray)
{
	int h, k, l;
	int npoint=xray.miller.size();
	int maxXBin, maxYBin, maxZBin;
	Real real, imag;

	Get3DVectorSize(lattice.density, maxXBin, maxYBin, maxZBin, "density");
	for (int a=0;a<maxXBin;a++)
	{
		for (int b=0;b<maxYBin;b++)
		{
			for (int c=0;c<maxZBin;c++)
			{
				for (int j=0;j<npoint;j++)
				{
					real=imag=0;
					h=int(xray.miller[j].Pos[X]-xray.minH);
					k=int(xray.miller[j].Pos[Y]-xray.minK);
					l=int(xray.miller[j].Pos[Z]-xray.minL);

					real+=xray.cosx[a][h]*xray.cosy[b][k]*xray.cosz[c][l];
					real-=xray.cosx[a][h]*xray.siny[b][k]*xray.sinz[c][l];
					real-=xray.sinx[a][h]*xray.cosy[b][k]*xray.sinz[c][l];
					real-=xray.sinx[a][h]*xray.siny[b][k]*xray.cosz[c][l];

					imag-=xray.sinx[a][h]*xray.siny[b][k]*xray.sinz[c][l];
					imag+=xray.sinx[a][h]*xray.cosy[b][k]*xray.cosz[c][l];
					imag+=xray.cosx[a][h]*xray.siny[b][k]*xray.cosz[c][l];
					imag+=xray.cosx[a][h]*xray.cosy[b][k]*xray.sinz[c][l];
					real*=lattice.density[a][b][c];
					imag*=lattice.density[a][b][c];
					/*
					if (lattice.density[a][b][c]!=0)
					{
						cout <<"a= "<<a<<" b= "<<b<<" c= "<<c<<endl;
						cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<endl;
						cout <<"density= "<<lattice.density[a][b][c]<<endl;
						cout <<"real= "<<real<<" imag= "<<imag<<endl;
						cout <<"cosx= "<<xray.cosx[a][h]<<endl;
						cout <<"cosy= "<<xray.cosy[b][k]<<endl;
						cout <<"cosz= "<<xray.cosz[c][l]<<endl;
						cout <<"sinx= "<<xray.sinx[a][h]<<endl;
						cout <<"siny= "<<xray.siny[b][k]<<endl;
						cout <<"sinz= "<<xray.sinz[c][l]<<endl;
						endProgram(__LINE__, __FILE__);
					}
					*/
					xray.real[j]+=real;
					xray.imag[j]+=imag;

				}
			}
		}
	}
}

void multiplyByParallelPipped(LatticeStruct &lattice, XRayStruct &xray)
{
	int npoint=xray.miller.size();
	Real trig;
	Real h, k, l;
	Real x, y, z;

	x=lattice.xCubeLength*0.5;
	y=lattice.yCubeLength*0.5;
	z=lattice.zCubeLength*0.5;
	for (int i=0;i<npoint;i++)
	{
		h=xray.miller[i].Pos[X];
		k=xray.miller[i].Pos[Y];
		l=xray.miller[i].Pos[Z];
		trig=1.0;
		if (h!=0) trig*=sin(2.0*pi*h*x)/(pi*h);
		else trig*=2.0*x;
		if (k!=0) trig*=sin(2.0*pi*k*y)/(pi*k);
		else trig*=2.0*y;
		if (l!=0) trig*=sin(2.0*pi*l*z)/(pi*l);
		else trig*=2.0*z;
		//cout <<"i= "<<i<<" real= "<<xray.real[i]<<" imag= "<<xray.imag[i]<<endl;
		xray.real[i]*=trig;
		xray.imag[i]*=trig;		
		//cout <<"i= "<<i<<" real= "<<xray.real[i]<<" imag= "<<xray.imag[i]<<endl;
	}
}

void calcScatteringFast(LatticeStruct &lattice, XRayStruct &xray)
{
	calcTrigForScattering(lattice, xray);
	calcRealImag(lattice, xray);	
}

void calcScattering(LatticeStruct &lattice, XRayStruct &xray)
{
	bool boolPrint=false;
	int npoint=xray.q.size();
	int MaxXBin, MaxYBin, MaxZBin;
	Real dotproduct, freal, fimag;
	Real x, y, z;

	fimag=0.0;
	cout <<"In calcScattering(LatticeStruct &lattice, XRayStruct &xray)"<<endl;	
	Get3DVectorSize(lattice.density, MaxXBin, MaxYBin, MaxZBin, "lattice.density");
	for (int i=0;i<MaxXBin;i++)
	{
		x=lattice.firstPos.Pos[X]+Real(i)*lattice.xCubeLength;
		for (int j=0;j<MaxYBin;j++)
		{
			y=lattice.firstPos.Pos[Y]+Real(j)*lattice.yCubeLength;
			for (int k=0;k<MaxZBin;k++)
			{
				z=lattice.firstPos.Pos[Z]+Real(k)*lattice.zCubeLength;
				freal=lattice.density[i][j][k];
				//freal=lattice.real[i][j][k];
				//fimag=lattice.imag[i][j][k];
				for (int point=0;point<npoint;point++)
				{
					if (npoint%1000==0) boolPrint=true;
					else boolPrint=false;
					dotproduct=x*xray.q[point].Pos[X];
					dotproduct+=y*xray.q[point].Pos[Y];
					dotproduct+=z*xray.q[point].Pos[Z];
					xray.real[point]+=freal*cos(dotproduct)-fimag*sin(dotproduct);
					xray.imag[point]+=fimag*cos(dotproduct)+freal*sin(dotproduct);
					if (boolPrint)
					{
						cout <<"point= "<<point<<endl;
						cout <<"x= "<<x<<" y= "<<y<<" z= "<<z<<endl;
						cout <<"qx= "<<xray.q[point].Pos[X]
							<<" qy= "<<xray.q[point].Pos[Y]
							<<" qz= "<<xray.q[point].Pos[Z]<<endl;
						cout <<"cos= "<<cos(dotproduct)<<" sin= "<<sin(dotproduct)<<endl;
						cout <<"freal= "<<freal<<" fimag= "<<fimag<<endl;
						cout <<"dotproduct= "<<dotproduct<<endl;
						cout <<"xray.real= "<<xray.real[point]<<" xray.imag= "<<xray.imag[point]<<endl;
					}
				}
			}
		}
	}
	calcIntensity(xray);
	calcPhase(xray);	
}

void calcScattering(vector<CubeStruct> &Cubes, XRayStruct &xray)
{
	int npoint=xray.q.size();
	int nreal=xray.real.size();
	int ncube=Cubes.size();
	Real dotproduct, freal, cubeSize=0.5;
	timeval start, end;
	if (npoint!=nreal)
	{
		cout <<"ERROR: npoint!=nreal."<<endl;
		cout <<"npoint= "<<npoint<<" nreal= "<<nreal<<endl;
		exit(EXIT_FAILURE);
	}
	gettimeofday(&start, NULL);
	zeroVector(xray.i);
	zeroVector(xray.real);
	zeroVector(xray.imag);
	for (int i=0;i<npoint;i++)
	{
		for (int j=0;j<ncube;j++)
		{
			dotproduct=Cubes[j].x*xray.q[i].Pos[X];
			dotproduct+=Cubes[j].y*xray.q[i].Pos[Y];
			dotproduct+=Cubes[j].z*xray.q[i].Pos[Z];
			if (isnan(dotproduct))
			{
				cout <<"ERROR: Dotproduct is nan"<<endl;
				exit(EXIT_FAILURE);
			}
			//xray.real[i]+=freal*cos(dotproduct)-fimag*sin(dotproduct);
			//xray.imag[i]+=fimag*cos(dotproduct)+freal*sin(dotproduct);
			freal=Cubes[j].density;
			freal*=sin(xray.q[i].Pos[X]*cubeSize*0.5);
			freal*=sin(xray.q[i].Pos[Y]*cubeSize*0.5);
			freal*=sin(xray.q[i].Pos[Z]*cubeSize*0.5);
			freal*=8.0/(xray.q[i].Pos[X]*xray.q[i].Pos[Y]*xray.q[i].Pos[Z]);
			xray.real[i]+=freal*cos(dotproduct);
			xray.imag[i]+=freal*sin(dotproduct);
		}	
	}
	gettimeofday(&end, NULL);
	calcIntensity(xray);
	calcPhase(xray);
}

void getMaxMillerIndexes(vector<PosStruct> &miller, int &maxH, int &maxK, int &maxL)
{
	int npoint=miller.size();

	if (npoint==0) error("No miller indexes", __LINE__, __FILE__);

	maxH=int(miller[0].Pos[X]);
	maxK=int(miller[0].Pos[Y]);
	maxL=int(miller[0].Pos[Z]);

	for (int i=0;i<npoint;i++)
	{
		if (int(miller[i].Pos[X])>maxH) maxH=int(miller[i].Pos[X]);
		if (int(miller[i].Pos[Y])>maxK) maxK=int(miller[i].Pos[Y]);
		if (int(miller[i].Pos[Z])>maxL) maxL=int(miller[i].Pos[Z]);
	}
}

void xray1Dto3D(XRayStruct &xray)
{
	int npoint=xray.i.size();
	int nmiller=xray.miller.size();
	int nreal=xray.real.size();
	int nimag=xray.imag.size();
	int nphase=xray.phase.size();
	int nq=xray.q.size();
	PosStruct tempPos;
	vector<PosStruct> vPos;
	vector< vector<PosStruct> > matrixPos;

	if (npoint==0 || npoint!=nmiller || npoint!=nreal || npoint!=nimag || npoint!=nphase || npoint!=nq)
	{
		string errorStr;
		errorStr="npoint= "+IntToStr(npoint)+" ";
		errorStr+="nmiller= "+IntToStr(nmiller)+" ";
		errorStr+="nreal= "+IntToStr(nreal)+" ";
		errorStr+="nimag= "+IntToStr(nimag)+" ";
		errorStr+="nphase= "+IntToStr(nphase)+" ";
		errorStr+="nq= "+IntToStr(nq);
		error(errorStr, __LINE__, __FILE__);
	}

	SafeAlloc(tempPos.Pos, 3, "tempPos");
	getMaxMillerIndexes(xray.miller, xray.maxH, xray.maxK, xray.maxL);
	Safe3DArrayAlloc(xray.i3d, xray.maxH+1, xray.maxK+1, xray.maxL+1, "xray.i3d");
	Safe3DArrayAlloc(xray.real3d, xray.maxH+1, xray.maxK+1, xray.maxL+1, "xray.real3d");
	Safe3DArrayAlloc(xray.imag3d, xray.maxH+1, xray.maxK+1, xray.maxL+1, "xray.imag3d");
	Safe3DArrayAlloc(xray.phase3d, xray.maxH+1, xray.maxK+1, xray.maxL+1, "xray.phase3d");
	SafeAlloc(vPos, tempPos, xray.maxL+1, "vPos");
	SafeAlloc(matrixPos, vPos, xray.maxK+1, "matrixPos");
	SafeAlloc(xray.q3d, matrixPos, xray.maxH+1, "xray.q3d");
	SafeAlloc(xray.miller3d, matrixPos, xray.maxH+1, "xray.miller3d");

	for (int i=0;i<npoint;i++)
	{
		int h,k,l;
		h=int(xray.miller[i].Pos[X]);
		k=int(xray.miller[i].Pos[Y]);
		l=int(xray.miller[i].Pos[Z]);
		xray.i3d[h][k][l]=xray.i[i];
		xray.real3d[h][k][l]=xray.real[i];
		xray.imag3d[h][k][l]=xray.imag[i];
		xray.phase3d[h][k][l]=xray.phase[i];
		xray.q3d[h][k][l]=xray.q[i];
		xray.miller3d[h][k][l]=xray.miller[i];
	}
}

void delete3DXRay(XRayStruct &xray)
{
	Delete3DArray(xray.i3d, xray.maxH+1, xray.maxK+1, xray.maxL+1);
	Delete3DArray(xray.real3d, xray.maxH+1, xray.maxK+1, xray.maxL+1);
	Delete3DArray(xray.imag3d, xray.maxH+1, xray.maxK+1, xray.maxL+1);
	Delete3DArray(xray.phase3d, xray.maxH+1, xray.maxK+1, xray.maxL+1);
}

void calcScatteringFFT(LatticeStruct &lattice, XRayStruct &xray)
{
	int max=0;
	Real ***RealDensity, ***ImagDensity;

	if (xray.maxH>max) max=xray.maxH;
	if (xray.maxK>max) max=xray.maxK;
	if (xray.maxL>max) max=xray.maxL;

	latticeToArray(lattice, RealDensity, ImagDensity);
	xray1Dto3D(xray);

	FFT3D(RealDensity, ImagDensity, xray.real3d, xray.imag3d, max);
}

void copyXRaySize(XRayStruct &xray1, XRayStruct &xray2)
{
	int npoint=xray1.q.size();
	xray2.q=xray1.q;

	SafeAlloc(xray2.i, npoint, "xray2.i");
	SafeAlloc(xray2.real, npoint, "xray2.real");
	SafeAlloc(xray2.imag, npoint, "xray2.imag");
	SafeAlloc(xray2.phase, npoint, "xray2.phase");
	SafeAlloc(xray2.F, npoint, "xray2.F");
}

void addXRays(XRayStruct &xray1, XRayStruct &xray2)
{
	int npoint=xray1.i.size();

	for (int i=0;i<npoint;i++)
	{
		xray2.real[i]+=xray1.real[i];
		xray2.imag[i]+=xray1.imag[i];
	}
	calcPhase(xray2);
	calcIntensity(xray2);
}

void calcScatteringCubeProtein(ProteinStruct &Protein, vector<CubeStruct> &Cubes, XRayStruct &xray)
{
	XRayStruct xRayCubes;

	copyXRaySize(xray, xRayCubes);
	calcScatteringBfactor(Protein, xray);
	calcScattering(Cubes, xRayCubes);
	addXRays(xRayCubes, xray);	
}

Real calcScatteringCubeProtein(ProteinStruct &Protein, vector<CubeStruct> &Cubes, XRayStruct &xray, XRayStruct &expXRay)
{
	calcScatteringCubeProtein(Protein, Cubes, xray, expXRay);
	return calcBestSqrDiff(xray.i, expXRay.i);
}

void combinePhaseAndIntensity(XRayStruct &xray, XRayStruct &expXRay, XRayStruct &combinedXRay)
{
	int npoint=xray.i.size();
	int ni=expXRay.i.size();
	Real sqrMagnitude;

	if (npoint!=ni)
	{
		string errorStr="Differnt number of calculated and experimental points. npoint= "+toStr(npoint)+" ni= "+toStr(ni);
		error(errorStr, __LINE__, __FILE__);
	}
	copyXRaySize(xray, combinedXRay);
	combinedXRay.miller=xray.miller;
	for (int i=0;i<npoint;i++)
	{
		sqrMagnitude=xray.imag[i]*xray.imag[i]+xray.real[i]*xray.real[i];
		combinedXRay.i[i]=expXRay.i[i];
		combinedXRay.real[i]=xray.real[i]*sqrt(expXRay.i[i]/sqrMagnitude);
		combinedXRay.imag[i]=xray.imag[i]*sqrt(expXRay.i[i]/sqrMagnitude);
	}
}

void combinePhaseAndIntensityForExperiment(XRayStruct &xray, XRayStruct &expXRay, XRayStruct &combinedXRay)
{
	int npoint=xray.i.size();
	int ni=expXRay.i.size();

	if (npoint!=ni)
	{
		string errorStr="Differnt number of calculated and experimental points. npoint= "+toStr(npoint)+" ni= "+toStr(ni);
		error(errorStr, __LINE__, __FILE__);
	}
	copyXRaySize(xray, combinedXRay);
	combinedXRay.miller=xray.miller;
	calcIntensity(xray);
	calcAmplitude(xray);
	calcAmplitude(expXRay);
	for (int i=0;i<npoint;i++)
	{
		combinedXRay.F[i]=2.0*expXRay.F[i]-xray.F[i];
		combinedXRay.i[i]=(2.0*expXRay.F[i]-xray.F[i])*(2.0*expXRay.F[i]-xray.F[i]);
		combinedXRay.F[i]=expXRay.F[i];
		combinedXRay.i[i]=expXRay.F[i]*expXRay.F[i];
		combinedXRay.real[i]=xray.real[i]*combinedXRay.F[i]/sqrt(xray.i[i]);
		combinedXRay.imag[i]=xray.imag[i]*combinedXRay.F[i]/sqrt(xray.i[i]);
	}
}

void normalizeDensity(Array3D &density)
{
	int MaxXBin, MaxYBin, MaxZBin;
	Real sum=0, scale=0.2;
	Real npoint;

	Get3DVectorSize(density, MaxXBin, MaxYBin, MaxZBin, "density");
	npoint=Real(MaxXBin*MaxYBin*MaxZBin);

	for (int i=0;i<MaxXBin;i++)
	{
		for (int j=0;j<MaxYBin;j++)
		{
			for (int k=0;k<MaxZBin;k++)
			{
				sum+=density[i][j][k];
			}
		}
	}
	
	for (int i=0;i<MaxXBin;i++)
	{
		for (int j=0;j<MaxYBin;j++)
		{
			for (int k=0;k<MaxZBin;k++)
			{
				density[i][j][k]*=(npoint/sum)*scale;
			}
		}
	}
}

Real calcDensity(XRayStruct &xray, Real x, Real y, Real z)
{
	int npoint=xray.real.size();
	Real real=0, imag=0, dotproduct;

	for (int point=0;point<npoint;point++)
	{
		dotproduct=xray.miller[point].Pos[X]*x;
		dotproduct+=xray.miller[point].Pos[Y]*y;
		dotproduct+=xray.miller[point].Pos[Z]*z;
		dotproduct*=2.0*pi;
		real+=xray.real[point]*cos(dotproduct)+xray.imag[point]*sin(dotproduct);
		imag+=xray.imag[point]*cos(dotproduct)-xray.real[point]*sin(dotproduct);
	}
	return sqrt(real*real+imag*imag);
}

void calcDensity(XRayStruct &xray, LatticeStruct &lattice)
{
	int MaxXBin, MaxYBin, MaxZBin;
	int npoint=xray.real.size();
	Real real, imag, dotproduct;
	Real x, y, z;

	Get3DVectorSize(lattice.density, MaxXBin, MaxYBin, MaxZBin, "density");
	for (int i=0;i<MaxXBin;i++)
	{
		x=lattice.firstPos.Pos[X]+Real(i)*lattice.xCubeLength;
		x=(Real(i)+0.5)/Real(MaxXBin);
		for (int j=0;j<MaxYBin;j++)
		{
			y=lattice.firstPos.Pos[Y]+Real(j)*lattice.yCubeLength;
			y=(Real(j)+0.5)/Real(MaxYBin);
			for (int k=0;k<MaxZBin;k++)
			{
				real=0;
				imag=0;
				z=lattice.firstPos.Pos[Z]+Real(k)*lattice.zCubeLength;
				z=(Real(k)+0.5)/Real(MaxZBin);
				for (int point=0;point<npoint;point++)
				{
					//dotproduct=x*xray.q[point].Pos[X];
					//dotproduct+=y*xray.q[point].Pos[Y];
					//dotproduct+=z*xray.q[point].Pos[Z];
					dotproduct=xray.miller[point].Pos[X]*x;
					dotproduct+=xray.miller[point].Pos[Y]*y;
					dotproduct+=xray.miller[point].Pos[Z]*z;
					dotproduct*=2.0*pi;
					real+=xray.real[point]*cos(dotproduct)+xray.imag[point]*sin(dotproduct);
					imag+=xray.imag[point]*cos(dotproduct)-xray.real[point]*sin(dotproduct);
					//cout <<"dotproduct= "<<dotproduct<<" cos= "<<cos(dotproduct)<<" sin= "<<sin(dotproduct)<<endl;
				}
				//cout <<"i= "<<i<<" j= "<<j<<" k= "<<k<<" real= "<<real<<" imag= "<<imag<<endl;
				lattice.density[i][j][k]=sqrt(real*real+imag*imag)/Real(npoint);
				//lattice.density[i][j][k]=sqrt(real*real+imag*imag);
				lattice.real[i][j][k]=real;
				lattice.imag[i][j][k]=imag;
			}
		}
	}	
	//normalizeDensity(lattice.density);
}

void calcDensity(XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice)
{
	XRayStruct combinedXRay;

	combinePhaseAndIntensity(xray, expXRay, combinedXRay);
	calcDensity(combinedXRay, lattice);
}

void calcExperimentalDensity(XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice)
{
	XRayStruct combinedXRay;

	combinePhaseAndIntensityForExperiment(xray, expXRay, combinedXRay);
	calcDensity(combinedXRay, lattice);
}

void calcDensity(ProteinStruct &Protein, XRayStruct &xray, LatticeStruct &lattice)
{
	cout <<"In calcDensity(Protein, xray, lattice)"<<endl;;
	calcScatteringFast(Protein, xray);
	cout <<"After calcScatteringFast"<<endl;
	calcDensity(xray, lattice);
}

Real calcProteinDensityMatch(ProteinStruct &Protein, XRayStruct &xray, LatticeStruct &calcLattice, LatticeStruct &expLattice)
{
	Real score;
	calcDensity(Protein, xray, calcLattice);
	//score=calcRMSD(calcLattice, expLattice);
	score=1.0-calcDensityCorrelation(calcLattice, expLattice);
	return score;
}

Real calcNativeDensityCorrelation(XRayStruct xray, LatticeStruct &lattice, ProteinStruct &Protein, XRayParamStruct &params)
{
	int Size;
	Real score, minScore;
	Matrix origins;
	ProteinStruct nativeProtein, tempProtein;
	LatticeStruct nativeLattice;
	XRayStruct nativeXRay;

	nativeLattice=lattice;
	nativeXRay=xray;
	nativeProtein=Protein;

	moveIntoUnitCell(Protein);
	getAlternateOrigins(Protein.spaceGroup, origins);
	Size=origins.size();
	readPdb(params.NativePdbFile, nativeProtein.Atoms);
	cout <<"NativePdbFile= "<<params.NativePdbFile<<endl;
	cout <<"natom= "<<nativeProtein.Atoms.size()<<endl;	
	calcSymmetryScattering(xray, Protein);
	calcIntensity(xray);
	calcSymmetryScattering(nativeXRay, nativeProtein);
	calcIntensity(nativeXRay);
	calcDensity(xray, lattice);
	calcDensity(nativeXRay, nativeLattice);
	minScore=1.0-calcDensityCorrelation(lattice, nativeLattice);
	for (int i=0;i<Size;i++)
	{
		tempProtein=nativeProtein;
		moveProteinFractional(tempProtein, origins[i][X], origins[i][Y], origins[i][Z]);
		correctForArbitraryAxis(tempProtein, Protein);
		moveIntoUnitCell(tempProtein);
		calcSymmetryScattering(nativeXRay, tempProtein);
		calcIntensity(nativeXRay);
		calcDensity(xray, lattice);
		calcDensity(nativeXRay, nativeLattice);
		score=1.0-calcDensityCorrelation(lattice, nativeLattice);
		if (score<minScore) minScore=score;
	}
	return minScore;
}

void calcPatterson(XRayStruct &xray)
{
	int Size=xray.patterson.real.size();
	int npoint=xray.q.size();
	Real dotproduct;
	for (int i=0;i<npoint;i++)
	{
		for (int j=0;j<Size;j++)
		{
			dotproduct=0;
			for (int k=0;k<3;k++)
			{
				dotproduct+=xray.q[i].Pos[k]*xray.patterson.u[j].Pos[k];
			}
			xray.patterson.real[j]+=xray.i[i]*cos(dotproduct);
			xray.patterson.imag[j]+=xray.i[i]*sin(dotproduct);
		}
	}
	for (int i=0;i<Size;i++)
	{
		//	cout <<"patterson.real["<<i<<"]= "<<xray.patterson.real[i]<<" imag= "<<xray.patterson.imag[i]<<endl;
	}
	//exit(EXIT_FAILURE);
}

Real calcPattersonOverlap(PattersonStruct &calc, PattersonStruct &experiment)
{
	int Size=calc.real.size();
	Real magCalc, magExperiment, sum=0;

	for (int i=0;i<Size;i++)
	{
		magCalc=calc.real[i]*calc.real[i]+calc.imag[i]*calc.imag[i];
		magExperiment=experiment.real[i]*experiment.real[i]+experiment.imag[i]*experiment.imag[i];
		sum+=sqrt(magCalc*magExperiment);
	}
	return sum;
}

Real calcPattersonOverlap(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay)
{
	Real overlap;
	calcScatteringComplex(Proteins, xray);
	calcPatterson(xray);
	overlap=calcPattersonOverlap(xray.patterson, expXRay.patterson);
	cout <<"overlap= "<<overlap<<endl;
	return overlap;
}

Real findMaxQ(vector<PosStruct> &q)
{
	int npoint=q.size();
	Real mag, max;

	max=calcMagnitude(q[0]);
	for (int i=1;i<npoint;i++)
	{
		mag=calcMagnitude(q[i]);
		if (mag>max) max=mag;
	}

	return max;
}

void calcExcludedScatteringSymmetry(ProteinStruct &Protein, XRayStruct &xray, LatticeStruct &lattice, XRayParamStruct &params)
{
	initializeLatticeForExcludedVolume(lattice, Protein, params);
	setExcludedVolumeDensity(lattice, Protein, params.BulkDensity);
	calcScatteringFast(lattice, xray);
}

void calcExcludedScattering(ProteinStruct &Protein, XRayStruct &xray, XRayParamStruct &params)
{
	int noperation=Protein.symmetryOperations.size();
	ProteinStruct tempProtein;
	LatticeStruct lattice;

	cout <<"In calcExcludedScattering"<<endl;
	cout <<"noperation= "<<noperation<<endl;
	zeroVector(xray.real);
	zeroVector(xray.imag);
	for (int i=0;i<noperation;i++)
	//for (int i=0;i<1;i++)
	{
		cout <<"i= "<<i<<endl;
		tempProtein=Protein;
		setSpaceGroupP1(tempProtein);
		applySymmetryOperation(tempProtein.Atoms, Protein.symmetryOperations[i]);
		calcExcludedScatteringSymmetry(tempProtein, xray, lattice, params);
	}
	multiplyByParallelPipped(lattice, xray);
	xray.excludedReal=xray.real;
	xray.excludedImag=xray.imag;
}

Real calcBFactorFix(XRayStruct &xray, XRayStruct &expXRay)
{
	int npoint1=xray.i.size();
	int npoint2=expXRay.i.size();
	int nq=xray.q.size();
	Real yIntercept=0, slope=0;
	Vector logRatio, q2;

	if (npoint1!=npoint2 || npoint1!=nq || npoint1==0)
	{
		string errorStr="npoint1= "+toStr(npoint1)+" npoint2= "+toStr(npoint2)+" nq= "+toStr(nq);
		error(errorStr, __LINE__, __FILE__);
	}

	SafeAlloc(logRatio, npoint1, "logRatio");
	SafeAlloc(q2, npoint1, "q2");
	//cout <<"q2\tlogRatio"<<endl;
	for (int i=0;i<npoint1;i++)
	{
		logRatio[i]=0.5*log(expXRay.i[i]/xray.i[i]);
		q2[i]=xray.qMag[i]*xray.qMag[i];
		//cout <<q2[i]<<"\t"<<logRatio[i]<<endl;
	}
	//cout <<endl;
	linearRegression(q2, logRatio, slope, yIntercept);
	slope*=-48.0*pi*pi;
	return slope;
}

void applyBFactor(XRayStruct &xray, Real bfactor)
{
	int npoint=xray.i.size();
	Real dwf;

	for (int i=0;i<npoint;i++)
	{
		dwf=calcDWF(xray.qMag[i], bfactor, xray.lookUp);
		xray.i[i]*=dwf*dwf;
	}
}

void fixBFactor(XRayStruct &xray, XRayStruct &expXRay)
{
	Real bfactor;
	//Real init_rfactor;
	//Real final_rfactor;

	//init_rfactor=RFactor(xray.i, expXRay.i);
	bfactor=calcBFactorFix(xray, expXRay);
	//cout <<"bfactor= "<<bfactor<<endl;
	applyBFactor(xray, bfactor);
	//final_rfactor=RFactor(xray.i, expXRay.i);
	//cout <<"init_rfactor= "<<init_rfactor<<" final_rfactor= "<<final_rfactor<<endl;
}

void average1DIntensity(XRayStruct &xray, Real binSize, Vector &intensity, Vector &q, Vector &num)
{
	int maxBin;
	int bin, nq=xray.q.size(), npoint=xray.i.size();
	Real max, mag;
	Vector total;

	if (nq!=npoint || npoint==0)
	{
		cout <<"ERROR: nq= "<<nq<<" npoint= "<<npoint<<endl;
		exit(EXIT_FAILURE);
	}

	max=findMaxQ(xray.q);
	maxBin=int(max/binSize+0.5);
	cout <<"maxQ= "<<max<<endl;
	SafeAlloc(intensity, maxBin+1, "intensity");
	SafeAlloc(q, maxBin+1, "q");
	SafeAlloc(total, maxBin+1, "total");
	SafeAlloc(num, maxBin+1, "num");
	for (int i=0;i<npoint;i++)
	{
		mag=calcMagnitude(xray.q[i]);
		bin=int(mag/binSize+0.5);
		total[bin]+=xray.i[i];
		num[bin]+=1.0;
		if (bin==29)
		{
			cout <<endl;
			cout <<"miller="<<endl;
			printVector(xray.miller[i].Pos);
			cout <<"i= "<<xray.i[i]<<endl;
			cout <<"total= "<<total[bin]<<" num= "<<num[bin]<<" bin= "<<bin<<endl;
			for (int j=0;j<NumAtomTypes;j++)
			{
				//cout <<"f= "<<xray.f[j][i]<<endl;
			}
			cout <<endl;
		}
	}

	for (int i=0;i<maxBin+1;i++)
	{
		if (num[i]!=0)
		{
			intensity[i]=total[i]/num[i];
			//cout <<"total= "<<total[i]<<" num= "<<num[i]<<" intensity= "<<intensity[i]<<endl;
		}
		q[i]=binSize*Real(i);
	}
}

void average1DIntensity(XRayStruct &xray, int maxBin, Vector &intensity, Vector &q)
{
	int bin, nq=xray.q.size(), npoint=xray.i.size();
	Real binSize, max, mag;
	Vector total, num;

	if (nq!=npoint || npoint==0)
	{
		cout <<"ERROR: nq= "<<nq<<" npoint= "<<npoint<<endl;
		exit(EXIT_FAILURE);
	}

	max=findMaxQ(xray.q);
	binSize=max/maxBin;
	cout <<"maxQ= "<<max<<endl;
	SafeAlloc(intensity, maxBin+1, "intensity");
	SafeAlloc(q, maxBin+1, "q");
	SafeAlloc(total, maxBin+1, "total");
	SafeAlloc(num, maxBin+1, "num");
	for (int i=0;i<npoint;i++)
	{
		mag=calcMagnitude(xray.q[i]);
		bin=int(mag/binSize+0.5);
		total[bin]+=xray.i[i];
		num[bin]+=1.0;
	}

	for (int i=0;i<maxBin;i++)
	{
		if (num[i]!=0)
		{
			intensity[i]=total[i]/num[i];
			cout <<"total= "<<total[i]<<" num= "<<num[i]<<" intensity= "<<intensity[i]<<endl;
		}
		q[i]=binSize*Real(i);
	}
}

void correctIntensity(XRayStruct &xray)
{
	int npoint=xray.i.size();
	int maxBin=xray.correction.size();
	int bin;
	Real mag;
	for (int i=0;i<npoint;i++)
	{
		mag=calcMagnitude(xray.q[i]);
		bin=int(mag/xray.correctionBinSize+0.5);
		if (bin<maxBin && xray.correction[bin]!=UNK_REAL)
		{
			xray.i[i]*=xray.correction[bin];
		}
	}
}

Real calcAveOccupancy(Vector &degOfFreedom, XRayParamStruct &params)
{
	int ndof=degOfFreedom.size();
	int start=params.NumCopies*6;
	Real sum=0;

	for (int i=start;i<ndof;i++)
	{
		sum+=degOfFreedom[i];
	}
	sum/=(ndof-start);
	return sum;
}

void correctOccupancyIntensity(XRayStruct &xray, Real aveOccupancy)
{
	int npoint=xray.i.size();
	int maxBin=xray.correction.size();
	int bin, occupancyBin;
	Real mag;

	occupancyBin=int(aveOccupancy/xray.OccupancyCorrection.occupancyBinSize+0.5);
	for (int i=0;i<npoint;i++)
	{
		mag=calcMagnitude(xray.q[i]);
		bin=int(mag/xray.OccupancyCorrection.correctionBinSize+0.5);
		if (bin<maxBin && xray.OccupancyCorrection.correction[occupancyBin][bin]!=UNK_REAL)
		{
			xray.i[i]*=xray.OccupancyCorrection.correction[occupancyBin][bin];
		}
	}
}

void correctOccupancyIntensity(XRayStruct &xray, Vector &degOfFreedom, XRayParamStruct &params)
{
	Real aveOccupancy;

	aveOccupancy=calcAveOccupancy(degOfFreedom, params);
	correctOccupancyIntensity(xray, aveOccupancy);
}

void correctOccupancyIntensity(XRayStruct &xray, ProteinStruct &Protein, XRayParamStruct &params)
{
	Real aveOccupancy;

	aveOccupancy=calcAveOccupancy(Protein);
	correctOccupancyIntensity(xray, aveOccupancy);
}

bool isSystematicAbsence(XRayStruct &xray, int index, ProteinStruct &Protein)
{
	//Checks if a given miller index will give zero intensity for a given 
	//space group.
	//Need to replace with something better.
	int ntest=10;
	Real d=0.001;
	AtomStruct tempAtom;
	ProteinStruct tempProtein;
	XRayStruct tempXRay;

	//tempXRay=xray;
	tempXRay.lookUp=xray.lookUp;
	tempXRay.maxH=xray.maxH;
	tempXRay.maxK=xray.maxK;
	tempXRay.maxL=xray.maxL;
	tempXRay.minH=xray.minH;
	tempXRay.minK=xray.minK;
	tempXRay.minL=xray.minL;
	tempXRay.i.clear();
	tempXRay.q.clear();
	tempXRay.qMag.clear();
	tempXRay.real.clear();
	tempXRay.imag.clear();
	tempXRay.miller.clear();
	tempXRay.f.clear();
	SafePushBack(tempXRay.i, Real(0), "i");
	SafePushBack(tempXRay.q, xray.q[index], "q");
	SafePushBack(tempXRay.qMag, xray.qMag[index], "qMag");
	SafePushBack(tempXRay.real, xray.real[index], "real");
	SafePushBack(tempXRay.imag, xray.imag[index], "imag");
	SafePushBack(tempXRay.miller, xray.miller[index], "miller");
	Safe2DAlloc(tempXRay.f, NumAtomTypes, 1, "f");
	for (int i=0;i<NumAtomTypes;i++)
	{
		tempXRay.f[i][0]=1.0;
	}
	tempAtom.Occupancy=1.0;
	tempProtein=Protein;
	tempProtein.Atoms.clear();
	SafePushBack(tempProtein.Atoms, tempAtom, "Atoms");

	for (int i=0;i<ntest;i++)
	{
		tempProtein.Atoms[0].x=randDouble(0, tempProtein.XBoxLength);
		tempProtein.Atoms[0].y=randDouble(0, tempProtein.YBoxLength);
		tempProtein.Atoms[0].z=randDouble(0, tempProtein.ZBoxLength);
		calcSymmetryScattering(tempXRay, tempProtein);
		if (tempXRay.i[0]>d) return false;
	}
	return true;
}

void removeSystematicAbsences(XRayStruct &xray, ProteinStruct &Protein)
{
	int ni=xray.i.size();
	int nmiller=xray.miller.size();
	int nq=xray.q.size();
	int nF=xray.F.size();
	int nreal=xray.real.size();
	int nimag=xray.imag.size();
	int nqmag=xray.qMag.size();
	int nerror=xray.Ferror.size();
	int ncentric=xray.centric.size();
	XRayStruct xrayOut;

	cout <<"In removeSystematicAbsences"<<endl;
	cout <<"ni= "<<xray.i.size()<<endl;
	cout <<"nmiller= "<<xray.miller.size()<<endl;
	cout <<"nq= "<<xray.q.size()<<endl;
	cout <<"nF= "<<xray.F.size()<<endl;
	cout <<"nreal= "<<xray.real.size()<<endl;
	cout <<"nimag= "<<xray.imag.size()<<endl;
	cout <<"nqmag= "<<xray.qMag.size()<<endl;
	cout <<"nerror= "<<xray.Ferror.size()<<endl;
	cout <<"ncentric= "<<xray.centric.size()<<endl;

	Safe2DAlloc(xrayOut.f, NumAtomTypes, 0, "f");
	for (int i=0;i<nmiller;i++)
	{
		if (!isSystematicAbsence(xray, i, Protein))
		{
			if (i<ni) SafePushBack(xrayOut.i, xray.i[i], "i");
			if (i<nmiller) SafePushBack(xrayOut.miller, xray.miller[i], "miller");
			if (i<nq) SafePushBack(xrayOut.q, xray.q[i], "q");
			if (i<nF) SafePushBack(xrayOut.F, xray.F[i], "F");
			if (i<nreal) SafePushBack(xrayOut.real, xray.real[i], "real");
			if (i<nimag) SafePushBack(xrayOut.imag, xray.imag[i], "imag");
			if (i<nqmag) SafePushBack(xrayOut.qMag, xray.qMag[i], "qMag");
			if (i<nerror) SafePushBack(xrayOut.Ferror, xray.Ferror[i], "Ferror");
			if (i<ncentric) SafePushBack(xrayOut.centric, bool(xray.centric[i]), "centric");
			for (int j=0;j<NumAtomTypes;j++)
			{
				SafePushBack(xrayOut.f[j], xray.f[j][i], "f");
			}
		}
	}
	xray.miller=xrayOut.miller;
	xray.i=xrayOut.i;
	xray.q=xrayOut.q;
	xray.F=xrayOut.F;
	xray.f=xrayOut.f;
	xray.real=xrayOut.real;
	xray.imag=xrayOut.imag;
	xray.qMag=xrayOut.qMag;
	xray.Ferror=xrayOut.Ferror;
	xray.centric=xrayOut.centric;
	cout <<"ni= "<<xray.i.size()<<endl;
	cout <<"nmiller= "<<xray.miller.size()<<endl;
	cout <<"nq= "<<xray.q.size()<<endl;
	cout <<"nF= "<<xray.F.size()<<endl;
	cout <<"nreal= "<<xray.real.size()<<endl;
	cout <<"nimag= "<<xray.imag.size()<<endl;
	cout <<"nqmag= "<<xray.qMag.size()<<endl;
	cout <<"nerror= "<<xray.Ferror.size()<<endl;
	cout <<"ncentric= "<<xray.centric.size()<<endl;
}

void removeNegativeIntensity(XRayStruct &xray)
{
	int ni=xray.i.size();
	int nmiller=xray.miller.size();
	int nq=xray.q.size();
	int nF=xray.F.size();
	XRayStruct xrayOut;

	for (int i=0;i<ni;i++)
	{
		cout <<"i= "<<i<<endl;
		if (xray.i[i]>0)
		{
			SafePushBack(xrayOut.i, xray.i[i], "i");
			if (i<nmiller) SafePushBack(xrayOut.miller, xray.miller[i], "miller");
			if (i<nq) SafePushBack(xrayOut.q, xray.q[i], "q");
			if (i<nF) SafePushBack(xrayOut.F, xray.F[i], "F");
		}
	}
	xray.miller=xrayOut.miller;
	xray.i=xrayOut.i;
	xray.q=xrayOut.q;
	xray.F=xrayOut.F;
}

void removeNegativeMiller(XRayStruct &xray)
{
	int nmiller=xray.miller.size();

	for (int i=0;i<nmiller;i++)
	{
		if (xray.miller[i].Pos[X]<0)
		{
			xray.miller[i].Pos[X]=-xray.miller[i].Pos[X];
			xray.miller[i].Pos[Y]=-xray.miller[i].Pos[Y];
			xray.miller[i].Pos[Z]=-xray.miller[i].Pos[Z];
		}
	}
}

void compareAverageIntensities(XRayStruct &xray, XRayStruct &expXRay)
{
	int maxBin=20;
	Vector q, intensity, expIntensity;
	XRayStruct xrayMatching, expXRayMatching;	
	findCorrespondingPoints(xray, expXRay, xrayMatching, expXRayMatching);
	average1DIntensity(xrayMatching, maxBin, intensity, q);
	average1DIntensity(expXRayMatching, maxBin, expIntensity, q);

	cout <<"q\tIobs\tIcalc"<<endl;
	for (int i=0;i<maxBin;i++)
	{
		cout <<q[i]<<"\t"<<expIntensity[i]<<"\t"<<intensity[i]<<endl;
	}
}

bool findCorrepondingPoint(PosStruct &expPos, XRayStruct &xray, int &xbin, int &ybin, int &zbin)
{
	int MaxXBin, MaxYBin, MaxZBin;
	Real d=0.1;
	xbin=int(expPos.Pos[X]);
	ybin=int(expPos.Pos[Y]);
	zbin=int(expPos.Pos[Z]);
	Get3DVectorSize(xray.miller3d, MaxXBin, MaxYBin, MaxZBin, "xray.miller3d");
	if (xbin>=0 && xbin<MaxXBin && ybin>=0 && ybin<MaxYBin && zbin>=0 && zbin<MaxZBin)
	{
		if (vectorEqual(xray.miller3d[xbin][ybin][zbin].Pos, expPos.Pos, d))
		{
			return true;
		}
	}
	return false;
}

void setMaxMiller(XRayStruct &xray)
{
	int nmiller=xray.miller.size();

	if (nmiller==0)
	{
		error("nmiller= 0", __LINE__, __FILE__);	
	}
	xray.maxH=int(xray.miller[0].Pos[X]);
	xray.maxK=int(xray.miller[0].Pos[Y]);
	xray.maxL=int(xray.miller[0].Pos[Z]);

	for (int i=0;i<nmiller;i++)
	{
		if (xray.miller[i].Pos[X]>xray.maxH) 
		{
			xray.maxH=int(xray.miller[i].Pos[X]);
		}
		if (xray.miller[i].Pos[Y]>xray.maxK) 
		{
			xray.maxK=int(xray.miller[i].Pos[Y]);
		}
		if (xray.miller[i].Pos[Z]>xray.maxL) 
		{
			xray.maxL=int(xray.miller[i].Pos[Z]);
		}
	}
}

void setMinMiller(XRayStruct &xray)
{
	int nmiller=xray.miller.size();

	if (nmiller==0)
	{
		error("nmiller= 0", __LINE__, __FILE__);	
	}
	xray.minH=int(xray.miller[0].Pos[X]);
	xray.minK=int(xray.miller[0].Pos[Y]);
	xray.minL=int(xray.miller[0].Pos[Z]);
	for (int i=0;i<nmiller;i++)
	{
		if (xray.miller[i].Pos[X]<xray.minH) 
		{
			xray.minH=int(xray.miller[i].Pos[X]);
		}
		if (xray.miller[i].Pos[Y]<xray.minK) 
		{
			xray.minK=int(xray.miller[i].Pos[Y]);
		}
		if (xray.miller[i].Pos[Z]<xray.minL) 
		{
			xray.minL=int(xray.miller[i].Pos[Z]);
		}
	}
}

void makeThreeDtoOneD(XRayStruct &xray, vector< vector< vector<int> > > &threeDtoOneD)
{
	int h, k, l;
	int nmiller=xray.miller.size();

	setMaxMiller(xray);
	setMinMiller(xray);
	Safe3DAlloc(threeDtoOneD, UNK_INT, xray.maxH-xray.minH+1, xray.maxK-xray.minK+1, xray.maxL-xray.minL+1, "threeDtoOneD");
	for (int i=0;i<nmiller;i++)
	{
		h=int(xray.miller[i].Pos[X]-xray.minH);
		k=int(xray.miller[i].Pos[Y]-xray.minK);
		l=int(xray.miller[i].Pos[Z]-xray.minL);
		if (h>=0 && k>=0 && l>=0) threeDtoOneD[h][k][l]=i;
		else
		{
			string errorStr="h= "+toStr(h)+" k= "+toStr(k);
			errorStr+=" l= "+toStr(l);
			error(errorStr, __LINE__, __FILE__);
		}
	}
}

void findCorrespondingPoints(XRayStruct &xrayIn, XRayStruct &expXRay, XRayStruct &xrayOut, XRayStruct &expXRayOut)
{
	int Size, index, nsym;
	int nerror, nerror2;
	int h, k, l;
	int ni, nq, nf, nF, nqmag, nreal, nimag, nphase;
	int ni2, nq2, nf2, nF2, nqmag2, nreal2, nimag2, nphase2, nmiller2;
	vector< vector< vector<int> > > threeDtoOneD;
	Vector v;

	//cout <<"In findCorrespondingPoints"<<endl;

	ni=xrayIn.i.size();
	nq=xrayIn.q.size();
	nF=xrayIn.F.size();
	nqmag=xrayIn.qMag.size();
	nreal=xrayIn.real.size();
	nimag=xrayIn.imag.size();
	nphase=xrayIn.phase.size();
	nerror=xrayIn.Ferror.size();
	nsym=xrayIn.symReal.size();
	Get2DVectorSize(xrayIn.f, Size, nf);	

	ni2=expXRay.i.size();
	nq2=expXRay.q.size();
	nF2=expXRay.F.size();
	nqmag2=expXRay.qMag.size();
	nreal2=expXRay.real.size();
	nimag2=expXRay.imag.size();
	nphase2=expXRay.phase.size();
	nmiller2=expXRay.miller.size();
	nerror2=expXRay.Ferror.size();
	Get2DVectorSize(expXRay.f, Size, nf2);	

	//cout <<"nerror= "<<nerror<<" nerror2= "<<nerror2<<endl;	
	//cout <<"ni= "<<ni<<endl;

	Safe2DAlloc(xrayOut.f, NumAtomTypes, 0, "f");
	Safe2DAlloc(expXRayOut.f, NumAtomTypes, 0, "f");
	Safe2DAlloc(xrayOut.symReal, nsym, 0, "symReal");
	Safe2DAlloc(xrayOut.symImag, nsym, 0, "symImag");
	makeThreeDtoOneD(xrayIn, threeDtoOneD);
	for (int i=0;i<nmiller2;i++)
	{
		h=int(expXRay.miller[i].Pos[X]-xrayIn.minH);
		k=int(expXRay.miller[i].Pos[Y]-xrayIn.minK);
		l=int(expXRay.miller[i].Pos[Z]-xrayIn.minL);
		//cout <<"h= "<<h<<" k= "<<k<<" l= "<<l<<endl;
		//cout <<"maxH= "<<xrayIn.maxH<<" maxK= "<<xrayIn.maxK<<" maxL= "<<xrayIn.maxL<<endl;
		//cout <<"minH= "<<xrayIn.minH<<" minK= "<<xrayIn.minK<<" minL= "<<xrayIn.minL<<endl;
		if (h<=xrayIn.maxH-xrayIn.minH && k<=xrayIn.maxK-xrayIn.minK && l<=xrayIn.maxL-xrayIn.minL && h>=0 && k>=0 && l>=0)
		{
			index=threeDtoOneD[h][k][l];
			//cout <<"index= "<<index<<endl;
			if (index!=UNK_INT)
			{
				if (index<ni) SafePushBack(xrayOut.i, xrayIn.i[index], "xrayOut.i");
				if (index<nq) SafePushBack(xrayOut.q, xrayIn.q[index], "xrayOut.q");
				if (index<nF) SafePushBack(xrayOut.F, xrayIn.F[index], "xrayOut.F");
				if (index<nqmag) SafePushBack(xrayOut.qMag, xrayIn.qMag[index], "xrayOut.qMag");
				if (index<nreal) SafePushBack(xrayOut.real, xrayIn.real[index], "xrayOut.real");
				if (index<nimag) SafePushBack(xrayOut.imag, xrayIn.imag[index], "xrayOut.imag");
				if (index<nphase) SafePushBack(xrayOut.phase, xrayIn.phase[index], "xrayOut.phase");
				if (index<nerror) SafePushBack(xrayOut.Ferror, xrayIn.Ferror[index], "xrayOut.Ferror");
				for (int j=0;j<NumAtomTypes;j++)
				{
					if (index<nf) SafePushBack(xrayOut.f[j], xrayIn.f[j][index], "xrayOut.f");	
				}
				for (int j=0;j<nsym;j++)
				{
					SafePushBack(xrayOut.symReal[j], xrayIn.symReal[j][index], "xrayOut.symReal");
					SafePushBack(xrayOut.symImag[j], xrayIn.symImag[j][index], "xrayOut.symReal");
				}
				SafePushBack(xrayOut.miller, xrayIn.miller[index], "xrayOut.miller");
				if (i<ni2) SafePushBack(expXRayOut.i, expXRay.i[i], "xrayOut.i");
				if (i<nq2) SafePushBack(expXRayOut.q, expXRay.q[i], "xrayOut.q");
				if (i<nF2) SafePushBack(expXRayOut.F, expXRay.F[i], "xrayOut.F");
				if (i<nqmag2) SafePushBack(expXRayOut.qMag, expXRay.qMag[i], "xrayOut.qMag");
				if (i<nreal2) SafePushBack(expXRayOut.real, expXRay.real[i], "xrayOut.real");
				if (i<nimag2) SafePushBack(expXRayOut.imag, expXRay.imag[i], "xrayOut.imag");
				if (i<nphase2) SafePushBack(expXRayOut.phase, expXRay.phase[i], "xrayOut.phase");
				if (i<nerror2) SafePushBack(expXRayOut.Ferror, expXRay.Ferror[i], "xrayOut.Ferror");
				for (int j=0;j<NumAtomTypes;j++)
				{
					if (i<nf2) SafePushBack(expXRayOut.f[j], expXRay.f[j][i], "xrayOut.f");	
				}
				SafePushBack(expXRayOut.miller, expXRay.miller[i], "xrayOut.miller");
			}
		}
	}
	xrayOut.lookUp=xrayIn.lookUp;
	xrayOut.maxH=xrayIn.maxH;
	xrayOut.maxK=xrayIn.maxK;
	xrayOut.maxL=xrayIn.maxL;

	nerror=xrayOut.Ferror.size();
	nerror2=expXRayOut.Ferror.size();

	//cout <<"nerror= "<<nerror<<" nerror2= "<<nerror2<<endl;	
	//cout <<"niOut= "<<xrayOut.i.size()<<endl;
//	endProgram(__LINE__, __FILE__);
}

bool isMatching(XRayStruct &xray1, XRayStruct &xray2)
{
	int nmiller1=xray1.miller.size();
	int nmiller2=xray2.miller.size();
//	timeval start, end;

	//gettimeofday(&start, NULL);
	if (nmiller1!=nmiller2) return false;
	for (int i=0;i<nmiller1;i++)
	{
		for (int j=0;j<3;j++)
		{
			if (xray1.miller[i].Pos[j]!=xray2.miller[i].Pos[j])
			{
				return false;
			}
		}
	}
	//gettimeofday(&end, NULL);
	//cout <<"isMatching took "<<calcTimeDiff(start, end)<<endl;
	return true;
}

void removeNonoverlapping(XRayStruct &xray, XRayStruct &expXRay)
{
	XRayStruct xrayMatching, xrayExpMatching;
	xrayMatching.score=xray.score;
	xrayMatching.stdDev=xray.stdDev;
	xrayMatching.scoreIndex=xray.scoreIndex;
	xrayMatching.correction=xray.correction;
	xrayMatching.correctionBinSize=xray.correctionBinSize;
	xrayMatching.OccupancyCorrection.occupancyBinSize=xray.OccupancyCorrection.occupancyBinSize;
	xrayMatching.OccupancyCorrection.correctionBinSize=xray.OccupancyCorrection.occupancyBinSize;
	xrayMatching.OccupancyCorrection.correction=xray.OccupancyCorrection.correction;
	xrayMatching.hBinSize=xray.hBinSize;
	xrayMatching.kBinSize=xray.kBinSize;
	xrayMatching.lBinSize=xray.lBinSize;
	xrayMatching.realContinuous=xray.realContinuous;
	xrayMatching.imagContinuous=xray.imagContinuous;
	xrayMatching.maxContinuousH=xray.maxContinuousH;
	xrayMatching.maxContinuousK=xray.maxContinuousK;
	xrayMatching.maxContinuousL=xray.maxContinuousL;
	xrayExpMatching.correction=expXRay.correction;
	xrayExpMatching.correctionBinSize=expXRay.correctionBinSize;
	xrayExpMatching.OccupancyCorrection.occupancyBinSize=expXRay.OccupancyCorrection.occupancyBinSize;
	xrayExpMatching.OccupancyCorrection.correctionBinSize=expXRay.OccupancyCorrection.occupancyBinSize;
	xrayExpMatching.OccupancyCorrection.correction=expXRay.OccupancyCorrection.correction;
	xrayExpMatching.hBinSize=expXRay.hBinSize;
	xrayExpMatching.kBinSize=expXRay.kBinSize;
	xrayExpMatching.lBinSize=expXRay.lBinSize;
	xrayExpMatching.realContinuous=expXRay.realContinuous;
	xrayExpMatching.imagContinuous=expXRay.imagContinuous;
	xrayExpMatching.maxContinuousH=expXRay.maxContinuousH;
	xrayExpMatching.maxContinuousK=expXRay.maxContinuousK;
	xrayExpMatching.maxContinuousL=expXRay.maxContinuousL;
	findCorrespondingPoints(xray, expXRay, xrayMatching, xrayExpMatching);
	xray=xrayMatching;
	expXRay=xrayExpMatching;
}

void combineMillerIndexes(XRayStruct &xray, int index, vector<bool> &combine, XRayStruct &xrayOut)
{
	int nmiller=xray.miller.size();
	int nq=xray.q.size();
	Real sum=0, count=0;

	for (int i=0;i<nmiller;i++)
	{
		if (abs(xray.miller[i].Pos[X])==abs(xray.miller[index].Pos[X]) && abs(xray.miller[i].Pos[Y])==abs(xray.miller[index].Pos[Y]) && abs(xray.miller[i].Pos[Z])==abs(xray.miller[index].Pos[Z]))
		{
			sum+=xray.i[i];
			count+=1.0;
			combine[i]=true;
		}
	}
	for (int i=0;i<3;i++)
	{
		xray.miller[index].Pos[i]=abs(xray.miller[index].Pos[i]);
	}
	SafePushBack(xrayOut.miller, xray.miller[index], "miller");
	if (nq>index) SafePushBack(xrayOut.q, xray.q[index], "miller");
	SafePushBack(xrayOut.i, sum/count, "i");
	SafePushBack(xrayOut.real, sqrt(sum/count), "real");
	SafePushBack(xrayOut.imag, Real(0), "imag");
	SafePushBack(xrayOut.phase, Real(0), "phase");

}

void combineMillerIndexes(XRayStruct &xray)
{
	vector<bool> combined;
	int nmiller=xray.miller.size();
	XRayStruct xrayOut;

	SafeAlloc(combined, nmiller, "combined");
	for (int i=0;i<nmiller;i++)
	{
		if (!combined[i]) 
		{
			combineMillerIndexes(xray, i, combined, xrayOut);
		}
	}
	xray.miller=xrayOut.miller;
	xray.i=xrayOut.i;
	xray.real=xrayOut.real;
	xray.imag=xrayOut.imag;
	xray.phase=xrayOut.phase;
	//xray.q=xrayOut.miller; //delete

	//printXRayJV("/home/jouko/project/mr/1PWA_test.txt", xray);

	nmiller=xray.miller.size();
}

Real calcRFactor(XRayStruct &xray, XRayStruct &expXRay, string xRayScoreType)
{
	Real score;
	XRayStruct xrayMatching, xrayExpMatching;

	if (!isMatching(xray, expXRay))
	{
		findCorrespondingPoints(xray, expXRay, xrayMatching, xrayExpMatching);
		score=quantifyMatch(xrayMatching, xrayExpMatching, xRayScoreType);
	}
	else
	{
		score=quantifyMatch(xray, expXRay, xRayScoreType);
	}
	return score;
}

Real calcRFactor(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params)
{
	Real score;
	XRayStruct xrayMatching, xrayExpMatching;

	if (!isMatching(xray, expXRay))
	{
		findCorrespondingPoints(xray, expXRay, xrayMatching, xrayExpMatching);
		score=quantifyWeightedMatch(xrayMatching, xrayExpMatching, params);
	}
	else
	{
		score=quantifyWeightedMatch(xray, expXRay, params);
	}

	return score;
}

Real calcExperimentalScaleFactor(XRayStruct &xray, XRayStruct &expXRay)
{
	XRayStruct xrayMatching, xrayExpMatching;
	findCorrespondingPoints(xray, expXRay, xrayMatching, xrayExpMatching);
	return calcScaleForSqrDiff(xrayMatching.i, xrayExpMatching.i);
}

void filterOutHighResolution(XRayStruct &xray, int hmax, int kmax, int lmax)
{
	int npoint=xray.miller.size();
	int nreal=xray.real.size();
	int nimag=xray.imag.size();
	int ni=xray.i.size();
	int nF=xray.F.size();
	int nphase=xray.phase.size();
	int nError=xray.Ferror.size();
	int nerror=xray.error.size();
	int ncentric=xray.centric.size();
	XRayStruct xrayOut;
	cout <<"In filter npoint= "<<npoint<<endl;

	cout <<"npoint= "<<npoint<<endl;
	cout <<"nreal= "<<nreal<<endl;
	cout <<"nimag= "<<nimag<<endl;
	cout <<"ni= "<<ni<<endl;
	cout <<"nF= "<<nF<<endl;
	cout <<"nphase= "<<nphase<<endl;
	cout <<"nerror= "<<nerror<<endl;
	cout <<"ncentric= "<<ncentric<<endl;

	for (int i=0;i<npoint;i++)
	{
		if (abs(xray.miller[i].Pos[X])<=Real(hmax) && abs(xray.miller[i].Pos[Y])<=Real(kmax) && abs(xray.miller[i].Pos[Z])<=Real(lmax))
		{
			SafePushBack(xrayOut.miller, xray.miller[i], "miller");
			if (i<nreal) SafePushBack(xrayOut.real, xray.real[i], "real");
			if (i<nimag) SafePushBack(xrayOut.imag, xray.imag[i], "imag");
			if (i<ni) SafePushBack(xrayOut.i, xray.i[i], "i");
			if (i<nF) SafePushBack(xrayOut.F, xray.F[i], "F");
			if (i<nphase) SafePushBack(xrayOut.phase, xray.phase[i], "phase");
			if (i<nError) SafePushBack(xrayOut.Ferror, xray.Ferror[i], "Ferror");
			if (i<nerror) SafePushBack(xrayOut.error, xray.error[i], "error");
			if (i<ncentric) SafePushBack(xrayOut.centric, bool(xray.centric[i]), "centric");
		}
	}
	
	xray=xrayOut;

	cout <<"npoint= "<<xray.miller.size()<<endl;
	cout <<"nreal= "<<xray.real.size()<<endl;
	cout <<"nimag= "<<xray.imag.size()<<endl;
	cout <<"ni= "<<xray.i.size()<<endl;
	cout <<"nF= "<<xray.F.size()<<endl;
	cout <<"nphase= "<<xray.phase.size()<<endl;
	cout <<"nerror= "<<xray.Ferror.size()<<endl;
	cout <<"ncentric= "<<xray.centric.size()<<endl;
	cout <<"Leaving filterOutHighResolution"<<endl;
}

void normalizeXRay(Vector &intensity)
{
	int npoint=intensity.size();
	Real max;

	max=findMax(intensity);
	for (int i=0;i<npoint;i++)
	{
		intensity[i]/=max;
	}
}

void printSym(XRayStruct &xray)
{
	int nprot, nsym, npoint;

	Get3DVectorSize(xray.complexSymTempReal, nprot, nsym, npoint);
	for (int i=0;i<nprot;i++)
	{
		for (int j=0;j<nsym;j++)
		{
			for (int k=0;k<npoint;k++)
			{
				cout <<"real= "<<xray.complexSymTempReal[i][j][k]
				<<" imag= "<<xray.complexSymTempImag[i][j][k]
				<<" i= "<<i<<" j= "<<j<<" k= "<<k
				<<" miller= "<<xray.miller[k].Pos[X]
				<<" "<<xray.miller[k].Pos[Y]
				<<" "<<xray.miller[k].Pos[Z]<<endl;
			}
		}
	}
}

void printXRayJV(string XRayOutput, XRayStruct &xray)
{
	ofstream file;
	int npoint=xray.i.size();
	int nmiller=xray.miller.size();
	int nq=xray.q.size();
	int nreal=xray.real.size();
	int nimag=xray.imag.size();
	int nphase=xray.phase.size();
	AddIndexToFile(XRayOutput);
	OpenFileForWriting(XRayOutput, file);

	cout <<"In printXRay npoint= "<<npoint<<endl;
	/*
	   for (int i=0;i<npoint;i++)
	   {
	   file <<setiosflags(ios::right)<<setw(4)<<xray.miller[i].Pos[X]
	   <<setiosflags(ios::right)<<setw(4)<<xray.miller[i].Pos[Y]
	   <<setiosflags(ios::right)<<setw(4)<<xray.miller[i].Pos[Z]
	   <<setiosflags(ios::right)<<setw(11)<<setprecision(3)
	   <<setiosflags(ios::fixed)<<xray.i[i]<<"\t";
	   file.unsetf(ios::fixed);
	   file <<resetiosflags(ios::right)<<xray.phase[i]<<xray.real[i]
	   <<"\t"<<xray.imag[i]<<"\t"<<xray.q[i].Pos[X]<<"\t"
	   <<xray.q[i].Pos[Y]<<"\t"<<xray.q[i].Pos[Z]<<endl;
	   }	
	 */	
	if (npoint!=nmiller)
	{
		cout <<"ERROR: npoint= "<<npoint<<" nmiller= "<<nmiller<<endl;
		exit(EXIT_FAILURE);
	}
	cout <<"npoint= "<<npoint<<" nmiller= "<<nmiller<<" nq= "<<nq<<" nreal= "<<nreal<<" nimag= "<<nimag<<" nphase= "<<nphase<<endl;
	for (int i=0;i<npoint;i++)
	{
		file <<xray.miller[i].Pos[X]<<"\t"
			<<xray.miller[i].Pos[Y]<<"\t"
			<<xray.miller[i].Pos[Z]<<"\t"
			<<xray.i[i]<<"\t"<<xray.phase[i]<<"\t"<<xray.real[i]<<"\t"
			<<xray.imag[i]<<"\t"<<xray.q[i].Pos[X]<<"\t"
			<<xray.q[i].Pos[Y]<<"\t"<<xray.q[i].Pos[Z]<<endl;
	}
	cout <<"Leaving printXRay"<<endl;	
}

void convertISigmaToFSigma(XRayStruct &xray)
{
	int ni=xray.i.size();
	int nerror=xray.error.size();

	if (ni==nerror)
	{
		SafeAlloc(xray.Ferror, ni, "Ferror");
		for (int i=0;i<ni;i++)
		{
			xray.Ferror[i]=xray.error[i]/(2.0*sqrt(xray.i[i]));
		}
	}

}

void printCifFile(string XRayOutput, XRayStruct &xray, ProteinStruct &Protein)
{
	ofstream file;
	int npoint=xray.i.size();
	int nmiller=xray.miller.size();
	int nerror=xray.error.size();
	int nError=xray.Ferror.size();
	int nphase=xray.phase.size();
	int width;
	cout <<"In printCifFile"<<endl;
	cout <<"XRayOutput= "<<XRayOutput<<endl;
	if (npoint!=nmiller || npoint==0)
	{
		string errorStr;
		errorStr="npoint= "+toStr(npoint);
		errorStr+=" nphase= "+toStr(nphase);
		errorStr+=" nmiller= "+toStr(nmiller);	
		error(errorStr, __LINE__, __FILE__);
	}
	SafeAlloc(xray.error, Real(1.0), npoint,  "error");
	convertISigmaToFSigma(xray);
	nerror=xray.error.size();
	nError=xray.Ferror.size();
	OpenFileForWriting(XRayOutput, file);

	file <<"data_rcalcsf"<<endl;
	file <<"#"<<endl;
	file <<"_cell.length_a       "<<Protein.XBoxLength<<endl
		<<"_cell.length_b       "<<Protein.YBoxLength<<endl
		<<"_cell.length_c       "<<Protein.ZBoxLength<<endl
		<<"_cell.angle_alpha   "<<Protein.alpha<<endl
		<<"_cell.angle_beta    "<<Protein.beta<<endl
		<<"_cell.angle_gamma   "<<Protein.gamma<<endl;

	file <<"#"<<endl;
	file <<"loop_"<<endl
		<<"_refln.crystal_id"<<endl
		<<"_refln.wavelength_id"<<endl
		<<"_refln.scale_group_code"<<endl
		<<"_refln.index_h"<<endl
		<<"_refln.index_k"<<endl
		<<"_refln.index_l"<<endl
		<<"_refln.status"<<endl
		//		<<"_refln.intensity_calc"<<endl;
		//	if (nerror>0) file <<"_refln.intensity_sigma"<<endl;
		<<"_refln.F_meas_au"<<endl;
	//if (nerror>0) file <<"_refln.intensity_sigma"<<endl;
	if (nError>0) file <<"_refln.F_meas_sigma_au"<<endl;
	//file <<"_refln.phase_calc"<<endl;

	//width=int(log(findMax(xray.i)))+2;
	width=int(log(findMax(xray.F)))+1;
	for (int i=0;i<npoint;i++)
	{
		//cout <<"i= "<<i<<" free= "<<xray.free[i]<<" x= "<<xray.miller[i].Pos[X]<<" y= "<<xray.miller[i].Pos[Y]<<" z= "<<xray.miller[i].Pos[Z]<<endl;
	}
	for (int i=0;i<npoint;i++)
	{
		if (xray.i[i]>0)
		{
			file <<"1 1 1"<<setiosflags(ios::right)<<setw(5)<<xray.miller[i].Pos[X]
				<<setiosflags(ios::right)<<setw(5)<<xray.miller[i].Pos[Y]
				<<setiosflags(ios::right)<<setw(5)<<xray.miller[i].Pos[Z];

			if (randDouble(0, 1.0)<0.048) file <<setiosflags(ios::right)<<" f  ";
			else file <<setiosflags(ios::left)<<" o  ";
			file <<setiosflags(ios::fixed)
				<<setiosflags(ios::right)<<setw(width)<<setprecision(3)
				<<setiosflags(ios::right)<<sqrt(xray.i[i]);
				//if (xray.free[i]) file <<"f ";
				//else file <<"o ";
			if (nerror>0) 
			{
				file <<" ";
				file <<setiosflags(ios::fixed)<<xray.error[i]<<" ";
			}
			if (nError>0) 
			{
				//file <<setiosflags(ios::fixed)<<xray.Ferror[i]<<" ";
			}
			if (nphase>i)
			{
				//file.unsetf(ios::fixed);
				//file <<resetiosflags(ios::right)<<xray.phase[i]*RAD_TO_DEGREE;
			}
			file.unsetf(ios::fixed);
			file <<endl;
		}
	}
	file <<"#END OF REFLECTIONS"<<endl;
}

void printXRay(string XRayOutput, XRayStruct &xray, ProteinStruct &Protein)
{
	string ext;

	ext=GetExtension(XRayOutput);

	if (ext=="txt") printXRayJV(XRayOutput,  xray);
	else if (ext=="cif") printCifFile(XRayOutput, xray, Protein);
	else
	{
		string errorStr="Unrecognized extension for x-ray file "+ext;
		errorStr+=" XRayOutput= "+XRayOutput;
		error(errorStr, __LINE__, __FILE__);
	}
}
