# include <iostream>
# include <vector>

using namespace std;


# include "../../LibraryFiles/Structures.h"
# include "../../LibraryFiles/VectorManip.h"
# include "../../LibraryFiles/IOUtils.h"
# include "../../LibraryFiles/ReadPdb.h"
# include "../../LibraryFiles/PrintPdb.h"
# include "../../LibraryFiles/MathUtils.h"
# include "../../LibraryFiles/AtomUtils.h"
# include "../../LibraryFiles/ASA.h"
# include "../../LibraryFiles/MinMax.h"
# include "LatticeStruct.h"
# include "LatticeUtils.h"
# include "Params.h"

LatticeStruct& LatticeStruct::operator = (const LatticeStruct &lattice)
{
	cubeSize=lattice.cubeSize;
	a=lattice.a;
	b=lattice.b;
	c=lattice.c;
	alpha=lattice.alpha;
	beta=lattice.beta;
	gamma=lattice.gamma;
	xCubeLength=lattice.xCubeLength;
	yCubeLength=lattice.yCubeLength;
	zCubeLength=lattice.zCubeLength;
	firstPos=lattice.firstPos;
	density=lattice.density;
	real=lattice.real;
	imag=lattice.imag;


	return *this;
}

void indexToPosition(LatticeStruct &lattice, int xbin, int ybin, int zbin, Real &x, Real &y, Real &z)
{
	x=lattice.firstPos.Pos[X]+Real(xbin)*lattice.xCubeLength;
	y=lattice.firstPos.Pos[Y]+Real(ybin)*lattice.yCubeLength;
	z=lattice.firstPos.Pos[Z]+Real(zbin)*lattice.zCubeLength;
}

void setExcludedVolumeDensity(LatticeStruct &lattice, ProteinStruct &Protein, Real bulkDensity)
{
	Real x, y, z;
	Real fracx, fracy, fracz;
	Real overlap, unitCellVolume;
	int MaxXBin, MaxYBin, MaxZBin;

	Get3DVectorSize(lattice.density, MaxXBin, MaxYBin, MaxZBin, "lattice.density");
	calcFracToCart(Protein.a, Protein.b, Protein.c, Protein.fracToCart);
	unitCellVolume=calcUnitCellVolume(Protein);
	cout <<"In setExcludedVolumeDensity"<<endl;
	for (int i=0;i<MaxXBin;i++)
	{
		for (int j=0;j<MaxYBin;j++)
		{
			for (int k=0;k<MaxZBin;k++)
			{
				indexToPosition(lattice, i, j, k, fracx, fracy, fracz);
				convertToCartesianCoordinates(fracx, fracy, fracz, Protein.fracToCart, x, y, z);
				overlap=isOverlap(Protein.Atoms, x, y, z);
				if (overlap==0) lattice.density[i][j][k]=-bulkDensity*unitCellVolume;
			}
		}
	}
	//string outputDensityFile="/home/joukov/project/mr/pdb_out/1AHO_Density.pdb";
	//printDensityFrac(outputDensityFile, lattice, Protein);
	//endProgram(__LINE__, __FILE__);	
	cout <<"Leaving setExcludedVolumeDensity"<<endl;
}

void initializeLattice(LatticeStruct &lattice, int XCubes, int YCubes, int ZCubes)
{
	Safe3DAlloc(lattice.density, XCubes, YCubes, ZCubes, "lattice.density");
}

void initializeLatticeForExcludedVolume(LatticeStruct &lattice, ProteinStruct &Protein, XRayParamStruct &params)
{
	int maxXBin, maxYBin, maxZBin;
	Real buffer=8.0;
	Real x1, y1, z1;
	Real x2, y2, z2;
	Real fracx1, fracy1, fracz1;
	Real fracx2, fracy2, fracz2;
	Real xmin, ymin, zmin;
	Real xmax, ymax, zmax;

	MinMaxProtein(xmin, ymin, zmin, xmax, ymax, zmax, Protein.Atoms);
	SafeAlloc(lattice.firstPos.Pos, 3, "firstPos");
	x1=xmin-buffer;
	y1=ymin-buffer;
	z1=zmin-buffer;
	x2=xmax+buffer;
	y2=ymax+buffer;
	z2=zmax+buffer;
	convertToFractionalCoordinates(x1, y1, z1, Protein.cartToFrac, fracx1, fracy1, fracz1);
	convertToFractionalCoordinates(x2, y2, z2, Protein.cartToFrac, fracx2, fracy2, fracz2);
	lattice.firstPos.Pos[X]=fracx1;
	lattice.firstPos.Pos[Y]=fracy1;
	lattice.firstPos.Pos[Z]=fracz1;
	lattice.xCubeLength=params.CubeSize/params.a;
	lattice.yCubeLength=params.CubeSize/params.b;
	lattice.zCubeLength=params.CubeSize/params.c;

	maxXBin=int((fracx2-fracx1)/lattice.xCubeLength+0.5);
	maxYBin=int((fracy2-fracy1)/lattice.yCubeLength+0.5);
	maxZBin=int((fracz2-fracz1)/lattice.zCubeLength+0.5);

	cout <<"maxXBin= "<<maxXBin<<" maxYBin= "<<maxYBin<<" maxZBin= "<<maxZBin<<endl;
	initializeLattice(lattice, maxXBin, maxYBin, maxZBin);
	cout <<"Leaving initializeLatticeForExcludedVolume"<<endl;
}
	
void latticeToArray(LatticeStruct &lattice, Real ***&RealDensity, Real ***&ImagDensity)
{
	int MaxXBin, MaxYBin, MaxZBin;

	Get3DVectorSize(lattice.density, MaxXBin, MaxYBin, MaxZBin, "lattice.density");
	Safe3DArrayAlloc(RealDensity, MaxXBin, MaxYBin, MaxZBin, "RealDensity");
	Safe3DArrayAlloc(ImagDensity, MaxXBin, MaxYBin, MaxZBin, "ImagDensity");

	for (int i=0;i<MaxXBin;i++)
	{
		for (int j=0;j<MaxYBin;j++)
		{
			for (int k=0;k<MaxZBin;k++)
			{
				RealDensity[i][j][k]=lattice.density[i][j][k];
				ImagDensity[i][j][k]=0;
			}
		}
	}
}

void fixBin(int maxBin, int &bin)
{
	int mod=0;

	if (bin>=maxBin)
	{
		mod=int(bin/maxBin);
		bin-=mod*maxBin;
	}
	else if (bin<0)
	{
		mod=int(abs(bin+1)/maxBin)+1;
		bin+=mod*maxBin;
	}
}

void fixBins(int maxXBin, int maxYBin, int maxZBin, int &xbin, int &ybin, int &zbin)
{
	fixBin(maxXBin, xbin);
	fixBin(maxYBin, ybin);
	fixBin(maxZBin, zbin);
}

void positionToIndex(Real xmin, Real ymin, Real zmin, Real cubeSize, Real x, Real y, Real z, int &xbin, int &ybin, int &zbin)
{
	xbin=int((x-xmin)/cubeSize+0.5);
	ybin=int((y-ymin)/cubeSize+0.5);
	zbin=int((z-zmin)/cubeSize+0.5);
}

void positionToIndex(Real xmin, Real ymin, Real zmin, Real xCubeLength, Real yCubeLength, Real zCubeLength, Real x, Real y, Real z, int &xbin, int &ybin, int &zbin)
{
	xbin=int((x-xmin)/xCubeLength+0.5);
	ybin=int((y-ymin)/yCubeLength+0.5);
	zbin=int((z-zmin)/zCubeLength+0.5);
}

void positionToIndex(LatticeStruct &lattice, Real x, Real y, Real z, int &xbin, int &ybin, int &zbin)
{
	int maxXBin, maxYBin, maxZBin;

	Get3DVectorSize(lattice.density, maxXBin, maxYBin, maxZBin, "lattice.density");
	positionToIndex(lattice.firstPos.Pos[X], lattice.firstPos.Pos[Y], lattice.firstPos.Pos[Z], lattice.xCubeLength, lattice.yCubeLength, lattice.zCubeLength, x, y, z, xbin, ybin, zbin);
	fixBins(maxXBin, maxYBin, maxZBin, xbin, ybin, zbin);
}

Real fixFractionalCoordinate(Real &x)
{
	if (x>1) x-=Real(int(x));
	if (x<0) x+=Real(int(abs(x)))+1.0;
	return x;
}

void fixFractionalCoordinates(Real &x, Real &y, Real &z)
{
	fixFractionalCoordinate(x);
	fixFractionalCoordinate(y);
	fixFractionalCoordinate(z);
}

void fractionalPositionToIndex(LatticeStruct &lattice, Real x, Real y, Real z, int &xbin, int &ybin, int &zbin)
{
	int maxXBin, maxYBin, maxZBin;

	Get3DVectorSize(lattice.density, maxXBin, maxYBin, maxZBin, "lattice.density");
	fixFractionalCoordinates(x, y, z);
	xbin=int(x*Real(maxXBin));
	ybin=int(y*Real(maxYBin));
	zbin=int(z*Real(maxZBin));

	fixBins(maxXBin, maxYBin, maxZBin, xbin, ybin, zbin);

}

void scaleLatticeDensity(Array3D &density, Real scale)
{
	int maxXBin, maxYBin, maxZBin;

	Get3DVectorSize(density, maxXBin, maxYBin, maxZBin, "density");

	for (int i=0;i<maxXBin;i++)
	{
		for (int j=0;j<maxYBin;j++)
		{
			for (int k=0;k<maxZBin;k++)
			{
				density[i][j][k]*=scale;
			}
		}
	}
}

Real calcLatticeElectrons(LatticeStruct &lattice)
{
	int maxXBin, maxYBin, maxZBin;
	Real electrons=0;

	Get3DVectorSize(lattice.density, maxXBin, maxYBin, maxZBin, "density");

	for (int i=0;i<maxXBin;i++)
	{
		for (int j=0;j<maxYBin;j++)
		{
			for (int k=0;k<maxZBin;k++)
			{
				electrons+=lattice.density[i][j][k];
			}
		}
	}
	electrons*=lattice.xCubeLength*lattice.yCubeLength*lattice.zCubeLength;

	return electrons;
}

void normalizeDensity(ProteinStruct &Protein, LatticeStruct &lattice)
{
	int nsym=Protein.symmetryOperations.size();
	Real mass, latticeMass, scale;

	mass=calcProteinElectrons(Protein.Atoms);
	mass*=Real(nsym);
	latticeMass=calcLatticeElectrons(lattice);
	scale=mass/latticeMass;
	scaleLatticeDensity(lattice.density, scale);
}

Real calcRMSD(LatticeStruct &lattice1, LatticeStruct &lattice2)
{
	int MaxXBin1, MaxYBin1, MaxZBin1;
	int MaxXBin2, MaxYBin2, MaxZBin2;
	Real sum=0, ncube, diff;

	Get3DVectorSize(lattice1.density, MaxXBin1, MaxYBin1, MaxZBin1, "den1");
	Get3DVectorSize(lattice2.density, MaxXBin2, MaxYBin2, MaxZBin2, "den2");

	if (MaxXBin1!=MaxXBin2 || MaxYBin1!=MaxYBin2 || MaxZBin1!=MaxZBin2)
	{
		cout <<"ERROR: The lattice sizes do not match"<<endl;
		exit(EXIT_FAILURE);
	}
	ncube=Real(MaxXBin1)*Real(MaxYBin1)*Real(MaxZBin1);
	for (int i=0;i<MaxXBin1;i++)
	{
		for (int j=0;j<MaxYBin1;j++)
		{
			for (int k=0;k<MaxZBin1;k++)
			{
				diff=lattice1.density[i][j][k]-lattice2.density[i][j][k];
				sum+=diff*diff;
			}
		}
	}
	return sqrt(sum/ncube);
}

Real calcDensityCorrelation(LatticeStruct &lattice1, LatticeStruct &lattice2)
{
	int MaxXBin1, MaxYBin1, MaxZBin1;
	int MaxXBin2, MaxYBin2, MaxZBin2;
	Vector x, y;

	Get3DVectorSize(lattice1.density, MaxXBin1, MaxYBin1, MaxZBin1, "den1");
	Get3DVectorSize(lattice2.density, MaxXBin2, MaxYBin2, MaxZBin2, "den2");

	if (MaxXBin1!=MaxXBin2 || MaxYBin1!=MaxYBin2 || MaxZBin1!=MaxZBin2)
	{
		cout <<"ERROR: The lattice sizes do not match"<<endl;
		exit(EXIT_FAILURE);
	}
	for (int i=0;i<MaxXBin1;i++)
	{
		for (int j=0;j<MaxYBin1;j++)
		{
			for (int k=0;k<MaxZBin1;k++)
			{
				SafePushBack(x, lattice1.density[i][j][k], "x");
				SafePushBack(y, lattice2.density[i][j][k], "y");
			}
		}
	}
	return correlation(x, y);
}

void readPdb(string Pdb, LatticeStruct &lattice)
{
	fstream file;
	bool HetAtom;
	string line, AtomName, duplicate, ResidueName, ChainName, SegID, ID;
	int AtomNumber, ResidueNum;
	int xbin, ybin, zbin;
	int MaxXBin, MaxYBin, MaxZBin;
	Real Occupancy, BFactor;
	Real x, y, z;
	Real xmin, ymin, zmin;
	Real xmax, ymax, zmax;
	cout <<"In readPdb"<<endl;	
	MinMaxPdb(Pdb, xmin, ymin, zmin, xmax, ymax, zmax);
	SafeAlloc(lattice.firstPos.Pos, 3, "pos");
	cout <<"After MinMaxPdb"<<endl;
	lattice.firstPos.Pos[X]=xmin;
	lattice.firstPos.Pos[Y]=ymin;
	lattice.firstPos.Pos[Z]=zmin;
	cout <<"Before positionToIndex"<<endl;
	positionToIndex(lattice, xmax, ymax, zmax, MaxXBin, MaxYBin, MaxZBin);
	cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
	Safe3DAlloc(lattice.density, MaxXBin+1, MaxYBin+1, MaxZBin+1, "lattice.density");
	OpenFile(Pdb, file, "pdb file");
	getline(file, line);

	while (true)
	{
		cout <<line<<endl;
		if (file.eof()) break;
		pdbLineToInfo(line, HetAtom, AtomNumber, AtomName, duplicate, ResidueName, ChainName, ResidueNum, x, y, z, Occupancy, BFactor, SegID, ID);
		positionToIndex(lattice, x, y, z, xbin, ybin, zbin);
		if (HetAtom) lattice.density[xbin][ybin][zbin]=BFactor;
		getline(file, line);
	}
	lattice.real=lattice.density;
}

void printDensity(string outputDensityFile, LatticeStruct &lattice)
{
	int MaxXBin, MaxYBin, MaxZBin;
	ofstream file;
	Real x, y, z;

	Get3DVectorSize(lattice.density, MaxXBin, MaxYBin, MaxZBin, "density");
	OpenFileForWriting(outputDensityFile, file);

	for (int i=0;i<MaxXBin;i++)
	{
		x=lattice.firstPos.Pos[X]+Real(i)*lattice.xCubeLength;
		for (int j=0;j<MaxYBin;j++)
		{
			y=lattice.firstPos.Pos[Y]+Real(j)*lattice.yCubeLength;
			for (int k=0;k<MaxZBin;k++)
			{
				z=lattice.firstPos.Pos[Z]+Real(k)*lattice.zCubeLength;
				printPdbLine(file, true, 1, "O", "HOH", 1, "A", x, y, z, 1.0, lattice.density[i][j][k], "");			
			}
		}
	}
}

void printDensityFrac(string outputDensityFile, LatticeStruct &lattice, ProteinStruct &Protein)
{
	int MaxXBin, MaxYBin, MaxZBin;
	ofstream file;
	Real x, y, z;
	Real x2, y2, z2;
	Real unitCellVolume=calcUnitCellVolume(Protein);

	Get3DVectorSize(lattice.density, MaxXBin, MaxYBin, MaxZBin, "density");
	OpenFileForWriting(outputDensityFile, file);

	for (int i=0;i<MaxXBin;i++)
	{
		x=lattice.firstPos.Pos[X]+Real(i)*lattice.xCubeLength;
		for (int j=0;j<MaxYBin;j++)
		{
			y=lattice.firstPos.Pos[Y]+Real(j)*lattice.yCubeLength;
			for (int k=0;k<MaxZBin;k++)
			{
				z=lattice.firstPos.Pos[Z]+Real(k)*lattice.zCubeLength;
				if (lattice.density[i][j][k]!=0)
				{
					convertToCartesianCoordinates(x, y, z, Protein.fracToCart, x2, y2, z2);
					printPdbLine(file, true, 1, "O", "HOH", 1, "A", x2, y2, z2, 1.0, lattice.density[i][j][k]/unitCellVolume, "");			
				}
			}
		}
	}
}

