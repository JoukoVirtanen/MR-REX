# include <iostream>
# include <vector>
# include <algorithm>

using namespace std;

# include "../../LibraryFiles/Structures.h"
# include "../../LibraryFiles/AtomUtils.h"
# include "../../LibraryFiles/ReadPdb.h"
# include "../../LibraryFiles/PrintPdb.h"
# include "../../LibraryFiles/rmsd.h"
# include "../../LibraryFiles/MathUtils.h"
# include "../../LibraryFiles/TMScore.h"

# include "RMSD.h"

void moveProteinFractional(ProteinStruct &Protein, Real fracX, Real fracY, Real fracZ)
{
	Real dx, dy, dz;

	dx=fracX*Protein.a.x+fracY*Protein.b.x+fracZ*Protein.c.x;
	dy=fracX*Protein.a.y+fracY*Protein.b.y+fracZ*Protein.c.y;
	dz=fracX*Protein.a.z+fracY*Protein.b.z+fracZ*Protein.c.z;
	moveAtoms(Protein.Atoms, dx, dy, dz);
}

void correctForArbitraryAxis(ProteinStruct &Protein1, ProteinStruct &Protein2)
{
	string spaceGroup=Protein1.spaceGroup;
	Real dx, dy, dz;
	Real aveX1, aveY1, aveZ1;
	Real aveX2, aveY2, aveZ2;
	Real fracX, fracY, fracZ;

	if (spaceGroup=="P3" || spaceGroup=="P31" || spaceGroup=="P32" || spaceGroup=="R3" || spaceGroup=="H3" || spaceGroup=="P4" || spaceGroup=="P41" || spaceGroup=="P42" || spaceGroup=="P43" || spaceGroup=="I4" || spaceGroup=="I41" || spaceGroup=="P6" || spaceGroup=="P61" || spaceGroup=="P65" || spaceGroup=="P62" || spaceGroup=="P64" || spaceGroup=="P63")
	{
		//Z-axis is arbitray
		calcCenter(Protein1.Atoms, aveX1, aveY1, aveZ1);
		calcCenter(Protein2.Atoms, aveX2, aveY2, aveZ2);
		dx=aveX2-aveX1;
		dy=aveY2-aveY1;
		dz=aveZ2-aveZ1;
		convertToFractionalCoordinates(dx, dy, dz, Protein1.cartToFrac, fracX, fracY, fracZ);

		dx=fracZ*Protein1.c.x;
		dy=fracZ*Protein1.c.y;
		dz=fracZ*Protein1.c.z;
		moveAtoms(Protein1.Atoms, dx, dy, dz);
	}
	else if (spaceGroup=="P2" || spaceGroup=="P21" || spaceGroup=="C2" || spaceGroup=="A2" || spaceGroup=="I2" || spaceGroup=="P1211" || spaceGroup=="C121" || spaceGroup=="I121" || spaceGroup=="P121")
	{
		//Y-axis is arbitrary
		calcCenter(Protein1.Atoms, aveX1, aveY1, aveZ1);
		calcCenter(Protein2.Atoms, aveX2, aveY2, aveZ2);
		dx=aveX2-aveX1;
		dy=aveY2-aveY1;
		dz=aveZ2-aveZ1;
		convertToFractionalCoordinates(dx, dy, dz, Protein1.cartToFrac, fracX, fracY, fracZ);
		//cout <<"fracX= "<<fracX<<endl;
		//cout <<"fracY= "<<fracY<<endl;
		//cout <<"fracZ= "<<fracZ<<endl;

		dx=fracY*Protein1.b.x;
		dy=fracY*Protein1.b.y;
		dz=fracY*Protein1.b.z;

		moveAtoms(Protein1.Atoms, dx, dy, dz);
/*
		calcCenter(Protein1.Atoms, aveX1, aveY1, aveZ1);
		calcCenter(Protein2.Atoms, aveX2, aveY2, aveZ2);
		dx=aveX2-aveX1;
		dy=aveY2-aveY1;
		dz=aveZ2-aveZ1;
		convertToFractionalCoordinates(dx, dy, dz, Protein1.cartToFrac, fracX, fracY, fracZ);
		cout <<"fracX= "<<fracX<<endl;
		cout <<"fracY= "<<fracY<<endl;
		cout <<"fracZ= "<<fracZ<<endl;
		endProgram(__LINE__, __FILE__);
*/
	}
	else if (spaceGroup=="P1" || spaceGroup=="A1")
	{
		calcCenter(Protein1.Atoms, aveX1, aveY1, aveZ1);
		calcCenter(Protein2.Atoms, aveX2, aveY2, aveZ2);
		dx=aveX2-aveX1;
		dy=aveY2-aveY1;
		dz=aveZ2-aveZ1;
		moveAtoms(Protein1.Atoms, dx, dy, dz);
	}

}

Real calcRMSDSymmetryOperation(ProteinStruct Protein1, ProteinStruct &Protein2, int index)
{
	string spaceGroup=Protein1.spaceGroup;
	int i, j, k;
	Real rmsd;
	Real aveX1, aveY1, aveZ1;
	Real aveX2, aveY2, aveZ2;
	Real dx, dy, dz;
	Real fracX, fracY, fracZ;

	applySymmetryOperation(Protein1.Atoms, Protein1.symmetryOperations[index]);
	calcCenter(Protein1.Atoms, aveX1, aveY1, aveZ1);
	calcCenter(Protein2.Atoms, aveX2, aveY2, aveZ2);

	dx=aveX2-aveX1;
	dy=aveY2-aveY1;
	dz=aveZ2-aveZ1;
	convertToFractionalCoordinates(dx, dy, dz, Protein1.cartToFrac, fracX, fracY, fracZ);
	

	i=int(abs(fracX)+0.5);
	j=int(abs(fracY)+0.5);
	k=int(abs(fracZ)+0.5);
	if (fracX<0) i=-i;
	if (fracY<0) j=-j;
	if (fracZ<0) k=-k;
	dx=Real(i)*Protein1.a.x+Real(j)*Protein1.b.x+Real(k)*Protein1.c.x;
	dy=Real(i)*Protein1.a.y+Real(j)*Protein1.b.y+Real(k)*Protein1.c.y;
	dz=Real(i)*Protein1.a.z+Real(j)*Protein1.b.z+Real(k)*Protein1.c.z;
	moveAtoms(Protein1.Atoms, dx, dy, dz);
	correctForArbitraryAxis(Protein1, Protein2);
	rmsd=RMSDNoAlign(Protein1.Atoms, Protein2.Atoms);
	//cout <<"fracX= "<<fracX<<" fracY= "<<fracY<<" fracZ= "<<fracZ<<" rmsd= "<<rmsd<<endl;
	return rmsd;
}

Real calcRMSDReference(ProteinStruct &Protein1, ProteinStruct &Protein2, Real x, Real y, Real z, int &index)
{
	string spaceGroup=Protein1.spaceGroup;
	int i, j, k;
	Real rmsd;
	Real aveX1, aveY1, aveZ1;
	Real aveX2, aveY2, aveZ2;
	Real dx, dy, dz;
	Real fracX, fracY, fracZ;

	dx=Real(x)*Protein1.a.x+Real(y)*Protein1.b.x+Real(z)*Protein1.c.x;
	dy=Real(x)*Protein1.a.y+Real(y)*Protein1.b.y+Real(z)*Protein1.c.y;
	dz=Real(x)*Protein1.a.z+Real(y)*Protein1.b.z+Real(z)*Protein1.c.z;

	moveAtoms(Protein1.Atoms, dx, dy, dz);

	applySymmetryOperation(Protein1.Atoms, Protein1.symmetryOperations[index]);
	calcCenter(Protein1.Atoms, aveX1, aveY1, aveZ1);
	calcCenter(Protein2.Atoms, aveX2, aveY2, aveZ2);

	dx=aveX2-aveX1;
	dy=aveY2-aveY1;
	dz=aveZ2-aveZ1;

	convertToFractionalCoordinates(dx, dy, dz, Protein1.cartToFrac, fracX, fracY, fracZ);

	i=int(abs(fracX)+0.5);
	j=int(abs(fracY)+0.5);
	k=int(abs(fracZ)+0.5);
	if (fracX<0) i=-i;
	if (fracY<0) j=-j;
	if (fracZ<0) k=-k;
	dx=Real(i)*Protein1.a.x+Real(j)*Protein1.b.x+Real(k)*Protein1.c.x;
	dy=Real(i)*Protein1.a.y+Real(j)*Protein1.b.y+Real(k)*Protein1.c.y;
	dz=Real(i)*Protein1.a.z+Real(j)*Protein1.b.z+Real(k)*Protein1.c.z;
	moveAtoms(Protein1.Atoms, dx, dy, dz);
	correctForArbitraryAxis(Protein1, Protein2);
	rmsd=RMSDNoAlign(Protein1.Atoms, Protein2.Atoms);
	//cout <<"rmsd= "<<rmsd<<endl;
	//if (rmsd<1.0) cout <<"index= "<<index<<endl;
	return rmsd;
}

Real calcRMSDSymmetryOperations(ProteinStruct Protein1, ProteinStruct &Protein2, int &bestIndex)
{
	int noperations=Protein1.symmetryOperations.size();
	Real bestRMSD, rmsd;
	bestRMSD=calcRMSDSymmetryOperation(Protein1, Protein2, 0);
	bestIndex=0;
	for (int i=1;i<noperations;i++)
	{
		rmsd=calcRMSDSymmetryOperation(Protein1, Protein2, i);
		if (rmsd<bestRMSD) 
		{
			bestRMSD=rmsd;
			bestIndex=i;
		}
	}
	return bestRMSD;
}

Real calcRMSD(ProteinStruct Protein1, ProteinStruct &Protein2, Real fracX, Real fracY, Real fracZ, int &bestIndex)
{
	Real rmsd;
	moveProteinFractional(Protein1, fracX, fracY, fracZ);
	rmsd=calcRMSDSymmetryOperations(Protein1, Protein2, bestIndex);
	return rmsd;
}

void getAlternateOrigins(string spaceGroup, Matrix &origin)
{
	//Source: http://www.ccp4.ac.uk/html/alternate_origins.html
	Vector tempOrigin;

	SafeAlloc(tempOrigin, 3, "tempOrigin");

	if (spaceGroup=="H32" || spaceGroup=="R32")
	{
		//I don't think H32 belongs here.
		Safe2DAlloc(origin, 6, 3, "origin");

		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=0;
		origin[1][1]=0;
		origin[1][2]=0.5;

		origin[2][0]=1.0/3.0;
		origin[2][1]=2.0/3.0;
		origin[2][2]=1.0/6.0;

		origin[3][0]=1.0/3.0;
		origin[3][1]=2.0/3.0;
		origin[3][2]=2.0/3.0;

		origin[4][0]=2.0/3.0;
		origin[4][1]=1.0/3.0;
		origin[4][2]=1.0/3.0;

		origin[5][0]=2.0/3.0;
		origin[5][1]=1.0/3.0;
		origin[5][2]=5.0/6.0;
	}
	else if (spaceGroup=="P1" || spaceGroup=="A1")
	{
		Safe2DAlloc(origin, 1, 3, "origin");

		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;
	}
	else if (spaceGroup=="P321" || spaceGroup=="P3121" || spaceGroup=="P3221" || spaceGroup=="P622" || spaceGroup=="P6122" || spaceGroup=="P6522" || spaceGroup=="P6222" || spaceGroup=="P6422" || spaceGroup=="P6322")
	{
		Safe2DAlloc(origin, 2, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=0;
		origin[1][1]=0;
		origin[1][2]=0.5;
	}
	else if (spaceGroup=="P422" || spaceGroup=="P4212" || spaceGroup=="P4122" || spaceGroup=="P41212" || spaceGroup=="P4222" || spaceGroup=="P42212" || spaceGroup=="P4322" || spaceGroup=="P43212" || spaceGroup=="I4122")
	{
		Safe2DAlloc(origin, 4, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=0;
		origin[1][1]=0;
		origin[1][2]=0.5;

		origin[2][0]=0.5;
		origin[2][1]=0.5;
		origin[2][2]=0;

		origin[3][0]=0.5;
		origin[3][1]=0.5;
		origin[3][2]=0.5;
	}
	else if (spaceGroup=="P2" || spaceGroup=="P21" || spaceGroup=="C2" || spaceGroup=="A2" || spaceGroup=="I2" || spaceGroup=="P1211" || spaceGroup=="C121" || spaceGroup=="I121" || spaceGroup=="P121")
	{
		//Check P1211 and C121
		Safe2DAlloc(origin, 4, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=0;
		origin[1][1]=0;
		origin[1][2]=0.5;
		
		origin[2][0]=0.5;
		origin[2][1]=0;
		origin[2][2]=0;

		origin[3][0]=0.5;
		origin[3][1]=0;
		origin[3][2]=0.5;
	}
	else if (spaceGroup=="P1211" || spaceGroup=="P222" || spaceGroup=="P2122" || spaceGroup=="P2212" || spaceGroup=="P2221" || spaceGroup=="P22121" || spaceGroup=="P21221" || spaceGroup=="P21212" || spaceGroup=="P212121" || spaceGroup=="C2221" || spaceGroup=="C222" || spaceGroup=="I222" || spaceGroup=="I212121" || spaceGroup=="F432" || spaceGroup=="F4132" || spaceGroup=="P-1" || spaceGroup=="P121/C1")
	{
		Safe2DAlloc(origin, 8, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=0;
		origin[1][1]=0;
		origin[1][2]=0.5;
		
		origin[2][0]=0;
		origin[2][1]=0.5;
		origin[2][2]=0;

		origin[3][0]=0;
		origin[3][1]=0.5;
		origin[3][2]=0.5;

		origin[4][0]=0.5;
		origin[4][1]=0;
		origin[4][2]=0;

		origin[5][0]=0.5;
		origin[5][1]=0;
		origin[5][2]=0.5;

		origin[6][0]=0.5;
		origin[6][1]=0.5;
		origin[6][2]=0;

		origin[7][0]=0.5;
		origin[7][1]=0.5;
		origin[7][2]=0.5;

	}	
	else if (spaceGroup=="F222" || spaceGroup=="F23")
	{
		Safe2DAlloc(origin, 16, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=0;
		origin[1][1]=0;
		origin[1][2]=0.5;
		
		origin[2][0]=0;
		origin[2][1]=0.5;
		origin[2][2]=0;

		origin[3][0]=0;
		origin[3][1]=0.5;
		origin[3][2]=0.5;

		origin[4][0]=0.25;
		origin[4][1]=0.25;
		origin[4][2]=0.25;

		origin[5][0]=0.25;
		origin[5][1]=0.25;
		origin[5][2]=0.75;

		origin[6][0]=0.25;
		origin[6][1]=0.75;
		origin[6][2]=0.25;

		origin[7][0]=0.25;
		origin[7][1]=0.75;
		origin[7][2]=0.75;

		origin[8][0]=0.5;
		origin[8][1]=0;
		origin[8][2]=0;

		origin[9][0]=0.5;
		origin[9][1]=0;
		origin[9][2]=0.5;

		origin[10][0]=0.5;
		origin[10][1]=0.5;
		origin[10][2]=0;

		origin[11][0]=0.5;
		origin[11][1]=0.5;
		origin[11][2]=0.5;

		origin[12][0]=0.75;
		origin[12][1]=0.25;
		origin[12][2]=0.25;

		origin[13][0]=0.75;
		origin[13][1]=0.25;
		origin[13][2]=0.75;

		origin[14][0]=0.75;
		origin[14][1]=0.75;
		origin[14][2]=0.25;

		origin[15][0]=0.75;
		origin[15][1]=0.75;
		origin[15][2]=0.75;

	}
	else if (spaceGroup=="P4" || spaceGroup=="P41" || spaceGroup=="P42" || spaceGroup=="P43" || spaceGroup=="I4" || spaceGroup=="I41")
	{
		//Z-axis is also arbitrary
		Safe2DAlloc(origin, 2, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=0.5;
		origin[1][1]=0.5;
		origin[1][2]=0;
	}
	else if (spaceGroup=="I422")
	{
		//Z-axis is also arbitrary
		Safe2DAlloc(origin, 2, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=0.5;
		origin[1][1]=0.5;
		origin[1][2]=0;
	}
	else if (spaceGroup=="P3" || spaceGroup=="P31" || spaceGroup=="P32" || spaceGroup=="R3" || spaceGroup=="H3")
	{
		//Z-axis is also arbitray
		Safe2DAlloc(origin, 3, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=1.0/3.0;
		origin[1][1]=2.0/3.0;
		origin[1][2]=0;

		origin[2][0]=2.0/3.0;
		origin[2][1]=1.0/3.0;
		origin[2][2]=0;

	}
	else if (spaceGroup=="P312" || spaceGroup=="P3112" || spaceGroup=="P3212")
	{
		Safe2DAlloc(origin, 6, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=0;
		origin[1][1]=0;
		origin[1][2]=0.5;

		origin[2][0]=1.0/3.0;
		origin[2][1]=2.0/3.0;
		origin[2][2]=0;

		origin[3][0]=1.0/3.0;
		origin[3][1]=2.0/3.0;
		origin[3][2]=0.5;

		origin[4][0]=2.0/3.0;
		origin[4][1]=1.0/3.0;
		origin[4][2]=0;

		origin[5][0]=2.0/3.0;
		origin[5][1]=1.0/3.0;
		origin[5][2]=0.5;

	}
	else if (spaceGroup=="P6" || spaceGroup=="P61" || spaceGroup=="P65" || spaceGroup=="P62" || spaceGroup=="P64" || spaceGroup=="P63")
	{
		//Z-axis is arbitrary
		Safe2DAlloc(origin, 1, 3, "origin");
	}
	else if (spaceGroup=="P23" || spaceGroup=="P213" || spaceGroup=="I23" || spaceGroup=="I213" || spaceGroup=="P432" || spaceGroup=="P4232" || spaceGroup=="P4332" || spaceGroup=="P4132" || spaceGroup=="I432" || spaceGroup=="I4132")
	{
		Safe2DAlloc(origin, 2, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;

		origin[1][0]=0.5;
		origin[1][1]=0.5;
		origin[1][2]=0.5;
		
	}
	else if (spaceGroup=="C121")
	{
		//No source for this.  Check
		Safe2DAlloc(origin, 1, 3, "origin");
		
		origin[0][0]=0;
		origin[0][1]=0;
		origin[0][2]=0;
	}
	else
	{
		error("Unrecognized space group= "+spaceGroup, __LINE__, __FILE__);
		for (Real i=-1.0;i<=1.1;i+=0.5)
		{
			for (Real j=-1.0;j<=1.1;j+=0.5)
			{
				for (Real k=-1.0;k<=1.1;k+=0.5)
				{
					tempOrigin[0]=i;
					tempOrigin[1]=j;
					tempOrigin[2]=k;
					SafePushBack(origin, tempOrigin, "origin");
				}
			}
		}
	}
}

Real calcRMSD(ProteinStruct &Protein1, ProteinStruct &Protein2, string superimposeType)
{
	int bestIndex, veryBestIndex, norigin;
	Real bestRMSD, rmsd;
	Real bestDx, bestDy, bestDz;
	Matrix origin;
	ProteinStruct tempProtein;

	tempProtein=Protein1;
	//superimposeType="RMSD";
	if (superimposeType=="RMSD") 
	{
		RMSD(Protein2.Atoms, tempProtein.Atoms);
	}
	else if (superimposeType=="TMSCORE")
	{
		for (int i=0;i<10;i++)
		{
			printAtomInfo(tempProtein.Atoms[i]);
		}
		calcTMScore(Protein2.Atoms, tempProtein.Atoms);
	}
	else if (superimposeType=="NONE")
	{
		tempProtein.Atoms=Protein2.Atoms;
	}
	else
	{
		string errorStr="Unknown superimposeType: "+superimposeType;
		errorStr+=" acceptable values are RMSD, TMSCORE, and NONE";
		error(errorStr, __LINE__, __FILE__);
	}
	bestRMSD=calcRMSD(Protein1, tempProtein, 0, 0, 0, bestIndex);
	veryBestIndex=bestIndex;
	bestDx=0;
	bestDy=0;
	bestDz=0;
	getAlternateOrigins(Protein1.spaceGroup, origin);
	norigin=origin.size();

	for (int i=0;i<norigin;i++)
	{
		rmsd=calcRMSD(Protein1, tempProtein, origin[i][0], origin[i][1], origin[i][2], bestIndex);
		if (rmsd<bestRMSD) 
		{
			bestRMSD=rmsd;
			veryBestIndex=bestIndex;
			bestDx=origin[i][0];
			bestDy=origin[i][1];
			bestDz=origin[i][2];
		}
	}
	calcRMSDReference(Protein1, tempProtein, bestDx, bestDy, bestDz, veryBestIndex);
	return bestRMSD;
}

Real calcRMSD_NoAlign(ProteinStruct &Protein1, ProteinStruct &Protein2)
{
	//int veryBestIndex;
	int bestIndex, norigin;
	Real bestRMSD, rmsd;
	//Real bestDx, bestDy, bestDz;
	Matrix origin;

	bestRMSD=calcRMSD(Protein1, Protein2, 0, 0, 0, bestIndex);
	//veryBestIndex=bestIndex;
	//bestDx=0;
	//bestDy=0;
	//bestDz=0;
	getAlternateOrigins(Protein1.spaceGroup, origin);
	norigin=origin.size();

	for (int i=1;i<norigin;i++)
	{
		rmsd=calcRMSD(Protein1, Protein2, origin[i][0], origin[i][1], origin[i][2], bestIndex);
		if (rmsd<bestRMSD) 
		{
			bestRMSD=rmsd;
			//veryBestIndex=bestIndex;
			//bestDx=origin[i][0];
			//bestDy=origin[i][1];
			//bestDz=origin[i][2];
		}
	}
	return bestRMSD;
}

Real calcRMSD(vector<ProteinStruct> &Proteins1, vector<ProteinStruct> &Proteins2, vector<int> &permutation, string superimposeType)
{
	int Size;
	Real rmsd, rmsd2;


	Size=permutation.size();
	//cout <<"Size= "<<Size<<endl;
	rmsd2=0;
	for (int i=0;i<Size;i++)
	{
		rmsd=calcRMSD(Proteins1[i], Proteins2[permutation[i]], superimposeType);
		rmsd2+=rmsd*rmsd;
	}
	rmsd/=Real(Size);
	rmsd=sqrt(rmsd2);

	return rmsd;
}

Real calcPointXtalRmsd(ProteinStruct &Protein1, ProteinStruct &Protein2)
{
	int ca;
	Real avex1, avey1, avez1;
	Real avex2, avey2, avez2;
	ProteinStruct tempProtein1, tempProtein2;

	tempProtein1=Protein1;
	tempProtein2=Protein2;
	calcCenter(Protein1.Atoms, avex1, avey1, avez1);
	calcCenter(Protein2.Atoms, avex2, avey2, avez2);
	tempProtein1.Atoms.clear();
	tempProtein2.Atoms.clear();
	ca=getAtomIndex(Protein1.Atoms, 1, "CA");
	SafePushBack(tempProtein1.Atoms, Protein1.Atoms[ca], "tempProtein1.Atoms");
	SafePushBack(tempProtein2.Atoms, Protein2.Atoms[ca], "tempProtein2.Atoms");

	tempProtein1.Atoms[0].x=avex1;
	tempProtein1.Atoms[0].y=avey1;
	tempProtein1.Atoms[0].z=avez1;

	tempProtein2.Atoms[0].x=avex2;
	tempProtein2.Atoms[0].y=avey2;
	tempProtein2.Atoms[0].z=avez2;
	return calcRMSD_NoAlign(tempProtein1, tempProtein2);
}

Real calcRMSD(Matrix &rmsdPair, vector<int> &permutation)
{
	int nprot=permutation.size();
	int Size1, Size2;
	Real rmsd, rmsd2=0;

	Get2DVectorSize(rmsdPair, Size1, Size2);

	if (nprot!=Size1 || nprot!=Size2 || nprot==0)
	{
		string errorStr="nprot= "+toStr(nprot)+" Size1= "+toStr(Size1);
		errorStr+=" Size2= "+toStr(Size2);
		error(errorStr, __LINE__, __FILE__);
	}

	for (int i=0;i<nprot;i++)
	{
		rmsd2+=rmsdPair[i][permutation[i]]*rmsdPair[i][permutation[i]];
	}
	rmsd2/=Real(nprot);
	rmsd=sqrt(rmsd2);

	return rmsd;
}

void getMatchingAtomsNoChain(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2, vector<AtomStruct> &AtomsOut)
{
        int natom1=Atoms1.size();
        int natom2=Atoms2.size();


        if (natom1==0 || natom2==0)
        {
                string errorStr="natom1= "+toStr(natom1);
                errorStr+=" natom2= "+toStr(natom2);
        }
        AtomsOut.clear();
	atomAlias(Atoms1);
	atomAlias(Atoms2);
        for (int i=0;i<natom1;i++)
        {
                for (int j=0;j<natom2;j++)
                {
                        //if (Atoms1[i].AtomName==Atoms2[j].AtomName && Atoms1[i].ResidueNum-firstResidue1==Atoms2[j].ResidueNum-firstResidue2)
                        if (Atoms1[i].AtomName==Atoms2[j].AtomName && Atoms1[i].ResidueNum==Atoms2[j].ResidueNum && Atoms1[i].ResidueName==Atoms2[j].ResidueName)
                        {
                                SafePushBack(AtomsOut, Atoms1[i], "AtomsOut");
                        }
                }
        }
}


Real calcRMSD(vector<ProteinStruct> &Proteins1, vector<ProteinStruct> &Proteins2, string superimposeType)
{
	int nprot=Proteins1.size();
	vector<int> permutation;
	Real bestRmsd, rmsd, min;
	Matrix rmsdPair;
	ProteinStruct tempProtein1, tempProtein2;

	SafeAlloc(permutation, nprot, "permutation");
	Safe2DAlloc(rmsdPair, nprot, nprot, "rmsdPair");
	for (int i=0;i<nprot;i++)
	{
		permutation[i]=i;
	}

	for (int i=0;i<nprot;i++)
	{
		for (int j=0;j<nprot;j++)
		{
			tempProtein1=Proteins1[i];
			tempProtein2=Proteins2[j];
			getMatchingAtomsNoChain(Proteins1[i].Atoms, Proteins2[j].Atoms, tempProtein1.Atoms);	
			getMatchingAtomsNoChain(Proteins2[j].Atoms, Proteins1[i].Atoms, tempProtein2.Atoms);	
			//rmsdPair[i][j]=RMSD_CA(tempProtein1.Atoms, tempProtein2.Atoms);
			rmsdPair[i][j]=calcRMSD(tempProtein1, tempProtein2, superimposeType);
		}
	}
	min=rmsdPair[0][0];
	for (int i=0;i<nprot;i++)
	{
		for (int j=0;j<nprot;j++)
		{
			if (min>rmsdPair[i][j]) min=rmsdPair[i][j];
		}
	}
	printMatrix(rmsdPair);
	bestRmsd=calcRMSD(rmsdPair, permutation);
	do
	{
		rmsd=calcRMSD(rmsdPair, permutation);
		if (rmsd<bestRmsd) bestRmsd=rmsd;
	} while (next_permutation(permutation.begin(), permutation.end()));
	cout <<"finalRMSD= "<<bestRmsd<<endl;
	return bestRmsd;
}
