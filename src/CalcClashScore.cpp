# include <iostream>
# include <vector>
# include <sys/time.h>
# include <algorithm>

using namespace std;


# include "../../LibraryFiles/TypeDef.h"
# include "../../LibraryFiles/Structures.h"
# include "../../LibraryFiles/normalize.h"
# include "../../LibraryFiles/AtomUtils.h"
# include "../../LibraryFiles/MinMax.h"
# include "../../LibraryFiles/MonteCarlo.h"
# include "../../LibraryFiles/MonteCarloTemplate.h"
# include "../../LibraryFiles/VectorManip.h"
# include "../../LibraryFiles/time.h"
# include "../../LibraryFiles/MemoryUsage.h"
# include "../../LibraryFiles/PrintPdb.h"
# include "../../LibraryFiles/IOUtils.h"
# include "../../LibraryFiles/StringUtils.h"
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

void calcClashScoreRecord(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2, Real &clashCC, Real &clashCS, Real &clashSS)
{
	int natom1=Atoms1.size();
	int natom2=Atoms2.size();
	Real cutOff=3.0;
	Real cutOff2=cutOff*cutOff;
	Real clash, dist2;
	Real dx, dy, dz;

	clashCC=0;
	clashCS=0;
	clashSS=0;
	for (int i=0;i<natom1;i++)
	{
		if (Atoms1[i].Occupancy!=0)
		{	
			for (int j=0;j<natom2;j++)
			{
				if (Atoms2[j].Occupancy!=0)
				{
					dx=Atoms1[i].x-Atoms2[j].x;
					dy=Atoms1[i].y-Atoms2[j].y;
					dz=Atoms1[i].z-Atoms2[j].z;
					dist2=dx*dx+dy*dy+dz*dz;
					if (dist2<cutOff2)
					{
						clash=cutOff2-dist2;
						if (Atoms1[i].core && Atoms2[j].core) clashCC+=clash;
						else if (!Atoms1[i].core && !Atoms2[j].core) clashSS+=clash;
						else clashCS+=clash;
						SafePushBack(Atoms1[i].clash, j, "clash");
						SafePushBack(Atoms1[j].clash, i, "clash");
					}
				}
			}
		}
	}
}

void calcClashScoreRecord(ProteinStruct &Protein, vector<ProteinStruct> &allProteins, int index, Real cutoff, Real &clashCC, Real &clashCS, Real &clashSS)
{
	int nprot=allProteins.size();
	Real r2;
	//Real buffer=3.0;
	Real dx, dy, dz;
	Real copyCC=0, copyCS=0, copySS=0;
	for (int i=index+1;i<nprot;i++)
	{
		dx=Protein.center[X]-allProteins[i].center[X];
		dy=Protein.center[Y]-allProteins[i].center[Y];
		dz=Protein.center[Z]-allProteins[i].center[Z];
		r2=dx*dx+dy*dy+dz*dz;
		if (r2<cutoff)
		{
			calcClashScoreRecord(Protein.Atoms, allProteins[i].Atoms, copyCC, copyCS, copySS);
			clashCC+=copyCC;
			clashCS+=copyCS;
			clashSS+=copySS;
		}
	}
}

Real calcClashScoreRecord(ProteinStruct &Protein, XRayParamStruct &params, Real &clashCC, Real &clashCS, Real &clashSS)
{
	Real r, cutoff, buffer=3.0;
	Real clashScore=0;
	Real proteinCC=0, proteinCS=0, proteinSS=0;
	vector<ProteinStruct> allProteins, Proteins;
	clashCC=0, clashCS=0, clashSS=0;

	//removeNonCA(Protein.Atoms);
	SafePushBack(Proteins, Protein, "Proteins");
	r=calcFurthestDistance(Protein.Atoms)+buffer;
	cutoff=4.0*r*r;
	makeUnitCellAndImagesFast3(Proteins, allProteins);
	calcClashScoreRecord(Proteins[0], allProteins, 0, cutoff, proteinCC, proteinCS, proteinSS);
	Protein=Proteins[0];
	clashCC+=proteinCC;
	clashCS+=proteinCS;
	clashSS+=proteinSS;
	clashScore=params.ClashWeightCoreCore*clashCC;
	clashScore+=params.ClashWeightCoreSurface*clashCS;
	clashScore+=params.ClashWeightSurfaceSurface*clashSS;
	return clashScore;
}

void calcClashScore(AtomStruct &Atom, vector<AtomStruct> &Atoms, Real &clashCC, Real &clashCS, Real &clashSS)
{
	int natom=Atoms.size();
	Real cutOff=3.0;
	Real cutOff2=cutOff*cutOff;
	Real clash, dist2;
	Real dx, dy, dz;

	clashCC=0;
	clashCS=0;
	clashSS=0;

	for (int i=0;i<natom;i++)
	{
		dx=Atom.x-Atoms[i].x;
		if (abs(dx)<cutOff)
		{
			dy=Atom.y-Atoms[i].y;
			if (abs(dy)<cutOff)
			{
				dz=Atom.z-Atoms[i].z;
				if (abs(dz)<cutOff)
				{
					dist2=dx*dx+dy*dy+dz*dz;
					if (dist2<cutOff2 && Atoms[i].Occupancy!=0 && Atom.Occupancy!=0)
					{
						clash=cutOff2-dist2;
						if (Atom.core && Atoms[i].core) clashCC+=clash;
						else if (!Atom.core && !Atoms[i].core) clashSS+=clash;
						else 
						{
							clashCS+=clash;
						}
					}
				}
			}
		}
	}
}

void calcClashScore(AtomStruct &Atom, int index, vector<AtomStruct> &Atoms, Vector &clashCC, Vector &clashCS, Vector &clashSS)
{
	int natom=Atoms.size();
	Real cutOff=3.0;
	Real cutOff2=cutOff*cutOff;
	Real clash, dist2;
	Real dx, dy, dz;

	for (int i=0;i<natom;i++)
	{
		dx=Atom.x-Atoms[i].x;
		if (abs(dx)<cutOff && Atoms[i].Occupancy!=0 && Atom.Occupancy!=0)
		//if (abs(dx)<cutOff)
		{
			dy=Atom.y-Atoms[i].y;
			if (abs(dy)<cutOff)
			{
				dz=Atom.z-Atoms[i].z;
				if (abs(dz)<cutOff)
				{
					dist2=dx*dx+dy*dy+dz*dz;
					if (dist2<cutOff2)
					{
						clash=cutOff2-dist2;
						if (Atom.core && Atoms[i].core) 
						{
							clashCC[i]+=clash;
							clashCC[index]+=clash;
						}
						else if (!Atom.core && !Atoms[i].core) 
						{
							clashSS[i]+=clash;
							clashSS[index]+=clash;
						}
						else 
						{
							clashCS[i]+=clash;
							clashCS[index]+=clash;
						}
					}
				}
			}
		}
	}
}

void calcClashScore(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2, Real &clashCC, Real &clashCS, Real &clashSS)
{
	int natom=Atoms1.size();
	Real atomCC=0, atomCS=0, atomSS=0;

	clashCC=0;
	clashCS=0;
	clashSS=0;
	for (int i=0;i<natom;i++)
	{
		calcClashScore(Atoms1[i], Atoms2, atomCC, atomCS, atomSS);
		clashCC+=atomCC;
		clashCS+=atomCS;
		clashSS+=atomSS;
	}
}

void calcClashScore(vector<AtomStruct> &Atoms1, ProteinStruct &Protein, Real atomCutoff, Real &clashCC, Real &clashCS, Real &clashSS)
{
	int natom=Atoms1.size();
	Real atomCC=0, atomCS=0, atomSS=0;
	Real dx, dy, dz, r2;

	clashCC=0;
	clashCS=0;
	clashSS=0;
	for (int i=0;i<natom;i++)
	{
		dx=Atoms1[i].x-Protein.center[X];
		dy=Atoms1[i].y-Protein.center[Y];
		dz=Atoms1[i].z-Protein.center[Z];
		r2=dx*dx+dy*dy+dz*dz;
		//if (r2<atomCutoff && Atoms1[i].Occupancy!=0)
		if (r2<atomCutoff)
		{
			calcClashScore(Atoms1[i], Protein.Atoms, atomCC, atomCS, atomSS);
			clashCC+=atomCC;
			clashCS+=atomCS;
			clashSS+=atomSS;
		}
	}
}

void calcClashScoreSlow(ProteinStruct &Protein, vector<ProteinStruct> &allProteins, int index, Real cutoff, Real &clashCC, Real &clashCS, Real &clashSS)
{
	int nprot=allProteins.size();
	Real r2;
	//Real buffer=3.0;
	Real dx, dy, dz;
	Real copyCC=0, copyCS=0, copySS=0;
	for (int i=index+1;i<nprot;i++)
	{
		dx=Protein.center[X]-allProteins[i].center[X];
		dy=Protein.center[Y]-allProteins[i].center[Y];
		dz=Protein.center[Z]-allProteins[i].center[Z];
		r2=dx*dx+dy*dy+dz*dz;
		if (r2<cutoff)
		{
			calcClashScore(Protein.Atoms, allProteins[i].Atoms, copyCC, copyCS, copySS);
			clashCC+=copyCC;
			clashCS+=copyCS;
			clashSS+=copySS;
		}
	}
}

void calcClashScore(ProteinStruct &Protein, vector<ProteinStruct> &allProteins, int index, Real cutoff, Real &clashCC, Real &clashCS, Real &clashSS)
{
	int nprot=allProteins.size();
	Real r2;
	Real buffer=3.0;
	Real furthestDistance, atomCutoff;
	Real dx, dy, dz;
	Real copyCC=0, copyCS=0, copySS=0;
	furthestDistance=calcFurthestDistance(Protein.Atoms)+2.0*buffer;
	atomCutoff=furthestDistance*furthestDistance;
	for (int i=index+1;i<nprot;i++)
	{
		dx=Protein.center[X]-allProteins[i].center[X];
		dy=Protein.center[Y]-allProteins[i].center[Y];
		dz=Protein.center[Z]-allProteins[i].center[Z];
		r2=dx*dx+dy*dy+dz*dz;
		if (r2<cutoff)
		{
			calcClashScore(Protein.Atoms, allProteins[i], atomCutoff, copyCC, copyCS, copySS);
			//calcClashScore(Protein.Atoms, allProteins[i].Atoms, copyCCSlow, copyCSSlow, copySSSlow);
			//cout <<"copyCC= "<<copyCC<<" copyCCSlow= "<<copyCCSlow<<endl;
			//cout <<"copyCS= "<<copyCS<<" copyCSSlow= "<<copyCSSlow<<endl;
			//cout <<"copySS= "<<copySS<<" copySSSlow= "<<copySSSlow<<endl;
			clashCC+=copyCC;
			clashCS+=copyCS;
			clashSS+=copySS;
			//cout <<"copyCC= "<<copyCC<<" copyCS= "<<copyCS<<" copySS= "<<copySS<<endl;
			//cout <<"clashCC= "<<clashCC<<" clashCS= "<<clashCS<<" clashSS= "<<clashSS<<endl;
		}
	}
}

void calcClashScoreByAtom(AtomStruct &Atom, int atomIndex, Vector &center, vector<ProteinStruct> &allProteins, int index, Real cutoff, Vector &clashCC, Vector &clashCS, Vector &clashSS)
{
	int nprot=allProteins.size();
	Real r2;
	//Real buffer=3.0;
	Real dx, dy, dz;

	for (int i=index+1;i<nprot;i++)
	{
		dx=center[X]-allProteins[i].center[X];
		dy=center[Y]-allProteins[i].center[Y];
		dz=center[Z]-allProteins[i].center[Z];
		r2=dx*dx+dy*dy+dz*dz;
		if (r2<cutoff)
		{
			calcClashScore(Atom, atomIndex, allProteins[i].Atoms, clashCC, clashCS, clashSS);
		}
	}
}

void calcClashScoreByAtom(ProteinStruct &Protein, vector<ProteinStruct> &allProteins, int index, Real cutoff, Vector &clashCC, Vector &clashCS, Vector &clashSS)
{
	int natom=Protein.Atoms.size();

	SafeAlloc(clashCC, natom, "clashCC");
	SafeAlloc(clashCS, natom, "calshCS");
	SafeAlloc(clashSS, natom, "clashSS");
	for (int i=0;i<natom;i++)
	{
		calcClashScoreByAtom(Protein.Atoms[i], i, Protein.center, allProteins, index, cutoff, clashCC, clashCS, clashSS);
	}	
}

Real calcClashScoreByAtom(ProteinStruct Protein, XRayParamStruct &params, Vector &clashCC, Vector &clashCS, Vector &clashSS, Vector &clashScoreAtom)
{
	int size;
	Real r, cutoff, buffer=3.0;
	Real clashScore=0;
	vector<ProteinStruct> allProteins, Proteins;

	removeNonCA(Protein.Atoms);
	SafePushBack(Proteins, Protein, "Proteins");
	r=calcFurthestDistance(Protein.Atoms)+buffer;
	cutoff=4.0*r*r;
	makeUnitCellAndImagesFast3(Proteins, allProteins);
	calcClashScoreByAtom(Proteins[0], allProteins, 0, cutoff, clashCC, clashCS, clashSS);
	size=clashCC.size();
	SafeAlloc(clashScoreAtom, size, "clashScoreAtom");
	for (int i=0;i<size;i++)
	{
		clashScoreAtom[i]=clashCC[i]*params.ClashWeightCoreCore;
		clashScoreAtom[i]+=clashCS[i]*params.ClashWeightCoreSurface;
		clashScoreAtom[i]+=clashSS[i]*params.ClashWeightSurfaceSurface;
		clashScore+=clashScoreAtom[i];
	}

	return clashScore;
}

Real calcClashScoreByAtom(ProteinStruct Protein, XRayParamStruct &params, Vector &clashScoreAtom)
{
	Real clashScore;
	Vector clashCC, clashCS, clashSS;
	
	clashScore=calcClashScoreByAtom(Protein, params, clashCC, clashCS, clashSS, clashScoreAtom);

	return clashScore;
}

Real calcClashScore(ProteinStruct Protein, XRayParamStruct &params, Real &clashCC, Real &clashCS, Real &clashSS)
{
	Real r, cutoff, buffer=3.0;
	Real clashScore=0;
	Real proteinCC=0, proteinCS=0, proteinSS=0;
	vector<ProteinStruct> allProteins, Proteins;
	//timeval start, end;
	clashCC=0, clashCS=0, clashSS=0;

	//gettimeofday(&start, NULL);
	removeNonCA(Protein.Atoms);
	//gettimeofday(&end, NULL);
	//cout <<"removeNonCA took "<<calcTimeDiff(start, end)<<endl;
	SafePushBack(Proteins, Protein, "Proteins");
	r=calcFurthestDistance(Protein.Atoms)+buffer;
	cutoff=4.0*r*r;
	//gettimeofday(&start, NULL);
	makeUnitCellAndImagesFast3(Proteins, allProteins);
	//makeUnitCellAndImages(Proteins, allProteins);
	//gettimeofday(&end, NULL);
	//cout <<"makeUnitCellAndImagesFast took "<<calcTimeDiff(start, end)<<endl;
	//gettimeofday(&start, NULL);
	calcClashScore(Proteins[0], allProteins, 0, cutoff, proteinCC, proteinCS, proteinSS);
	//gettimeofday(&end, NULL);
	//cout <<"calcClashScore took "<<calcTimeDiff(start, end)<<endl;
	clashCC+=proteinCC;
	clashCS+=proteinCS;
	clashSS+=proteinSS;
	clashScore=params.ClashWeightCoreCore*clashCC;
	clashScore+=params.ClashWeightCoreSurface*clashCS;
	clashScore+=params.ClashWeightSurfaceSurface*clashSS;
	return clashScore;
}

Real calcClashScore(ProteinStruct &Protein, XRayParamStruct &params)
{
	Real clashCC=0, clashCS=0, clashSS=0;
	return calcClashScore(Protein, params, clashCC, clashCS, clashSS);
}

void printClashScore(ProteinStruct &Protein, XRayParamStruct &params)
{
	Real totalClash=0;
	Real clashCC=0, clashCS=0, clashSS=0;
	totalClash=calcClashScore(Protein, params, clashCC, clashCS, clashSS);
	cout <<"clashCC= "<<clashCC<<" clashCS= "<<clashCS<<" clashSS= "<<clashSS<<" totalClash= "<<totalClash<<endl;
}

Real calcClashScore(vector<ProteinStruct> &Proteins, vector<ProteinStruct> &allProteins, Real &clashCC, Real &clashCS, Real &clashSS, XRayParamStruct &params)
{
	int nprot=Proteins.size();
	Real r, cutoff, clashScore=0, buffer=3.0;
	Real proteinCC=0, proteinCS=0, proteinSS=0;

	r=calcFurthestDistance(Proteins[0].Atoms)+buffer;
	cutoff=4.0*r*r;
	for (int i=0;i<nprot;i++) 
	{
		calcClashScore(Proteins[i], allProteins, i, cutoff, proteinCC, proteinCS, proteinSS);
		clashCC+=proteinCC;
		clashCS+=proteinCS;
		clashSS+=proteinSS;
	}
	clashScore=params.ClashWeightCoreCore*clashCC;
	clashScore+=params.ClashWeightCoreSurface*clashCS;
	clashScore+=params.ClashWeightSurfaceSurface*clashSS;
	return clashScore;

}

Real calcClashScore(Vector &degOfFreedom, ProteinStruct Protein, XRayParamStruct &params, Real &clashCC, Real &clashCS, Real &clashSS)
{
	int ndof=degOfFreedom.size();
	int nprot=ndof/6;
	Real clashScore=0;
	Vector subVector;
	vector<ProteinStruct> Proteins, allProteins;
	timeval start, end;

	clashCC=0, clashCS=0, clashSS=0;
	gettimeofday(&start, NULL);
	removeNonCA(Protein.Atoms);
	gettimeofday(&end, NULL);
	if (params.Verbose) cout <<"removeNonCA took "<<calcTimeDiff(start, end)<<endl;
	gettimeofday(&start, NULL);
	SafeAlloc(Proteins, Protein, nprot, "Proteins");
	gettimeofday(&end, NULL);
	if (params.Verbose) cout <<"SafeAlloc took "<<calcTimeDiff(start, end)<<endl;
	gettimeofday(&start, NULL);
	subVector=getSubVector(degOfFreedom, 0, 6*params.NumCopies-1);
	placeProteins(subVector, Protein, Proteins);
	gettimeofday(&end, NULL);
	if (params.Verbose) cout <<"placeProteins took "<<calcTimeDiff(start, end)<<endl;
	gettimeofday(&start, NULL);
	//makeUnitCellAndImages(Proteins, allProteins);
	makeUnitCellAndImagesFast3(Proteins, allProteins);
	gettimeofday(&end, NULL);
	if (params.Verbose) cout <<"makeUnitCellAndImages took "<<calcTimeDiff(start, end)<<endl;
	gettimeofday(&start, NULL);
	clashScore=calcClashScore(Proteins, allProteins, clashCC, clashCS, clashSS, params);
	gettimeofday(&end, NULL);
	if (params.Verbose) cout <<"calcClashScore took "<<calcTimeDiff(start, end)<<endl;
	return clashScore;
}



Real calcClashScore(Vector &degOfFreedom, ProteinStruct &Protein, XRayParamStruct &params)
{
	Real clashCC=0, clashCS=0, clashSS=0;
	return calcClashScore(degOfFreedom, Protein, params, clashCC, clashCS, clashSS);
}

void makeUnitCellAndImagesFast(vector<ProteinStruct> &Proteins)
{
	int nprot=Proteins.size();
	int noperations, index;
	int natom=Proteins[0].Atoms.size();
	Real dx, dy, dz;
	Real centerx, centery, centerz;
	VectorStruct a, b, c;

	if (nprot==0)
	{
		error("nprot= 0", __LINE__, __FILE__);
	}
	
	noperations=Proteins[0].symmetryOperations.size();
	a=Proteins[0].a;
	b=Proteins[0].b;
	c=Proteins[0].c;
	calcCenter(Proteins[0].Atoms, Proteins[0].center[X], Proteins[0].center[Y], Proteins[0].center[Z]);
	for (int i=1;i<noperations;i++)
	{
		applySymmetryOperation(Proteins[0].Atoms, Proteins[0].symmetryOperations[i], Proteins[i].Atoms);
		applySymmetryOperation(Proteins[0].center, Proteins[0].symmetryOperations[i], Proteins[i].center);
	}

	index=noperations;
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
					for (int i=0;i<noperations;i++)
					{
						centerx=Proteins[i].center[X]+dx;
						centery=Proteins[i].center[Y]+dy;
						centerz=Proteins[i].center[Z]+dz;
						for (int j=0;j<natom;j++)
						{
							Proteins[index].Atoms[j].x=Proteins[i].Atoms[j].x+dx;
							Proteins[index].Atoms[j].y=Proteins[i].Atoms[j].y+dy;
							Proteins[index].Atoms[j].z=Proteins[i].Atoms[j].z+dz;
						}
						Proteins[index].center[X]=centerx;
						Proteins[index].center[Y]=centery;
						Proteins[index].center[Z]=centerz;
						index++;
					}
				}
			}
		}
	}
}

Real calcClashScoreFast(Vector &degOfFreedom, vector<ProteinStruct> &Proteins, XRayParamStruct &params, Real &clashCC, Real &clashCS, Real &clashSS)
{
	Real clashScore=0;
	Real r, buffer=3.0, cutoff;
	ProteinStruct tempProtein;
	//timeval start, end;

	clashCC=0, clashCS=0, clashSS=0;
	//gettimeofday(&start, NULL);
	r=calcFurthestDistance(Proteins[0].Atoms)+buffer;
	//gettimeofday(&end, NULL);
	//cout <<"calcFurthestDistance took "<<calcTimeDiff(start, end)<<endl;
	cutoff=4.0*r*r;
	//gettimeofday(&start, NULL);
	tempProtein=Proteins[0];
	//gettimeofday(&end, NULL);
	//cout <<"Copy protein took "<<calcTimeDiff(start, end)<<endl;
	//gettimeofday(&start, NULL);
	placeProtein(Proteins[0], degOfFreedom);
	//gettimeofday(&end, NULL);
	//cout <<"placeProtein took "<<calcTimeDiff(start, end)<<endl;
	//gettimeofday(&start, NULL);
	makeUnitCellAndImagesFast(Proteins);
	//gettimeofday(&end, NULL);
	//cout <<"makeUnitCellAndImagesFast took "<<calcTimeDiff(start, end)<<endl;
	//gettimeofday(&start, NULL);
	calcClashScore(Proteins[0], Proteins, 0, cutoff, clashCC, clashCS, clashSS);
	//gettimeofday(&end, NULL);
	//cout <<"calcClashScore took "<<calcTimeDiff(start, end)<<endl;
	//gettimeofday(&start, NULL);
	Proteins[0]=tempProtein;
	//gettimeofday(&end, NULL);
	//cout <<"Copy protein took "<<calcTimeDiff(start, end)<<endl;
	clashScore=params.ClashWeightCoreCore*clashCC;
	clashScore+=params.ClashWeightCoreSurface*clashCS;
	clashScore+=params.ClashWeightSurfaceSurface*clashSS;
	//endProgram(__LINE__, __FILE__);
	return clashScore;
}

Real calcClashScoreFast(Vector &degOfFreedom, vector<ProteinStruct> &Proteins, XRayParamStruct &params)
{
	Real clashCC=0, clashCS=0, clashSS=0;
	return calcClashScoreFast(degOfFreedom, Proteins, params, clashCC, clashCS, clashSS);
} 

void printClashScore(Vector &degOfFreedom, ProteinStruct &Protein, XRayParamStruct &params)
{
	Real totalClash;
	Real clashCC=0, clashCS=0, clashSS=0;
	totalClash=calcClashScore(degOfFreedom, Protein, params, clashCC, clashCS, clashSS);
	cout <<"clashCC= "<<clashCC<<endl;
	cout <<"clashCS= "<<clashCS<<endl;
	cout <<"clashSS= "<<clashSS<<endl;
	cout <<"totalClash= "<<totalClash<<endl;
}

Real calcClashScore(vector<ProteinStruct> Proteins, XRayParamStruct &params)
{
	int nprot=Proteins.size();
	Real clashScore=0;
	Real clashCC=0, clashCS=0, clashSS=0;
	vector<ProteinStruct> allProteins;

	for (int i=0;i<nprot;i++) removeNonCA(Proteins[i].Atoms);
	makeUnitCellAndImages(Proteins, allProteins);
	clashScore=calcClashScore(Proteins, allProteins, clashCC, clashCS, clashSS, params);
	return clashScore;
}

Real calcClashScore(Vector &degOfFreedoms, vector<ProteinStruct> Proteins, XRayParamStruct &params)
{
	placeProteins(Proteins, degOfFreedoms);
	return calcClashScore(Proteins, params); 
}
