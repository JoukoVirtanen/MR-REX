# include <vector>
# include <iostream>
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
# include "../../LibraryFiles/PrintPdb.h"
# include "../../LibraryFiles/MinMax.h"
# include "XRayStruct.h"
# include "molecularReplacement.h"
# include "XRay.h"
# include "LatticeStruct.h"
# include "LatticeUtils.h"
# include "SpaceGroup.h"

void setSpaceGroupA1(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;

	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.translationVector[X]=Protein.b.x*0.5+Protein.c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=Protein.b.y*0.5+Protein.c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=Protein.b.z*0.5+Protein.c.z*0.5;
	
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	printSymmetryOperations(Protein.symmetryOperations);
}

void setSpaceGroupA2(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=Protein.b.x*0.5+Protein.c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=Protein.b.y*0.5+Protein.c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=Protein.b.z*0.5+Protein.c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupB112(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=Protein.a.x*0.5+Protein.c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=Protein.a.y*0.5+Protein.c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=Protein.a.z*0.5+Protein.c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=Protein.a.x*0.5+Protein.c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=Protein.a.y*0.5+Protein.c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=Protein.a.z*0.5+Protein.c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupB2(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=Protein.a.x*0.5+Protein.c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=Protein.a.y*0.5+Protein.c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=Protein.a.z*0.5+Protein.c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=Protein.a.x*0.5+Protein.c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=Protein.a.y*0.5+Protein.c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=Protein.a.z*0.5+Protein.c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupB2212(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=0;
	tempSymmetryOperation.translationVector[Y]=0;
	tempSymmetryOperation.translationVector[Z]=0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupC121(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupC1211(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupC222(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");
	
	//1555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//2555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//3555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//4555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//5555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//6555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//7555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//8555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupC2221(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");
	
	//1555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//2555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//3555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//4555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=0;
	tempSymmetryOperation.translationVector[Y]=0;
	tempSymmetryOperation.translationVector[Z]=0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//5555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//6555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//7555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//8555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

/*
void setSpaceGroupC2221(ProteinStruct &Protein)
{
	//This is wrong.  This is just for a test.
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");
	
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}
*/

void setSpaceGroupC21(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");
	
	//1555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//2555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//3555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//4555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

/*
void setSpaceGroupC4212(ProteinStruct &Protein)
{
	//This space group makes no sense
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");
	
	//1555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//2555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[X][Y]=0.0;
	tempSymmetryOperation.rotationMat[X][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=0.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=-1.41417;
	tempSymmetryOperation.rotationMat[Z][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=0;
	tempSymmetryOperation.translationVector[Y]=0;
	tempSymmetryOperation.translationVector[Z]=0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//3555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=-1.41417;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=-2.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0.0;
	tempSymmetryOperation.rotationMat[Z][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//4555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=1.41417;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=2.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=-1.41417;
	tempSymmetryOperation.rotationMat[Z][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//5555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=-1.41417;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=-2.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=1.41417;
	tempSymmetryOperation.rotationMat[Z][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//6555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=1.41417;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=2.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0.0;
	tempSymmetryOperation.rotationMat[Z][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//7555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=0.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=1.41417;
	tempSymmetryOperation.rotationMat[Z][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=0.0;
	tempSymmetryOperation.translationVector[Y]=0.0;
	tempSymmetryOperation.translationVector[Z]=0.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//8555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=0.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0.0;
	tempSymmetryOperation.rotationMat[Z][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=0.0;
	tempSymmetryOperation.translationVector[Y]=0.0;
	tempSymmetryOperation.translationVector[Z]=0.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//9555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=0.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0.0;
	tempSymmetryOperation.rotationMat[Z][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//10555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=0.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=-1.41417;
	tempSymmetryOperation.rotationMat[Z][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//11555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=-1.41417;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=-2.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0.0;
	tempSymmetryOperation.rotationMat[Z][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=0.0;
	tempSymmetryOperation.translationVector[Y]=0.0;
	tempSymmetryOperation.translationVector[Z]=0.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//12555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=1.41417;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=2.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=-1.41417;
	tempSymmetryOperation.rotationMat[Z][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=0.0;
	tempSymmetryOperation.translationVector[Y]=0.0;
	tempSymmetryOperation.translationVector[Z]=0.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//13555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=-1.41417;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=-2.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=1.41417;
	tempSymmetryOperation.rotationMat[Z][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=0.0;
	tempSymmetryOperation.translationVector[Y]=0.0;
	tempSymmetryOperation.translationVector[Z]=0.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//14555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=1.41417;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=2.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0.0;
	tempSymmetryOperation.rotationMat[Z][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=0.0;
	tempSymmetryOperation.translationVector[Y]=0.0;
	tempSymmetryOperation.translationVector[Z]=0.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//15555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=0.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=1.41417;
	tempSymmetryOperation.rotationMat[Z][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//16555
	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=0.0;
	tempSymmetryOperation.rotationMat[Z][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Y][X]=0.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=0.0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0.0;
	tempSymmetryOperation.rotationMat[Z][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}
*/

void setSpaceGroupF23(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupF4132(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

/*
void setSpaceGroupF422(ProteinStruct &Protein)
{
	//Inconsistant space group
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}
*/

void setSpaceGroupF222(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupF432(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupH3(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.666667+b.x*0.333333+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.666667+b.y*0.333333+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.666667+b.z*0.333333+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.666667+b.x*0.333333+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.666667+b.y*0.333333+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.666667+b.z*0.333333+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.666667+b.x*0.333333+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.666667+b.y*0.333333+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.666667+b.z*0.333333+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.333333+b.x*0.666667+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.333333+b.y*0.666667+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.333333+b.z*0.666667+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.333333+b.x*0.666667+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.333333+b.y*0.666667+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.333333+b.z*0.666667+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.333333+b.x*0.666667+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.333333+b.y*0.666667+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.333333+b.z*0.666667+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupH32(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	//SYM1
	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM2
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM3
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM4
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM5
	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM6
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM7
	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.666667+b.x*0.333333+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=a.y*0.666667+b.y*0.333333+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=a.z*0.666667+b.z*0.333333+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM8
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.666667+b.x*0.333333+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=a.y*0.666667+b.y*0.333333+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=a.z*0.666667+b.z*0.333333+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM9
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.666667+b.x*0.333333+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=a.y*0.666667+b.y*0.333333+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=a.z*0.666667+b.z*0.333333+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM10
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.666667+b.x*0.333333+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=a.y*0.666667+b.y*0.333333+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=a.z*0.666667+b.z*0.333333+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM11
	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.666667+b.x*0.333333+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=a.y*0.666667+b.y*0.333333+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=a.z*0.666667+b.z*0.333333+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM12
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.666667+b.x*0.333333+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=a.y*0.666667+b.y*0.333333+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=a.z*0.666667+b.z*0.333333+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM13
	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.333333+b.x*0.666667+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=a.y*0.333333+b.y*0.666667+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=a.z*0.333333+b.z*0.666667+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM14
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.333333+b.x*0.666667+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=a.y*0.333333+b.y*0.666667+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=a.z*0.333333+b.z*0.666667+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM15
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.333333+b.x*0.666667+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=a.y*0.333333+b.y*0.666667+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=a.z*0.333333+b.z*0.666667+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM16
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.333333+b.x*0.666667+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=a.y*0.333333+b.y*0.666667+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=a.z*0.333333+b.z*0.666667+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM17
	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.333333+b.x*0.666667+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=a.y*0.333333+b.y*0.666667+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=a.z*0.333333+b.z*0.666667+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//SYM18
	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.333333+b.x*0.666667+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=a.y*0.333333+b.y*0.666667+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=a.z*0.333333+b.z*0.666667+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//cout <<"a= "<<vectorToStr(a)<<endl;
	//cout <<"b= "<<vectorToStr(b)<<endl;
	//cout <<"c= "<<vectorToStr(c)<<endl;
	//printSymmetryOperations(Protein.symmetryOperations);
	//exit(EXIT_FAILURE);
}

void setSpaceGroupI121(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI1211(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI21(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI212121(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI213(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI222(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI23(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI4(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI41(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI41_A(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI411(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI4122(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI4132(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI_42D(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-0.002443;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0.001048;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.00035;
	tempSymmetryOperation.rotationMat[0][1]=0.998522;
	tempSymmetryOperation.rotationMat[0][2]=0.000698;

	tempSymmetryOperation.rotationMat[1][0]=-1.00148;
	tempSymmetryOperation.rotationMat[1][1]=0.00035;
	tempSymmetryOperation.rotationMat[1][2]=-0.001748;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.00035;
	tempSymmetryOperation.rotationMat[0][1]=-0.998522;
	tempSymmetryOperation.rotationMat[0][2]=0.001745;

	tempSymmetryOperation.rotationMat[1][0]=1.00148;
	tempSymmetryOperation.rotationMat[1][1]=-0.00035;
	tempSymmetryOperation.rotationMat[1][2]=0.0007;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0.000698;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=-0.001048;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=-0.000698;
	tempSymmetryOperation.rotationMat[0][2]=0.002444;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.00035;
	tempSymmetryOperation.rotationMat[0][1]=-0.998522;
	tempSymmetryOperation.rotationMat[0][2]=-0.000699;

	tempSymmetryOperation.rotationMat[1][0]=-1.00148;
	tempSymmetryOperation.rotationMat[1][1]=0.00035;
	tempSymmetryOperation.rotationMat[1][2]=-0.0007;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.00035;
	tempSymmetryOperation.rotationMat[0][1]=0.998522;
	tempSymmetryOperation.rotationMat[0][2]=-0.001745;

	tempSymmetryOperation.rotationMat[1][0]=1.00148;
	tempSymmetryOperation.rotationMat[1][1]=-0.00035;
	tempSymmetryOperation.rotationMat[1][2]=0.001748;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-0.002443;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0.001048;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.00035;
	tempSymmetryOperation.rotationMat[0][1]=0.998522;
	tempSymmetryOperation.rotationMat[0][2]=0.000698;

	tempSymmetryOperation.rotationMat[1][0]=-1.00148;
	tempSymmetryOperation.rotationMat[1][1]=0.00035;
	tempSymmetryOperation.rotationMat[1][2]=-0.001748;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.00035;
	tempSymmetryOperation.rotationMat[0][1]=-0.998522;
	tempSymmetryOperation.rotationMat[0][2]=0.001745;

	tempSymmetryOperation.rotationMat[1][0]=1.00148;
	tempSymmetryOperation.rotationMat[1][1]=-0.00035;
	tempSymmetryOperation.rotationMat[1][2]=0.0007;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0.000698;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=-0.001048;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=-0.000698;
	tempSymmetryOperation.rotationMat[0][2]=0.002444;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.00035;
	tempSymmetryOperation.rotationMat[0][1]=-0.998522;
	tempSymmetryOperation.rotationMat[0][2]=-0.000699;

	tempSymmetryOperation.rotationMat[1][0]=-1.00148;
	tempSymmetryOperation.rotationMat[1][1]=0.00035;
	tempSymmetryOperation.rotationMat[1][2]=-0.0007;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.00035;
	tempSymmetryOperation.rotationMat[0][1]=0.998522;
	tempSymmetryOperation.rotationMat[0][2]=-0.001745;

	tempSymmetryOperation.rotationMat[1][0]=1.00148;
	tempSymmetryOperation.rotationMat[1][1]=-0.00035;
	tempSymmetryOperation.rotationMat[1][2]=0.001748;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI_4C2(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI422(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupI432(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP1(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	
	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;

	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupP_1(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP1121(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	//cout <<"a= "<<vectorToStr(a)<<endl;
	//cout <<"b= "<<vectorToStr(b)<<endl;
	//cout <<"c= "<<vectorToStr(c)<<endl;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP121(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP121_N1(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP121_C1(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP1211(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP21221(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP213(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP21312(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP21212(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupP21212A(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP212121(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[X][X]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;
	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupP222(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP2221(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP3(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP3112(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP312(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP3212(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP3221(ProteinStruct &Protein)
{
	//Check to make sure this is correct
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	//1555	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//2555
	tempSymmetryOperation.rotationMat[X][X]=-0.5;
	tempSymmetryOperation.rotationMat[X][Y]=-0.866025;
	tempSymmetryOperation.rotationMat[X][Z]=0;

	tempSymmetryOperation.rotationMat[Y][X]=0.866025;
	tempSymmetryOperation.rotationMat[Y][Y]=-0.5;
	tempSymmetryOperation.rotationMat[Y][Z]=0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0;
	tempSymmetryOperation.rotationMat[Z][Y]=0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;

	tempSymmetryOperation.translationVector[X]=2.0*c.x/3.0;
	tempSymmetryOperation.translationVector[Y]=2.0*c.y/3.0;
	tempSymmetryOperation.translationVector[Z]=2.0*c.z/3.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//3555
	tempSymmetryOperation.rotationMat[X][X]=-0.5;
	tempSymmetryOperation.rotationMat[X][Y]=0.866025;
	tempSymmetryOperation.rotationMat[X][Z]=0;

	tempSymmetryOperation.rotationMat[Y][X]=-0.866025;
	tempSymmetryOperation.rotationMat[Y][Y]=-0.5;
	tempSymmetryOperation.rotationMat[Y][Z]=0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0;
	tempSymmetryOperation.rotationMat[Z][Y]=0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;

	tempSymmetryOperation.translationVector[X]=c.x/3.0;
	tempSymmetryOperation.translationVector[Y]=c.y/3.0;
	tempSymmetryOperation.translationVector[Z]=c.z/3.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//4555
	tempSymmetryOperation.rotationMat[X][X]=-0.5;
	tempSymmetryOperation.rotationMat[X][Y]=0.866025;
	tempSymmetryOperation.rotationMat[X][Z]=0;

	tempSymmetryOperation.rotationMat[Y][X]=0.866025;
	tempSymmetryOperation.rotationMat[Y][Y]=0.5;
	tempSymmetryOperation.rotationMat[Y][Z]=0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0;
	tempSymmetryOperation.rotationMat[Z][Y]=0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;

	tempSymmetryOperation.translationVector[X]=0;
	tempSymmetryOperation.translationVector[Y]=0;
	tempSymmetryOperation.translationVector[Z]=0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//5555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[X][Y]=0;
	tempSymmetryOperation.rotationMat[X][Z]=0;

	tempSymmetryOperation.rotationMat[Y][X]=0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0;
	tempSymmetryOperation.rotationMat[Z][Y]=0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;

	tempSymmetryOperation.translationVector[X]=c.x/3.0;
	tempSymmetryOperation.translationVector[Y]=c.y/3.0;
	tempSymmetryOperation.translationVector[Z]=c.z/3.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//6555
	tempSymmetryOperation.rotationMat[X][X]=-0.5;
	tempSymmetryOperation.rotationMat[X][Y]=-0.866025;
	tempSymmetryOperation.rotationMat[X][Z]=0;

	tempSymmetryOperation.rotationMat[Y][X]=-0.866025;
	tempSymmetryOperation.rotationMat[Y][Y]=0.5;
	tempSymmetryOperation.rotationMat[Y][Z]=0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0;
	tempSymmetryOperation.rotationMat[Z][Y]=0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;

	tempSymmetryOperation.translationVector[X]=2.0*c.x/3.0;
	tempSymmetryOperation.translationVector[Y]=2.0*c.y/3.0;
	tempSymmetryOperation.translationVector[Z]=2.0*c.z/3.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	
}

void setSpaceGroupP3121(ProteinStruct &Protein)
{
	//Check to make sure this is correct
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	//1555	
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[Y][Y]=1.0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	//2555
	tempSymmetryOperation.rotationMat[X][X]=-0.5;
	tempSymmetryOperation.rotationMat[X][Y]=-0.866025;
	tempSymmetryOperation.rotationMat[X][Z]=0;

	tempSymmetryOperation.rotationMat[Y][X]=0.866025;
	tempSymmetryOperation.rotationMat[Y][Y]=-0.5;
	tempSymmetryOperation.rotationMat[Y][Z]=0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0;
	tempSymmetryOperation.rotationMat[Z][Y]=0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;

	tempSymmetryOperation.translationVector[X]=c.x/3.0;
	tempSymmetryOperation.translationVector[Y]=c.y/3.0;
	tempSymmetryOperation.translationVector[Z]=c.z/3.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//3555
	tempSymmetryOperation.rotationMat[X][X]=-0.5;
	tempSymmetryOperation.rotationMat[X][Y]=0.866025;
	tempSymmetryOperation.rotationMat[X][Z]=0;

	tempSymmetryOperation.rotationMat[Y][X]=-0.866025;
	tempSymmetryOperation.rotationMat[Y][Y]=-0.5;
	tempSymmetryOperation.rotationMat[Y][Z]=0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0;
	tempSymmetryOperation.rotationMat[Z][Y]=0;
	tempSymmetryOperation.rotationMat[Z][Z]=1.0;

	tempSymmetryOperation.translationVector[X]=2.0*c.x/3.0;
	tempSymmetryOperation.translationVector[Y]=2.0*c.y/3.0;
	tempSymmetryOperation.translationVector[Z]=2.0*c.z/3.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//4555
	tempSymmetryOperation.rotationMat[X][X]=-0.5;
	tempSymmetryOperation.rotationMat[X][Y]=0.866025;
	tempSymmetryOperation.rotationMat[X][Z]=0;

	tempSymmetryOperation.rotationMat[Y][X]=0.866025;
	tempSymmetryOperation.rotationMat[Y][Y]=0.5;
	tempSymmetryOperation.rotationMat[Y][Z]=0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0;
	tempSymmetryOperation.rotationMat[Z][Y]=0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;

	tempSymmetryOperation.translationVector[X]=0;
	tempSymmetryOperation.translationVector[Y]=0;
	tempSymmetryOperation.translationVector[Z]=0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//5555
	tempSymmetryOperation.rotationMat[X][X]=1.0;
	tempSymmetryOperation.rotationMat[X][Y]=0;
	tempSymmetryOperation.rotationMat[X][Z]=0;

	tempSymmetryOperation.rotationMat[Y][X]=0;
	tempSymmetryOperation.rotationMat[Y][Y]=-1.0;
	tempSymmetryOperation.rotationMat[Y][Z]=0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0;
	tempSymmetryOperation.rotationMat[Z][Y]=0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;

	tempSymmetryOperation.translationVector[X]=2.0*c.x/3.0;
	tempSymmetryOperation.translationVector[Y]=2.0*c.y/3.0;
	tempSymmetryOperation.translationVector[Z]=2.0*c.z/3.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
	
	//6555
	tempSymmetryOperation.rotationMat[X][X]=-0.5;
	tempSymmetryOperation.rotationMat[X][Y]=-0.866025;
	tempSymmetryOperation.rotationMat[X][Z]=0;

	tempSymmetryOperation.rotationMat[Y][X]=-0.866025;
	tempSymmetryOperation.rotationMat[Y][Y]=0.5;
	tempSymmetryOperation.rotationMat[Y][Z]=0;
	
	tempSymmetryOperation.rotationMat[Z][X]=0;
	tempSymmetryOperation.rotationMat[Z][Y]=0;
	tempSymmetryOperation.rotationMat[Z][Z]=-1.0;

	tempSymmetryOperation.translationVector[X]=c.x/3.0;
	tempSymmetryOperation.translationVector[Y]=c.y/3.0;
	tempSymmetryOperation.translationVector[Z]=c.z/3.0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);
}

void setSpaceGroupP41(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP41212(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP4122(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP4132(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP4212(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP42212(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP4222(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP43(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP432(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP43212(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0.5+b.x*0.5+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=a.y*0.5+b.y*0.5+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=a.z*0.5+b.z*0.5+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP4322(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP4332(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0.5+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0.5+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0.5+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.5+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.5+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.5+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=-1.0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=1.0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.75+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.75+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.75+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.75+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.75+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.75+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.75+b.x*0.25+c.x*0.75;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.75+b.y*0.25+c.y*0.75;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.75+b.z*0.25+c.z*0.75;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=-1.0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=-1.0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=0;

	tempSymmetryOperation.translationVector[X]=+a.x*0.25+b.x*0.25+c.x*0.25;
	tempSymmetryOperation.translationVector[Y]=+a.y*0.25+b.y*0.25+c.y*0.25;
	tempSymmetryOperation.translationVector[Z]=+a.z*0.25+b.z*0.25+c.z*0.25;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP6(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP61(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.833333;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.833333;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.833333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.166667;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.166667;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.166667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP6122(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.833333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.833333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.833333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.166667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.166667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.166667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.833333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.833333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.833333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.166667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.166667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.166667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP22121(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0.5+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0.5+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0.5+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP62(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP622(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP6222(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP63(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP31(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP32(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP422(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP321(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP4(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP6322(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP64(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP42(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=-1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0;
	tempSymmetryOperation.rotationMat[0][1]=1.0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-1.0;
	tempSymmetryOperation.rotationMat[1][1]=0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP6422(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP65(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	//Uncomment
	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.666667;
	//tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0;
	//tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0;
	//tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.166667;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.166667;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.166667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=a.x*0+b.x*0+c.x*0.833333;
	tempSymmetryOperation.translationVector[Y]=a.y*0+b.y*0+c.y*0.833333;
	tempSymmetryOperation.translationVector[Z]=a.z*0+b.z*0+c.z*0.833333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupR32(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.329691;
	tempSymmetryOperation.rotationMat[0][1]=0.234083;
	tempSymmetryOperation.rotationMat[0][2]=0.914609;

	tempSymmetryOperation.rotationMat[1][0]=0.944089;
	tempSymmetryOperation.rotationMat[1][1]=-0.081745;
	tempSymmetryOperation.rotationMat[1][2]=-0.319396;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0.968774;
	tempSymmetryOperation.rotationMat[2][2]=-0.247945;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.329691;
	tempSymmetryOperation.rotationMat[0][1]=0.944089;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.234083;
	tempSymmetryOperation.rotationMat[1][1]=-0.081745;
	tempSymmetryOperation.rotationMat[1][2]=0.968774;

	tempSymmetryOperation.rotationMat[2][0]=0.914609;
	tempSymmetryOperation.rotationMat[2][1]=-0.319396;
	tempSymmetryOperation.rotationMat[2][2]=-0.247945;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.329691;
	tempSymmetryOperation.rotationMat[0][1]=-0.944089;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.944089;
	tempSymmetryOperation.rotationMat[1][1]=0.329691;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-0.247945;
	tempSymmetryOperation.rotationMat[1][2]=-0.968774;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=-0.968774;
	tempSymmetryOperation.rotationMat[2][2]=0.247945;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.329691;
	tempSymmetryOperation.rotationMat[0][1]=-0.234083;
	tempSymmetryOperation.rotationMat[0][2]=-0.914609;

	tempSymmetryOperation.rotationMat[1][0]=-0.234083;
	tempSymmetryOperation.rotationMat[1][1]=-0.918255;
	tempSymmetryOperation.rotationMat[1][2]=0.319396;

	tempSymmetryOperation.rotationMat[2][0]=-0.914609;
	tempSymmetryOperation.rotationMat[2][1]=0.319396;
	tempSymmetryOperation.rotationMat[2][2]=0.247945;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupP6522(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.166667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.166667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.166667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.833333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.833333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.833333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.666667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.666667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.666667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=-1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.333333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.333333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.333333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=-0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.166667;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.166667;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.166667;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.5;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.5;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.5;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=0.5;
	tempSymmetryOperation.rotationMat[0][1]=0.866025;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0.866025;
	tempSymmetryOperation.rotationMat[1][1]=-0.5;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=-1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0.833333;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0.833333;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0.833333;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void setSpaceGroupR3(ProteinStruct &Protein)
{
	SymmetryOperationStruct tempSymmetryOperation;
	VectorStruct a, b, c;

	a=Protein.a;
	b=Protein.b;
	c=Protein.c;

	Protein.symmetryOperations.clear();

	SafeAlloc(tempSymmetryOperation.translationVector, 3, "translation");
	Safe2DAlloc(tempSymmetryOperation.rotationMat, 3, 3, "rotationMat");

	tempSymmetryOperation.rotationMat[0][0]=1.0;
	tempSymmetryOperation.rotationMat[0][1]=0;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=0;
	tempSymmetryOperation.rotationMat[1][1]=1.0;
	tempSymmetryOperation.rotationMat[1][2]=0;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0;
	tempSymmetryOperation.rotationMat[2][2]=1.0;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.372502;
	tempSymmetryOperation.rotationMat[0][1]=-0.550908;
	tempSymmetryOperation.rotationMat[0][2]=0.746822;

	tempSymmetryOperation.rotationMat[1][0]=0.928031;
	tempSymmetryOperation.rotationMat[1][1]=-0.221128;
	tempSymmetryOperation.rotationMat[1][2]=0.299766;

	tempSymmetryOperation.rotationMat[2][0]=0;
	tempSymmetryOperation.rotationMat[2][1]=0.804738;
	tempSymmetryOperation.rotationMat[2][2]=0.59363;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

	tempSymmetryOperation.rotationMat[0][0]=-0.372502;
	tempSymmetryOperation.rotationMat[0][1]=0.928031;
	tempSymmetryOperation.rotationMat[0][2]=0;

	tempSymmetryOperation.rotationMat[1][0]=-0.550908;
	tempSymmetryOperation.rotationMat[1][1]=-0.221128;
	tempSymmetryOperation.rotationMat[1][2]=0.804738;

	tempSymmetryOperation.rotationMat[2][0]=0.746822;
	tempSymmetryOperation.rotationMat[2][1]=0.299766;
	tempSymmetryOperation.rotationMat[2][2]=0.59363;

	tempSymmetryOperation.translationVector[X]=+a.x*0+b.x*0+c.x*0;
	tempSymmetryOperation.translationVector[Y]=+a.y*0+b.y*0+c.y*0;
	tempSymmetryOperation.translationVector[Z]=+a.z*0+b.z*0+c.z*0;
	SafePushBack(Protein.symmetryOperations, tempSymmetryOperation);

}

void printSpaceGroups()
{
	cout <<"The valid space groups are: "<<endl;
	cout <<"A1" <<endl;
	cout <<"A2" <<endl; 
	cout <<"B112" <<endl;
	cout <<"B2" <<endl;
	cout <<"B2212" <<endl;
	cout <<"C121" <<endl;
	cout <<"C1211" <<endl;
	cout <<"C21" <<endl;
	cout <<"C222" <<endl;
	cout <<"C2221" <<endl;
	cout <<"F4132" <<endl;
	cout <<"F222" <<endl;
	cout <<"F23" <<endl;
	cout <<"F432" <<endl;
	cout <<"H3" <<endl;
	cout <<"H32" <<endl;
	cout <<"I121" <<endl;
	cout <<"I1211" <<endl;
	cout <<"I21" <<endl;
	cout <<"I212121" <<endl;
	cout <<"I213" <<endl;
	cout <<"I222" <<endl;
	cout <<"I23" <<endl;
	cout <<"I-42D" <<endl;
	cout <<"I-4C2" <<endl;
	cout <<"I4" <<endl;
	cout <<"I41" <<endl;
	cout <<"I41/A" <<endl;
	cout <<"I411" <<endl;
	cout <<"I4122" <<endl;
	cout <<"I4132" <<endl;
	cout <<"I422" <<endl;
	cout <<"I432" <<endl;
	cout <<"P1" <<endl;
	cout <<"P-1" <<endl;
	cout <<"P1121" <<endl;
	cout <<"P121" <<endl;
	cout <<"P121/C1" <<endl;
	cout <<"P121/N1" <<endl;
	cout <<"P1211" <<endl;
	cout <<"P21212" <<endl;
	cout <<"P212121" <<endl;
	cout <<"P21212A" <<endl;
	cout <<"P21221" <<endl;
	cout <<"P21312" <<endl;
	cout <<"P213" <<endl;
	cout <<"P222" <<endl;
	cout <<"P2221" <<endl;
	cout <<"P3" <<endl;
	cout <<"P31" <<endl;
	cout <<"P3112" <<endl;
	cout <<"P3121" <<endl;
	cout <<"P312" <<endl;
	cout <<"P32" <<endl;
	cout <<"P321" <<endl;
	cout <<"P3212" <<endl;
	cout <<"P3221" <<endl;
	cout <<"P4" <<endl;
	cout <<"P41212" <<endl;
	cout <<"P4132" <<endl;
	cout <<"P42" <<endl;
	cout <<"P4212" <<endl;
	cout <<"P422" <<endl;
	cout <<"P22121" <<endl;
	cout <<"P4222" <<endl;
	cout <<"P43" <<endl;
	cout <<"P4332" <<endl;
	cout <<"P6" <<endl;
	cout <<"P61" <<endl;
	cout <<"P6122" <<endl;
	cout <<"P62" <<endl;
	cout <<"P622" <<endl;
	cout <<"P6222" <<endl;
	cout <<"P63" <<endl;
	cout <<"P6322" <<endl;
	cout <<"P64" <<endl;
	cout <<"P6422" <<endl;
	cout <<"P65" <<endl;
	cout <<"P6522" <<endl;
	cout <<"R32" <<endl;
	cout <<"R3" <<endl;	
}

int spaceGroupStrToSpaceGroupInt(string spaceGroup)
{
	if (spaceGroup=="A1") return A1;
	else if (spaceGroup=="A2") return A2;
	else if (spaceGroup=="B112") return B112;
	else if (spaceGroup=="B2") return B2;
	else if (spaceGroup=="B2212") return B2212;
	else if (spaceGroup=="C121") return C121;
	else if (spaceGroup=="C1211") return C1211;
	else if (spaceGroup=="C21") return C21;
	else if (spaceGroup=="C222") return C222;
	else if (spaceGroup=="C2221") return C2221;
	else if (spaceGroup=="F4132") return F4132;
	else if (spaceGroup=="F222") return F222;
	else if (spaceGroup=="F23") return F23;
	else if (spaceGroup=="F432") return F432;
	else if (spaceGroup=="H3") return H3;
	else if (spaceGroup=="H32") return H32;
	else if (spaceGroup=="I121") return I121;
	else if (spaceGroup=="I1211") return I1211;
	else if (spaceGroup=="I21") return I21;
	else if (spaceGroup=="I212121") return I212121;
	else if (spaceGroup=="I213") return I213;
	else if (spaceGroup=="I222") return I222;
	else if (spaceGroup=="I23") return I23;
	else if (spaceGroup=="I-42D") return I_42D;
	else if (spaceGroup=="I-4C2") return I_4C2;
	else if (spaceGroup=="I4") return I4;
	else if (spaceGroup=="I41") return I41;
	else if (spaceGroup=="I41/A") return I41_A;
	else if (spaceGroup=="I411") return I411;
	else if (spaceGroup=="I4122") return I4122;
	else if (spaceGroup=="I4132") return I4132;
	else if (spaceGroup=="I422") return I422;
	else if (spaceGroup=="I432") return I432;
	else if (spaceGroup=="P1") return P1;
	else if (spaceGroup=="P-1") return P_1;
	else if (spaceGroup=="P1121") return P1121;
	else if (spaceGroup=="P121") return P121;
	else if (spaceGroup=="P121/C1") return P121_C1;
	else if (spaceGroup=="P121/N1") return P121_N1;
	else if (spaceGroup=="P1211") return P1211;
	else if (spaceGroup=="P21212") return P21212;
	else if (spaceGroup=="P212121") return P212121;
	else if (spaceGroup=="P21212A") return P21212A;
	else if (spaceGroup=="P21221") return P21221;
	else if (spaceGroup=="P21312") return P21312;
	else if (spaceGroup=="P213") return P213;
	else if (spaceGroup=="P222") return P222;
	else if (spaceGroup=="P2221") return P2221;
	else if (spaceGroup=="P3") return P3;
	else if (spaceGroup=="P31") return P31;
	else if (spaceGroup=="P3112") return P3112;
	else if (spaceGroup=="P3121") return P3121;
	else if (spaceGroup=="P312") return P312;
	else if (spaceGroup=="P32") return P32;
	else if (spaceGroup=="P321") return P321;
	else if (spaceGroup=="P3212") return P3212;
	else if (spaceGroup=="P3221") return P3221;
	else if (spaceGroup=="P4") return P4;
	else if (spaceGroup=="P41") return P41;
	else if (spaceGroup=="P41212") return P41212;
	else if (spaceGroup=="P4122") return P4122;
	else if (spaceGroup=="P4132") return P4132;
	else if (spaceGroup=="P42") return P42;
	else if (spaceGroup=="P4212") return P4212;
	else if (spaceGroup=="P422") return P422;
	else if (spaceGroup=="P22121") return P22121;
	else if (spaceGroup=="P42212") return P42212;
	else if (spaceGroup=="P4222") return P4222;
	else if (spaceGroup=="P43") return P43;
	else if (spaceGroup=="P432") return P432;
	else if (spaceGroup=="P4322") return P4322;
	else if (spaceGroup=="P4332") return P4332;
	else if (spaceGroup=="P43212") return P43212;
	else if (spaceGroup=="P6") return P6;
	else if (spaceGroup=="P61") return P61;
	else if (spaceGroup=="P6122") return P6122;
	else if (spaceGroup=="P62") return P62;
	else if (spaceGroup=="P622") return P622;
	else if (spaceGroup=="P6222") return P6222;
	else if (spaceGroup=="P63") return P63;
	else if (spaceGroup=="P6322") return P6322;
	else if (spaceGroup=="P64") return P64;
	else if (spaceGroup=="P6422") return P6422;
	else if (spaceGroup=="P65") return P65;
	else if (spaceGroup=="P6522") return P6522;
	else if (spaceGroup=="R32") return R32;
	else if (spaceGroup=="R3") return R3;
	else
	{
		string errorStr="Unknown space group "+spaceGroup;
		error(errorStr, __LINE__, __FILE__);
	}
	return UNK_INT;
}

void getSpaceGroups(vector<string> &spaceGroups)
{
	string str="UNK";
	SafeAlloc(spaceGroups, str, 87, "spaceGroups");
	spaceGroups[A1]="A1";
	spaceGroups[A2]="A2"; 
	spaceGroups[B112]="B112";
	spaceGroups[B2]="B2";
	spaceGroups[B2212]="B2212";
	spaceGroups[C121]="C121";
	spaceGroups[C1211]="C1211";
	spaceGroups[C21]="C21";
	spaceGroups[C222]="C222";
	spaceGroups[C2221]="C2221";
	spaceGroups[F4132]="F4132";
	spaceGroups[F222]="F222";
	spaceGroups[F23]="F23";
	spaceGroups[F432]="F432";
	spaceGroups[H3]="H3";
	spaceGroups[H32]="H32";
	spaceGroups[I121]="I121";
	spaceGroups[I1211]="I1211";
	spaceGroups[I21]="I21";
	spaceGroups[I212121]="I212121";
	spaceGroups[I213]="I213";
	spaceGroups[I222]="I222";
	spaceGroups[I23]="I23";
	spaceGroups[I_42D]="I-42D";
	spaceGroups[I_4C2]="I-4C2";
	spaceGroups[I4]="I4";
	spaceGroups[I41]="I41";
	spaceGroups[I41_A]="I41/A";
	spaceGroups[I411]="I411";
	spaceGroups[I4122]="I4122";
	spaceGroups[I4132]="I4132";
	spaceGroups[I422]="I422";
	spaceGroups[I432]="I432";
	spaceGroups[P1]="P1";
	spaceGroups[P_1]="P-1";
	spaceGroups[P1121]="P1121";
	spaceGroups[P121]="P121";
	spaceGroups[P121_C1]="P121/C1";
	spaceGroups[P121_N1]="P121/N1";
	spaceGroups[P1211]="P1211";
	spaceGroups[P21212]="P21212";
	spaceGroups[P212121]="P212121";
	spaceGroups[P21212A]="P21212A";
	spaceGroups[P21221]="P21221";
	spaceGroups[P21312]="P21312";
	spaceGroups[P213]="P213";
	spaceGroups[P222]="P222";
	spaceGroups[P2221]="P2221";
	spaceGroups[P3]="P3";
	spaceGroups[P31]="P31";
	spaceGroups[P3112]="P3112";
	spaceGroups[P3121]="P3121";
	spaceGroups[P312]="P312";
	spaceGroups[P32]="P32";
	spaceGroups[P321]="P321";
	spaceGroups[P3212]="P3212";
	spaceGroups[P3221]="P3221";
	spaceGroups[P4]="P4";
	spaceGroups[P41]="P41";
	spaceGroups[P41212]="P41212";
	spaceGroups[P4122]="P4122";
	spaceGroups[P4132]="P4132";
	spaceGroups[P42]="P42";
	spaceGroups[P4212]="P4212";
	spaceGroups[P422]="P422";
	spaceGroups[P22121]="P22121";
	spaceGroups[P42212]="P42212";
	spaceGroups[P4222]="P4222";
	spaceGroups[P43]="P43";
	spaceGroups[P432]="P432";
	spaceGroups[P4322]="P4322";
	spaceGroups[P4332]="P4332";
	spaceGroups[P43212]="P43212";
	spaceGroups[P6]="P6";
	spaceGroups[P61]="P61";
	spaceGroups[P6122]="P6122";
	spaceGroups[P62]="P62";
	spaceGroups[P622]="P622";
	spaceGroups[P6222]="P6222";
	spaceGroups[P63]="P63";
	spaceGroups[P6322]="P6322";
	spaceGroups[P64]="P64";
	spaceGroups[P6422]="P6422";
	spaceGroups[P65]="P65";
	spaceGroups[P6522]="P6522";
	spaceGroups[R32]="R32";
	spaceGroups[R3]="R3";	
}

void setSpaceGroup(ProteinStruct &Protein, string spaceGroup)
{
	if (spaceGroup=="A1") setSpaceGroupA1(Protein);
	else if (spaceGroup=="A2") setSpaceGroupA2(Protein); 
	else if (spaceGroup=="B112") setSpaceGroupB112(Protein);
	else if (spaceGroup=="B2") setSpaceGroupB2(Protein);
	else if (spaceGroup=="B2212") setSpaceGroupB2212(Protein);
	else if (spaceGroup=="C121") setSpaceGroupC121(Protein);
	else if (spaceGroup=="C1211") setSpaceGroupC1211(Protein);
	else if (spaceGroup=="C21") setSpaceGroupC21(Protein);
	else if (spaceGroup=="C222") setSpaceGroupC222(Protein);
	else if (spaceGroup=="C2221") setSpaceGroupC2221(Protein);
	//else if (spaceGroup=="C4212") setSpaceGroupC4212(Protein); //This space group makes no sense
	else if (spaceGroup=="F4132") setSpaceGroupF4132(Protein);
	//else if (spaceGroup=="F422") setSpaceGroupF422(Protein); //Inconsistant space group
	else if (spaceGroup=="F222") setSpaceGroupF222(Protein);
	else if (spaceGroup=="F23") setSpaceGroupF23(Protein);
	else if (spaceGroup=="F432") setSpaceGroupF432(Protein);
	else if (spaceGroup=="H3") setSpaceGroupH3(Protein);
	else if (spaceGroup=="H32") setSpaceGroupH32(Protein);
	else if (spaceGroup=="I121") setSpaceGroupI121(Protein);
	else if (spaceGroup=="I1211") setSpaceGroupI1211(Protein);
	else if (spaceGroup=="I21") setSpaceGroupI21(Protein);
	else if (spaceGroup=="I212121") setSpaceGroupI212121(Protein);
	else if (spaceGroup=="I213") setSpaceGroupI213(Protein);
	else if (spaceGroup=="I222") setSpaceGroupI222(Protein);
	else if (spaceGroup=="I23") setSpaceGroupI23(Protein);
	else if (spaceGroup=="I-42D") setSpaceGroupI_42D(Protein);
	else if (spaceGroup=="I-4C2") setSpaceGroupI_4C2(Protein);
	else if (spaceGroup=="I4") setSpaceGroupI4(Protein);
	else if (spaceGroup=="I41") setSpaceGroupI41(Protein);
	else if (spaceGroup=="I41/A") setSpaceGroupI41_A(Protein);
	else if (spaceGroup=="I411") setSpaceGroupI411(Protein);
	else if (spaceGroup=="I4122") setSpaceGroupI4122(Protein);
	else if (spaceGroup=="I4132") setSpaceGroupI4132(Protein);
	else if (spaceGroup=="I422") setSpaceGroupI422(Protein);
	else if (spaceGroup=="I432") setSpaceGroupI432(Protein);
	else if (spaceGroup=="P1") setSpaceGroupP1(Protein);
	else if (spaceGroup=="P-1") setSpaceGroupP_1(Protein);
	else if (spaceGroup=="P1121") setSpaceGroupP1121(Protein);
	else if (spaceGroup=="P121") setSpaceGroupP121(Protein);
	else if (spaceGroup=="P121/C1") setSpaceGroupP121_C1(Protein);
	else if (spaceGroup=="P121/N1") setSpaceGroupP121_N1(Protein);
	else if (spaceGroup=="P1211") setSpaceGroupP1211(Protein);
	else if (spaceGroup=="P21212") setSpaceGroupP21212(Protein);
	else if (spaceGroup=="P212121") setSpaceGroupP212121(Protein);
	else if (spaceGroup=="P21212A") setSpaceGroupP21212A(Protein);
	else if (spaceGroup=="P21221") setSpaceGroupP21221(Protein);
	else if (spaceGroup=="P21312") setSpaceGroupP21312(Protein);
	else if (spaceGroup=="P213") setSpaceGroupP213(Protein);
	else if (spaceGroup=="P222") setSpaceGroupP222(Protein);
	else if (spaceGroup=="P2221") setSpaceGroupP2221(Protein);
	else if (spaceGroup=="P3") setSpaceGroupP3(Protein);
	else if (spaceGroup=="P31") setSpaceGroupP31(Protein);
	else if (spaceGroup=="P3112") setSpaceGroupP3112(Protein);
	else if (spaceGroup=="P3121") setSpaceGroupP3121(Protein);
	else if (spaceGroup=="P312") setSpaceGroupP312(Protein);
	else if (spaceGroup=="P32") setSpaceGroupP32(Protein);
	else if (spaceGroup=="P321") setSpaceGroupP321(Protein);
	else if (spaceGroup=="P3212") setSpaceGroupP3212(Protein);
	else if (spaceGroup=="P3221") setSpaceGroupP3221(Protein);
	else if (spaceGroup=="P4") setSpaceGroupP4(Protein);
	else if (spaceGroup=="P41") setSpaceGroupP41(Protein);
	else if (spaceGroup=="P41212") setSpaceGroupP41212(Protein);
	else if (spaceGroup=="P4122") setSpaceGroupP4122(Protein);
	else if (spaceGroup=="P4132") setSpaceGroupP4132(Protein);
	else if (spaceGroup=="P42") setSpaceGroupP42(Protein);
	else if (spaceGroup=="P4212") setSpaceGroupP4212(Protein);
	else if (spaceGroup=="P422") setSpaceGroupP422(Protein);
	else if (spaceGroup=="P22121") setSpaceGroupP22121(Protein);
	else if (spaceGroup=="P42212") setSpaceGroupP42212(Protein);
	else if (spaceGroup=="P4222") setSpaceGroupP4222(Protein);
	else if (spaceGroup=="P43") setSpaceGroupP43(Protein);
	else if (spaceGroup=="P432") setSpaceGroupP432(Protein);
	else if (spaceGroup=="P4322") setSpaceGroupP4322(Protein);
	else if (spaceGroup=="P4332") setSpaceGroupP4332(Protein);
	else if (spaceGroup=="P43212") setSpaceGroupP43212(Protein);
	else if (spaceGroup=="P6") setSpaceGroupP6(Protein);
	else if (spaceGroup=="P61") setSpaceGroupP61(Protein);
	else if (spaceGroup=="P6122") setSpaceGroupP6122(Protein);
	else if (spaceGroup=="P62") setSpaceGroupP62(Protein);
	else if (spaceGroup=="P622") setSpaceGroupP622(Protein);
	else if (spaceGroup=="P6222") setSpaceGroupP6222(Protein);
	else if (spaceGroup=="P63") setSpaceGroupP63(Protein);
	else if (spaceGroup=="P6322") setSpaceGroupP6322(Protein);
	else if (spaceGroup=="P64") setSpaceGroupP64(Protein);
	else if (spaceGroup=="P6422") setSpaceGroupP6422(Protein);
	else if (spaceGroup=="P65") setSpaceGroupP65(Protein);
	else if (spaceGroup=="P6522") setSpaceGroupP6522(Protein);
	else if (spaceGroup=="R32") setSpaceGroupR32(Protein);
	else if (spaceGroup=="R3") setSpaceGroupR3(Protein);	
	else
	{
		printSpaceGroups();
		error("Unknown space group "+spaceGroup, __LINE__, __FILE__);
	}

	Protein.spaceGroup=spaceGroup;
	Protein.intSpaceGroup=spaceGroupStrToSpaceGroupInt(spaceGroup);
}
