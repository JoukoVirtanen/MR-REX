#ifndef _LatticeUtils_include_
#define _LatticeUtils_include_

#include "LatticeStruct.h"
#include "Params.h"

void indexToPosition(LatticeStruct &lattice, int xbin, int ybin, int zbin, Real &x, Real &y, Real &z);
void setExcludedVolumeDensity(LatticeStruct &lattice, ProteinStruct &Protein, Real bulkDensity);
void initializeLattice(LatticeStruct &lattice, int XCubes, int YCubes, int ZCubes);
void initializeLatticeForExcludedVolume(LatticeStruct &lattice, ProteinStruct &Protein, XRayParamStruct &params);
void latticeToArray(LatticeStruct &lattice, Real ***&RealDensity, Real ***&ImagDensity);
void fixBin(int maxBin, int &bin);
void fixBins(int maxXBin, int maxYBin, int maxZBin, int &xbin, int &ybin, int &zbin);
void positionToIndex(Real xmin, Real ymin, Real zmin, Real cubeSize, Real x, Real y, Real z, int &xbin, int &ybin, int &zbin);
void positionToIndex(Real xmin, Real ymin, Real zmin, Real xCubeLength, Real yCubeLength, Real zCubeLength, Real x, Real y, Real z, int &xbin, int &ybin, int &zbin);
void positionToIndex(LatticeStruct &lattice, Real x, Real y, Real z, int &xbin, int &ybin, int &zbin);
Real fixFractionalCoordinate(Real &x);
void fixFractionalCoordinates(Real &x, Real &y, Real &z);
void fractionalPositionToIndex(LatticeStruct &lattice, Real x, Real y, Real z, int &xbin, int &ybin, int &zbin);
void scaleLatticeDensity(Array3D &density, Real scale);
Real calcLatticeElectrons(LatticeStruct &lattice);
void normalizeDensity(ProteinStruct &Protein, LatticeStruct &lattice);
Real calcRMSD(LatticeStruct &lattice1, LatticeStruct &lattice2);
Real calcDensityCorrelation(LatticeStruct &lattice1, LatticeStruct &lattice2);
void readPdb(string Pdb, LatticeStruct &lattice);
void printDensity(string outputDensityFile, LatticeStruct &lattice);
void printDensityFrac(string outputDensityFile, LatticeStruct &lattice, ProteinStruct &Protein);
#endif
