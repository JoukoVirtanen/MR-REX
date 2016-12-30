#ifndef _LatticeStruct_included_
#define _LatticeStruct_included_

# include "../../LibraryFiles/Structures.h"

struct LatticeStruct
{
	Real cubeSize, a, b, c, alpha, beta, gamma;
	Real xCubeLength, yCubeLength, zCubeLength;
	PosStruct firstPos;
	Array3D density, real, imag;	

	LatticeStruct& operator = (const LatticeStruct &lattice);
};
#endif
