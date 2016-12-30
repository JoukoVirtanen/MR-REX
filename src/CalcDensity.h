#ifndef _CalcDensity_included_
#define _CalcDensity_included_

void calcDensity(string xrayDataFile, XRayParamStruct &params)
{
	LatticeStruct lattice;
	XRayStruct expXRay;

	cout <<"About to enter initializeLattice"<<endl;
	initializeLattice(lattice, params);
	cout <<"About to enter readExperimentalDataFile"<<endl;
	readExperimentalDataFile(xrayDataFile, expXRay);
	millerToQ(expXRay, lattice);
	cout <<"About to calcDensity"<<endl;
	calcDensity(expXRay, lattice);
	cout <<"About to printDensity"<<endl;
	printDensity(params.DensityOutput, lattice);
}

#endif
