# include <vector>
# include <iostream>
# include <algorithm>

using namespace std;

# include "/home/jouko/project/LibraryFiles/Structures.h"
# include "/home/jouko/project/LibraryFiles/MathUtils.h"
# include "/home/jouko/project/LibraryFiles/ReadPdb.h"
# include "/home/jouko/project/LibraryFiles/AtomUtils.h"

# include "Params.h"
# include "XRay.h"
# include "LatticeStruct.h"
# include "Initialize.h"
# include "molecularReplacement.h"
# include "CalcDensity.h"
# include "ReadExperimentalDataFile.h"
# include "ReadParameters.h"

int main(int argc, char *argv[])
{
	string xrayDataFile, paramsFile;
	XRayParamStruct params;
	time_t seconds1, seconds2;

	cout <<"STARTING"<<endl;
	if (argc!=3)
	{
		cout <<"Usage: "<<argv[0]<<" inputPdbFile paramsFile"<<endl;
		exit(EXIT_FAILURE);
	}

	xrayDataFile=argv[1];
	paramsFile=argv[2];

	SetDefaultParams(params);
	ReadParameterFile(paramsFile, params);

	seconds1=time(NULL);
	cout <<"About to enter calcScattering"<<endl;
	calcDensity(xrayDataFile, params);
	seconds2=time(NULL);
	cout <<"Run time: "<<seconds2-seconds1<<" seconds"<<endl;
	cout <<"DONE"<<endl;

	return EXIT_SUCCESS;
}
