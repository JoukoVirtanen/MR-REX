# include <vector>
# include <iostream>
# include <algorithm>

using namespace std;

# include "../../LibraryFiles/Structures.h"
# include "../../LibraryFiles/MathUtils.h"
# include "../../LibraryFiles/ReadPdb.h"
# include "../../LibraryFiles/AtomUtils.h"

# include "Params.h"
# include "XRay.h"
# include "LatticeStruct.h"
# include "Initialize.h"
# include "molecularReplacement.h"
# include "ReadExperimentalDataFile.h"
# include "ReadParameters.h"

int main(int argc, char *argv[])
{
	string inputPdbFile, paramsFile;
	XRayParamStruct params;

	time_t seconds1, seconds2;

	if (argc!=3)
	{
		cout <<"Usage: "<<argv[0]<<" inputPdbFile paramsFile"<<endl;
		exit(EXIT_FAILURE);
	}

	inputPdbFile=argv[1];
	paramsFile=argv[2];

	cout <<"Before SetDefaultParams"<<endl;
	setDefaultParams(params);
	cout <<"Before ReadParameterFile"<<endl;
	ReadParameterFile(paramsFile, params);

	cout <<"Before molecularReplacement"<<endl;
	seconds1=time(NULL);
	molecularReplacement(inputPdbFile, params);
	seconds2=time(NULL);
	cout <<"Run time: "<<seconds2-seconds1<<" seconds"<<endl;
	cout <<"DONE"<<endl;

	return EXIT_SUCCESS;
}
