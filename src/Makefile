CC=g++
#CFLAGS=-c -O5 -march=native -mtune=native -ffast-math -Wall
CFLAGS=-c -g -Wall
LDFLAGS=
lib=/home/joukov/project/LibraryFiles
mr=/home/joukov/project/mr/src
SOURCES=$(mr)/FlexibleEditing.cpp $(lib)/test_functions.cpp $(mr)/CalcClashScore.cpp $(mr)/Params.cpp $(lib)/ReadPrm.cpp $(lib)/ReadRtf.cpp $(lib)/TMScore.cpp $(lib)/rmsd.cpp $(mr)/RMSD.cpp $(lib)/ReadDcd.cpp $(lib)/MonteCarlo.cpp $(lib)/ASA.cpp $(mr)/XRay.cpp $(lib)/ReadGro.cpp $(lib)/normalize.cpp $(lib)/FormFactors.cpp $(lib)/IOUtils.cpp $(lib)/MathUtils.cpp $(lib)/AtomUtils.cpp $(mr)/LatticeUtils.cpp $(lib)/time.cpp $(lib)/PrintPdb.cpp $(lib)/ReadPdb.cpp $(lib)/GetAtomType.cpp $(lib)/GetAtomType2.cpp $(lib)/ResidueID.cpp $(lib)/AssignAtomIDs.cpp $(mr)/Initialize.cpp $(mr)/ReadExperimentalDataFile.cpp $(mr)/ReadParameters.cpp $(mr)/molecularReplacement.cpp $(lib)/Structures.cpp $(lib)/ReadPsf.cpp $(lib)/MemoryUsage.cpp $(lib)/error.cpp $(lib)/GetNextFrame.cpp $(mr)/SpaceGroup.cpp $(mr)/MaximumLikelihood.cpp $(mr)/mr.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mr_test

SOURCES2=$(mr)/EmpiricalMaximumLikelihood.cpp $(lib)/test_functions.cpp $(mr)/FlexibleEditing.cpp $(mr)/CalcClashScore.cpp $(mr)/Params.cpp $(lib)/ReadPrm.cpp $(lib)/ReadRtf.cpp $(lib)/TMScore.cpp $(lib)/rmsd.cpp $(mr)/RMSD.cpp $(lib)/ReadDcd.cpp $(lib)/ASA.cpp $(mr)/XRay.cpp $(lib)/ReadGro.cpp $(lib)/normalize.cpp $(lib)/FormFactors.cpp $(lib)/IOUtils.cpp $(lib)/MathUtils.cpp $(lib)/AtomUtils.cpp $(mr)/LatticeUtils.cpp $(lib)/time.cpp $(lib)/PrintPdb.cpp $(lib)/ReadPdb.cpp $(lib)/GetAtomType.cpp $(lib)/GetAtomType2.cpp $(lib)/ResidueID.cpp $(lib)/AssignAtomIDs.cpp $(mr)/Initialize.cpp $(mr)/ReadExperimentalDataFile.cpp $(mr)/ReadParameters.cpp $(mr)/molecularReplacement.cpp $(lib)/Structures.cpp $(lib)/ReadPsf.cpp $(lib)/MemoryUsage.cpp $(lib)/error.cpp $(lib)/GetNextFrame.cpp $(mr)/SpaceGroup.cpp $(mr)/MaximumLikelihood.cpp $(lib)/MonteCarlo.cpp $(mr)/CalcXRay.cpp
OBJECTS2=$(SOURCES2:.cpp=.o)
EXECUTABLE2=CalcXRay

mr: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

CalcXRay: $(SOURCES2) $(EXECUTABLE2)

$(EXECUTABLE2): $(OBJECTS2)
	$(CC) $(LDFLAGS) $(OBJECTS2) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
