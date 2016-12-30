#ifndef _molecularReplacement_include_
#define _molecularReplacement_include_

# include "../../LibraryFiles/TypeDef.h"
# include "../../LibraryFiles/minimize.h"
# include "../../LibraryFiles/MonteCarlo.h"
# include "XRayStruct.h"
# include "XRay.h"
# include "LatticeStruct.h"
void iterativeElectronDensityModification(XRayStruct &expXRay, XRayStruct &xray, LatticeStruct &lattice, XRayParamStruct &params);
Real quantifyMatch(XRayStruct &xray, XRayStruct &expXRay, string xRayScoreType);
Real calcLogSqrDiff(Vector &calci, Vector &expi);
Real quantifyWeightedMatch(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
Real quantifyMatch(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
Real phase(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
Real phase(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
Real phase(Vector degOfFreedom, MRStruct args);
Real calcFitness(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
Real calcFitness(Vector degOfFreedom, MRStruct args);
Real calcFitness(Vector &degOfFreedom, vector<ProteinStruct> Proteins, XRayStruct &XRay, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
void calcScattering(Vector &degOfFreedom, ProteinStruct Protein, XRayStruct &xray, LatticeStruct &lattice, XRayParamStruct &params);
void initializeXRayProtein(ProteinStruct &Protein, XRayStruct &xray);
Real calcFitnessChange(Vector &degOfFreedom, ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params, int pick);
Real calcFitnessTranslation(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &dPos);
Real calcFitnessTranslation(Vector dPos, MRStruct args);
void degOfFreedomToRotMat(Matrix &rotationMatrix, Real rx, Real ry, Real rz);
void degOfFreedomToRotMat(vector<ProteinStruct> &Proteins, Vector &degOfFreedom);
void degOfFreedomToRotMatAndTranslation(ProteinStruct &Protein, Vector &degOfFreedom);
void rotateProteins(vector<ProteinStruct> Proteins1, Vector &degOfFreedom, vector<ProteinStruct> &Proteins2);
void placeProtein(ProteinStruct &Protein1, Vector &degOfFreedom, ProteinStruct &Protein2);
void placeProtein(ProteinStruct &Protein, Vector &degOfFreedom);
void placeProteins(vector<ProteinStruct> &Proteins, Vector &degOfFreedom);
void placeProteins(Vector &degOfFreedom, ProteinStruct &Protein, vector<ProteinStruct> &Proteins);
void placeProteins(vector<ProteinStruct> &Proteins1, Vector &degOfFreedom, vector<ProteinStruct> &Proteins2);
void placeProteins(vector<ProteinStruct> Proteins, Vector &degOfFreedom, ProteinStruct &Complex, XRayStruct xray);
void placeProteins(vector<ProteinStruct> Proteins, Vector &degOfFreedom, ProteinStruct &Complex);
void makeUnitCell(ProteinStruct &Protein, vector<ProteinStruct> &allProteins);
void makeImages(vector<ProteinStruct> &allProteins);
void makeUnitCellAndImages(ProteinStruct &Protein, vector<ProteinStruct> &allProteins);
void makeUnitCellAndImages(vector<ProteinStruct> &Proteins, vector<ProteinStruct> &allProteins);
void makeUnitCellAndImagesFast(vector<ProteinStruct> &Proteins, vector<ProteinStruct> &allProteins);
void makeUnitCellAndImagesFast2(vector<ProteinStruct> &Proteins, vector<ProteinStruct> &allProteins);
void makeUnitCellAndImagesFast3(vector<ProteinStruct> &Proteins, vector<ProteinStruct> &allProteins);
Real calcRotationTranslationFitness(Vector &degOfFreedom, vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
Real calcRotationTranslationFitness(Vector degOfFreedom, MRStruct args);
Real penalizeOutOfBounds(Vector &degOfFreedom, Real x, Real y, Real z);
Real calcRotationTranslationFitness2(Vector degOfFreedom, vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
Real calcRotationTranslationFitness2(Vector degOfFreedom, MRStruct &args);
Real calcBFactorFitness(Vector bFactorScale, MRStruct args);
Real calcFormFactorFitness(Vector formFactorScale, MRStruct args);
Real MRRotateTranslate(ProteinStruct Protein, Vector degOfFreedom, XRayStruct xray, XRayStruct expXRay, LatticeStruct lattice, XRayParamStruct params);
Real MRRotationalGridSearch(Vector &bestGrid, Real rinc, int npoint, ProteinStruct Protein, XRayStruct xray, XRayStruct &expXRay, LatticeStruct lattice, XRayParamStruct &params, vector<Array3D> &degOfFreedoms, Array3D &scores);
Real MRRotationalGridSearch(Vector &bestGrid, ProteinStruct Protein, XRayStruct xray, XRayStruct &expXRay, LatticeStruct lattice, XRayParamStruct &params, vector<Array3D> &degOfFreedoms, Array3D &scores);
Real minimizeRotation(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom, Real maxHours);
Real minimizeRotation(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Real maxHours);
Real minimizeRotationOnly(ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom, Real maxHours);
bool minimizeRotationUsingEstimatedScore(ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom);
void minimizeRotationFast(ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom);
Real calcAnalyticalRotationalGrad(Vector &degOfFreedom, ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &rotGrad);
Real minimizeRotationAnalytical(vector<Real> &coefficient, MRStruct args, Real dScoreCriteria, Real maxHours);
Real minimizeRotationAnalytical(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
Real minimizeRotationTranslationAnalytical(Vector &degOfFreedom, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Real maxHours);
bool isLocalMinima(Array3D &scores, int i, int j, int k);
void findLocalMinima(vector< Array3D > &degOfFreedoms, Array3D &scores, Matrix &localMinima, Vector &localMinimaScores);
void getBestResults(Matrix &localMinima, Vector &localMinimaScores);
Real MRRotate(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params, Matrix &localMinima, Vector &degOfFreedom);
Real MRTranslationGridSearch(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &dPos);
Real minimizeTranslation(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom);
Real minimizeTranslationOnly(vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom);
Real minimizeRotationTranslation(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom);
Real MRTranslate(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &dPos);
void makeRandomRotationTranslation(ProteinStruct &Protein);
Real MRTranslate(ProteinStruct Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &degOfFreedom, Vector &dPos);
Real MRTranslate(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Matrix &degOfFreedoms, Vector &bestDegOfFreedom, Vector &bestDPos);
Real MRRotateTranslate(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
void printDofs(Matrix &degOfFreedoms);
void initializeDegOfFreedoms(Matrix &degOfFreedoms, string degOfFreedomFile);
void fixDegOfFreedoms(Vector &degOfFreedoms, Real x, Real y, Real z);
void fixDegOfFreedoms(Matrix &degOfFreedoms, Real x, Real y, Real z);
int pickRandom(MoveStruct &move);
void fixSpaceGroupPick(MoveStruct &move, string spaceGroup, int &pick);
int pickRandomDOF(MoveStruct &move, string spaceGroup, int nprot);
int pickRandomDOF(MoveStruct &move, string spaceGroup);
void fixMoveSize(int pick, Real xlength, Real ylength, Real zlength, Real &moveSize, Vector &coefficient);
void fixDOF(Real xlength, Real ylength, Real zlength, Vector &coefficient);
bool checkIfAllZero(Vector &coefficient, int nprot);
void fixDOF(Vector &coefficient, int nprot);
void fixMoveSizeOccupancy(Real &moveSize, Vector &coefficient);
void fixMoveSize(int pick, Real xlength, Real ylength, Real zlength, Real &moveSize, Vector &coefficient, int nprot);
void makeRandomChange(Vector &coefficient, MoveStruct &move, int &pick, Real xlength, Real ylength, Real zlength, string spaceGroup, int nprot);
Real getRx(Matrix &rotMat, Real ry);
Real getRz(Matrix &rotMat, Real ry);
void rotationMatrixToCoefficients(Matrix &rotMat, Real &rx, Real &ry, Real &rz);
void calcNCSTranslation(Vector &degOfFreedomNCS, VectorStruct &v, Matrix &rotMat, Real &dx, Real &dy, Real &dz);
void calcDOFFromNCS(Vector &degOfFreedomNCS, Vector &degOfFreedom);
void calcDOFFromNCS(Matrix &degOfFreedomNCS, Matrix &degOfFreedom);
void calcDOFFromNCS(Matrix &degOfFreedomNCS);
void fixMoveSizeNCS(int pick, Real xlength, Real ylength, Real zlength, Real &moveSize);
void makeRandomChangeNCS(Vector &degOfFreedom, Vector &degOfFreedomNCS, MoveStruct &move, int &pick, Real xlength, Real ylength, Real zlength, string spaceGroup);
Real calcScore(Vector &degOfFreedoms, vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
void calcScores(Matrix &degOfFreedoms, Vector &score, Matrix &scores, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
void calcStandardDeviations(Matrix &scores, Vector &stdDev);
Real monteCarloNCS(Vector &degOfFreedomNCS, vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params, Real temp, int nsteps, MoveStruct &move);
Real monteCarlo(Vector &degOfFreedom, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params, Real temp, int nsteps, MoveStruct &move);
void calcScoresWithClash(Matrix &degOfFreedoms, Vector &score, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
void printReplica(Vector &degOfFreedom, vector<ProteinStruct> Proteins, string outputFile, XRayParamStruct &params);
void printReplicas(Matrix &degOfFreedoms, vector<ProteinStruct> &Proteins, XRayParamStruct &params, string replicaStructuresOutput);
void roundCoordinates(vector<AtomStruct> &Atoms);
void calcScores(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params, Vector &scores, vector<string> &scoreNames);
void printScores(XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
void printScores(Vector &degOfFreedoms, vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
void printScores(Matrix &degOfFreedoms, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
void printXRay(Vector &degOfFreedoms, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayParamStruct &params, string xrayOutFile);
void printXRay(Matrix degOfFreedoms, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayParamStruct &params);
void optimizeAll(Matrix &degOfFreedoms, Vector &score, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
void initializeMove(int ndof, MoveStruct &move, ProteinStruct &Protein, int index, int nrep, XRayParamStruct &params);
void initializeMoveNCS(int ndof, MoveStruct &move, ProteinStruct &Protein, int index, int nrep, XRayParamStruct &params);
void initializeMoves(int nrep, int ndof, vector<MoveStruct> &moves, ProteinStruct &Protein, XRayParamStruct &params);
bool checkForGap(Vector &scores, Real stdDevBest);
Real calcRMSD(vector<ProteinStruct> &Proteins1, Vector &degOfFreedom, vector<ProteinStruct> Proteins2);
bool checkSimilarity(Matrix &degOfFreedom, vector<ProteinStruct> Proteins);
bool checkForConvergence(Vector &scores, Matrix &degOfFreedom, vector<ProteinStruct> &Proteins);
void calcRMSDMatrix(Matrix &degOfFreedoms, ProteinStruct &Protein, Matrix &rmsd);
void calcScatteringForAllModels(Vector &degOfFreedom, XRayStruct &xray, LatticeStruct &lattice, vector<ProteinStruct> &Proteins, XRayParamStruct &params, Matrix &real, Matrix &imag);
void calcScatteringForAllModels(Matrix &degOfFreedoms, XRayStruct &xray, LatticeStruct &lattice, vector<ProteinStruct> &Proteins, XRayParamStruct &params, Matrix &real, Matrix &imag);
Real getScoreForPermutation(Matrix &real, Matrix &imag, XRayStruct &xray, XRayStruct &expXRay, vector<int> &permutation, Matrix &rmsds, XRayParamStruct &params);
void initializePermutation(vector<int> &permutation, int nmodel, int nprot);
void calcScoresForPermutations(Matrix &real, Matrix &imag, XRayStruct &xray, XRayStruct expXRay, int nprot, Vector &scores, vector< vector<int> > &permutations, Matrix &rmsds, int nuse, XRayParamStruct &params);
void getNewDegOfFreedoms(Vector &scores, vector< vector<int> > &permutations, Matrix &degOfFreedoms);
void calcPermutationsOfProteins(Matrix &degOfFreedoms, XRayStruct xray, XRayStruct &expXRay, LatticeStruct &lattice, vector<ProteinStruct> &Proteins, XRayParamStruct &params);
void shuffleReplicas(Matrix &degOfFreedoms);
void printContinuationFile(string continuationFile, Matrix &degOfFreedoms);
void outputDegreesOfFreedom(Matrix &degOfFreedoms, int cycle, XRayParamStruct &params);
void outputDegreesOfFreedomHeader(XRayParamStruct &params);
Real calcOrientation(ProteinStruct &Protein1, ProteinStruct &Protein2);
void roundOccupancies(vector<AtomStruct> &Atoms);
void outputTrajectorySnapShot(Vector &degOfFreedom, vector<ProteinStruct> &Proteins, ProteinStruct &Native, XRayStruct &xray, XRayStruct &expXRay, string name, ofstream &file, XRayParamStruct &params);
void outputTrajectorySnapShot(Matrix &degOfFreedoms, int cycle, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
void printCorrectConformation(vector<ProteinStruct> Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
void printTrajectoryHeader(vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
void deleteMoves(vector<MoveStruct> &moves);
void cleanMemory(vector<MoveStruct> &moves);
Real replicaExchange(Matrix &degOfFreedoms, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params, Real minTemp, Real maxTemp, int nrep, int nsteps, int ncycles);
void initializeDegOfFreedomNoNCS(Matrix &degOfFreedoms, ProteinStruct &Protein, XRayParamStruct &params);
void initializeDegOfFreedomNCS(Matrix &degOfFreedoms, ProteinStruct &Protein, XRayParamStruct &params);
void initializeDegOfFreedom(Matrix &degOfFreedoms, ProteinStruct &Protein, XRayParamStruct &params);
void addOccupanciesToDegOfFreedom(Matrix &degOfFreedoms, int NumOccupancySegments);
void turnOnOptimizeOccupancies(Matrix &degOfFreedoms, vector<ProteinStruct> &Proteins, XRayStruct &xray, XRayParamStruct &params, XRayParamStruct &copyParams);
Real MRReplicaExchange(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
Real optimizeBFactorScale(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
Real optimizeFormFactorScale(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, XRayParamStruct &params);
void molecularReplacement(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay, LatticeStruct &lattice, XRayParamStruct &params);
void molecularReplacement(string inputPdbFile, XRayParamStruct &params);
#endif
