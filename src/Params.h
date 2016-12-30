#ifndef _XRayParams_included_
#define _XRayParams_included_

using namespace std;

#include "../../LibraryFiles/TypeDef.h"

#if _MSC_VER > 1000
#pragma once
#endif

struct XRayParamStruct
{
	bool BulkSolventCorrection, CalcRotationPossibleScore;
	bool EndIfConverged, FilterOutHighResolution, MinimizeAfterEachStep;
	bool MinimizeEverySolution, OptimizeBFactorScale;
	bool OptimizeFormFactorScale, OptimizeOccupancies;
	bool OptimizeSegmentOccupanciesAfterMR, Phase, PrintOutput;
	bool PrintTrajectory, RemoveProtein, RemoveUnknownAtoms, RemoveWaters;
	bool ReplicaExchange, RigidBodyOptimization, RmsdDependantELL;
	bool RotateTranslate, SegmentsFromFile, SetTempsBasedOnScores;
	bool ShuffleDegreesOfFreedom, UseFastTranslationRotation;
	bool UseDynamicTemperature, UseFastTranslation, UseVeryFastTranslation;
	bool UseExpLookUp, UseMillerIndexesFromExperiment, UseNCS, Verbose;
	string CenteredAndRotatedOutputPdb, ClashScoreUse, CorrectionFile;
	string CubeFile, ContinuationFile, DcdFilePaths, DegOfFreedomFile;
	string DegOfFreedomTrajectoryFile, DensityOutput, ExcludedVolumeMethod;
	string ExcludedVolumeRadii;
	string ErrorHistogramFile, ErrorPathsFile;
	string FastTranslationRotationInterpolation;
	string LogFiles, MinimizeRotationTranslation, NativePdbFile;
	string OccupancyCorrectionFile, OccupanciesPdbDir, PdbFile, PdbFiles;
	string PdbPathsFile, ProbabilitiesFile, ReplicaStructuresOutput;
	string ReplicaXRayOutput, ScoreFileBase, SegmentsFile, SpaceGroup;
	string TrajectoryOutputFile, TranslationRotationMethod, XRayInput;
	string XRayOutput, XRayScoreType;
	int ContinuousHValues, ContinuousKValues, ContinuousLValues;
	int MaxTrajectoryStructures, NumCopies, NumCycles, NumOccupancySegments;
	int NumQxValues, NumQyValues, NumQzValues, NumReplicas, MRIterations;
	int PrintTrajectoryInterval, PrintTrajectoryMax, RandomSeed;
	int RotationGridSearchPoints, ShuffleInterval, StepsPerCycle;
	int TranslationGridSearchPoints, XCubes, YCubes, ZCubes;
	Real AsaCutoff, a, b, c, alpha, beta, gamma, BFactorScale, BulkDensity;
	Real CalphaRadius, ClashWeightCoreCore, ClashWeightCoreSurface;
	Real ClashWeightSurfaceSurface, CubeSize, DiffOfSquaresWeight;
	Real ExcludedVolumeRadiusHydrogen, ExcludedVolumeRadiusCarbon;
	Real ExcludedVolumeRadiusNitrogen, ExcludedVolumeRadiusOxygen;
	Real ExcludedVolumeRadiusSulfur, LogSquareDiffWeight, MaxTemp;
	Real MaxTempFrac, MinTemp, MaxReplicaExchangeTime;
	Real MaxReplicaExchangeClashTime;
	Real MaxReplicaExchangeOptimizeOccupancyTime, MaxAngleChange;
	Real MinAngleChange, MinTranslationFraction, MaxTranslationFraction;
	Real MinTempFrac, MLRF0Weight, MaxQx, MaxQy, MaxQz, OccupancyWeight;
	Real PearsonWeight, ProductWeight, RFactorWeight, RotationWeight, scale;
	Real XCubeLength, YCubeLength, ZCubeLength;
};

void setDefaultParams(XRayParamStruct &params);

#endif
