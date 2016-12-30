# include <iostream>
# include <fstream>
# include <string>
# include <cstring>
# include <sstream>
# include <time.h>
# include <vector>
# include <iomanip>
# include <complex>
# include <numeric>

using namespace std;

# include "../../LibraryFiles/error.h"
# include "../../LibraryFiles/Constants.h"
# include "Params.h"

void setDefaultParams(XRayParamStruct &params)
{
	params.AsaCutoff=0.2;
	params.a=UNK_REAL;
	params.b=UNK_REAL;
	params.c=UNK_REAL;
	params.alpha=UNK_REAL;
	params.beta=UNK_REAL;
	params.gamma=UNK_REAL;
	params.BulkSolventCorrection=false;
	params.CalcRotationPossibleScore=true;
	params.EndIfConverged=false;
	params.FilterOutHighResolution=false;
	params.MinimizeAfterEachStep=false;
	params.MinimizeEverySolution=false;
	params.OptimizeBFactorScale=true;
	params.OptimizeFormFactorScale=true;
	params.OptimizeOccupancies=false;
	params.OptimizeSegmentOccupanciesAfterMR=false;
	params.Phase=false;
	params.PrintOutput=false;
	params.PrintTrajectory=false;
	params.RemoveProtein=false;
	params.RemoveUnknownAtoms=false;
	params.RemoveWaters=false;
	params.ReplicaExchange=false;
	params.RigidBodyOptimization=false;
	params.RmsdDependantELL=false;
	params.RotateTranslate=false;
	params.SegmentsFromFile=false;
	params.SetTempsBasedOnScores=true;
	params.ShuffleDegreesOfFreedom=false;
	params.UseDynamicTemperature=false;
	params.UseExpLookUp=false;
	params.UseFastTranslationRotation=false;
	params.UseFastTranslation=false;
	params.UseMillerIndexesFromExperiment=false;
	params.UseNCS=false;
	params.UseVeryFastTranslation=true;
	params.Verbose=false;
	params.CenteredAndRotatedOutputPdb="";
	params.ClashScoreUse="AtEnd";
	params.ContinuationFile="";
	params.CorrectionFile="";
	params.CubeFile="";
	params.DcdFilePaths="";
	params.DegOfFreedomFile="";
	params.DegOfFreedomTrajectoryFile="";
	params.DensityOutput="";
	params.ErrorHistogramFile="";
	params.ErrorPathsFile="";
	params.ExcludedVolumeMethod="Spheres";
	params.FastTranslationRotationInterpolation="Linear";
	params.LogFiles="";
	params.MinimizeRotationTranslation="Analytical";
	params.NativePdbFile="";
	params.OccupancyCorrectionFile="";
	params.OccupanciesPdbDir="";
	params.PdbFile="";
	params.PdbPathsFile="";
	params.ProbabilitiesFile="";
	params.ScoreFileBase="";
	params.SegmentsFile="";
	params.SpaceGroup="";
	params.TrajectoryOutputFile="";
	params.TranslationRotationMethod="RFactor";
	params.XRayInput="";
	params.XRayOutput="";
	params.XRayScoreType="DiffOfSquares";
	params.ContinuousHValues=100;
	params.ContinuousKValues=100;
	params.ContinuousLValues=100;
	params.MaxTrajectoryStructures=1000;
	params.NumCopies=1;
	params.NumCycles=30000;
	params.NumOccupancySegments=10;
	params.NumQxValues=16;
	params.NumQyValues=16;
	params.NumQzValues=16;
	params.NumReplicas=300;
	params.MRIterations=1;
	params.PrintTrajectoryInterval=100;
	params.PrintTrajectoryMax=100;
	params.RandomSeed=0;
	params.ReplicaStructuresOutput="";
	params.ReplicaXRayOutput="";
	params.RotationGridSearchPoints=31;
	params.ShuffleInterval=1000;
	params.StepsPerCycle=20;
	params.TranslationGridSearchPoints=31;
	params.XCubes=100;
	params.YCubes=100;
	params.ZCubes=100;
	params.BFactorScale=1.0;
	params.BulkDensity=0.334;
	params.CalphaRadius=4.0;
	params.ClashWeightCoreCore=1000.0;
	params.ClashWeightCoreSurface=10.0;
	params.ClashWeightSurfaceSurface=1.0;
	params.CubeSize=0.5;
	params.DiffOfSquaresWeight=0.1;
	params.ExcludedVolumeRadii="Default";
	params.ExcludedVolumeRadiusHydrogen=0.902255;
	params.ExcludedVolumeRadiusCarbon=1.43693;
	params.ExcludedVolumeRadiusNitrogen=1.24527;
	params.ExcludedVolumeRadiusOxygen=1.22099;
	params.ExcludedVolumeRadiusSulfur=2.19596;
	params.scale=1.0;
	params.LogSquareDiffWeight=0;
	params.MaxAngleChange=0.12;
	params.MaxReplicaExchangeTime=4.0;
	params.MaxReplicaExchangeClashTime=0.5;
	params.MaxReplicaExchangeOptimizeOccupancyTime=0.0;
	params.MaxTemp=0.3;
	params.MaxTempFrac=0.25;
	params.MaxTranslationFraction=0.2;
	params.MaxQx=3.0;
	params.MaxQy=3.0;
	params.MaxQz=3.0;
	params.MinAngleChange=0.05;
	params.MinTemp=0.02;
	params.MinTempFrac=0.005;
	params.MinTranslationFraction=0.1;
	params.MLRF0Weight=0;
	params.OccupancyWeight=0;
	params.PearsonWeight=1.0;
	params.ProductWeight=0;
	params.RFactorWeight=1.0;
	params.RotationWeight=0;
	params.XCubeLength=0.5;
	params.YCubeLength=0.5;
	params.ZCubeLength=0.5;
	cout <<"Leaving SetDefaultPar"<<endl;
}

