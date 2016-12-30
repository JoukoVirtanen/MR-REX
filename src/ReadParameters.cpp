# include <vector>
# include <iostream>
# include <algorithm>

using namespace std;

# include "../../LibraryFiles/IOUtils.h"
# include "../../LibraryFiles/StringUtils.h"
# include "Params.h"

void ReadParameterFile(string ParameterFile, XRayParamStruct &params)
{
	int NumLines;
	string par, parameter;
	vector<string> lines, str;

	readLines(ParameterFile, lines, "parameter file");
	ProcessLinesForReadingParameters(lines);
	NumLines=lines.size();

	for (int i=0;i<NumLines;i++)
	{
		Tokenize(lines[i], " ", str);
		par=str[0];
		parameter=str[1];
		if (par=="AsaCutoff") params.AsaCutoff=toReal(parameter);
		else if (par=="a") params.a=StrToFloat(parameter);
		else if (par=="b") params.b=StrToFloat(parameter);
		else if (par=="c") params.c=StrToFloat(parameter);
		else if (par=="alpha") params.alpha=StrToFloat(parameter);
		else if (par=="beta") params.beta=StrToFloat(parameter);
		else if (par=="gamma") params.gamma=StrToFloat(parameter);
		else if (par=="BFactorScale") params.BFactorScale=StrToFloat(parameter);
		else if (par=="BulkDensity") params.BulkDensity=StrToFloat(parameter);
		else if (par=="BulkSolventCorrection") params.BulkSolventCorrection=StrToBool(parameter);
		else if (par=="CalcRotationPossibleScore") params.CalcRotationPossibleScore=StrToBool(parameter);
		else if (par=="CalphaRadius") params.CalphaRadius=toReal(parameter);
		else if (par=="CenteredAndRotatedOutputPdb") params.CenteredAndRotatedOutputPdb=parameter;
		else if (par=="ClashScoreUse") params.ClashScoreUse=parameter;
		else if (par=="ClashWeightCoreCore") params.ClashWeightCoreCore=toReal(parameter);
		else if (par=="ClashWeightCoreSurface") params.ClashWeightCoreSurface=toReal(parameter);
		else if (par=="ClashWeightSurfaceSurface") params.ClashWeightSurfaceSurface=toReal(parameter);
		else if (par=="ContinuationFile") params.ContinuationFile=parameter;
		else if (par=="ContinuousHValues") params.ContinuousHValues=StrToInt(parameter);
		else if (par=="ContinuousKValues") params.ContinuousKValues=StrToInt(parameter);
		else if (par=="ContinuousLValues") params.ContinuousLValues=StrToInt(parameter);
		else if (par=="CorrectionFile") params.CorrectionFile=parameter;
		else if (par=="CubeFile") params.CubeFile=parameter;
		else if (par=="CubeSize") params.CubeSize=StrToFloat(parameter);
		else if (par=="DcdFilePaths") params.DcdFilePaths=parameter;
		else if (par=="DiffOfSquaresWeight") params.DiffOfSquaresWeight=toReal(parameter);
		else if (par=="DegOfFreedomFile") params.DegOfFreedomFile=parameter;
		else if (par=="DegOfFreedomTrajectoryFile") params.DegOfFreedomTrajectoryFile=parameter;
		else if (par=="DensityOutput") params.DensityOutput=parameter;
		else if (par=="EndIfConverged") params.EndIfConverged=StrToBool(parameter);
		else if (par=="ErrorHistogramFile") params.ErrorHistogramFile=parameter;
		else if (par=="ErrorPathsFile") params.ErrorPathsFile=parameter;
		else if (par=="ExcludedVolumeMethod") params.ExcludedVolumeMethod=parameter;
		else if (par=="ExcludedVolumeRadii") params.ExcludedVolumeRadii=parameter;
		else if (par=="ExcludedVolumeRadiusCarbon") params.ExcludedVolumeRadiusCarbon=StrToFloat(parameter);
		else if (par=="ExcludedVolumeRadiusHydrogen") params.ExcludedVolumeRadiusHydrogen=StrToFloat(parameter);
		else if (par=="ExcludedVolumeRadiusNitrogen") params.ExcludedVolumeRadiusNitrogen=StrToFloat(parameter);
		else if (par=="ExcludedVolumeRadiusOxygen") params.ExcludedVolumeRadiusOxygen=StrToFloat(parameter);
		else if (par=="ExcludedVolumeRadiusSulfur") params.ExcludedVolumeRadiusSulfur=StrToFloat(parameter);
		else if (par=="FastTranslationRotationInterpolation") params.FastTranslationRotationInterpolation=parameter;
		else if (par=="FilterOutHighResolution") params.FilterOutHighResolution=StrToBool(parameter);
		else if (par=="LogFiles") params.LogFiles=parameter;
		else if (par=="LogSquareDiffWeight") params.LogSquareDiffWeight=toReal(parameter);
		else if (par=="MinimizeAfterEachStep") params.MinimizeAfterEachStep=StrToBool(parameter);
		else if (par=="MinimizeEverySolution") params.MinimizeEverySolution=StrToBool(parameter);
		else if (par=="MinimizeRotationTranslation") params.MinimizeRotationTranslation=parameter;
		else if (par=="MRIterations") params.MRIterations=StrToInt(parameter);
		else if (par=="MaxAngleChange") params.MaxAngleChange=toReal(parameter);
		else if (par=="MaxReplicaExchangeTime") params.MaxReplicaExchangeTime=toReal(parameter);
		else if (par=="MaxReplicaExchangeClashTime") params.MaxReplicaExchangeClashTime=toReal(parameter);
		else if (par=="MaxReplicaExchangeOptimizeOccupancyTime") params.MaxReplicaExchangeOptimizeOccupancyTime=toReal(parameter);
		else if (par=="MaxQx") params.MaxQx=StrToFloat(parameter);
		else if (par=="MaxQy") params.MaxQy=StrToFloat(parameter);
		else if (par=="MaxQz") params.MaxQz=StrToFloat(parameter);
		else if (par=="MaxTemp") params.MaxTemp=StrToFloat(parameter);
		else if (par=="MaxTempFrac") params.MaxTempFrac=toReal(parameter);
		else if (par=="MaxTranslationFraction") params.MaxTranslationFraction=toReal(parameter);
		else if (par=="MinAngleChange") params.MinAngleChange=toReal(parameter);
		else if (par=="MinTemp") params.MinTemp=StrToFloat(parameter);
		else if (par=="MinTempFrac") params.MinTempFrac=toReal(parameter);
		else if (par=="MinTranslationFraction") params.MinTranslationFraction=toReal(parameter);
		else if (par=="MaxTrajectoryStructures") params.MaxTrajectoryStructures=toInt(parameter);
		else if (par=="NativePdbFile") params.NativePdbFile=parameter;
		else if (par=="NumCopies") params.NumCopies=StrToInt(parameter);
		else if (par=="NumCycles") params.NumCycles=StrToInt(parameter);
		else if (par=="NumOccupancySegments") params.NumOccupancySegments=StrToInt(parameter);
		else if (par=="NumQxValues") params.NumQxValues=StrToInt(parameter);
		else if (par=="NumQyValues") params.NumQyValues=StrToInt(parameter);
		else if (par=="NumQzValues") params.NumQzValues=StrToInt(parameter);
		else if (par=="NumReplicas") params.NumReplicas=StrToInt(parameter);
		else if (par=="OccupanciesPdbDir") params.OccupanciesPdbDir=parameter;
		else if (par=="OccupancyCorrectionFile") params.OccupancyCorrectionFile=parameter;
		else if (par=="OccupancyWeight") params.OccupancyWeight=toReal(parameter);
		else if (par=="OptimizeBFactorScale") params.OptimizeBFactorScale=StrToBool(parameter);
		else if (par=="OptimizeFormFactorScale") params.OptimizeFormFactorScale=StrToBool(parameter);
		else if (par=="OptimizeOccupancies") params.OptimizeOccupancies=StrToBool(parameter);
		else if (par=="OptimizeSegmentOccupanciesAfterMR") params.OptimizeSegmentOccupanciesAfterMR=StrToBool(parameter);
		else if (par=="PdbFile") params.PdbFile=parameter;
		else if (par=="PdbPathsFile") params.PdbPathsFile=parameter;
		else if (par=="PearsonWeight") params.PearsonWeight=toReal(parameter);
		else if (par=="Phase") params.Phase=StrToBool(parameter);
		else if (par=="PrintOutput") params.PrintOutput=StrToBool(parameter);
		else if (par=="PrintTrajectory") params.PrintTrajectory=StrToBool(parameter);
		else if (par=="PrintTrajectoryInterval") params.PrintTrajectoryInterval=toInt(parameter);
		else if (par=="PrintTrajectoryMax") params.PrintTrajectoryMax=toInt(parameter);
		else if (par=="ProbabilitiesFile") params.ProbabilitiesFile=parameter;
		else if (par=="ProductWeight") params.ProductWeight=toReal(parameter);
		else if (par=="MLRF0Weight") params.MLRF0Weight=toReal(parameter);
		else if (par=="RandomSeed") params.RandomSeed=StrToInt(parameter);
		else if (par=="ReplicaStructuresOutput") params.ReplicaStructuresOutput=parameter;
		else if (par=="ReplicaXRayOutput") params.ReplicaXRayOutput=parameter;
		else if (par=="RemoveProtein") params.RemoveProtein=StrToBool(parameter);
		else if (par=="RemoveUnknownAtoms") params.RemoveUnknownAtoms=StrToBool(parameter);
		else if (par=="RemoveWaters") params.RemoveWaters=StrToBool(parameter);
		else if (par=="ReplicaExchange") params.ReplicaExchange=StrToBool(parameter);
		else if (par=="RFactorWeight") params.RFactorWeight=toReal(parameter);
		else if (par=="RigidBodyOptimization") params.RigidBodyOptimization=StrToBool(parameter);
		else if (par=="RmsdDependantELL") params.RmsdDependantELL=StrToBool(parameter);
		else if (par=="RotationGridSearchPoints") params.RotationGridSearchPoints=StrToInt(parameter);
		else if (par=="RotationWeight") params.RotationWeight=toReal(parameter);
		else if (par=="RotateTranslate") params.RotateTranslate=StrToBool(parameter);
		else if (par=="ScoreFileBase") params.ScoreFileBase=parameter;
		else if (par=="SegmentsFromFile") params.SegmentsFromFile=StrToBool(parameter);
		else if (par=="SegmentsFile") params.SegmentsFile=parameter;
		else if (par=="SetTempsBasedOnScores") params.SetTempsBasedOnScores=StrToBool(parameter);
		else if (par=="ShuffleDegreesOfFreedom") params.ShuffleDegreesOfFreedom=StrToBool(parameter);
		else if (par=="ShuffleInterval") params.ShuffleInterval=toInt(parameter);
		else if (par=="SpaceGroup") params.SpaceGroup=parameter;
		else if (par=="StepsPerCycle") params.StepsPerCycle=StrToInt(parameter);
		else if (par=="TrajectoryOutputFile") params.TrajectoryOutputFile=parameter;
		else if (par=="TranslationGridSearchPoints") params.TranslationGridSearchPoints=StrToInt(parameter);
		else if (par=="TranslationRotationMethod") params.TranslationRotationMethod=parameter;
		else if (par=="UseDynamicTemperature") params.UseDynamicTemperature=StrToBool(parameter);
		else if (par=="UseExpLookUp") params.UseExpLookUp=StrToBool(parameter);
		else if (par=="UseFastTranslationRotation") params.UseFastTranslationRotation=StrToBool(parameter);
		else if (par=="UseFastTranslation") params.UseFastTranslation=StrToBool(parameter);
		else if (par=="UseMillerIndexesFromExperiment") params.UseMillerIndexesFromExperiment=StrToBool(parameter);
		else if (par=="UseNCS") params.UseNCS=StrToBool(parameter);
		else if (par=="UseVeryFastTranslation") params.UseVeryFastTranslation=StrToBool(parameter);
		else if (par=="Verbose") params.Verbose=StrToBool(parameter);
		else if (par=="XCubeLength") params.XCubeLength=StrToFloat(parameter);
		else if (par=="XCubes") params.XCubes=StrToInt(parameter);
		else if (par=="XRayInput") params.XRayInput=parameter;
		else if (par=="XRayOutput") params.XRayOutput=parameter;
		else if (par=="XRayScoreType") params.XRayScoreType=parameter;
		else if (par=="YCubeLength") params.YCubeLength=StrToFloat(parameter);
		else if (par=="YCubes") params.YCubes=StrToInt(parameter);
		else if (par=="ZCubeLength") params.ZCubeLength=StrToFloat(parameter);
		else if (par=="ZCubes") params.ZCubes=StrToInt(parameter);
		else if (par=="scale") params.scale=StrToFloat(parameter);
		else
		{
			cout <<"ERROR: Unrecognized parameter "<<par<<endl;
			exit(EXIT_FAILURE);
		}
	}
}
