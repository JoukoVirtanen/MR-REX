import sys
sys.path.append('/home/jouko/Scripts/python/lib/')
import IOUtils

class ParamsClass:
	def __init__(self):
		self.a='UNK_REAL'
		self.b='UNK_REAL'
		self.c='UNK_REAL'
		self.alpha='UNK_REAL'
		self.beta='UNK_REAL'
		self.gamma='UNK_REAL'
		self.BulkSolventCorrection='false'
		self.FilterOutHighResolution='false'
		self.MinimizeAfterEachStep='false'
		self.MinimizeEverySolution='false'
		self.OptimizeBFactorScale='true'
		self.OptimizeFormFactorScale='true'
		self.Phase='false'
		self.PrintOutput='false'
		self.RemoveProtein='false'
		self.RemoveUnknownAtoms='false'
		self.removeWaters='false'
		self.ReplicaExchange='false'
		self.RigidBodyOptimization='false'
		self.RotateTranslate='false'
		self.SetTempsBasedOnScores='true'
		self.UseExpLookUp='false'
		self.UseFastRotationTranslation='false'
		self.UseFastTranslation='false'
		self.UseMillerIndexesFromExperiment='false'
		self.ClashScoreUse='AtEnd'
		self.CubeFile=''
		self.DcdFilePaths=''
		self.DegOfFreedomFile=''
		self.DensityOutput=''
		self.FastTranslationRotationInterpolation='Hessian'
		self.PdbFile=''
		self.SpaceGroup=''
		self.TranslationRotationMethod='RFactor'
		self.XRayInput=''
		self.XRayOutput=''
		self.XRayScoreType='DiffOfSquares'
		self.ContinuousHValues='50'
		self.ContinuousKValues='50'
		self.ContinuousLValues='50'
		self.MaxTrajectoryStructures='1000'
		self.NumCopies='1'
		self.NumCycles='100'
		self.NumQxValues='16'
		self.NumQyValues='16'
		self.NumQzValues='16'
		self.NumReplicas='70'
		self.MRIterations='1'
		self.RandomSeed='0'
		self.ReplicaStructuresOutput=''
		self.RotationGridSearchPoints='31'
		self.StepsPerCycle='100'
		self.TranslationGridSearchPoints='31'
		self.XCubes='100'
		self.YCubes='100'
		self.ZCubes='100'
		self.BFactorScale='1.0'
		self.BulkDensity='0.334'
		self.ClashWeightCoreCore='1000.0'
		self.ClashWeightCoreSurface='10.0'
		self.ClashWeightSurfaceSurface='1.0'
		self.CubeSize='0.5'
		self.ExcludedVolumeRadii='Default'
		self.ExcludedVolumeRadiusHydrogen='0.902255'
		self.ExcludedVolumeRadiusCarbon='1.43693'
		self.ExcludedVolumeRadiusNitrogen='1.24527'
		self.ExcludedVolumeRadiusOxygen='1.22099'
		self.ExcludedVolumeRadiusSulfur='2.19596'
		self.scale='1.0'
		self.MaxReplicaExchangeTime='3.0'
		self.MaxTemp='1.0e14'
		self.MaxQx='3.0'
		self.MaxQy='3.0'
		self.MaxQz='3.0'
		self.MinTemp='1.0e7'
		self.XCubeLength='0.5'
		self.YCubeLength='0.5'
		self.ZCubeLength='0.5'

def lines_to_params(lines, params):
	for i in range(len(lines)):
		lines[i]=lines[i].strip()
		if (len(lines[i].split(' '))==2):
			par=lines[i].split(' ')[0]
			parameter=lines[i].split(' ')[1]
			if (par=='a'): params.a=parameter
			elif (par=='a'): params.a=parameter
			elif (par=='b'): params.b=parameter
			elif (par=='c'): params.c=parameter
			elif (par=='alpha'): params.alpha=parameter
			elif (par=='beta'): params.beta=parameter
			elif (par=='gamma'): params.gamma=parameter
			elif (par=='BulkSolventCorrection'): params.BulkSolventCorrection=parameter
			elif (par=='FilterOutHighResolution'): params.FilterOutHighResolution=parameter
			elif (par=='MinimizeAfterEachStep'): params.MinimizeAfterEachStep=parameter
			elif (par=='MinimizeEverySolution'): params.MinimizeEverySolution=parameter
			elif (par=='OptimizeBFactorScale'): params.OptimizeBFactorScale=parameter
			elif (par=='OptimizeFormFactorScale'): params.OptimizeFormFactorScale=parameter
			elif (par=='Phase'): params.Phase=parameter
			elif (par=='PrintOutput'): params.PrintOutput=parameter
			elif (par=='RemoveProtein'): params.RemoveProtein=parameter
			elif (par=='RemoveUnknownAtoms'): params.RemoveUnknownAtoms=parameter
			elif (par=='removeWaters'): params.removeWaters=parameter
			elif (par=='ReplicaExchange'): params.ReplicaExchange=parameter
			elif (par=='RigidBodyOptimization'): params.RigidBodyOptimization=parameter
			elif (par=='RotateTranslate'): params.RotateTranslate=parameter
			elif (par=='SetTempsBasedOnScores'): params.SetTempsBasedOnScores=parameter
			elif (par=='UseExpLookUp'): params.UseExpLookUp=parameter
			elif (par=='UseFastRotationTranslation'): params.UseFastRotationTranslation=parameter
			elif (par=='UseFastTranslation'): params.UseFastTranslation=parameter
			elif (par=='UseMillerIndexesFromExperiment'): params.UseMillerIndexesFromExperiment=parameter
			elif (par=='ClashScoreUse'): params.ClashScoreUse=parameter
			elif (par=='CubeFile'): params.CubeFile=parameter
			elif (par=='DcdFilePaths'): params.DcdFilePaths=parameter
			elif (par=='DegOfFreedomFile'): params.DegOfFreedomFile=parameter
			elif (par=='DensityOutput'): params.DensityOutput=parameter
			elif (par=='FastTranslationRotationInterpolation'): params.FastTranslationRotationInterpolation=parameter
			elif (par=='PdbFile'): params.PdbFile=parameter
			elif (par=='SpaceGroup'): params.SpaceGroup=parameter
			elif (par=='TranslationRotationMethod'): params.TranslationRotationMethod=parameter
			elif (par=='XRayInput'): params.XRayInput=parameter
			elif (par=='XRayOutput'): params.XRayOutput=parameter
			elif (par=='XRayScoreType'): params.XRayScoreType=parameter
			elif (par=='ContinuousHValues'): params.ContinuousHValues=parameter
			elif (par=='ContinuousKValues'): params.ContinuousKValues=parameter
			elif (par=='ContinuousLValues'): params.ContinuousLValues=parameter
			elif (par=='MaxTrajectoryStructures'): params.MaxTrajectoryStructures=parameter
			elif (par=='NumCopies'): params.NumCopies=parameter
			elif (par=='NumCycles'): params.NumCycles=parameter
			elif (par=='NumQxValues'): params.NumQxValues=parameter
			elif (par=='NumQyValues'): params.NumQyValues=parameter
			elif (par=='NumQzValues'): params.NumQzValues=parameter
			elif (par=='NumReplicas'): params.NumReplicas=parameter
			elif (par=='MRIterations'): params.MRIterations=parameter
			elif (par=='RandomSeed'): params.RandomSeed=parameter
			elif (par=='ReplicaStructuresOutput'): params.ReplicaStructuresOutput=parameter
			elif (par=='RotationGridSearchPoints'): params.RotationGridSearchPoints=parameter
			elif (par=='StepsPerCycle'): params.StepsPerCycle=parameter
			elif (par=='TranslationGridSearchPoints'): params.TranslationGridSearchPoints=parameter
			elif (par=='XCubes'): params.XCubes=parameter
			elif (par=='YCubes'): params.YCubes=parameter
			elif (par=='ZCubes'): params.ZCubes=parameter
			elif (par=='BFactorScale'): params.BFactorScale=parameter
			elif (par=='BulkDensity'): params.BulkDensity=parameter
			elif (par=='ClashWeightCoreCore'): params.ClashWeightCoreCore=parameter
			elif (par=='ClashWeightCoreSurface'): params.ClashWeightCoreSurface=parameter
			elif (par=='ClashWeightSurfaceSurface'): params.ClashWeightSurfaceSurface=parameter
			elif (par=='CubeSize'): params.CubeSize=parameter
			elif (par=='ExcludedVolumeRadii'): params.ExcludedVolumeRadii=parameter
			elif (par=='ExcludedVolumeRadiusHydrogen'): params.ExcludedVolumeRadiusHydrogen=parameter
			elif (par=='ExcludedVolumeRadiusCarbon'): params.ExcludedVolumeRadiusCarbon=parameter
			elif (par=='ExcludedVolumeRadiusNitrogen'): params.ExcludedVolumeRadiusNitrogen=parameter
			elif (par=='ExcludedVolumeRadiusOxygen'): params.ExcludedVolumeRadiusOxygen=parameter
			elif (par=='ExcludedVolumeRadiusSulfur'): params.ExcludedVolumeRadiusSulfur=parameter
			elif (par=='scale'): params.scale=parameter
			elif (par=='MaxReplicaExchangeTime'): params.MaxReplicaExchangeTime=parameter
			elif (par=='MaxTemp'): params.MaxTemp=parameter
			elif (par=='MaxQx'): params.MaxQx=parameter
			elif (par=='MaxQy'): params.MaxQy=parameter
			elif (par=='MaxQz'): params.MaxQz=parameter
			elif (par=='MinTemp'): params.MinTemp=parameter
			elif (par=='XCubeLength'): params.XCubeLength=parameter
			elif (par=='YCubeLength'): params.YCubeLength=parameter
			elif (par=='ZCubeLength'): params.ZCubeLength=parameter
			else: 
				print 'ERROR unkown parameter ' + par
				sys.exit()

def read_params_file(ParamsFile, params):
	lines=IOUtils.readlines(ParamsFile)
	lines_to_params(lines, params)

def write_params_file(ParamsFile, params):
	IOUtils.delete_file(ParamsFile)
	file=open(ParamsFile, 'a')
	file.write('a ' + params.a + '\n')
	file.write('b ' + params.b + '\n')
	file.write('c ' + params.c + '\n')
	file.write('alpha ' + params.alpha + '\n')
	file.write('beta ' + params.beta + '\n')
	file.write('gamma ' + params.gamma + '\n')
	file.write('BulkSolventCorrection ' + params.BulkSolventCorrection + '\n')
	file.write('FilterOutHighResolution ' + params.FilterOutHighResolution + '\n')
	file.write('MinimizeAfterEachStep ' + params.MinimizeAfterEachStep + '\n')
	file.write('MinimizeEverySolution ' + params.MinimizeEverySolution + '\n')
	file.write('OptimizeBFactorScale ' + params.OptimizeBFactorScale + '\n')
	file.write('OptimizeFormFactorScale ' + params.OptimizeFormFactorScale + '\n')
	file.write('Phase ' + params.Phase + '\n')
	file.write('PrintOutput ' + params.PrintOutput + '\n')
	file.write('RemoveProtein ' + params.RemoveProtein + '\n')
	file.write('RemoveUnknownAtoms ' + params.RemoveUnknownAtoms + '\n')
	file.write('removeWaters ' + params.removeWaters + '\n')
	file.write('ReplicaExchange ' + params.ReplicaExchange + '\n')
	file.write('RigidBodyOptimization ' + params.RigidBodyOptimization + '\n')
	file.write('RotateTranslate ' + params.RotateTranslate + '\n')
	file.write('SetTempsBasedOnScores ' + params.SetTempsBasedOnScores + '\n')
	file.write('UseExpLookUp ' + params.UseExpLookUp + '\n')
	file.write('UseFastRotationTranslation ' + params.UseFastRotationTranslation + '\n')
	file.write('UseFastTranslation ' + params.UseFastTranslation + '\n')
	file.write('UseMillerIndexesFromExperiment ' + params.UseMillerIndexesFromExperiment + '\n')
	file.write('ClashScoreUse ' + params.ClashScoreUse + '\n')
	file.write('CubeFile ' + params.CubeFile + '\n')
	file.write('DcdFilePaths ' + params.DcdFilePaths + '\n')
	file.write('DegOfFreedomFile ' + params.DegOfFreedomFile + '\n')
	file.write('DensityOutput ' + params.DensityOutput + '\n')
	file.write('FastTranslationRotationInterpolation ' + params.FastTranslationRotationInterpolation + '\n')
	file.write('PdbFile ' + params.PdbFile + '\n')
	file.write('SpaceGroup ' + params.SpaceGroup + '\n')
	file.write('TranslationRotationMethod ' + params.TranslationRotationMethod + '\n')
	file.write('XRayInput ' + params.XRayInput + '\n')
	file.write('XRayOutput ' + params.XRayOutput + '\n')
	file.write('XRayScoreType ' + params.XRayScoreType + '\n')
	file.write('ContinuousHValues ' + params.ContinuousHValues + '\n')
	file.write('ContinuousKValues ' + params.ContinuousKValues + '\n')
	file.write('ContinuousLValues ' + params.ContinuousLValues + '\n')
	file.write('MaxTrajectoryStructures ' + params.MaxTrajectoryStructures + '\n')
	file.write('NumCopies ' + params.NumCopies + '\n')
	file.write('NumCycles ' + params.NumCycles + '\n')
	file.write('NumQxValues ' + params.NumQxValues + '\n')
	file.write('NumQyValues ' + params.NumQyValues + '\n')
	file.write('NumQzValues ' + params.NumQzValues + '\n')
	file.write('NumReplicas ' + params.NumReplicas + '\n')
	file.write('MRIterations ' + params.MRIterations + '\n')
	file.write('RandomSeed ' + params.RandomSeed + '\n')
	file.write('ReplicaStructuresOutput ' + params.ReplicaStructuresOutput + '\n')
	file.write('RotationGridSearchPoints ' + params.RotationGridSearchPoints + '\n')
	file.write('StepsPerCycle ' + params.StepsPerCycle + '\n')
	file.write('TranslationGridSearchPoints ' + params.TranslationGridSearchPoints + '\n')
	file.write('XCubes ' + params.XCubes + '\n')
	file.write('YCubes ' + params.YCubes + '\n')
	file.write('ZCubes ' + params.ZCubes + '\n')
	file.write('BFactorScale ' + params.BFactorScale + '\n')
	file.write('BulkDensity ' + params.BulkDensity + '\n')
	file.write('ClashWeightCoreCore ' + params.ClashWeightCoreCore + '\n')
	file.write('ClashWeightCoreSurface ' + params.ClashWeightCoreSurface + '\n')
	file.write('ClashWeightSurfaceSurface ' + params.ClashWeightSurfaceSurface + '\n')
	file.write('CubeSize ' + params.CubeSize + '\n')
	file.write('ExcludedVolumeRadii ' + params.ExcludedVolumeRadii + '\n')
	file.write('ExcludedVolumeRadiusHydrogen ' + params.ExcludedVolumeRadiusHydrogen + '\n')
	file.write('ExcludedVolumeRadiusCarbon ' + params.ExcludedVolumeRadiusCarbon + '\n')
	file.write('ExcludedVolumeRadiusNitrogen ' + params.ExcludedVolumeRadiusNitrogen + '\n')
	file.write('ExcludedVolumeRadiusOxygen ' + params.ExcludedVolumeRadiusOxygen + '\n')
	file.write('ExcludedVolumeRadiusSulfur ' + params.ExcludedVolumeRadiusSulfur + '\n')
	file.write('scale ' + params.scale + '\n')
	file.write('MaxReplicaExchangeTime ' + params.MaxReplicaExchangeTime + '\n')
	file.write('MaxTemp ' + params.MaxTemp + '\n')
	file.write('MaxQx ' + params.MaxQx + '\n')
	file.write('MaxQy ' + params.MaxQy + '\n')
	file.write('MaxQz ' + params.MaxQz + '\n')
	file.write('MinTemp ' + params.MinTemp + '\n')
	file.write('XCubeLength ' + params.XCubeLength + '\n')
	file.write('YCubeLength ' + params.YCubeLength + '\n')
	file.write('ZCubeLength ' + params.ZCubeLength + '\n')

def copy_params(params1, params2):
	params2.a=params1.a
	params2.b=params1.b
	params2.c=params1.c
	params2.alpha=params1.alpha
	params2.beta=params1.beta
	params2.gamma=params1.gamma
	params2.BulkSolventCorrection=params1.BulkSolventCorrection
	params2.FilterOutHighResolution=params1.FilterOutHighResolution
	params2.MinimizeAfterEachStep=params1.MinimizeAfterEachStep
	params2.MinimizeEverySolution=params1.MinimizeEverySolution
	params2.OptimizeBFactorScale=params1.OptimizeBFactorScale
	params2.OptimizeFormFactorScale=params1.OptimizeFormFactorScale
	params2.Phase=params1.Phase
	params2.PrintOutput=params1.PrintOutput
	params2.RemoveProtein=params1.RemoveProtein
	params2.RemoveUnknownAtoms=params1.RemoveUnknownAtoms
	params2.removeWaters=params1.removeWaters
	params2.ReplicaExchange=params1.ReplicaExchange
	params2.RigidBodyOptimization=params1.RigidBodyOptimization
	params2.RotateTranslate=params1.RotateTranslate
	params2.SetTempsBasedOnScores=params1.SetTempsBasedOnScores
	params2.UseExpLookUp=params1.UseExpLookUp
	params2.UseFastRotationTranslation=params1.UseFastRotationTranslation
	params2.UseFastTranslation=params1.UseFastTranslation
	params2.UseMillerIndexesFromExperiment=params1.UseMillerIndexesFromExperiment
	params2.ClashScoreUse=params1.ClashScoreUse
	params2.CubeFile=params1.CubeFile
	params2.DcdFilePaths=params1.DcdFilePaths
	params2.DegOfFreedomFile=params1.DegOfFreedomFile
	params2.DensityOutput=params1.DensityOutput
	params2.FastTranslationRotationInterpolation=params1.FastTranslationRotationInterpolation
	params2.PdbFile=params1.PdbFile
	params2.SpaceGroup=params1.SpaceGroup
	params2.TranslationRotationMethod=params1.TranslationRotationMethod
	params2.XRayInput=params1.XRayInput
	params2.XRayOutput=params1.XRayOutput
	params2.XRayScoreType=params1.XRayScoreType
	params2.ContinuousHValues=params1.ContinuousHValues
	params2.ContinuousKValues=params1.ContinuousKValues
	params2.ContinuousLValues=params1.ContinuousLValues
	params2.MaxTrajectoryStructures=params1.MaxTrajectoryStructures
	params2.NumCopies=params1.NumCopies
	params2.NumCycles=params1.NumCycles
	params2.NumQxValues=params1.NumQxValues
	params2.NumQyValues=params1.NumQyValues
	params2.NumQzValues=params1.NumQzValues
	params2.NumReplicas=params1.NumReplicas
	params2.MRIterations=params1.MRIterations
	params2.RandomSeed=params1.RandomSeed
	params2.ReplicaStructuresOutput=params1.ReplicaStructuresOutput
	params2.RotationGridSearchPoints=params1.RotationGridSearchPoints
	params2.StepsPerCycle=params1.StepsPerCycle
	params2.TranslationGridSearchPoints=params1.TranslationGridSearchPoints
	params2.XCubes=params1.XCubes
	params2.YCubes=params1.YCubes
	params2.ZCubes=params1.ZCubes
	params2.BFactorScale=params1.BFactorScale
	params2.BulkDensity=params1.BulkDensity
	params2.ClashWeightCoreCore=params1.ClashWeightCoreCore
	params2.ClashWeightCoreSurface=params1.ClashWeightCoreSurface
	params2.ClashWeightSurfaceSurface=params1.ClashWeightSurfaceSurface
	params2.CubeSize=params1.CubeSize
	params2.ExcludedVolumeRadii=params1.ExcludedVolumeRadii
	params2.ExcludedVolumeRadiusHydrogen=params1.ExcludedVolumeRadiusHydrogen
	params2.ExcludedVolumeRadiusCarbon=params1.ExcludedVolumeRadiusCarbon
	params2.ExcludedVolumeRadiusNitrogen=params1.ExcludedVolumeRadiusNitrogen
	params2.ExcludedVolumeRadiusOxygen=params1.ExcludedVolumeRadiusOxygen
	params2.ExcludedVolumeRadiusSulfur=params1.ExcludedVolumeRadiusSulfur
	params2.scale=params1.scale
	params2.MaxReplicaExchangeTime=params1.MaxReplicaExchangeTime
	params2.MaxTemp=params1.MaxTemp
	params2.MaxQx=params1.MaxQx
	params2.MaxQy=params1.MaxQy
	params2.MaxQz=params1.MaxQz
	params2.MinTemp=params1.MinTemp
	params2.XCubeLength=params1.XCubeLength
	params2.YCubeLength=params1.YCubeLength
	params2.ZCubeLength=params1.ZCubeLength

