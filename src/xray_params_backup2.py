import sys
sys.path.append('/home/jouko/Scripts/python/lib/')
import IOUtils

class ParamsClass:
	def __init__(self):
		self.OptimizeBFactorScale='true'
		self.OptimizeFormFactorScale='true'
		self.Phase='false'
		self.PrintOutput='false'
		self.RemoveUnknownAtoms='false'
		self.removeWaters='false'
		self.RotateTranslate='false'
		self.UseExpLookUp='false'
		self.CubeFile=''
		self.DensityOutput=''
		self.PdbFile=''
		self.TranslationRotationMethod='RFactor'
		self.XRayInput=''
		self.XRayOutput=''
		self.XRayScoreType='DiffOfSquares'
		self.NumQxValues='16'
		self.NumQyValues='16'
		self.NumQzValues='16'
		self.MRIterations='1'
		self.XCubes='100'
		self.YCubes='100'
		self.ZCubes='100'
		self.BulkDensity='0.334'
		self.CubeSize='0.5'
		self.ExcludedVolumeRadii='Default'
		self.ExcludedVolumeRadiusHydrogen='0.902255'
		self.ExcludedVolumeRadiusCarbon='1.43693'
		self.ExcludedVolumeRadiusNitrogen='1.24527'
		self.ExcludedVolumeRadiusOxygen='1.22099'
		self.ExcludedVolumeRadiusSulfur='2.19596'
		self.scale='1.0'
		self.MaxQx='3.0'
		self.MaxQy='3.0'
		self.MaxQz='3.0'
		self.XCubeLength='0.5'
		self.YCubeLength='0.5'
		self.ZCubeLength='0.5'

def lines_to_params(lines, params):
	for i in range(len(lines)):
		lines[i]=lines[i].strip()
		if (len(lines[i].split(' '))==2):
			par=lines[i].split(' ')[0]
			parameter=lines[i].split(' ')[1]
			if (par=='OptimizeBFactorScale'): params.OptimizeBFactorScale=parameter
			elif (par=='OptimizeBFactorScale'): params.OptimizeBFactorScale=parameter
			elif (par=='OptimizeFormFactorScale'): params.OptimizeFormFactorScale=parameter
			elif (par=='Phase'): params.Phase=parameter
			elif (par=='PrintOutput'): params.PrintOutput=parameter
			elif (par=='RemoveUnknownAtoms'): params.RemoveUnknownAtoms=parameter
			elif (par=='removeWaters'): params.removeWaters=parameter
			elif (par=='RotateTranslate'): params.RotateTranslate=parameter
			elif (par=='UseExpLookUp'): params.UseExpLookUp=parameter
			elif (par=='CubeFile'): params.CubeFile=parameter
			elif (par=='DensityOutput'): params.DensityOutput=parameter
			elif (par=='PdbFile'): params.PdbFile=parameter
			elif (par=='TranslationRotationMethod'): params.TranslationRotationMethod=parameter
			elif (par=='XRayInput'): params.XRayInput=parameter
			elif (par=='XRayOutput'): params.XRayOutput=parameter
			elif (par=='XRayScoreType'): params.XRayScoreType=parameter
			elif (par=='NumQxValues'): params.NumQxValues=parameter
			elif (par=='NumQyValues'): params.NumQyValues=parameter
			elif (par=='NumQzValues'): params.NumQzValues=parameter
			elif (par=='MRIterations'): params.MRIterations=parameter
			elif (par=='XCubes'): params.XCubes=parameter
			elif (par=='YCubes'): params.YCubes=parameter
			elif (par=='ZCubes'): params.ZCubes=parameter
			elif (par=='BulkDensity'): params.BulkDensity=parameter
			elif (par=='CubeSize'): params.CubeSize=parameter
			elif (par=='ExcludedVolumeRadii'): params.ExcludedVolumeRadii=parameter
			elif (par=='ExcludedVolumeRadiusHydrogen'): params.ExcludedVolumeRadiusHydrogen=parameter
			elif (par=='ExcludedVolumeRadiusCarbon'): params.ExcludedVolumeRadiusCarbon=parameter
			elif (par=='ExcludedVolumeRadiusNitrogen'): params.ExcludedVolumeRadiusNitrogen=parameter
			elif (par=='ExcludedVolumeRadiusOxygen'): params.ExcludedVolumeRadiusOxygen=parameter
			elif (par=='ExcludedVolumeRadiusSulfur'): params.ExcludedVolumeRadiusSulfur=parameter
			elif (par=='scale'): params.scale=parameter
			elif (par=='MaxQx'): params.MaxQx=parameter
			elif (par=='MaxQy'): params.MaxQy=parameter
			elif (par=='MaxQz'): params.MaxQz=parameter
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
	file.write('OptimizeBFactorScale ' + params.OptimizeBFactorScale + '\n')
	file.write('OptimizeFormFactorScale ' + params.OptimizeFormFactorScale + '\n')
	file.write('Phase ' + params.Phase + '\n')
	file.write('PrintOutput ' + params.PrintOutput + '\n')
	file.write('RemoveUnknownAtoms ' + params.RemoveUnknownAtoms + '\n')
	file.write('removeWaters ' + params.removeWaters + '\n')
	file.write('RotateTranslate ' + params.RotateTranslate + '\n')
	file.write('UseExpLookUp ' + params.UseExpLookUp + '\n')
	file.write('CubeFile ' + params.CubeFile + '\n')
	file.write('DensityOutput ' + params.DensityOutput + '\n')
	file.write('PdbFile ' + params.PdbFile + '\n')
	file.write('TranslationRotationMethod ' + params.TranslationRotationMethod + '\n')
	file.write('XRayInput ' + params.XRayInput + '\n')
	file.write('XRayOutput ' + params.XRayOutput + '\n')
	file.write('XRayScoreType ' + params.XRayScoreType + '\n')
	file.write('NumQxValues ' + params.NumQxValues + '\n')
	file.write('NumQyValues ' + params.NumQyValues + '\n')
	file.write('NumQzValues ' + params.NumQzValues + '\n')
	file.write('MRIterations ' + params.MRIterations + '\n')
	file.write('XCubes ' + params.XCubes + '\n')
	file.write('YCubes ' + params.YCubes + '\n')
	file.write('ZCubes ' + params.ZCubes + '\n')
	file.write('BulkDensity ' + params.BulkDensity + '\n')
	file.write('CubeSize ' + params.CubeSize + '\n')
	file.write('ExcludedVolumeRadii ' + params.ExcludedVolumeRadii + '\n')
	file.write('ExcludedVolumeRadiusHydrogen ' + params.ExcludedVolumeRadiusHydrogen + '\n')
	file.write('ExcludedVolumeRadiusCarbon ' + params.ExcludedVolumeRadiusCarbon + '\n')
	file.write('ExcludedVolumeRadiusNitrogen ' + params.ExcludedVolumeRadiusNitrogen + '\n')
	file.write('ExcludedVolumeRadiusOxygen ' + params.ExcludedVolumeRadiusOxygen + '\n')
	file.write('ExcludedVolumeRadiusSulfur ' + params.ExcludedVolumeRadiusSulfur + '\n')
	file.write('scale ' + params.scale + '\n')
	file.write('MaxQx ' + params.MaxQx + '\n')
	file.write('MaxQy ' + params.MaxQy + '\n')
	file.write('MaxQz ' + params.MaxQz + '\n')
	file.write('XCubeLength ' + params.XCubeLength + '\n')
	file.write('YCubeLength ' + params.YCubeLength + '\n')
	file.write('ZCubeLength ' + params.ZCubeLength + '\n')

def copy_params(params1, params2):
	params2.OptimizeBFactorScale=params1.OptimizeBFactorScale
	params2.OptimizeFormFactorScale=params1.OptimizeFormFactorScale
	params2.Phase=params1.Phase
	params2.PrintOutput=params1.PrintOutput
	params2.RemoveUnknownAtoms=params1.RemoveUnknownAtoms
	params2.removeWaters=params1.removeWaters
	params2.RotateTranslate=params1.RotateTranslate
	params2.UseExpLookUp=params1.UseExpLookUp
	params2.CubeFile=params1.CubeFile
	params2.DensityOutput=params1.DensityOutput
	params2.PdbFile=params1.PdbFile
	params2.TranslationRotationMethod=params1.TranslationRotationMethod
	params2.XRayInput=params1.XRayInput
	params2.XRayOutput=params1.XRayOutput
	params2.XRayScoreType=params1.XRayScoreType
	params2.NumQxValues=params1.NumQxValues
	params2.NumQyValues=params1.NumQyValues
	params2.NumQzValues=params1.NumQzValues
	params2.MRIterations=params1.MRIterations
	params2.XCubes=params1.XCubes
	params2.YCubes=params1.YCubes
	params2.ZCubes=params1.ZCubes
	params2.BulkDensity=params1.BulkDensity
	params2.CubeSize=params1.CubeSize
	params2.ExcludedVolumeRadii=params1.ExcludedVolumeRadii
	params2.ExcludedVolumeRadiusHydrogen=params1.ExcludedVolumeRadiusHydrogen
	params2.ExcludedVolumeRadiusCarbon=params1.ExcludedVolumeRadiusCarbon
	params2.ExcludedVolumeRadiusNitrogen=params1.ExcludedVolumeRadiusNitrogen
	params2.ExcludedVolumeRadiusOxygen=params1.ExcludedVolumeRadiusOxygen
	params2.ExcludedVolumeRadiusSulfur=params1.ExcludedVolumeRadiusSulfur
	params2.scale=params1.scale
	params2.MaxQx=params1.MaxQx
	params2.MaxQy=params1.MaxQy
	params2.MaxQz=params1.MaxQz
	params2.XCubeLength=params1.XCubeLength
	params2.YCubeLength=params1.YCubeLength
	params2.ZCubeLength=params1.ZCubeLength

