# include <iostream>
# include <vector>
# include <algorithm>

using namespace std;

# include "../../LibraryFiles/error.h"
# include "../../LibraryFiles/IOUtils.h"
# include "../../LibraryFiles/StringUtils.h"
# include "../../LibraryFiles/VectorManip.h"
# include "XRay.h"

void readExperimentalDataFileJV(string experimentalDataFile, XRayStruct &xray)
{
	fstream file;
	string line;
	Real i, phase, real, imag;
	PosStruct tempPos;

	cout <<"In readExperimentalDataFileJV"<<endl;

	tempPos.Pos.resize(3);
	OpenFile(experimentalDataFile, file, "X-Ray file");
	getline(file, line);

	while (true)
	{
		if (file.eof()) break;
		tempPos.Pos[X]=StrToFloat(GetWord2(line, 1));
		tempPos.Pos[Y]=StrToFloat(GetWord2(line, 2));
		tempPos.Pos[Z]=StrToFloat(GetWord2(line, 3));
		i=StrToFloat(GetWord2(line, 4));
		phase=StrToFloat(GetWord2(line, 5));
		real=StrToFloat(GetWord2(line, 6));
		imag=StrToFloat(GetWord2(line, 7));
		SafePushBack(xray.miller, tempPos, "xray.miller");
		SafePushBack(xray.i, i, "xray.i");
		SafePushBack(xray.phase, phase, "xray.phase");
		SafePushBack(xray.real, real, "xray.real");
		SafePushBack(xray.imag, imag, "xray.imag");
		tempPos.Pos[X]=toReal(GetWord2(line, 8));
		tempPos.Pos[Y]=toReal(GetWord2(line, 9));
		tempPos.Pos[Z]=toReal(GetWord2(line, 10));
		SafePushBack(xray.q, tempPos, "xray.q");
		getline(file, line);
	}
	calcAmplitude(xray);
	cout <<"xray.i.size()= "<<xray.i.size()<<" xray.real.size()= "
	<<xray.real.size()<<" xray.imag.size()= "<<xray.imag.size()
	<<" xray.miller.size()= "<<xray.miller.size()<<" xray.q.size()= "
	<<xray.q.size()<<endl;
}

void readCifFile(string experimentalDataFile, XRayStruct &xray)
{
	fstream file;
	bool inLoop, inReflectionData, isValid;
	string line;
	vector<string> str, str2;
	int nline, ncolumn=0, loopLine, Size, Size2;
	int index_h_column, index_k_column, index_l_column;
	int F_meas_au_column, F_meas_sigma_au_column;
	int intensity_meas_column, intensity_sigma_column;
	int status_column;
	int nintensity=0, nreal=0, nimag=0, nmiller=0;
	Real i, real;
	PosStruct tempMiller;

	cout <<"In readCifFile"<<endl;

	loopLine=0;
	inLoop=false;
	inReflectionData=false;
	intensity_meas_column=UNK_INT;
	intensity_sigma_column=UNK_INT;
	F_meas_au_column=UNK_INT;
	F_meas_sigma_au_column=UNK_INT;
	index_h_column=UNK_INT;
	index_k_column=UNK_INT;
	index_l_column=UNK_INT;
	status_column=UNK_INT;


	xray.i.clear();
	xray.real.clear();
	xray.imag.clear();
	xray.miller.clear();
	xray.q.clear();

	tempMiller.Pos.resize(3);
	OpenFile(experimentalDataFile, file, "X-Ray cif file");
	getline(file, line);
	nline=0;
	while (true)
	{
		if (file.eof()) break;
		Trim(line);
		//cout <<line<<endl;
		if (line.substr(0,14)=="_cell.length_a")
		{
			xray.a=toReal(GetWord2(line, 2));
		}
		if (line.substr(0,14)=="_cell.length_b")
		{
			xray.b=toReal(GetWord2(line, 2));
		}
		if (line.substr(0,14)=="_cell.length_c")
		{
			xray.c=toReal(GetWord2(line, 2));
		}
		if (line.substr(0,17)=="_cell.angle_alpha")
		{
			xray.alpha=toReal(GetWord2(line, 2));
		}
		if (line.substr(0,16)=="_cell.angle_beta")
		{
			xray.beta=toReal(GetWord2(line, 2));
		}
		if (line.substr(0,17)=="_cell.angle_gamma")
		{
			xray.gamma=toReal(GetWord2(line, 2));
		}

		if (line.substr(0,1)!="_") 
		{
			if (inLoop) ncolumn=nline-loopLine-1;
			inLoop=false;
		}
		if (line.substr(0,5)=="loop_") 
		{
			inLoop=true;
			loopLine=nline;
		}
		if (inLoop)
		{
			cout <<line<<endl;
		}
		if (inLoop)
		{
			if (line=="_refln.index_h") 
			{
				index_h_column=nline-loopLine-1;
				cout <<"index_h_column= "<<index_h_column<<endl;
			}
			if (line=="_refln.index_k") 
			{
				index_k_column=nline-loopLine-1;
				cout <<"index_k_column= "<<index_k_column<<endl;
			}
			if (line=="_refln.index_l") 
			{
				index_l_column=nline-loopLine-1;
				cout <<"index_l_column= "<<index_l_column<<endl;
			}
			if (line=="_refln.intensity_meas" || line=="_refln.intensity_calc" || line=="_refln.F_squared_meas" || line=="_refln.pdbx_I_plus")
			{
				intensity_meas_column=nline-loopLine-1;
				cout <<"intensity_meas_column= "<<intensity_meas_column<<endl;
			}
			if (line=="_refln.intensity_sigma")
			{
				intensity_sigma_column=nline-loopLine-1;
				cout <<"intensity_sigma_column= "<<intensity_sigma_column<<endl;
			}
			if (line=="_refln.F_meas_au" || line=="_refln.pdbx_F_plus")
			{
				F_meas_au_column=nline-loopLine-1;
				cout <<"F_meas_au_column= "<<F_meas_au_column<<endl;
			}
			if (line=="_refln.F_meas_sigma_au")
			{
				F_meas_sigma_au_column=nline-loopLine-1;
				cout <<"F_meas_sigma_au_column= "<<F_meas_sigma_au_column<<endl;
			}
			if (line=="_refln.status")
			{
				status_column=nline-loopLine-1;
			}
		}
		//if (F_meas_au_column!=UNK_INT && intensity_meas_column!=UNK_INT) F_meas_au_column=UNK_INT;
		if (inReflectionData && line.length()==0) inReflectionData=true;
		//cout <<"index_h_column= "<<index_h_column<<" index_k_column= "<<index_k_column<<" index_l_column= "<<index_l_column<<" inReflectionData= "<<inReflectionData<<endl;
		if (!inLoop && index_h_column!=UNK_INT && index_k_column!=UNK_INT && index_l_column!=UNK_INT)
		{
			//cout <<"inLoop= "<<inLoop<<" index_h_column= "<<index_h_column<<" index_k_column= "<<index_k_column<<" index_l_column= "
			//<<index_l_column<<" F_meas_au_column= "<<F_meas_au_column<<endl;
			if (line.length()>0) inReflectionData=true;
		}
		//cout <<"inReflectionData= "<<inReflectionData<<endl;
		if (inReflectionData && line.substr(0,1)=="#") break;
		if (inReflectionData)
		{
			Trim(line);
			Tokenize2(line, " ", str);
			Size=str.size();
			while (Size<ncolumn)
			{
				getline(file, line);
				Tokenize2(line, " ", str2);
				Size2=str2.size();
				for (int j=0;j<Size2;j++)
				{
					SafePushBack(str, str2[j], "str");
				}
				Size=str.size();
			}
			for (int j=0;j<Size;j++)
			{
				//cout <<str[j]<<"\t"<<endl;
			}
			//cout <<endl;
			if (Size==ncolumn)
			{
				//cout <<"h= "<<str[index_h_column]<<endl;
				//cout <<"k= "<<str[index_k_column]<<endl;
				//cout <<"l= "<<str[index_l_column]<<endl;
				//cout <<"isReal(h)= "<<isReal(str[index_h_column])<<endl;
				//cout <<"isReal(k)= "<<isReal(str[index_k_column])<<endl;
				//cout <<"isReal(l)= "<<isReal(str[index_l_column])<<endl;
				isValid=true;
				if (!isReal(str[index_h_column])) isValid=false;
				if (!isReal(str[index_k_column])) isValid=false;
				if (!isReal(str[index_l_column])) isValid=false;
				if (F_meas_au_column!=UNK_INT && !isReal(str[F_meas_au_column])) isValid=false;
				if (F_meas_sigma_au_column!=UNK_INT && !isReal(str[F_meas_sigma_au_column])) isValid=false;
				if (intensity_meas_column!=UNK_INT && !isReal(str[intensity_meas_column])) isValid=false;
				if (isValid)
				{
					
					tempMiller.Pos[X]=StrToReal(str[index_h_column]);
					tempMiller.Pos[Y]=StrToReal(str[index_k_column]);
					tempMiller.Pos[Z]=StrToReal(str[index_l_column]);
					SafePushBack(xray.miller, tempMiller, "xray.miller");
					if (F_meas_au_column!=UNK_INT)
					{
						//cout <<"In F_meas_au_column"<<endl;
						if (isReal(str[F_meas_au_column]))
						{
							real=StrToReal(str[F_meas_au_column]);
							SafePushBack(xray.real, real, "xray.real");
							SafePushBack(xray.imag, Real(0), "xray.imag");
							//cout <<"real= "<<real<<endl;
							//cout <<"Pushed back to real and imag"<<endl;
						}
					}
					if (F_meas_sigma_au_column!=UNK_INT)
					{
						if (isReal(str[F_meas_sigma_au_column]))
						{
							SafePushBack(xray.Ferror, StrToReal(str[F_meas_sigma_au_column]), "Ferror");
						}
					}
					if (intensity_sigma_column!=UNK_INT)
					{
						if (isReal(str[intensity_sigma_column]))
						{
							SafePushBack(xray.error, StrToReal(str[intensity_sigma_column]), "error");
						}
					}
					if (intensity_meas_column!=UNK_INT && F_meas_au_column==UNK_INT)
					{
						if (isReal(str[intensity_meas_column]))
						{
							i=StrToReal(str[intensity_meas_column]);
							real=sqrt(i);
							SafePushBack(xray.i, i, "xray.i");
							SafePushBack(xray.real, real, "xray.real");
							SafePushBack(xray.imag, Real(0), "xray.imag");
						}
					}
					if (status_column!=UNK_INT)
					{
						if (str[status_column]=="f" || str[status_column]=="1")
						{
							SafePushBack(xray.free, true, "free");
						}
						else SafePushBack(xray.free, false, "free");
					}
				}
			}
		}
		nmiller=xray.miller.size();
		nintensity=xray.i.size();
		nreal=xray.real.size();
		nimag=xray.imag.size();
		if (nmiller!=nreal || nmiller!=nimag)
		{
			string errorStr;
			errorStr="nmiller= "+IntToStr(nmiller)+" nreal= "+IntToStr(nreal)+" nimag= "+IntToStr(nimag);
			error(errorStr, __LINE__, __FILE__);
		}
		getline(file, line);
		nline++;
	}
	if (nintensity!=nmiller)
	{
		calcIntensity(xray);
		calcPhase(xray);
	}
	calcAmplitude(xray);
	cout <<"xray.i.size()= "<<xray.i.size()<<" xray.real.size()= "
	<<xray.real.size()<<" xray.imag.size()= "<<xray.imag.size()
	<<" xray.miller.size()= "<<xray.miller.size()<<" xray.q.size()= "
	<<xray.q.size()<<endl;
}

void addErrors(XRayStruct &xray)
{
	int ni, nerror;

	ni=xray.i.size();
	nerror=xray.Ferror.size();

	if (nerror!=ni) SafeAlloc(xray.Ferror, ni, "fError");

	for (int i=0;i<ni;i++)
	{
		if (xray.Ferror[i]==0) xray.Ferror[i]=xray.F[i]*0.05;
	}
}

void readExperimentalDataFile(string experimentalDataFile, XRayStruct &xray)
{
	string ext;
	int nerror, nError;

	ext=GetExtension(experimentalDataFile);

	if (ext=="txt") readExperimentalDataFileJV(experimentalDataFile,  xray);
	else if (ext=="cif") readCifFile(experimentalDataFile, xray);
	else
	{
		string errorStr="Unrecognized extension for x-ray file "+ext;
		error(errorStr, __LINE__, __FILE__);
	}
	nerror=xray.error.size();
	nError=xray.Ferror.size();
	if (nerror==0 && nError==0) addErrors(xray);
	cout <<"In readExperimentalDataFile miller.size()= "<<xray.miller.size()<<endl;
}
