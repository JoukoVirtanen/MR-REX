# include <vector>
# include <iostream>

using namespace std;

# include "../../LibraryFiles/MathUtils.h"
# include "../../LibraryFiles/normalize.h"
# include "../../LibraryFiles/MinMax.h"
# include "../../LibraryFiles/time.h"
# include "XRayStruct.h"
# include "XRay.h"

Real calcWoolfson(Real fObs, Real fCalc, Real fVar)
{
	//Woolfson, M. M. (1956). Acta Cryst. 9, 804-810.
	Real w;

	w=sqrt(2.0/(pi*fVar))*cosh(fObs*fCalc/fVar)/exp((fObs*fObs+fCalc*fCalc)/(2.0*fVar));
	return w;
}

Real calcLogWoolfson(Real fObs, Real fCalc, Real fVar)
{
	//Woolfson, M. M. (1956). Acta Cryst. 9, 804-810.
	Real w, x, cutoff=15.0;

	x=fObs*fCalc/fVar;
	w=log(sqrt(2.0/(pi*fVar)));
	if (x<cutoff) w+=log(cosh(x));
	else w+=x-log(2);
	w-=(fObs*fObs+fCalc*fCalc)/(2.0*fVar);
	return w;
}

Real calcRice(Real fObs, Real fCalc, Real fVar)
{
	//Sim, G. A. (1959). Acta Cryst. 12, 813-815
	//The distributio of phase angles for structures...
	Real r, i0;

	i0=calcBesselI0(2.0*fObs*fCalc/fVar);
	r=(2.0*fObs/fVar)*i0/exp((fObs*fObs+fCalc*fCalc)/fVar);
	return r;
}

Real calcLogRice(Real fObs, Real fCalc, Real fVar)
{
	//Sim, G. A. (1959). Acta Cryst. 12, 813-815
	//The distributio of phase angles for structures...
	Real r, x;

	x=2.0*fObs*fCalc/fVar;
	r=log(2.0);
	r+=log(fObs);
	r-=log(fVar);
	r+=calcLogBesselI0(x);
	r-=(fObs*fObs+fCalc*fCalc)/fVar;

	return r;
}

Real calcMLRF0_NoScale(vector<Real> &fObs, vector<Real> &fCalc, vector<Real> &fError, vector<bool> &centric)
{
	//J. Appl. Cryst. (2007). 40, 658-674
	int nobs=fObs.size();
	int ncalc=fCalc.size();
	int nerror=fError.size();
	int ncentric=centric.size();
	Real mlrf0=0, fVar;
	Real sumFObs2, sumFCalc2;

	if (nobs==0 || nobs!=ncalc || nobs!=nerror || nobs!=ncentric)
	{
		string errorStr="nobs= "+IntToStr(nobs)+" ncalc= ";
		errorStr+=IntToStr(ncalc)+" nerror= "+IntToStr(nerror);
		errorStr+=" ncentric= "+IntToStr(ncentric);
		error(errorStr, __LINE__, __FILE__);
	}
	
	sumFObs2=calcSumSqr(fObs);
	cout <<"sumFObs2= "<<sumFObs2<<endl;
	sumFCalc2=calcSumSqr(fCalc);
	cout <<"sumFCalc2= "<<sumFCalc2<<endl;

	for (int i=0;i<nobs;i++)
	{
		cout <<"i= "<<i<<endl;
		fVar=sumFObs2+fCalc[i]*fCalc[i]-sumFCalc2+fError[i];
		if (centric[i]) mlrf0+=log(calcWoolfson(fObs[i], 0, fVar));
		else mlrf0+=log(calcRice(fObs[i], 0, fVar));
	}
	return mlrf0;
}

Real calcMLRF0(vector<Real> &fObs, vector<Real> &fCalc, vector<Real> &fError, vector<bool> &centric)
{
	//Get scale that minimizes difference of squares.  Change later?
	Real scale, mlrf0;
	Vector fCalcScaled;
	cout <<"Before calcScaleForSqrDiff"<<endl;
	scale=calcScaleForSqrDiff(fCalc, fObs);
	cout <<"scale= "<<scale<<endl;
	multiply(scale, fCalc, fCalcScaled);
	cout <<"Before calcMLRF0_NoScale"<<endl;
	mlrf0=calcMLRF0_NoScale(fObs, fCalcScaled, fError, centric);
	cout <<"mlrf0"<<endl;
	return mlrf0;
}

Real calcMLRF0(XRayStruct &xray, XRayStruct &expXRay)
{
	Real mlrf0;
	calcAmplitude(expXRay);
	calcAmplitude(xray);
	mlrf0=calcMLRF0(xray.F, expXRay.F, expXRay.Ferror, xray.centric);
	cout <<"mlrf0= "<<mlrf0<<endl;
	return mlrf0;
}

Real calcMLTF_NoScale(Vector &fObs, Vector &fCalc, Vector &fMove, Vector &fFix, Vector &fError, vector<bool> &centric)
{
	//J. Appl. Cryst. (2007). 40, 658-674
	int nobs=fObs.size();
	int ncalc=fCalc.size();
	int nmove=fMove.size();
	int nerror=fError.size();
	int ncentric=centric.size();
	Real mltf=0, fVar;
	timeval start, end;
	cout <<"In calcMLTF"<<endl;
	if (nobs==0 || nobs!=ncalc || nobs!=nmove || nobs!=nerror || nobs!=ncentric)
	{
		string errorStr="nobs= "+IntToStr(nobs)+" ncalc= ";
		errorStr+=IntToStr(ncalc)+" nmove= "+IntToStr(nmove);
		errorStr+=" nerror= "+IntToStr(nerror);
		errorStr+=" ncentric= "+IntToStr(ncentric);
		error(errorStr, __LINE__, __FILE__);
	}

	gettimeofday(&start, NULL);
	for (int i=0;i<nobs;i++)
	{
		fVar=fCalc[i]*fCalc[i]-fFix[i]*fFix[i]-fMove[i]*fMove[i]+fError[i]*fError[i];
		if (fVar==0 || isnan(fVar)) fVar=fObs[i]*0.1;
		//cout <<"fObs= "<<fObs[i]<<" fCalc= "<<fCalc[i]<<" fVar= "<<fVar<<endl;
		if (centric[i]) mltf+=calcLogWoolfson(fObs[i], fCalc[i], fVar);
		else mltf+=calcLogRice(fObs[i], fCalc[i], fVar);	
	}
	gettimeofday(&end, NULL);
	cout <<"calcMLTF took "<<calcTimeDiff(start, end)<<endl;
	mltf=-mltf;
	return mltf;
}

Real calcMLTF(Vector &fObs, Vector fCalc, Vector &fMove, Vector &fFix, Vector &fError, vector<bool> &centric)
{
	int npoint=fCalc.size();
	Real scale=calcScaleForSqrDiff(fCalc, fObs);

	for (int i=0;i<npoint;i++)
	{
		fCalc[i]*=scale;
	}
	fMove=fCalc;
	return calcMLTF_NoScale(fObs, fCalc, fMove, fFix, fError, centric);
}

Real calcMLTF(XRayStruct &xray, XRayStruct &expXRay, int nmove)
{
	int npoint=xray.F.size();
	Real mltf;
	//calcMove(xray, nmove);
	//calcFixed(xray, nmove);
	SafeAlloc(xray.Ffix, npoint, "xray.Ffix");
	mltf=calcMLTF(expXRay.F, xray.F, xray.Fmove, xray.Ffix, expXRay.Ferror, expXRay.centric);

	return mltf;
}

void calcMeanIntensity(XRayStruct &xray, Vector &meanI)
{
	//Given the amplitudes of symmetry mates calculates the
	//average intensity over a random walk of phases.
	//Does so by evaluating the following integral.
	//inegral[ (sigma(r_i*cos(theta_i))^2 + (sigma(r_i*sin(theta_i))^2 ] dtheta
	int nprot, nsym, npoint;

	Get3DVectorSize(xray.complexSymReal, nprot, nsym, npoint, "complexSym");


	SafeAlloc(meanI, npoint, "rootMeanI");
	for (int i=0;i<npoint;i++)
	{
		meanI[i]=0;
		for (int j=0;j<nsym;j++)
		{
			for (int k=0;k<nprot;k++)
			{
				meanI[i]+=xray.complexSymReal[k][j][i]*xray.complexSymReal[k][j][i];
				meanI[i]+=xray.complexSymImag[k][j][i]*xray.complexSymImag[k][j][i];
			}
		}
	}
}

void calcRootMeanIntensity(XRayStruct &xray, Vector &rootMeanI)
{
	//Given the amplitudes of symmetry mates calculates the
	//average intensity over a random walk of phases.
	//Does so by evaluating the following integral.
	//inegral[ (sigma(r_i*cos(theta_i))^2 + (sigma(r_i*sin(theta_i))^2 ] dtheta
	int npoint;

	calcMeanIntensity(xray, rootMeanI);
	npoint=rootMeanI.size();
	for (int i=0;i<npoint;i++)
	{
		rootMeanI[i]=sqrt(rootMeanI[i]);
	}
}

void calcAverageSquaredIntensity(XRayStruct &xray, Vector &aveSqrIntensity)
{
	//Calculates the average squared intensity over a random walk of phases.
	//Does so by evaluating the following integral.
	//inegral[ (sigma(r_i*cos(theta_i))^2 + (sigma(r_i*sin(theta_i))^2 ]^2 dtheta
	int nsym=xray.symReal.size();
	int npoint=0;
	Real sum;
	Real a2, a4, b2;
	//Real a2, a4, b2, b4;

	if (nsym>0)
	{
		npoint=xray.symReal[0].size();
	}
	else
	{
		string errorStr="No info on amplitudes of symmetry mates.";
		error(errorStr, __LINE__, __FILE__);
	}

	SafeAlloc(aveSqrIntensity, npoint, "aveSqrIntensity");
	for (int i=0;i<npoint;i++)
	{
		sum=0;
		for (int j=0;j<nsym;j++)
		{
			a2=xray.symReal[j][i]*xray.symReal[j][i];
			a2+=xray.symImag[j][i]*xray.symImag[j][i];
			a4=a2*a2;
			//sum+=3.0*a4/8.0;
			//sum+=a4/8.0;
			sum+=a4*0.5;
			for (int k=j;k<nsym;k++)
			{
				b2=xray.symReal[j][i]*xray.symReal[j][i];
				b2+=xray.symImag[j][i]*xray.symImag[j][i];
				//b4=a2*a2;
				//sum+=3.0*(real2a*real2b+imag2a*imag2b)/2.0;
				//sum+=(real2a*real2b+imag2a*imag2b)/2.0;
				sum+=2.5*a2*b2;
			}
		}
		aveSqrIntensity[i]=sum;
	}
}

void calcAverageAndStdDevIntensity(XRayStruct &xray, Vector &rootMeanIntensity, Vector &stdDevIntensity)
{
	int npoint;
	int index=22;
	Real variance;
	Vector rootMeanI, aveSqrIntensity;

	cout <<"Before calcRootMeanIntensity(xray, rootMeanI)"<<endl;
	calcRootMeanIntensity(xray, rootMeanI);
	cout <<"Before calcAverageSquaredIntensity(xray, aveSqrIntensity)"<<endl;
	calcAverageSquaredIntensity(xray, aveSqrIntensity);	
	cout <<"After calcAverageSquaredIntensity(xray, aveSqrIntensity)"<<endl;
	npoint=rootMeanI.size();
	SafeAlloc(stdDevIntensity, npoint);
	for (int i=0;i<npoint;i++)
	{
		variance=aveSqrIntensity[i]-rootMeanI[i]*rootMeanI[i];
		stdDevIntensity[i]=sqrt(variance);
	}
	cout <<"rootMeanI["<<index<<"]= "<<rootMeanI[index]<<" stdDevIntensity= "<<stdDevIntensity[index]<<endl;
}

Real calcRandomIntensity(XRayStruct &xray, int index)
{
	int nsym, nprot, npoint;
	Real real=0, imag=0, dp;

	Get3DVectorSize(xray.complexSymReal, nprot, nsym, npoint, "complexSymReal");
	for (int i=0;i<nprot;i++)
	{
		for (int j=0;j<nsym;j++)
		{
			dp=randDouble(0, 2.0*pi);
			real+=cos(dp)*xray.complexSymReal[i][j][index];
			real-=sin(dp)*xray.complexSymImag[i][j][index];
			imag+=sin(dp)*xray.complexSymReal[i][j][index];
			imag+=cos(dp)*xray.complexSymImag[i][j][index];
		}
	}
	//endProgram(__LINE__, __FILE__);
	//return sqrt(real*real+imag*imag);
	return real*real+imag*imag;
}

void calcIntensityHistogram(XRayStruct &xray, Vector &histogram)
{
	int max=20000000;
	int index=22;
	int Size;
	Real intensity, binSize=100.0, x=0;
	Vector intensities;

	for (int i=0;i<max;i++)
	{
		intensity=calcRandomIntensity(xray, index);
		SafePushBack(intensities, intensity, "intensities");
	}
	calcHistogram(intensities, histogram, binSize);
	Size=histogram.size();
	int nmax=10000;
	int end=min(nmax, Size);
	for (int i=0;i<end;i++)
	{
		cout <<x<<"\t"<<histogram[i]<<endl;
		x+=binSize;
	}
	for (int i=Size;i<nmax;i++)
	{
		cout <<x<<"\t"<<0<<endl;
		x+=binSize;

	}
}

void calcIntensityHistogram(ProteinStruct &Protein, XRayStruct &xray, Vector &histogram)
{
	calcSymmetryScattering(xray, Protein, 0);
	calcIntensityHistogram(xray, histogram);
}

void calcMinMaxIntensity(XRayStruct &xray, int index, Real &min, Real &max, Vector &amplitudes)
{
	int nprot=xray.complexSymReal.size();
	int nsym=xray.complexSymReal[0].size();
	int nmax, count=0;
	Real amplitude, largest;

	for (int i=0;i<nsym;i++)
	{
		for (int j=0;j<nprot;j++)
		{
			amplitude=xray.complexSymReal[j][i][index]*xray.complexSymReal[j][i][index];
			amplitude+=xray.complexSymImag[j][i][index]*xray.complexSymImag[j][i][index];
			amplitude=sqrt(amplitude);
			amplitudes[count]=amplitude;
			count++;
		}
	}
	nmax=findMax(amplitudes, largest);
	max=largest;
	min=largest;
	for (int i=0;i<nsym;i++)
	{
		if (i!=nmax)
		{
			max+=amplitudes[i];
			min-=amplitudes[i];
		}
	}
	if (min<0) min=0;
}

Real calcRotationPossibleScore(XRayStruct &xray, XRayStruct &expXRay)
{
	int npoint;
	//int npoint1, npoint2;
	int nprot=xray.complexSymReal.size();
	int nsym;
	Real score=0, min, max, scale;
	Vector amplitudes, meanI;

	//calcAmplitude(xray);
	//calcAmplitude(expXRay);
	npoint=xray.miller.size();
	//npoint1=xray.F.size();	
	//npoint2=expXRay.F.size();	

	if (nprot==0)
	{
		error("nprot= 0", __LINE__, __FILE__);
	}
/*
	if (npoint1!=npoint2)
	{
		string errorStr="npoint1= "+toStr(npoint1)+" npoint2= ";
		errorStr+=toStr(npoint2);
		error(errorStr, __LINE__, __FILE__);
	}
*/
	//calcRootMeanIntensity(xray, rootMeanI);
	calcMeanIntensity(xray, meanI);
	scale=calcScaleForSqrDiff(meanI, expXRay.i);
	scale=sqrt(scale);
	//scale=calcScaleForSqrDiff(xray.F, expXRay.F);
	nsym=xray.complexSymReal[0].size();
	SafeAlloc(amplitudes, nprot*nsym, "amplitudes");
	for (int i=0;i<npoint;i++)
	{
		calcMinMaxIntensity(xray, i, min, max, amplitudes);
		min*=scale;
		max*=scale;
		//cout <<min<<"\t"<<max<<"\t"<<expXRay.i[i]<<endl;
		if (expXRay.i[i]>max*max)
		{
			score+=(max*max-expXRay.i[i])*(max*max-expXRay.i[i]);
		}
		if (expXRay.i[i]<min*min)
		{
			score+=(min*min-expXRay.i[i])*(min*max-expXRay.i[i]);
		}
/*
		if (expXRay.i[i]>max)
		{
			score+=(max-expXRay.F[i])*(max-expXRay.F[i]);
		}
		if (expXRay.i[i]<min)
		{
			score+=(min-expXRay.F[i])*(min-expXRay.F[i]);
		}
*/
	}
	return score;
}

Real calcRotationPossibleScore(ProteinStruct &Protein, XRayStruct &xray, XRayStruct &expXRay)
{
	removeNonoverlapping(xray, expXRay);	
	calcSymmetryScattering(xray, Protein, 0);
	return calcRotationPossibleScore(xray, expXRay);

}
