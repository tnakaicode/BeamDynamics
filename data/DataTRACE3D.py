import sys
sys.path.append('../lib/')

from numpy import *
from matplotlib.pyplot import *
from numpy.linalg import inv
from .ConversionFunctions import *
from .ellipse import *

Path0 = 'sigma_trace3d/'

Sinjection = matrix([
[2.8834 , 0.000 , 0.000 , 0.0000 , 0.000 , 0.000],
[3.4922 ,-0.154 , 0.000 , 0.0000 , 0.000 , 0.000],
[3.7727 , 0.000 , 0.000 , 0.0000 , 0.000 , 0.000],
[5.5126 , 0.000 , 0.000 ,-0.9000 , 0.000 , 0.000],
[5.8750 , 0.000 , 0.000 , 0.0000 , 0.000 , 0.000],
[1.2536 , 0.000 , 0.000 , 0.0000 , 0.000 , 0.984]],float)

SinjectionLC = matrix([
[2.9989 , 0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.0000],
[3.4066 ,-0.2140 , 0.0000 , 0.0000 , 0.0000 , 0.0000],
[3.5736 , 0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.0000],
[5.8419 , 0.0000 , 0.0000 ,-0.9000 , 0.0000 , 0.0000],
[5.5170 , 0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.0000],
[1.1617 , 0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.9790]],float)

Bend90 = matrix([
[6.2429 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[3.7663 , 0.906 , 0.000 , 0.000 , 0.000 , 0.000],
[13.8997, 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[13.4632, 0.000 , 0.000 , 0.982 , 0.000 , 0.000],
[18.8279, 0.000 , 0.000 ,-0.869 ,-0.932 , 0.000],
[1.2944 , 0.000 , 0.000 ,-0.921 ,-0.960 , 0.989]],float)



# Final Modified Sigma Matrices
#==================================================================== 0000 A
S0000 = matrix([
[5.9546 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[3.7470 , 0.895 , 0.000 , 0.000 , 0.000 , 0.000],
[4.9859 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[5.4023 , 0.000 , 0.000 , 0.942 , 0.000 , 0.000],
[24.8496, 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[1.2935 , 0.000 , 0.000 , 0.000 , 0.000 , 0.999]],float)

S0000LC = matrix([
[5.3088 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[3.4066 , 0.835 , 0.000 , 0.000 , 0.000 , 0.000],
[5.6937 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[5.8423 , 0.000 , 0.000 , 0.962 , 0.000 , 0.000],
[22.7092, 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[1.1616 , 0.000 , 0.000 , 0.000 , 0.000 , 0.999]],float)

#==================================================================== 1600 A
S1600 =  matrix([
[5.7018 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[3.9816 , 0.899 , 0.000 , 0.000 , 0.000 , 0.000],
[5.2000 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[5.4841 , 0.000 , 0.000 , 0.894 , 0.000 , 0.000],
[25.1032, 0.000 , 0.000 , 0.188 , 0.473 , 0.000],
[1.2938 , 0.000 , 0.000 , 0.215 , 0.493 , 0.999]],float)

#S1600NG=  matrix([
#[6.0420 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
#[3.7500 , 0.898 , 0.000 , 0.000 , 0.000 , 0.000],
#[5.1966 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
#[5.9436 , 0.000 , 0.000 , 0.912 , 0.000 , 0.000],
#[25.1022, 0.000 , 0.000 , 0.191 , 0.450 , 0.000],
#[1.2937 , 0.000 , 0.000 , 0.217 , 0.470 , 0.999]],float)

S1600NG=  matrix([
[6.3816 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[3.7599 , 0.910 , 0.000 , 0.000 , 0.000 , 0.000],
[5.6911 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[5.9285 , 0.000 , 0.000 , 0.915 , 0.000 , 0.000],
[26.2963, 0.000 , 0.000 , 0.170 , 0.448 , 0.000],
[1.2944 , 0.000 , 0.000 , 0.199 , 0.472 , 0.999]],float)

S1600LC=  matrix([
[5.0726 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[3.6450 , 0.843 , 0.000 , 0.000 , 0.000 , 0.000],
[5.8895 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[5.6828 , 0.000 , 0.000 , 0.921 , 0.000 , 0.000],
[22.9414, 0.000 , 0.000 , 0.137 , 0.400 , 0.000],
[1.1616 , 0.000 , 0.000 , 0.171 , 0.429 , 0.998]],float)

S1600NGH=  matrix([
[6.0681 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[4.5034 , 0.794 , 0.000 , 0.000 , 0.000 , 0.000],
[5.2285 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[5.4055 , 0.000 , 0.000 , 0.945 , 0.000 , 0.000],
[25.1140, 0.146 , 0.594 , 0.000 , 0.000 , 0.000],
[1.2937 , 0.186 , 0.621 , 0.000 , 0.000 , 0.998]],float)

#==================================================================== 3200 A

S3120 =  matrix([
[5.5371 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[4.5848 , 0.920 , 0.000 , 0.000 , 0.000 , 0.000],
[5.4908 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[6.0884 , 0.000 , 0.000 , 0.789 , 0.000 , 0.000],
[25.2144, 0.000 , 0.000 , 0.357 , 0.787 , 0.000],
[1.2941 , 0.000 , 0.000 , 0.409 , 0.813 , 0.998]],float)

S3120 =  matrix([
[5.5371 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[4.5848 , 0.920 , 0.000 , 0.000 , 0.000 , 0.000],
[5.4908 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[6.0884 , 0.000 , 0.000 , 0.789 , 0.000 , 0.000],
[25.2144, 0.000 , 0.000 , 0.357 , 0.787 , 0.000],
[1.2941 , 0.000 , 0.000 , 0.409 , 0.813 , 0.998]],float)

S3120NG =  matrix([
[6.4777 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[3.7630 , 0.913 , 0.000 , 0.000 , 0.000 , 0.000],
[6.0167 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[7.1834 , 0.000 , 0.000 , 0.856 , 0.000 , 0.000],
[26.5075, 0.000 , 0.000 , 0.330 , 0.722 , 0.000],
[1.2946 , 0.000 , 0.000 , 0.384 , 0.757 , 0.998]],float)

#==================================================================== 4450 A

S4450 =  matrix([
[6.0062 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[6.5468 , 0.967 , 0.000 , 0.000 , 0.000 , 0.000],
[6.2312 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[6.7172 , 0.000 , 0.000 , 0.598 , 0.000 , 0.000],
[26.1901, 0.000 , 0.000 , 0.546 , 0.962 , 0.000],
[1.2951 , 0.000 , 0.000 , 0.662 , 0.961 , 0.995]],float)

S4450NG =  matrix([
[6.0062 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[6.5468 , 0.967 , 0.000 , 0.000 , 0.000 , 0.000],
[6.2312 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000],
[6.7172 , 0.000 , 0.000 , 0.598 , 0.000 , 0.000],
[26.1901, 0.000 , 0.000 , 0.546 , 0.962 , 0.000],
[1.2951 , 0.000 , 0.000 , 0.662 , 0.961 , 0.995]],float)


SigmaInj,dS=ConvertT3D(Sinjection); SigmaInjLC,dS=ConvertT3D(SinjectionLC)
Sigma0000,dS = ConvertT3D(S0000); Sigma0000LC,dS = ConvertT3D(S0000LC)
Sigma1600,dS = ConvertT3D(S1600); Sigma1600LC,dS = ConvertT3D(S1600LC); Sigma1600NG,dS = ConvertT3D(S1600NG)
Sigma3120,dS = ConvertT3D(S3120); Sigma3120NG,dS = ConvertT3D(S3120NG);
Sigma4450,dS = ConvertT3D(S4450); Sigma4450NG,dS = ConvertT3D(S4450NG);
SigmaBend90,dS = ConvertT3D(Bend90);

savetxt(Path0+'SigmaInjection.dat',SigmaInj); savetxt(Path0+'SigmaInjectionLC.dat',SigmaInjLC)

savetxt(Path0+'Trace3DSigma_I_0.dat',Sigma0000);
savetxt(Path0+'Trace3DSigma_I_0LC.dat',Sigma0000LC)

savetxt(Path0+'Trace3DSigma_I_1600.dat',Sigma1600);
savetxt(Path0+'Trace3DSigma_I_1600LC.dat',Sigma1600LC)
savetxt(Path0+'Trace3DSigma_I_1600NG.dat',Sigma1600NG)

savetxt(Path0+'Trace3DSigma_I_3120.dat',Sigma3120);
savetxt(Path0+'Trace3DSigma_I_3120NG.dat',Sigma3120NG);

savetxt(Path0+'Trace3DSigma_I_4450.dat',Sigma4450);
savetxt(Path0+'Trace3DSigma_I_4450NG.dat',Sigma4450NG);


savetxt(Path0+'Trace3DSigmaBend90.dat',SigmaBend90);


if True: 
	figure(1,figsize=(8,8));
	E0 = ellipse(Sigma0000)
	E1 = ellipse(Sigma0000LC)
	M = E0.MismatchFactor(E1,Type=1)

	subplot(2,2,1)
	E0.PlotXX1()
	E1.PlotXX1()
	text(0,0,'M=%0.4f' % M[1],va='center',ha='center',color='r',size=16)
	legend((r'1.000 mA','0.001 mA'),loc=2)

	subplot(2,2,2)
	E0.PlotYY1()
	E1.PlotYY1()
	text(0,0,'M=%0.4f' % M[1],va='center',ha='center',color='r',size=16)

	subplot(2,2,3)
	E0.PlotZZ1()
	E1.PlotZZ1()
	text(0,0,'M=%0.4f' % M[2],va='center',ha='center',color='r',size=16)

	subplot(2,2,4)
	E0.PlotXY()
	E1.PlotXY()
	text(0,0,'M=%0.4f' % M[3],va='center',ha='center',color='r',size=16)

	suptitle(r'Space Charge Effects: 1mA versus 1$\mu$A',size=16)


show()
