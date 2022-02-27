#FFLAGS=-O3 -g -ffixed-line-length-132
FFLAGS=-O3 -g -Mextend
#CFLAGS = -lm -O3 -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
CFLAGS = -O3 -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -I/opt/local/include
#CFLAGS = -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk/usr/lib/crt1.o /usr/lib/crt1.o

F77= pgfortran
CC = cc
#CC = pgcc

.c.o: 
	${CC}  $(CFLAGS) -c $*.c

.f.o: 
	${F77}  $(FFLAGS) -c $*.f

CatMergeMagType: CatMergeMagType.o Sphere.o DateTime.o CatUtil.o
	${CC} CatMergeMagType.c Sphere.c DateTime.c CatUtil.c -o CatMergeMagType

CatMerge: CatMerge.o Sphere.o DateTime.o CatUtil.o
	${CC} CatMerge.c Sphere.c DateTime.c CatUtil.c -o CatMerge

GetGeolFormation: GetGeolFormation.o Sphere.o DateTime.o CatUtil.o
	${CC} GetGeolFormation.c Sphere.c DateTime.c CatUtil.c -o GetGeolFormation

PointSelect: PointSelect.o Sphere.o DateTime.o CatUtil.o
	${CC} PointSelect.c Sphere.c DateTime.c CatUtil.c -o PointSelect

CatSelect: CatSelect.o Sphere.o DateTime.o CatUtil.o
	${CC} CatSelect.c Sphere.c DateTime.c CatUtil.c -o CatSelect

Cat2Mw: Cat2Mw.o Sphere.o DateTime.o CatUtil.o
	${CC} Cat2Mw.c Sphere.c DateTime.c CatUtil.c -o Cat2Mw

CatTbl2Prime: CatTbl2Prime.o Sphere.o DateTime.o CatUtil.o
	${CC} CatTbl2Prime.o Sphere.o DateTime.o CatUtil.o -o CatTbl2Prime

Cat2Prime: Cat2Prime.o Sphere.o DateTime.o CatUtil.o
	${CC} Cat2Prime.c Sphere.c DateTime.c CatUtil.c -o Cat2Prime

CMTMerge: CMTMerge.o Sphere.o DateTime.o CatUtil.o
	${CC} CMTMerge.o Sphere.o DateTime.o CatUtil.o -o CMTMerge

ISCMechMerge: ISCMechMerge.o Sphere.o DateTime.o CatUtil.o
	${CC} ISCMechMerge.o Sphere.o DateTime.o CatUtil.o -o ISCMechMerge

SocalMech: SocalMech.o Sphere.o DateTime.o CatUtil.o
	${CC} SocalMech.o Sphere.o DateTime.o CatUtil.o -o SocalMech

CatInterface: CatInterface.o Sphere.o DateTime.o CatUtil.o
	${CC} CatInterface.o Sphere.o DateTime.o CatUtil.o -o CatInterface

CMTMergeTbl: CMTMergeTbl.o Sphere.o DateTime.o CatUtil.o
	${CC} CMTMergeTbl.o Sphere.o DateTime.o CatUtil.o -o CMTMergeTbl

CatSort: CatSort.o Sphere.o DateTime.o CatUtil.o
	${CC} CatSort.c -o CatSort DateTime.o CatUtil.o

CatSmooth: CatSmooth.o DateTime.o CatUtil.o Sphere.o
	${CC} CatSmooth.o -o CatSmooth DateTime.o CatUtil.o Sphere.o

CatDecluster: CatDecluster.o Sphere.o DateTime.o CatUtil.o
	${CC} CatDecluster.c -o CatDecluster Sphere.o DateTime.o CatUtil.o

ASCE-matchdisag: ASCE-matchdisag.o Sphere.o DateTime.o CatUtil.o
	${CC} ASCE-matchdisag.c -o ASCE-matchdisag Sphere.o DateTime.o CatUtil.o

RupTime: RupTime.o Sphere.o 
	${CC} RupTime.c -o RupTime Sphere.o 

IRIS: IRIS.o Sphere.o DateTime.o CatUtil.o
	${CC} IRIS.c -o IRIS Sphere.o DateTime.o CatUtil.o

Cat2Breq: Cat2Breq.o Sphere.o DateTime.o CatUtil.o Station.o PStimes.o
	${F77} Cat2Breq.c -o Cat2Breq Sphere.o DateTime.o CatUtil.o PStimes.o iaspei-tau/src/libiasp.a libget.a Station.o -Mnomain

Kiknet2Sac: Kiknet2Sac.o Seis.o CatUtil.o DateTime.o dcops.o  cfft.o
	${CC} -lm -O3 Kiknet2Sac.o -o Kiknet2Sac Seis.o CatUtil.o DateTime.o dcops.o cfft.o

SacRatio: SacRatio.o Seis.o CatUtil.o DateTime.o dcops.o  cfft.o
	${CC} -lm -O3 SacRatio.o -o SacRatio Seis.o CatUtil.o DateTime.o dcops.o cfft.o

DisToLines: DisToLines.o Sphere.o DateTime.o CatUtil.o
	${CC} DisToLines.c Sphere.c DateTime.c CatUtil.c -o DisToLines

DisToLineCDF: DisToLineCDF.o Sphere.o DateTime.o CatUtil.o CDFUtil.o
	${CC} DisToLineCDF.o Sphere.o DateTime.o CatUtil.o CDFUtil.o -L/opt/lib -lnetcdf -o DisToLineCDF

libCDFUtil: CDFUtil.o
	libtool -static -o libCDFUtil.a CDFUtil.o

JSHIS2pfdha: JSHIS2pfdha.o Sphere.o DateTime.o CatUtil.o
	${CC} JSHIS2pfdha.c Sphere.c DateTime.c CatUtil.c -o JSHIS2pfdha

cont2grid_v1.2: cont2grid_v1.2.c Sphere.c spline.c
	${CC} cont2grid_v1.2.c -o cont2grid_v1.2 Sphere.c spline.c

cont2grid_v1.3: cont2grid_v1.3.c Sphere.c spline.c
	${CC} cont2grid_v1.3.c -o cont2grid_v1.3 Sphere.c spline.c

cont2grid_v1.4: cont2grid_v1.4.c Sphere.c spline.c
	${CC} cont2grid_v1.4.c -o cont2grid_v1.4 Sphere.c spline.c

cont2grid_v1.5: cont2grid_v1.5.c Sphere.c spline.c
	${CC} cont2grid_v1.5.c -o cont2grid_v1.5 Sphere.c spline.c

Resample_Invall-v1.0: Resample_Invall-v1.0.c 
	${CC} Resample_Invall-v1.0.c -o Resample_Invall-v1.0 

