#include <iostream>
#include "iso2d.h"
#include "iso2d_error.h"
#include <time.h>
#include <cstring>
#include<cstdio>
#include <cstdlib>
#if 1

int GenMesh2D(int nbpt, int nbelm, double *bpt, int* belm, 
			int *npt, double **gpt, int *nelm, int** elm)
{
	ISO2D::DTIso2D generator;
	clock_t startGenTime;
	clock_t endGenTime = 0.0, prevClock = 0.0;
	ISO2D::POINT minW,maxW;
	ISO2D::POINT minN = {-1., -1.}, maxN = {1., 1.};
	prevClock = startGenTime = clock();
	char fr2file[256], ba2file[256], pl2file[256], vtkfile[256];

	memset(ba2file, 0, sizeof(ba2file));
	strcpy(ba2file, "test");
	strcat(ba2file, ".ba2");

	memset(pl2file, 0, sizeof(pl2file));
	strcpy(pl2file, "test");
	strcat(pl2file, ".pl2");

	memset(vtkfile, 0, sizeof(vtkfile));
	strcpy(vtkfile, "test");
	strcat(vtkfile, ".vtk");

	generator.crtEnv();
	//	printf( "1-001-1\n" );
	if (generator.readFr2(nbpt, nbelm, bpt, belm) == false)
	//if(generator.readFr2("test0.fr2") == false)
	{		
		return 1;
	}
	
#ifdef _3D_DECOM_
	generator.readBa3(ba2file);
#else
	generator.readBa2(ba2file);
#endif

	generator.calcBox(&minW, &maxW);
	generator.scaleFactor(minW, maxW, minN, maxN);
	generator.scaleBoxPnt();
	generator.setupInitTri();

#if 0 /* no scale */	
	generator.scaGeom();
	generator.scaBkGrnd();
#endif
	generator.calcDens();
	printf("Initialization time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();

	generator.bndPntInst();
	printf("Boundary insertion time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();

	generator.recoverBnds();
	
	generator.clrOuterEles();
	printf("Boundary recover time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();
	generator.writevtk("test1.vtk");
	while (generator.innerPntInst()) ;
	printf("Field point insertion time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();
	generator.writevtk("test2.vtk");
	generator.smooth();
	printf("Smooth time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();

	generator.rmvEmpNods();
	generator.rmvEmpEles();	
	//	generator.updateBndParent();
#if 0
	generator.rescaGeom();
#endif /* no scale */
	//	printf( "1-006\n" );
	generator.writePl2(npt, gpt, nelm, elm);
	//generator.writevtk(vtkfile);
	generator.output();
	printf("Create Points:%d. Accepted Points:%d.\n", ISO2D::g_nCreatePnts, ISO2D::g_nAcceptPnts);

	printf("Postprocessing (including output) time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();

	printf("Total time: %fs.\n", ((double)(clock()-startGenTime))/CLOCKS_PER_SEC);

	return 0;
}
#else
int main(int argc, char **argv)
{
	ISO2D::DTIso2D generator;
	clock_t startGenTime;
	clock_t endGenTime = 0.0, prevClock = 0.0;
	ISO2D::POINT minW,maxW;
	ISO2D::POINT minN = {-1., -1.}, maxN = {1., 1.};
	prevClock = startGenTime = clock();
	char fr2file[256], ba2file[256], pl2file[256], vtkfile[256];
	
	memset(fr2file, 0, sizeof(fr2file));
	strcpy(fr2file, argv[1]);
	strcat(fr2file, ".fr2");
	
	memset(ba2file, 0, sizeof(ba2file));
	strcpy(ba2file, argv[1]);
	strcat(ba2file, ".ba2");
	
	memset(pl2file, 0, sizeof(pl2file));
	strcpy(pl2file, argv[1]);
	strcat(pl2file, ".pl2");
	
	memset(vtkfile, 0, sizeof(vtkfile));
	strcpy(vtkfile, argv[1]);
	strcat(vtkfile, ".vtk");


	generator.crtEnv();
	//	printf( "1-001-1\n" );
	if (generator.readFr2(fr2file) == false)
	{		
		return 1;
	}
#ifdef _3D_DECOM_
	generator.readBa3(ba2file);
#else
	generator.readBa2(ba2file);
#endif
	
	generator.calcBox(&minW, &maxW);
	generator.scaleFactor(minW, maxW, minN, maxN);
	generator.scaleBoxPnt();
	generator.setupInitTri();
	
#if 0 /* no scale */	
	generator.scaGeom();
	generator.scaBkGrnd();
#endif
	generator.calcDens();
	printf("Initialization time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();

	generator.bndPntInst();
	printf("Boundary insertion time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();

	generator.recoverBnds();
	generator.clrOuterEles();
	printf("Boundary recover time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();

	while (generator.innerPntInst()) ;
	printf("Field point insertion time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();

	generator.smooth();
	printf("Smooth time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();

	generator.rmvEmpNods();
	generator.rmvEmpEles();	
//	generator.updateBndParent();
#if 0
	generator.rescaGeom();
#endif /* no scale */
	//	printf( "1-006\n" );
	generator.writePl2(pl2file);
	generator.writevtk(vtkfile);
	generator.output();
	printf("Create Points:%d. Accepted Points:%d.\n", ISO2D::g_nCreatePnts, ISO2D::g_nAcceptPnts);
	
	printf("Postprocessing (including output) time: %fs.\n", ((double)(clock()-prevClock))/CLOCKS_PER_SEC);
	prevClock = clock();

	printf("Total time: %fs.\n", ((double)(clock()-startGenTime))/CLOCKS_PER_SEC);
	
	return 0;
}
#endif
