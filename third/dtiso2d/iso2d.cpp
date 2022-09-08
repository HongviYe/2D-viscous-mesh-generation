/* *************************************************************
 * Delaunay triangulation in 2-dimensions
 *      with automatic point creation
 *
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 04, 27
 * 
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 * **************************************************************/
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <fstream>
#include <time.h>
#include <iostream>
#include <cstring>
#include "../geom/geom_func.h"
#include <cstdlib>
#include "iso2d.h"
#include "iso2d_error.h"
#pragma optimize("",off)
using namespace std;
extern double insphere(double *pa, double *pb, double *pc, double *pd, double *pe);
extern double orient2d(double *pa, double *pb, double *pc);
extern double orient3d(double *pa, double *pb, double *pc, double *pd);

extern "C" {
	extern double exactinit();
	extern double incircle(double*, double*, double*, double*);
}
namespace ISO2D
{

POINT g_cors[] = 
{
	{-4.0, -4.0}, {-4.0, 4.0}, {4.0, -4.0}, {4.0, 4.0}
};
/*
 * global variables
 * 0.75 for yao
 */
REAL g_alpha = 0.70228615;

INTEGER g_nCreatePnts = 0;
INTEGER g_nAcceptPnts = 0;
/*
 * global functions
 */
static REAL Min(REAL a, REAL b)
{
	return a<b ? a : b;
}

static REAL Area2(POINT a, POINT b, POINT c)
{
	return GEOM_FUNC::orient2d(b, c, a);// (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
}
                                                                                                                                               
/* check whether c is on the left side of ab
 */
static bool Left(POINT a, POINT b, POINT c)
{
    return Area2(a, b, c) > 0;
}
                                                                                                                                               
/* check whether a,b,c is collinear
 */
static bool Collinear(POINT a, POINT b, POINT c)
{
    return Area2(a, b, c) == 0;
}

/* Exclusive or: T iff exactly one argument is true.
 */
static bool Xor(bool x, bool y)
{
    /*The arguments are negated to ensure that they are 0/1 values.*/
    return !x ^ !y;
}
                                                                                 
/* check whether ab intersectes with cd properly
 */
static bool IntersectProp(POINT a, POINT b, POINT c, POINT d)
{
    /*Eliminate improper cases. */
    if (Collinear(a, b, c) || Collinear(a, b, d) || Collinear(c, d, a) || Collinear(c, d, b))
        return false;
    return Xor(Left(a, b, c), Left(a, b, d)) && Xor(Left(c, d, a), Left(c, d, b));
}
                                  
DTIso2D::DTIso2D()
{
	int i = 0;
	m_pElems = NULL;
	m_pNodes = NULL;
	m_pBnds = NULL;
	m_nElems = 0;
	m_nNodes = 0;
	m_nBnds = 0;
	m_nAllocElems = 0;
	m_nAllocNodes = 0;
	m_nAllocBnds = 0;
	m_scale = 1.0;
	for (i = 0; i < DIM; i++)
	{
		m_cenW[i] = m_cenN[i] = 0.0;
		m_minW[i] = m_maxW[i] = 0.0;
		m_minN[i] = m_maxN[i] = 0.0;
	}

	m_nBkGrndElem = 0;
	m_nBkGrndNode = 0;
	m_nPntSrcNum = 0;
	m_nLneSrcNum = 0;
	m_nTriSrcNum = 0;
	m_pntSources = NULL;
	m_lineSources = NULL;
	m_triSources = NULL;
		
	m_nLocInd = 0;
	m_nCurInd = 0;
	GEOM_FUNC::exactinit();
	
	g_cors[0][0] = g_cors[0][1] = g_cors[1][0] = g_cors[2][1] = -4.0;
	g_cors[1][1] = g_cors[2][0] = g_cors[3][0] = g_cors[3][1] = 4.0;
}

DTIso2D::~DTIso2D()
{
	if (m_pElems)
	{
		free(m_pElems);
		m_pElems = NULL;
	}
	if (m_pNodes)
	{
		free(m_pNodes);
		m_pNodes = NULL;
	}
	if (m_pBnds)
	{
		free(m_pBnds);
		m_pBnds = NULL;
	} //  [12/16/2005]
}

bool DTIso2D::crtEnv()
{
	m_pElems = (Elem*)malloc(sizeof(Elem) * INIT_ALLOC_ELE_NUM);
	if (NULL == m_pElems)
	{
			printf( "Not enough memory...\n" );
			exit(1);
			goto CLEAR;
	}		
	m_pNodes = (Node*)malloc(sizeof(Node) * INIT_ALLOC_NOD_NUM);
	if (NULL == m_pNodes)
	{
			printf( "Not enough memory...\n" );
			exit(1);
			goto CLEAR;
	}		
	m_pBnds = (Bnd*)malloc(sizeof(Bnd) * INIT_ALLOC_BND_NUM);
	if (NULL == m_pBnds)
	{
			printf( "Not enough memory...\n" );
			exit(1);
			goto CLEAR;
	}		
	
	m_nAllocElems = INIT_ALLOC_ELE_NUM;
	m_nAllocNodes = INIT_ALLOC_NOD_NUM;
	m_nAllocBnds = INIT_ALLOC_BND_NUM;
	
	/* setup initial triangulation */
//	setupInitTri();
	
	return true;
CLEAR:
	if (m_pElems)
	{
		free(m_pElems);
		m_pElems = NULL;
	}
	if (m_pNodes)
	{
		free(m_pNodes);
		m_pNodes = NULL;
	}
	if (m_pBnds)
	{
		free(m_pBnds);
		m_pBnds = NULL;
	}
	return false;
}

bool DTIso2D::readFr2(const char* fname)
{
	FILE* fp = NULL;
	INTEGER iTok = 0;
	INTEGER i;
	if (fname)
	{
		fp = fopen(fname, "r");
		if (!fp)
		{
			printf("cannot read file %s\n", fname);
			return false;
		}
		fscanf(fp, "%d%d%d%d%d%d%d", &m_nBnds, &iTok, &iTok, &iTok, &iTok, &iTok, &iTok);
		if (m_nBnds==0)
		{
			fclose(fp);
			return false;
		}
		if (m_nBnds > INIT_ALLOC_BND_NUM)
		{
			printf("Error: Not enough memory for input boundaries!\n");
			fclose(fp);
			return false;
		}
		if (m_nBnds + INIT_NOD_NUM > INIT_ALLOC_NOD_NUM)
		{
			printf("Error: Not enough memory for input nodes!\n");
			fclose(fp);
			return false;
		}
		for (i = 0; i < m_nBnds; i++)
		{
	#ifdef _DOUBLE
			fscanf(fp, "%d%lf%lf", &iTok, &((m_pNodes)[i + INIT_NOD_NUM].pt)[0],
				&((m_pNodes)[i + INIT_NOD_NUM].pt)[1]);
	#else
			fscanf(fp, "%d%f%f", &iTok, &((m_pNodes)[i + INIT_NOD_NUM].pt)[0],
				&((m_pNodes)[i + INIT_NOD_NUM].pt)[1]);
	#endif
		} 
		m_nNodes = m_nBnds + INIT_NOD_NUM;
		for (i = 0; i < m_nBnds; i++)
		{
			fscanf(fp, "%d%d%d%d%d", &iTok, &(m_pBnds)[i].beg, &(m_pBnds)[i].end,
				&(m_pBnds)[i].curve, &(m_pBnds)[i].loop);
			//fscanf(fp, "%d%d%d%d", &iTok, &(m_pBnds)[i].beg, &(m_pBnds)[i].end,&(m_pBnds)[i].curve);
			//m_pBnds[i].curve = 0;
// 			m_pBnds[i].loop = 0;
// 			m_pBnds[i].ele = NULL_ELEM;
		}
		fclose(fp);
		return true;
	}
	return false;
}

bool DTIso2D::readFr2(int nbpt, int nbelm, double *bpt, int* belm)
{
	FILE* fp = NULL;
	INTEGER iTok = 0;
	INTEGER i;

	m_nBnds = nbpt;

	if (m_nBnds > INIT_ALLOC_BND_NUM)
	{
		printf("Error: Not enough memory for input boundaries!\n");
		return false;
	}
	if (m_nBnds + INIT_NOD_NUM > INIT_ALLOC_NOD_NUM)
	{
			printf("Error: Not enough memory for input nodes!\n");
		return false;
	}
	for (i = 0; i < m_nBnds; i++)
	{
		((m_pNodes)[i + INIT_NOD_NUM].pt)[0] = bpt[2*i +0];
		((m_pNodes)[i + INIT_NOD_NUM].pt)[1] = bpt[2*i +1];
	} 
	m_nNodes = m_nBnds + INIT_NOD_NUM;

	for (i = 0; i < m_nBnds; i++)
	{
		(m_pBnds)[i].beg = belm[2*i+0] + 1;
		(m_pBnds)[i].end = belm[2*i+1] + 1;
		(m_pBnds)[i].curve = 1;
		(m_pBnds)[i].loop = 1;
	}
	
	return true;
}

bool DTIso2D::readBa2(const char* fname)
{
	FILE* fp = NULL;
	INTEGER iTok = 0;
	INTEGER i,j,k;
	char cTok[100];
	REAL area;
	if (fname)
	{
		fp = fopen(fname, "r");
		if (!fp)
		{
			printf("Warning: cannot read file %s\n", fname);
			return false;
		}
		fgets(cTok,100,fp);
		fscanf(fp, "%d%d\n", &m_nBkGrndNode, &m_nBkGrndElem);
		m_pBkGrndNode = (Node *)malloc(sizeof(Node)*m_nBkGrndNode);
		for (i=0; i<m_nBkGrndNode; i++)
		{
			//read node information
			fscanf(fp, "%d", &iTok);
			for (j=0; j<DIM; j++)
			{
#ifdef _DOUBLE
				fscanf(fp, "%lf", &(m_pBkGrndNode[i].pt[j]));
#else
				fscanf(fp, "%f", &(m_pBkGrndNode[i].pt[j]));
#endif
			}
#ifdef _DOUBLE
			fscanf(fp, "%lf", &(m_pBkGrndNode[i].density));
#else
			fscanf(fp, "%f", &(m_pBkGrndNode[i].density));
#endif
			fgets(cTok,100,fp);
		}
		m_pBkGrndElem = (Elem *)malloc(sizeof(Elem)*m_nBkGrndElem);
		for (i=0; i<m_nBkGrndElem; i++)
		{
			fscanf(fp, "%d", &iTok);
			for (j=0; j<=DIM; j++)
			{
				fscanf(fp, "%d", &(m_pBkGrndElem[i].form[j]));
			}
			area = Area2(m_pBkGrndNode[m_pBkGrndElem[i].form[0]-1].pt,
				m_pBkGrndNode[m_pBkGrndElem[i].form[1]-1].pt,
				m_pBkGrndNode[m_pBkGrndElem[i].form[2]-1].pt);
			if (area < 0)
			{
				k = m_pBkGrndElem[i].form[0];
				m_pBkGrndElem[i].form[0] = m_pBkGrndElem[i].form[DIM];
				m_pBkGrndElem[i].form[DIM] = k;
			}
			fgets(cTok,100,fp);
		}
		fgets(cTok,100,fp);
		fscanf(fp, "%d%d%d\n",&m_nPntSrcNum,&m_nLneSrcNum,&iTok);
		fgets(cTok,100,fp);
		m_pntSources = (PointSource *)malloc(sizeof(PointSource)*m_nPntSrcNum);
		for (i=0; i<m_nPntSrcNum; i++)
		{
			fgets(cTok,100,fp);
			for (j=0; j<DIM; j++)
			{
#ifdef _DOUBLE
				fscanf(fp, "%lf", &m_pntSources[i].pt[j]);
#else
				fscanf(fp, "%f", &m_pntSources[i].pt[j]);
#endif
			}
#ifdef _DOUBLE
			fscanf(fp, "%lf%lf%lf\n", &m_pntSources[i].rIntensity,&m_pntSources[i].rInnerRad,&m_pntSources[i].rOuterRad);
#else
			fscanf(fp, "%f%f%f\n", &m_pntSources[i].rIntensity,&m_pntSources[i].rInnerRad,&m_pntSources[i].rOuterRad);
#endif
			assert(m_pntSources[i].rOuterRad-m_pntSources[i].rInnerRad > EPS_ZERO_SQ);
		}
		m_lineSources = (LineSource *)malloc(sizeof(LineSource)*m_nLneSrcNum);
		fgets(cTok,100,fp);
		for (i=0; i<m_nLneSrcNum; i++)
		{
			fgets(cTok,100,fp);
			for (j=0; j<2; j++)
			{
				for (k=0; k<DIM; k++)
				{
#ifdef _DOUBLE
					fscanf(fp, "%lf", &m_lineSources[i].points[j].pt[k]);
#else
					fscanf(fp, "%f", &m_lineSources[i].points[j].pt[k]);
#endif
				}
#ifdef _DOUBLE
				fscanf(fp, "%lf%lf%lf\n", &m_lineSources[i].points[j].rIntensity, &m_lineSources[i].points[j].rInnerRad,
					&m_lineSources[i].points[j].rOuterRad);
#else
				fscanf(fp, "%f%f%f\n", &m_lineSources[i].points[j].rIntensity, &m_lineSources[i].points[j].rInnerRad,
					&m_lineSources[i].points[j].rOuterRad);
#endif
				assert(m_lineSources[i].points[j].rOuterRad-m_lineSources[i].points[j].rInnerRad > EPS_ZERO_SQ);
			}
			assert(sqrt(distSqua(m_lineSources[i].points[0].pt,m_lineSources[i].points[1].pt)) > EPS_ZERO_SQ);
		}
		fclose(fp);
		return true;
	}
	return false;
}

bool DTIso2D::readBa3(const char* fname)
{
	FILE* fp = NULL;
	INTEGER iTok = 0;
	INTEGER i,j,k;
	char cTok[100];
	REAL x, y, z;

	if (fname)
	{
		fp = fopen(fname, "r");
		if (!fp)
		{
			printf("cannot read file %s\n", fname);
			return false;
		}
		fgets(cTok,100,fp);
		fscanf(fp, "%d%d%d%d%d\n", &m_nBkGrndNode, &m_nBkGrndElem, 
			&m_nPntSrcNum, &m_nLneSrcNum, &m_nTriSrcNum);
		if (m_nBkGrndNode > 0)
			m_pBkGrndNode = (Node *)malloc(sizeof(Node)*m_nBkGrndNode);
		else
			m_pBkGrndNode = NULL;

		for (i=0; i<m_nBkGrndNode; i++)
		{
			//read node information
			fscanf(fp, "%d", &iTok);
			for (j=0; j<3; j++)
			{
#ifdef _DOUBLE
				fscanf(fp, "%lf", &(m_pBkGrndNode[i].pt[j]));
#else
				fscanf(fp, "%f", &(m_pBkGrndNode[i].pt[j]));
#endif
			}
			fscanf(fp, "%lf%lf%lf%lf", &x, &y, &z, &m_pBkGrndNode[i].density);
			fscanf(fp, "%lf%lf%lf%lf", &x, &y, &z, &m_pBkGrndNode[i].density);
			fscanf(fp, "%lf%lf%lf%lf", &x, &y, &z, &m_pBkGrndNode[i].density);
// 			fgets(cTok,100,fp);
// 			fgets(cTok,100,fp);
// 			fgets(cTok,100,fp);
// 			fgets(cTok,100,fp);
		}
		if (m_nBkGrndElem > 0)
			m_pBkGrndElem = (Elem *)malloc(sizeof(Elem)*m_nBkGrndElem);
		else
			m_pBkGrndElem = NULL;

		fgets(cTok, 100, fp);
		for (i=0; i<m_nBkGrndElem; i++)
		{
//			fscanf(fp, "%d", &iTok);
//			for (j=0; j<=DIM; j++)
//			{
//				fscanf(fp, "%d", &(m_pBkGrndElem[i].form[j]));
//			}
//			area = ::Area2(m_pBkGrndNode[m_pBkGrndElem[i].form[0]-1].pt,
//				m_pBkGrndNode[m_pBkGrndElem[i].form[1]-1].pt,
//				m_pBkGrndNode[m_pBkGrndElem[i].form[2]-1].pt);
//			if (area < 0)
//			{
//				k = m_pBkGrndElem[i].form[0];
//				m_pBkGrndElem[i].form[0] = m_pBkGrndElem[i].form[DIM];
//				m_pBkGrndElem[i].form[DIM] = k;
//			} // skip the element information
			fgets(cTok,100,fp);
		}
		fgets(cTok,100,fp);
		if (m_nPntSrcNum > 0)
			m_pntSources = (PointSource *)malloc(sizeof(PointSource)*m_nPntSrcNum);
		for (i=0; i<m_nPntSrcNum; i++)
		{
			fgets(cTok,100,fp);
			for (j=0; j<3; j++)
			{
#ifdef _DOUBLE
				fscanf(fp, "%lf", &m_pntSources[i].pt[j]);
#else
				fscanf(fp, "%f", &m_pntSources[i].pt[j]);
#endif
			}
#ifdef _DOUBLE
			fscanf(fp, "%lf%lf%lf\n", &m_pntSources[i].rIntensity,&m_pntSources[i].rInnerRad,&m_pntSources[i].rOuterRad);
#else
			fscanf(fp, "%f%f%f\n", &m_pntSources[i].rIntensity,&m_pntSources[i].rInnerRad,&m_pntSources[i].rOuterRad);
#endif
			assert(m_pntSources[i].rOuterRad-m_pntSources[i].rInnerRad > EPS_ZERO_SQ);
		}
		if (m_nLneSrcNum > 0)
			m_lineSources = (LineSource *)malloc(sizeof(LineSource)*m_nLneSrcNum);
		fgets(cTok,100,fp);
		for (i=0; i<m_nLneSrcNum; i++)
		{
			fgets(cTok,100,fp);
			for (j=0; j<2; j++)
			{
				for (k=0; k<3; k++)
				{
#ifdef _DOUBLE
					fscanf(fp, "%lf", &m_lineSources[i].points[j].pt[k]);
#else
					fscanf(fp, "%f", &m_lineSources[i].points[j].pt[k]);
#endif
				}
#ifdef _DOUBLE
				fscanf(fp, "%lf%lf%lf\n", &m_lineSources[i].points[j].rIntensity, &m_lineSources[i].points[j].rInnerRad,
					&m_lineSources[i].points[j].rOuterRad);
#else
				fscanf(fp, "%f%f%f\n", &m_lineSources[i].points[j].rIntensity, &m_lineSources[i].points[j].rInnerRad,
					&m_lineSources[i].points[j].rOuterRad);
#endif
				assert(m_lineSources[i].points[j].rOuterRad-m_lineSources[i].points[j].rInnerRad > EPS_ZERO_SQ);
			}
		}
		if (m_nTriSrcNum > 0)
			m_triSources = (TriangleSource *)malloc(sizeof(TriangleSource)*m_nTriSrcNum);
		fgets(cTok,100,fp);
		for (i=0; i<m_nTriSrcNum; i++)
		{
			fgets(cTok,100,fp);
			for (j=0; j<3; j++)
			{
				for (k=0; k<3; k++)
				{
#ifdef _DOUBLE
					fscanf(fp, "%lf", &m_triSources[i].points[j].pt[k]);
#else
					fscanf(fp, "%f", &m_triSources[i].points[j].pt[k]);
#endif
				}
#ifdef _DOUBLE
				fscanf(fp, "%lf%lf%lf\n", &m_triSources[i].points[j].rIntensity, &m_triSources[i].points[j].rInnerRad,
					&m_triSources[i].points[j].rOuterRad);
#else
				fscanf(fp, "%f%f%f\n", &m_triSources[i].points[j].rIntensity, &m_triSources[i].points[j].rInnerRad,
					&m_triSources[i].points[j].rOuterRad);
#endif
				assert(m_triSources[i].points[j].rOuterRad-m_triSources[i].points[j].rInnerRad > EPS_ZERO_SQ);
			}
		}
		fclose(fp);
		return true;
	}
	return false;
}

/*
 * calculate outer box
 */
bool DTIso2D::calcBox(POINT *minW, POINT* maxW)
{
	INTEGER i;
	int j;
	if (minW && maxW && m_nNodes > INIT_NOD_NUM)
	{
		for (j = 0; j < DIM; j++)
			(*minW)[j] = (*maxW)[j] = (m_pNodes[INIT_NOD_NUM].pt)[j];
	
		for (i = INIT_NOD_NUM; i < m_nNodes; i++)
		{
			for (j = 0; j < DIM; j++)
			{
				if ((m_pNodes[i].pt)[j] > (*maxW)[j])
					(*maxW)[j] = (m_pNodes[i].pt)[j];
				if ((m_pNodes[i].pt)[j] < (*minW)[j])
					(*minW)[j] = (m_pNodes[i].pt)[j];
			}
		}
		return true;
	}
	return false;
}
	

/*
 *	scale the geometry
 */
bool DTIso2D::scaGeom()
{
	int i, j;
	POINT wp;
	for (i = INIT_NOD_NUM; i < m_nNodes; i++)
	{
		for (j = 0; j < DIM; j++)
			wp[j] = (m_pNodes[i].pt)[j];
		wToN(wp, &(m_pNodes[i].pt));
		m_pNodes[i].density *= m_scale;
	}
	return true;
}

/*
 *	scale the geometry back
 */
bool DTIso2D::rescaGeom()
{
	int i, j;
	POINT wp;
	for (i = 0; i < m_nNodes; i++)
	{
		for (j = 0; j < DIM; j++)
			wp[j] = (m_pNodes[i].pt)[j];
		nToW(wp, &(m_pNodes[i].pt));
		m_pNodes[i].density /= m_scale;
	}
	for (i=0; i < m_nElems; i++)
	{
		for (j = 0; j < DIM; j++)
		{
			wp[j] = (m_pElems[i].cen)[j];
		}
		nToW(wp, &(m_pElems[i].cen));
		m_pElems[i].rad /= m_scale*m_scale;
	}
	
	return true;
}

/*
 *	scale the background mesh and source control
 */
bool DTIso2D::scaBkGrnd()
{
	int i, j, k;
#ifdef _3D_DECOM_
	POINT3D wp;
#else
	POINT wp;
#endif

#ifdef _3D_DECOM_
	for (i=0; i<m_nBkGrndNode; i++)
	{
		m_pBkGrndNode[i].density *= m_scale;
	}
#else
	for (i=0; i<m_nBkGrndNode; i++)
	{
		for (j=0; j<DIM; j++)
			wp[j] = m_pBkGrndNode[i].pt[j];
		wToN(wp, &(m_pBkGrndNode[i].pt));
		m_pBkGrndNode[i].density *= m_scale;
	}
#endif


	for (i=0; i<m_nPntSrcNum; i++)
	{
#ifdef _3D_DECOM_
		for (j=0; j < 3; j++)
		{
			wp[j] = m_pntSources[i].pt[j];
		}
		wToN3D(wp, &(m_pntSources[i].pt));
#else
		for (j=0; j < DIM; j++)
		{
			wp[j] = m_pntSources[i].pt[j];
		}
		wToN(wp, &(m_pntSources[i].pt));
#endif
		m_pntSources[i].rInnerRad *= m_scale;
		m_pntSources[i].rOuterRad *= m_scale;
		m_pntSources[i].rIntensity *= m_scale;
	}
	for (i=0; i<m_nLneSrcNum; i++)
	{
		for (k=0; k < 2; k++)
		{
#ifdef _3D_DECOM_
			for (j=0; j < 3; j++)
			{
				wp[j] = m_lineSources[i].points[k].pt[j];
			}
			wToN3D(wp, &(m_lineSources[i].points[k].pt));
#else
			for (j=0; j < DIM; j++)
			{
				wp[j] = m_lineSources[i].points[k].pt[j];
			}
			wToN(wp, &(m_lineSources[i].points[k].pt));
#endif
			m_lineSources[i].points[k].rInnerRad *= m_scale;
			m_lineSources[i].points[k].rOuterRad *= m_scale;
			m_lineSources[i].points[k].rIntensity *= m_scale;
		}
	}
#ifdef _3D_DECOM_
	for (i=0; i<m_nTriSrcNum; i++)
	{
		for (k=0; k < 3; k++)
		{
			for (j=0; j < 3; j++)
			{
				wp[j] = m_triSources[i].points[k].pt[j];
			}
			wToN3D(wp, &(m_triSources[i].points[k].pt));
			m_triSources[i].points[k].rInnerRad *= m_scale;
			m_triSources[i].points[k].rOuterRad *= m_scale;
			m_triSources[i].points[k].rIntensity *= m_scale;
		}
	}
#endif
	return true;
}

/*
 * calculate density values 
 */
bool DTIso2D::calcDens()
{
	int i;
	INTEGER iNod;
	REAL st = 1.0e12;
	INTEGER ctUnmatch = 0;

	for (i = INIT_NOD_NUM; i < m_nNodes; i++)
	{
		m_pNodes[i].density = 0.0;
	}
	for (i = 0; i < m_nBnds; i++)
	{
		Node* pNd1 = &(m_pNodes[m_pBnds[i].beg + INIT_NOD_NUM-1]); 
		Node* pNd2 = &(m_pNodes[m_pBnds[i].end + INIT_NOD_NUM-1]); 
		REAL dt = sqrt(distSqua(pNd1->pt, pNd2->pt));
		pNd1->density += dt;
		pNd2->density += dt;
	}
	for (iNod = INIT_NOD_NUM; iNod < m_nNodes; iNod++)
	{
		m_pNodes[iNod].density /= 2.0;

		st = spacFrmSrc(m_pNodes[iNod].pt);
		/* ��ȡȫ���ܶ� */
		if (m_nBkGrndNode > 0 && m_pBkGrndNode)
			st = Min(st, m_pBkGrndNode[0].density);
	                             
#ifdef _DEBUG
		if (st/m_pNodes[iNod].density >= 2.0 || st/m_pNodes[iNod].density <= 0.5)
		{/* unmatch in the mesh density */
			ctUnmatch++;
		}
#endif /* _DEBUG */

		m_pNodes[iNod].density = Min(m_pNodes[iNod].density, st);
	}

#ifdef _DEBUG
	printf("%d of totally %d nodes have unmatched spacing value.\n", ctUnmatch, m_nNodes);
#endif
	
	return true;
}

 /* ---------------------------------
	* for normalizing the coordinates |
  * ---------------------------------*/
	
/*
 * scale the coordinates range [minW, maxW] to a range of [minN, maxN],
 * and return according center & scale value
 */
bool DTIso2D::scaleFactor(POINT minW, POINT maxW, POINT minN, POINT maxN)
{
	int i, j;
	REAL ss[DIM], ssmin;
	ss[0] = (maxN[0] - minN[0]) / (maxW[0] - minW[0]);
	ssmin = ss[0];
	j = 0;
	for (i = 1; i < DIM; i++)
	{
		ss[i] = (maxN[i] - minN[i]) / (maxW[i] - minW[i]);
		if (ss[i] < ssmin)
		{
			ssmin = ss[i];
			j = i;
		}	
	}

	m_scale = ssmin;

	for (i = 0; i < DIM; i++)
	{
		m_cenN[i] = 0.5 * (minN[i] + maxN[i]);
		m_cenW[i] = 0.5 * (minW[i] + maxW[i]);
		m_minW[i] = minW[i];
		m_maxW[i] = maxW[i];
		m_minN[i] = minN[i];
		m_maxN[i] = maxN[i];
	}
	
	return true;
}
						
/* 
 * world cordinates to normalization coordinates 
 */
bool DTIso2D::wToN(POINT wp, POINT *np)
{
	int i;
	if (np)
	{
		for (i = 0; i < DIM; i++)
			(*np)[i] = m_cenN[i] + m_scale * (wp[i] - m_cenW[i]);
		return true;
	}
	return false;
}
#ifdef _3D_DECOM_
bool DTIso2D::wToN3D(POINT3D wp, POINT3D *np)
{
	int i;
	if (np)
	{
		for (i = 0; i < 2; i++) 
			(*np)[i] = m_cenN[i] + m_scale * (wp[i] - m_cenW[i]);
		(*np)[2] = m_scale*wp[2]; // add this [11/26/2005]
		return true;
	}
	return false;
}
#endif		
/* 
 * normalization cordinates to world coordinates 
 */
bool DTIso2D::nToW(POINT np, POINT *wp)
{
	int i;
	if (wp && m_scale != 0.0)
	{
		for (i = 0; i < DIM; i++)
			(*wp)[i] = m_cenW[i] + (np[i] - m_cenN[i]) / m_scale;
		return true;
	}
	return false;
}

int DTIso2D::scaleBoxPnt()
{
	INTEGER i;
	POINT np, wp;
	for (i = 0; i < INIT_NOD_NUM; i++)
	{
		np[0] = g_cors[i][0];
		np[1] = g_cors[i][1];
		nToW(np, &wp);
		g_cors[i][0] = wp[0];
		g_cors[i][1] = wp[1];
	}
	
	
	return 1;
}

int DTIso2D::setupInitTri()
{
	int i, j, k;
	FORM_PNTS pnts;
	REAL d_a;
	if (m_pElems && m_pNodes)
	{
		for (i = 0; i < INIT_NOD_NUM; i++)
		{
			for (j = 0; j < DIM; j++)
				((m_pNodes)[i].pt)[j] = g_cors[i][j];
		}
		for (i = 0; i < INIT_TRI_NUM; i++)
		{
			for (j = 0; j <= DIM; j++)
			{
				(m_pElems[i].form)[j] = 
					g_form[i][j];
				(m_pElems[i].neig)[j] = g_neig[i][j];
			}
			for (j = 0; j <= DIM; j++)
				for (k = 0; k < DIM; k++)
					pnts[j][k] = g_cors[g_form[i][j]][k];
			calcElePar(pnts, &d_a, &(m_pElems[i].cen), &(m_pElems[i].rad));
			m_pElems[i].iReserved = 0;
		}
	}
	m_nElems = INIT_TRI_NUM;
//	m_nNodes = INIT_NOD_NUM;
	return 1;
}

int DTIso2D::findFirstEle(POINT pnt, INTEGER *ele)
{
	//loop from the last elements
	INTEGER iElem = m_nElems - 1; /*current element where tree search is performed*/
	INTEGER iSrch = iElem; /*for tree search*/
	Elem* pElem = NULL, *pSrch = NULL;
	INTEGER iNxt = 0; /*index referring to the current index in the neighboring array*/
	INTEGER iSel = NULL_ELEM; /*index referring to the selected element for the next searching*/
	const REAL MAX_DT = std::numeric_limits<double>::max();
	REAL min_dt = MAX_DT;
	REAL dt, cri;
	int nCase = -2; /*-2: unknown; -1: error; 0: degeneracy; 1: successful*/
	INTEGER i;
	for (i = m_nElems-1; i >= 0; i--)
	{
		iElem = i;
		iSrch = i;
		if (isDelEle(iSrch))
		{
			continue;
		}
		iNxt = 0;
		min_dt = MAX_DT;
		pElem = &(m_pElems)[iSrch];
		do
		{
		 	pSrch = &(m_pElems)[iSrch];
			if (!isTstEle(iSrch) && !isDelEle(iSrch))
			{
				addTstEle(iSrch); //?

				//use predict.cxx

				auto p1=  m_pNodes[pSrch->form[0]].pt;
				auto p2 = m_pNodes[pSrch->form[1]].pt;
				auto p3 = m_pNodes[pSrch->form[2]].pt;
				//old form may cause Floating point accuracy problem
				dt = distSqua(pSrch->cen, pnt);
				//cri = dt - pSrch->rad;

				// new form
				cri =-incircle(p1, p2, p3, pnt);

				if (cri < 0)
				{/*incircle cirteria is broken*/
					*ele = iSrch;
					nCase = 1;
				}		
				else if (cri > 0)/*incircle criteria is kept, perform tree search*/
					nCase = -3; //flag 
				else /* four points in the same circle */
					nCase = 0;	
			}

			if (nCase != 0 && nCase != 1)
			{
				if (nCase == -3)
				{
	//				addTstEle(iSrch);//?
					if (dt < min_dt)
					{
						iSel = iSrch;
						min_dt = dt;
					}
					nCase = -2; //restore the flag
				}
				if (iNxt > DIM)
				{
					iNxt = 0;
					if (iSel == NULL_ELEM || iElem == iSel) //Still the same, break to for loop
					{
						nCase = -2;
						break;
					}
					iElem = iSel;
					pElem = &(m_pElems)[iElem];
					iSrch = iSel;
					min_dt = MAX_DT;
	#ifdef _ERROR_CHK
				/*	for (i = 0; i <DIM; i++)
						if (!isTstEle((pElem->neig)[i]]) && !isDelEle(pElem->neig)[i]))
							break;
					if (i >= DIM)
					{
						printf("cannot find first triangle!\n");
						//go to error handle
					}*/
	#endif //_ERROR_CHK
					continue;
				}
				while (iNxt <= DIM && (iSrch = (pElem->neig)[iNxt++]) == NULL_NEIG);
			}
		}
		while (nCase == -2 && iSrch != -1);
		if (nCase == 1)
			break;
	}

	//clear test flags
	clrTstEles();

	return nCase;
}

/*
 * boundary point insertion
 */
/*
 * boundary point insertion
 */
int DTIso2D::bndPntInst()
{
	int iSucc = 0, i, iFail = 0;
	INTEGER iNod;
	POINT pnt;
	REAL disturbDist = 0.0;

	for (i = 0; i < m_nBnds; i++)
		m_lstInstBndNods.push_back(i + INIT_NOD_NUM);
	
 	while (iSucc != m_nBnds)
	{
		assert(!m_lstInstBndNods.empty());
		std::list<INTEGER>::iterator it_fir = 
			m_lstInstBndNods.begin();
		iNod = *it_fir;
		for (i = 0; i < DIM; i++)
			pnt[i] = (m_pNodes[iNod].pt)[i];

		POINT deb;
		nToW(pnt,&deb);
		if (addBndPnt(pnt, iNod) == 1)
		{
			m_lstInstBndNods.erase(it_fir);	
			iSucc++;
			iFail = 0;
		}
		else
		{
#ifdef _VERBOS
			cout<<"Node "<<iNod<<": Failed"<<endl;
#endif
			iFail++;
			if (iFail >= m_lstInstBndNods.size())
			{
#ifdef _VERBOS
				cout<<"Node "<<iNod<<": Disturbed"<<endl;
#endif
				if (!isDisturbed(iNod))
				{
					for (i=0; i<DIM; i++)
					{
						pnt[i] = m_pNodes[iNod].pt[i];
					}
					addDistInfo(iNod, pnt);
				}
				//disturb the point
				disturbDist = EPS_DISTURB * m_pNodes[iNod].density;
				for (i = 0; i < DIM; i++)
				{
					(m_pNodes[iNod].pt)[i] = (m_pNodes[iNod].pt)[i] + 
						disturbDist*(g_cors[0][i] - g_cors[INIT_NOD_NUM - 1][i]);
				}
			}
			else //delay the insertion of iNod
			{
				m_lstInstBndNods.erase(it_fir);
				m_lstInstBndNods.push_back(iNod);
			}
		}
#ifdef _ERROR_CHK
//		char middle[20];
//		sprintf(middle,"Mid%d.dt2",iNod);
		//		writeDt2(middle);
#endif
	}
	recDistPnts();
	return 1;
}

/*
 * add boundary points
 */
#pragma optimize("",off)
int DTIso2D::addBndPnt(POINT pnt, INTEGER iNod)
{
	int i, k, m;
	INTEGER ele, iElem, iSrch;
	Elem *pElem = NULL, *pSrch = NULL;
	INTEGER nnew[DIM+1];
	POINT pnew[DIM+1];
	POINT cen;
	REAL d_a, rad;
	INTEGER iEmp; /*empty position for new created element*/
	INTEGER iTakEle;
	int nfc = 0, npc = 0;
	int nBegin = 0;
#ifdef _VERBOS
	cout<<"Inserting boundary node "<<iNod<<". Current elements num: "<<m_nElems<<endl;
#endif
//	char s[50];
//	gets(s);
	int nCase = findFirstEle(pnt, &ele);
	if (nCase == 1)//find it
	{
		//perform tree search
//		addTstEle(ele);
		addDeleted(ele);
		addTreeSearch(ele);
		while (!isTreeSearchEmpty())
		{
			iElem = pickTreeSearch();
			pElem = &(m_pElems)[iElem];
			for (i = 0; i <= DIM; i++)
			{
				iSrch = (pElem->neig)[i];
				pSrch = &(m_pElems)[iSrch];
				if (iSrch == NULL_NEIG)
				{
					nBegin = 1;
				}
				else if (/*!isTstEle(iSrch) &&*/ !isDelEle(iSrch))
				{
//					addTstEle(iSrch);
					
					REAL dt = distSqua(pSrch->cen, pnt); /*distance square between two points*/
					REAL cri = dt - pSrch->rad;
					if (cri < 0)
					{/*incircle cirteria is broken*/
					//	printf("volume is near zero!!\n");
						addDeleted(iSrch);
						addTreeSearch(iSrch);
					}		
					else if (cri > EPS_ZERO_SQ)/*incircle criteria is kept, perform tree search*/
					{
						nBegin = 1;
					}//else if (cri > EPS_ZERO_SQ)
					else
					{/* four points in the same circle */
						nCase = 0;
						goto RECOVERED;
					}
				}/*if (!isTested(pSrch) && isValid(pSrch))*/
				if (nBegin)
				{
					/*find the nonshared node index*/
//					for (j = 0; j <= DIM; j++)
//					{
//						no = (pElem->form)[j];
//						if (no != (pSrch->form)[0] && no != (pSrch->form)[1] && no != (pSrch->form)[2])
//							break;
//					}
//#ifdef _ERROR_CHK
//					if (j > DIM)
//					{
//						printf("can not find uncommon point!!\n");
//						//go to error handle
//					}
//					
//#endif //_ERROR_CHK				
					for (k = 0; k < DIM; k++)
					{
						nnew[k] = (pElem->form)[(i+k+1)%(DIM+1)];//[(j+k+1)%(DIM+1)];
						for (m = 0; m < DIM; m++)
							pnew[k][m] = (m_pNodes[nnew[k]].pt)[m];
					}
					nnew[k] = iNod; /*add a node*/
					for (m = 0; m < DIM; m++)
						pnew[DIM][m] = pnt[m];
					calcElePar(pnew, &d_a, &cen, &rad);
//#ifdef _ERROR_CHK
					if (d_a < 0.0)
					{
					//	this->writevtk("123.vtk");
						printf("negative volume\n!!");
						nCase = 3;
						goto RECOVERED;
						//go to error handle
					}
					if (d_a < EPS_ZERO_SQ)
					{
						printf("volume is near zero!!\n");
						nCase = 3;
						goto RECOVERED;
						//go to error handle
					}
//#endif //_ERROR_CHK
					nfc = nfc + 1;
					
					iEmp = crtEle(nnew, cen, rad); /*create a new element & return its location*/
					((m_pElems)[iEmp].neig)[DIM]/*[0]*/ = iSrch;
//					writeDt2("Temp.dt2");
        			
   					/*
   					 * check two nodes of the created element. If they are taken by a new created element,
   					 * then update the neighboring information correspondingly
   					 */
   					for (k = 0; k < DIM; k++)
   					{
   						iTakEle = getNodTakEle(nnew[k]);
   						if (iTakEle >= 0)
   						{//taken
   							((m_pElems)[iEmp].neig)[(k+1)%DIM] = iTakEle;
   							if (((m_pElems)[iTakEle].form)[0] == nnew[k])
   			 					((m_pElems)[iTakEle].neig)[1] = iEmp;
   							else
   			 					((m_pElems)[iTakEle].neig)[0]/*[2]*/ = iEmp;
   						}
   						else
   						{
   							addNodTakEle(nnew[k], iEmp);
							npc = npc + 1;
   						}
					}
					nBegin = 0;
				}
			}/*for (i = 0; i < DIM; i++)*/
		}/*while (!isTreeSearchEmpty())*/
	}/*if (nCase == 1)*/
	else
	{

	this->writevtk("current.vtk");
	int nCase = findFirstEle(pnt, &ele);
		cout<<"***Error: Cannot find the first element for Node "<<iNod<<endl;
		nCase = 3;
		goto RECOVERED;
	}
//#ifdef _ERROR_CHK
	if (nfc != (DIM - 1) * npc - 4 * (DIM - 2))
	{
		printf("the number of faces & the number of nodes is imbalant!!");
		nCase = 3;
		goto RECOVERED;
		//go to error handle
	}
//#endif //_ERROR_CHK

	/*update neighboring info. of elements near the cavity*/
	updateNeigInfo();
	updateNewEleLoc();
	goto FINSHED;
RECOVERED:
	/*error happens, try to recover*/
	undoInsertion();
	clrTreeSearch();
FINSHED:
	clrAddEles();
	clrDelEles();
	//clear test flags
//	clrTstEles();
	//clear taken flags
	clrTakNods();
	if (nCase == 1) //successful
	{
	}
	return nCase;
}
#pragma optimize("",on)
/* --------------------------------------
 * function for tree searching          |
 * -------------------------------------*/
INTEGER DTIso2D::pickTreeSearch()
{
	INTEGER iEle = NULL_ELEM;
	if (!m_stkTreeSear.empty())
	{
		iEle = m_stkTreeSear.top();
		m_stkTreeSear.pop();
	}
	return iEle;
}

int DTIso2D::addTreeSearch(INTEGER ele)
{
	m_stkTreeSear.push(ele);
	return m_stkTreeSear.size();
}

bool DTIso2D::isTreeSearchEmpty()
{
	return m_stkTreeSear.empty();
}

/* -----------------------------------------
 * function for deleted elements          |
 * ---------------------------------------*/
int DTIso2D::addDeleted(INTEGER iEle)
{
	m_pElems[iEle].rad = -m_pElems[iEle].rad;
	m_vecDelEles.push_back(iEle);
	return m_vecDelEles.size();
}

/* ------------------------------------------------------
 * function for location management of new elements     |
 * -----------------------------------------------------*/
INTEGER DTIso2D::getNewEleLoc()
{
	INTEGER i;
	INTEGER iLoc;
	if (m_nLocInd <= 0)
	{
		return m_nElems;// + 1;
	}
	else
	{
		i = m_nLocInd-1;
		do {
			int j,size;
			bool bFind = false;
			size = m_vecDelEles.size();
			for (j=0; j<size; j++)
			{
				if (m_vecDelEles[j] == m_vecEmpLocs[i])
				{
					bFind = true;
					break;
				}
			}
			if (!bFind)  //Ok, deleted long ago, so reuse it
			{
				iLoc = m_vecEmpLocs[i];
				m_vecEmpLocs.erase(m_vecEmpLocs.begin()+i);
				m_nLocInd--;
				return iLoc;
			}
			i--;
		} while(i>=0);
		return m_nElems;
	}
}

/* ------------------------------------------------------
 * function for management of node-taken elements       |
 * -----------------------------------------------------*/
INTEGER DTIso2D::getNodTakEle(INTEGER iNod)
{
	int size = m_vecTakNods.size();
	for (int i=0; i<size; i++)
	{
		if (m_vecTakNods[i] == iNod)
		{
			m_vecTakNods.erase(i+m_vecTakNods.begin());
			return (m_pNodes)[iNod].iReserved;
		}
	}
	return -1;
}

bool DTIso2D::addNodTakEle(INTEGER iNod, INTEGER iEle)
{
	(m_pNodes)[iNod].iReserved = iEle;
	m_vecTakNods.push_back(iNod);
	return true;
}

/*
 *update neighboring info. of elements near the cavity
 */
bool DTIso2D::updateNeigInfo()
{
	bool bChk = true;
	int i, j;
	INTEGER iNeig, iNgNg;
	for (i = 0; i < m_vecAddEles.size(); i++)
	{
		iNeig = ((m_pElems)[m_vecAddEles[i]].neig)[DIM];
		if (iNeig >= 0) //a valid neighbour
		{
			for (j = 0; j <= DIM; j++)
			{
				if (m_pElems[iNeig].form[j] != m_pElems[m_vecAddEles[i]].form[0] &&
					m_pElems[iNeig].form[j] != m_pElems[m_vecAddEles[i]].form[1]) 					
				{ // Now, this new element should be the jth neighbour of Element iNeig
					iNgNg = ((m_pElems)[iNeig].neig)[j];
					if (iNgNg >= 0 && isDelEle(iNgNg)) // This should be always true if iNgNg>=0
					{
						((m_pElems)[iNeig].neig)[j] = m_vecAddEles[i];
						break;
					}
				}
			}
#ifdef _ERROR_CHK
			if (j > DIM)
			{
				bChk = false;
				printf("cannot set neighboring info. of elements near the cavity correctly!\n");
				break;
			}
#endif //_ERROR_CHK
		}
	}
	return bChk;
}
	
/*
 * calcuate element parameters
 * d_a: double area
 * cen: center of circumcenter
 * rad: radius of circumcenter
 */
bool DTIso2D::calcElePar(POINT pnt[], REAL *d_a, POINT *cen, REAL *rad)
{
	REAL a, b, d, e, d1, d2, d3, c, f, te, t1, t2;
	a = 2.0 * (pnt[2][0] - pnt[1][0]);
	b = 2.0 * (pnt[2][1] - pnt[1][1]);
	d = 2.0 * (pnt[0][0] - pnt[1][0]);
	e = 2.0 * (pnt[0][1] - pnt[1][1]);
	d1 = pnt[0][0] * pnt[0][0] + pnt[0][1] * pnt[0][1];
	d2 = pnt[1][0] * pnt[1][0] + pnt[1][1] * pnt[1][1];
	d3 = pnt[2][0] * pnt[2][0] + pnt[2][1] * pnt[2][1];
	c = d3 - d2;
	f = d1 - d2;
	te = a * e - b * d;
	t1 = c * e - f * b;
	t2 = a * f - c * d;
	*d_a = GEOM_FUNC::orient2d(&(pnt[0][0]),&(pnt[1][0]),&(pnt[2][0]));
	(*cen)[0] = t1 / te;
	(*cen)[1] = t2 / te;
	*rad = distSqua(*cen, pnt[0]);
	return true;	
}

/*
 * create a new element
 */
INTEGER DTIso2D::crtEle(INTEGER form[], POINT cen, REAL rad)
{
	INTEGER iLoc = getNewEleLoc();
	int i;
	for (i = 0; i <= DIM; i++)
	{
		((m_pElems)[iLoc].form)[i] = form[i];
		((m_pElems)[iLoc].neig)[i] = NULL_NEIG;
	}
	for (i = 0; i < DIM; i++)
	{
		(m_pElems[iLoc].cen)[i] = cen[i];
	}
	m_pElems[iLoc].rad = rad;
	m_pElems[iLoc].iReserved = 0;
	if (iLoc == m_nElems)
		m_nElems++;
	m_vecAddEles.push_back(iLoc);
	return iLoc;
} 

/*
 * function for elements property
 */
bool DTIso2D::isDelEle(INTEGER iEle)
{
	assert(iEle >= 0 && iEle < m_nElems);
	return (m_pElems)[iEle].rad <= 0.0;
}

/*
 * check if an elememt is deleted
 */
bool DTIso2D::isTstEle(INTEGER iEle)
{
	assert(iEle >= 0 && iEle < m_nElems);
	return (m_pElems)[iEle].iReserved == 1;
}

/*
 * add an element to tested vector
 */
bool DTIso2D::addTstEle(INTEGER iEle)
{
	assert(iEle >= 0 && iEle < m_nElems);
	m_vecTstEles.push_back(iEle);
	enabEleTst(iEle, true);
	return true;
}

/*
 * clear flags of all tested elements
 */
bool DTIso2D::clrTstEles()
{
	int i;
	for (i = 0; i < m_nElems; i++)
		enabEleTst(i, false);
	m_vecTstEles.clear();
	return true;
}

/*
 * enable test flag of an element
 */
bool DTIso2D::enabEleTst(INTEGER iEle, bool flag)
{
	assert(iEle >= 0 && iEle < m_nElems);
	(m_pElems)[iEle].iReserved = flag ? 1 : 0;
	return flag;
}

/* ---------------------------------
 * for disturbance point info.     |
 * ---------------------------------*/
bool DTIso2D::isDisturbed(INTEGER iNod)
{
	std::list<DistInfo*>::iterator it = m_lstDistInfo.begin();
	while (it != m_lstDistInfo.end())
	{
		DistInfo* pDI = *it;
		if (pDI->iNod == iNod)
			return true;

		++it;
	}
	return false;
}

bool DTIso2D::addDistInfo(INTEGER iNod, POINT pnt)
{
	int i;
	DistInfo* pDI = new DistInfo;
	if (pDI)
	{
		pDI->iNod = iNod;
		for (i = 0; i < DIM; i++)
			(pDI->old_pt)[i] = pnt[i];
		m_lstDistInfo.push_back(pDI);
		return true;
	}
	return false;
}

/*
 * calculate square of distance between two points
 */
REAL distSqua(POINT p1, POINT p2)
{
	REAL s = 0.0;
	int i;
	for (i = 0; i < DIM; i++)
		s += (p2[i] - p1[i]) * (p2[i] - p1[i]);
	return s;
}
#ifdef _3D_DECOM_
REAL distSqua_3D(POINT3D p1, POINT3D p2)
{
	REAL s = 0.0;
	int i;
	for (i = 0; i < 3; i++)
		s += (p2[i] - p1[i]) * (p2[i] - p1[i]);
	return s;
}
#endif

int DTIso2D::writeDt2(char *fname)
{
	INTEGER i;
	ofstream out;
	out.open(fname, ios::out);

	out<<"Node Count: "<<m_nNodes<<endl;
	out<<"Boundary Count: "<<m_nBnds<<endl;
	out<<"Element Count: "<<m_nElems<<endl;

	out<<"Deleted Element Count: "<<m_vecDelEles.size()<<endl;
	out<<"Added Element Count: "<<m_vecAddEles.size()<<endl;
	out<<"Tested Element Count: "<<m_vecTstEles.size()<<endl;

	out<<"Taken Nodes Count: "<<m_vecTakNods.size()<<endl;

	out<<"Node Information:"<<endl;
	for (i=0; i<m_nNodes; i++)
	{
		out<<i<<"\t"<<m_pNodes[i].pt[0]<<"\t"<<m_pNodes[i].pt[1]<<"\t"<<m_pNodes[i].density<<"\t"<<m_pNodes[i].iReserved<<endl;
	}

	out<<"Boundary Information:"<<endl;
	for (i=0; i<m_nBnds; i++)
	{
		out<<i<<"\t"<<m_pBnds[i].beg<<"\t"<<m_pBnds[i].end<<"\t"<<m_pBnds[i].ele<<"\t"<<m_pBnds[i].curve<<"\t"<<m_pBnds[i].loop<<endl;
	}

	out<<"Element Information:"<<endl;
	for (i=0; i<m_nElems; i++)
	{
		out<<i<<"\t"<<m_pElems[i].form[0]<<"\t"<<m_pElems[i].form[1]<<"\t"<<m_pElems[i].form[2]
			<<"\t"<<m_pElems[i].cen[0]<<"\t"<<m_pElems[i].cen[1]<<"\t"<<m_pElems[i].rad
			<<"\t"<<m_pElems[i].neig[0]<<"\t"<<m_pElems[i].neig[1]<<"\t"<<m_pElems[i].neig[2]
			<<"\t"<<m_pElems[i].iReserved<<endl;
	}
	
	out.close();
	return 0;
}

#if 0
int DTIso2D::writePl2(char *fname)
{
	INTEGER i;
	ofstream out;
	out.open(fname, ios::out);	

	out<<m_nElems<<"\t"<<m_nNodes<<"\t"<<m_nBnds<<endl;
	for (i=0; i<m_nNodes; i++)
	{
		out<<i+1<<"\t"<<m_pNodes[i].pt[0]<<"\t"<<m_pNodes[i].pt[1]<<endl;
	}

	for (i=0; i<m_nElems; i++)
	{
		out<<i+1<<"\t"<<m_pElems[i].form[0]+1<<"\t"<<m_pElems[i].form[1]+1<<
			"\t"<<m_pElems[i].form[2]+1<<"\t"<<0<<"\t"<<0<<endl;
	}

	for (i=0; i<m_nBnds; i++)
	{
		out<<i+1<<"\t"<<m_pBnds[i].beg+1<<"\t"<<m_pBnds[i].end+1<<"\t"<<m_pBnds[i].ele+1<<"\t"<<0<<"\t"<<0<<endl;
	}
	out.close();
	return 0;
}
#else
int DTIso2D::writePl2(char *fname)
{
	INTEGER i;
	FILE *fp = fopen(fname, "w");
	if (!fp)
	{
		printf("Cannot open file %s.", fname);
		return 2;
	}
	
	fprintf(fp, "%d %d %d\n", m_nElems, m_nNodes, m_nBnds);
	for (i=0; i<m_nNodes; i++)
	{
		fprintf(fp, "%d %f %f\n", i+1, m_pNodes[i].pt[0], m_pNodes[i].pt[1]);
	}
	
	for (i=0; i<m_nElems; i++)
	{
		fprintf(fp, "%d %d %d %d 0 0\n", i+1, m_pElems[i].form[0]+1,
			m_pElems[i].form[1]+1, m_pElems[i].form[2]+1);
	}
	
	for (i=0; i<m_nBnds; i++)
	{
		fprintf(fp, "%d %d %d %d 0 0\n", i+1, m_pBnds[i].beg+1, 
			m_pBnds[i].end+1, m_pBnds[i].ele+1);
	}
	fclose(fp);

	return 0;
}
#endif

int DTIso2D::writePl2(int *npt, double **gpt, int *nelm, int** elm)
{
	int i;

	*npt = m_nNodes;
	*nelm = m_nElems;

	*gpt = new double [2*m_nNodes];
	*elm = new int[3*m_nElems];

	for (i=0; i<m_nNodes; i++)
	{
		(*gpt)[i*2+0] = m_pNodes[i].pt[0];
		(*gpt)[i*2+1] = m_pNodes[i].pt[1];
	}

	for (i=0; i<m_nElems; i++)
	{
		(*elm)[i*3+0] = m_pElems[i].form[0];
		(*elm)[i*3+1] = m_pElems[i].form[1];
		(*elm)[i*3+2] = m_pElems[i].form[2];
	}

	return 0;
}

int DTIso2D::writevtk(char *filename)
{
	int i, j, npt, nelm, sidx=0;
	FILE *fout = NULL;
	
	npt = m_nNodes;
	nelm = m_nElems;
	
	fout = fopen(filename, "w");
	if (fout == NULL)
	{
		printf("Can not open file %s\n", filename);
		return 0;
	}
	
	//vtk
	fprintf(fout, "# vtk DataFile Version 2.0\n");
	fprintf(fout, "boundary layer mesh\n");
	fprintf(fout, "ASCII\n");
	fprintf(fout, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fout, "POINTS %d double\n", m_nNodes);
	
	for (i=0; i<npt; i++)
	{
		fprintf(fout, "%f %f 0.0\n", m_pNodes[i].pt[0], m_pNodes[i].pt[1]);
	}
	
	fprintf(fout, "CELLS %d %d\n", m_nElems, m_nElems*4);
		
	for (i=0; i<nelm; i++)
	{
		fprintf(fout, "3 %d %d %d\n", m_pElems[i].form[0], 
			m_pElems[i].form[1], m_pElems[i].form[2]);
	}
	
	fprintf(fout, "CELL_TYPES %d\n", m_nElems);
		
	for(i=0; i< m_nElems; i++)
	{
		fprintf(fout, "%d\n", 5);
	}

	fclose(fout);
	return 0;
}

int DTIso2D::writeSpa(char *fname)
{
	INTEGER i;
	ofstream out;
	out.open(fname, ios::out);	

	out<<m_nNodes<<endl;
	for (i=0; i<m_nNodes; i++)
	{
		out<<i+1<<"\t"<<m_pNodes[i].density<<endl;
	}
	out.close();
	return 0;
}

int DTIso2D::writeNgb(char *fname)
{
	INTEGER i;
	int j;
	std::ofstream out;
	out.open(fname, ios::out);

	out<<m_nElems<<endl;
	for (i=0; i<m_nElems; i++)
	{
		out<<i+1;
		for (j=0; j<=DIM; j++)
		{
			if (m_pElems[i].neig[j] != NULL_NEIG)
			{
				out<<"\t"<<m_pElems[i].neig[j]+1;
			}
			else
			{
				out<<"\t"<<NULL_NEIG;
			}
		}
		out<<endl;
	}
	return 0;
}
/*
 * clear the vector after all a point is added
 */
bool DTIso2D::clrAddEles()
{
	m_vecAddEles.clear();
	return true;
}

/*
 * update the empty location vector after a point is added
 */
bool DTIso2D::updateNewEleLoc()
{
	int i;
	int size;
	size = m_vecDelEles.size();
	for (i=0; i<size; i++)
	{
		m_vecEmpLocs.push_back(m_vecDelEles[i]);
	}
	m_nLocInd = m_vecEmpLocs.size();
	return true;
}

/*
 * clear the vector after a point is added
 */
bool DTIso2D::clrTakNods()
{
	m_vecTakNods.clear();
	return true;
}

/*
 * undo current insertion operation if some errors happen
 */
bool DTIso2D::undoInsertion()
{
	int i,size;
	INTEGER iNod;
	size = m_vecDelEles.size();
	for (i=0; i<size; i++)
	{
		iNod = m_vecDelEles[i];
		m_pElems[iNod].rad = -m_pElems[iNod].rad;
	}
	size = m_vecAddEles.size();
	for (i=0; i<size; i++)
	{
		iNod = m_vecAddEles[i];
		m_pElems[iNod].rad = -m_pElems[iNod].rad;
		m_vecEmpLocs.push_back(iNod);
	}
	m_nLocInd = m_vecEmpLocs.size();
	return true;
}

/*
 * clear the vector after a point is added
 */
bool DTIso2D::clrDelEles()
{
	m_vecDelEles.clear();
	return true;
}

/*
 * clear the tree search stack
 */ 
bool DTIso2D::clrTreeSearch()
{
	while (!m_stkTreeSear.empty())
	{
		m_stkTreeSear.pop();
	}
	return true;
}

/*
 * recover all boundarys
 */
bool DTIso2D::recoverBnds()
{
	INTEGER iElem, iSrch;
	INTEGER i,j;
	int k;
	for (i=0; i<m_nBnds; i++)
	{
	}
	//determine which element every bound belongs to 
	for (i=0; i<m_nBnds; i++)
	{
		//first, determine which element this bound node belongs to
		bool bFind = false;
		for (j=0; j<m_nElems; j++)
		{
			if (isDelEle(j))
			{
				continue;
			}
			for (k=0; k<=DIM; k++)
			{
				if (m_pElems[j].form[k] == m_pBnds[i].beg+INIT_NOD_NUM-1)
				{
					iElem = j;
					bFind = true;
					break;
				}
			}
			if (bFind)
			{
				break;
			}
		}
		if (!bFind)
		{
			cout<<"***Error: Cannot find element for boundary "<<i<<endl;
			continue;
		}
		// then, check whether this boundary needs swaping
		bFind = false;
		iSrch = iElem;
		do {
			int index;
			for (j=0; j<=DIM; j++)
			{
				if (m_pElems[iSrch].form[j] == m_pBnds[i].beg+INIT_NOD_NUM-1)
				{
					index = j;
					break;
				}
			}
			if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[i].end+INIT_NOD_NUM-1)
			{
				//Find this bound edge
				bFind = true;
				break;
			}
			iSrch = m_pElems[iSrch].neig[(index+1)%(DIM+1)];
			if (iSrch == NULL_NEIG)
			{
				break;
			}
		} while(iSrch != iElem);
		if (!bFind)
		{
			//Try the other direction
			iSrch = iElem;
			do {
				int index;
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iSrch].form[j] == m_pBnds[i].beg+INIT_NOD_NUM-1)
					{
						index = j;
						break;
					}
				}
				if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[i].end+INIT_NOD_NUM-1)
				{
					//Find this bound edge
					bFind = true;
					break;
				}
				iSrch = m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
				if (iSrch == NULL_NEIG)
				{
					break;
				}
			} while(iSrch != iElem);
		}
		if (bFind == false)
		{
			//Search the first edge that needs swapping
			INTEGER iCur, iPrev;
			int index = -1;
			bool bSwapFind = false;
			iSrch = iElem;
			do {
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iSrch].form[j] == m_pBnds[i].beg+INIT_NOD_NUM-1)
					{
						index = j;
						break;
					}
				}
				if (index != -1) 
				{
// 					if (Collinear(m_pNodes[m_pBnds[i].beg+INIT_NOD_NUM-1].pt,
// 						m_pNodes[m_pBnds[i].end+INIT_NOD_NUM-1].pt,
// 						m_pNodes[m_pElems[iSrch].form[(index+1)%(DIM+1)]].pt) || 
// 						Collinear(m_pNodes[m_pBnds[i].beg+INIT_NOD_NUM-1].pt,
// 						m_pNodes[m_pBnds[i].end+INIT_NOD_NUM-1].pt,
// 						m_pNodes[m_pElems[iSrch].form[(index+2)%(DIM+1)]].pt))
// 					{
// 							cout<<"***Error: Collinear case found during recovering edge "<<i<<endl;
// 					}
					if (IntersectProp(m_pNodes[m_pBnds[i].beg+INIT_NOD_NUM-1].pt,
						m_pNodes[m_pBnds[i].end+INIT_NOD_NUM-1].pt,
						m_pNodes[m_pElems[iSrch].form[(index+1)%(DIM+1)]].pt,
						m_pNodes[m_pElems[iSrch].form[(index+2)%(DIM+1)]].pt))
					{
// 						SwapEdge edge;
// 						edge.iElemLeft = iSrch;
// 						edge.iElemRight = m_pElems[iSrch].neig[index];
// 						edge.iDiag[0] = m_pElems[iSrch].form[(index+1)%(DIM+1)];
// 						edge.iDiag[DIM-1] = m_pElems[iSrch].form[(index+2)%(DIM+1)];
// 						addSwapEdge(&edge);
// 						iCur = m_pElems[iSrch].neig[index];
// 						iPrev = iSrch;
// 						bSwapFind = true;
// 						break;
						addSwapEdge(iSrch, m_pElems[iSrch].neig[index], 
							m_pElems[iSrch].form[(index+1)%(DIM+1)], m_pElems[iSrch].form[(index+2)%(DIM+1)]);
						iCur = m_pElems[iSrch].neig[index];
						iPrev = iSrch;
						bSwapFind = true;
						break;
					}	
				}
				iSrch = m_pElems[iSrch].neig[(index+1)%(DIM+1)];
				if (iSrch == NULL_NEIG)
				{
					break;
				}
			} while(iSrch != iElem);
			if (!bSwapFind)
			{
				//Try the other direction
				iSrch = iElem;
				do {
					for (j=0; j<=DIM; j++)
					{
						if (m_pElems[iSrch].form[j] == m_pBnds[i].beg+INIT_NOD_NUM-1)
						{
							index = j;
							break;
						}
					}
					if (index != -1) 
					{
// 						if (Collinear(m_pNodes[m_pBnds[i].beg+INIT_NOD_NUM-1].pt,
// 							m_pNodes[m_pBnds[i].end+INIT_NOD_NUM-1].pt,
// 							m_pNodes[m_pElems[iSrch].form[(index+1)%(DIM+1)]].pt) || 
// 							Collinear(m_pNodes[m_pBnds[i].beg+INIT_NOD_NUM-1].pt,
// 							m_pNodes[m_pBnds[i].end+INIT_NOD_NUM-1].pt,
// 							m_pNodes[m_pElems[iSrch].form[(index+2)%(DIM+1)]].pt))
// 						{
// 							cout<<"***Error: Collinear case found during recovering edge "<<i<<endl;
// 						}
						if (IntersectProp(m_pNodes[m_pBnds[i].beg+INIT_NOD_NUM-1].pt,
							m_pNodes[m_pBnds[i].end+INIT_NOD_NUM-1].pt,
							m_pNodes[m_pElems[iSrch].form[(index+1)%(DIM+1)]].pt,
							m_pNodes[m_pElems[iSrch].form[(index+2)%(DIM+1)]].pt))
						{
// 							SwapEdge edge;
// 							edge.iElemLeft = iSrch;
// 							edge.iElemRight = m_pElems[iSrch].neig[index];
// 							edge.iDiag[0] = m_pElems[iSrch].form[(index+1)%(DIM+1)];
// 							edge.iDiag[DIM-1] = m_pElems[iSrch].form[(index+2)%(DIM+1)];
// 							addSwapEdge(&edge);
// 							iCur = m_pElems[iSrch].neig[index];
// 							iPrev = iSrch;
// 							bSwapFind = true;
// 							break;
							addSwapEdge(iSrch, m_pElems[iSrch].neig[index], 
								m_pElems[iSrch].form[(index+1)%(DIM+1)], m_pElems[iSrch].form[(index+2)%(DIM+1)]);
							iCur = m_pElems[iSrch].neig[index];
							iPrev = iSrch;
							bSwapFind = true;
							break;
						}	
					}
					iSrch = m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
					if (iSrch == NULL_NEIG)
					{
						break;
					}
				} while(iSrch != iElem);
			}
			if (!bSwapFind)
			{
				cout<<"***Error: Cannot find the first edge need swapping during recovering edge "<<i<<endl;
			}
			//Search all other edges that need swapping
			while (1)
			{
				int index;
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iCur].neig[j] == iPrev)
					{
						index = j;
						break;
					}
				}
				if (Collinear(m_pNodes[m_pBnds[i].beg+INIT_NOD_NUM-1].pt,
					m_pNodes[m_pBnds[i].end+INIT_NOD_NUM-1].pt,
					m_pNodes[m_pElems[iCur].form[index]].pt))
				{
					if (m_pElems[iCur].form[index] != m_pBnds[i].end+INIT_NOD_NUM-1)
					{
						cout<<"***Error: Collinear case found during recovering edge "<<i<<endl;
					}
					break;
				}
				if (IntersectProp(m_pNodes[m_pBnds[i].beg+INIT_NOD_NUM-1].pt,
					m_pNodes[m_pBnds[i].end+INIT_NOD_NUM-1].pt,
					m_pNodes[m_pElems[iCur].form[index]].pt,
					m_pNodes[m_pElems[iCur].form[(index+1)%(DIM+1)]].pt))
				{
// 					SwapEdge edge;
// 					iPrev = iCur;
// 					iCur = m_pElems[iCur].neig[(index+DIM)%(DIM+1)];
// 					edge.iElemLeft = iPrev;
// 					edge.iElemRight = iCur;
// 					edge.iDiag[0] = m_pElems[iPrev].form[index];
// 					edge.iDiag[DIM-1] = m_pElems[iPrev].form[(index+1)%(DIM+1)];
// 					addSwapEdge(&edge);
					iPrev = iCur;
					iCur = m_pElems[iCur].neig[(index+DIM)%(DIM+1)];
					addSwapEdge(iPrev, iCur, m_pElems[iPrev].form[index], m_pElems[iPrev].form[(index+1)%(DIM+1)]);
				}
				else if (IntersectProp(m_pNodes[m_pBnds[i].beg+INIT_NOD_NUM-1].pt,
					m_pNodes[m_pBnds[i].end+INIT_NOD_NUM-1].pt,
					m_pNodes[m_pElems[iCur].form[index]].pt,
					m_pNodes[m_pElems[iCur].form[(index+DIM)%(DIM+1)]].pt))
				{
// 					SwapEdge edge;
// 					iPrev = iCur;
// 					iCur = m_pElems[iCur].neig[(index+1)%(DIM+1)];
// 					edge.iElemLeft = iPrev;
// 					edge.iElemRight = iCur;
// 					edge.iDiag[0] = m_pElems[iPrev].form[(index+DIM)%(DIM+1)];
// 					edge.iDiag[DIM-1] = m_pElems[iPrev].form[index];
// 					addSwapEdge(&edge);
					iPrev = iCur;
					iCur = m_pElems[iCur].neig[(index+1)%(DIM+1)];
					addSwapEdge(iPrev, iCur, m_pElems[iPrev].form[(index+DIM)%(DIM+1)], m_pElems[iPrev].form[index]);
				}
				else
				{
					cout<<"***Error: Cannot find next edge need swapping during recovering edge "<<i<<endl;
					exit(0);
				}
			}
			while (!isSwapingFinished())
			{
				bool isswap=swapEdge(i);
				//if (!isswap) {
				//	moveAlong();
				//}
			}
		}
	}	

	return true;
}

/*
 * check whether swapping has finished
 */
bool DTIso2D::isSwapingFinished()
{
	return m_lstSwapEdge.empty();
}

/*
 * swap the current edge
 * -------            -------
 * |\    |            |    /|
 * | \   |            |   / |
 * |  \  |  ------->  |  /  | 
 * |   \ |            | /   | 
 * |    \|            |/    |
 * -------            -------
 */
bool DTIso2D::swapEdge(INTEGER iBnd)
{
	INTEGER iLeft,iRight;
	INTEGER iNeigLeftBeg,iNeigLeftEnd,iNeigRightBeg,iNeigRightEnd;
	INTEGER oldDiag[DIM],newDiag[DIM];	
	int i,j;
	SwapEdge *prevEdge, *nextEdge, *curEdge;

	assert(m_nCurInd>=0 && m_nCurInd<m_lstSwapEdge.size());	

	curEdge = curSwapEdge();
	iLeft = curEdge->iElemLeft;
	iRight = curEdge->iElemRight;
#ifdef _VERBOS
	cout<<"Swap Elem "<<iLeft<<" and Elem "<<iRight<<endl;
#endif
	for (i=0; i<DIM; i++)
	{
		oldDiag[i] = curEdge->iDiag[i];
	}

	/* Find the new diagonal, just rotate the old diagonal clockwise */
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iLeft].form[i] == oldDiag[0])
		{
			newDiag[0] = m_pElems[iLeft].form[(i+DIM)%(DIM+1)];
			break;
		}
	}
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iRight].form[i] == oldDiag[DIM-1])
		{
			newDiag[DIM-1] = m_pElems[iRight].form[(i+DIM)%(DIM+1)];
			break;
		}
	}
	
	if (!IntersectProp(m_pNodes[oldDiag[0]].pt,m_pNodes[oldDiag[DIM-1]].pt,
		m_pNodes[newDiag[0]].pt,m_pNodes[newDiag[DIM-1]].pt))
	{
		//for (j = 0; j < DIM; j++)
		//{
		//	curEdge->iDiag[j] = newDiag[j];
		//}

		moveAlong();

		//Unable to swap these two diagonals
		return false;
	}
	/* Only two outer neighbours should have a new neighbour */
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iLeft].form[i] == oldDiag[DIM-1])
		{
			iNeigLeftEnd = m_pElems[iLeft].neig[i];
			iNeigLeftBeg = m_pElems[iLeft].neig[(i+DIM)%(DIM+1)];
			if (iNeigLeftEnd != NULL_NEIG)
			{
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iNeigLeftEnd].neig[j] == iLeft)
					{
						m_pElems[iNeigLeftEnd].neig[j] = iRight;
						break;
					}
				}
			}
			break;
		}
	}
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iRight].form[i] == oldDiag[0])
		{
			iNeigRightBeg = m_pElems[iRight].neig[i];
			iNeigRightEnd = m_pElems[iRight].neig[(i+DIM)%(DIM+1)];
			if (iNeigRightBeg != NULL_NEIG)
			{
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iNeigRightBeg].neig[j] == iRight)
					{
						m_pElems[iNeigRightBeg].neig[j] = iLeft;
						break;
					}
				}
			}
			break;
		}
	}
	/* Form the new elements*/
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iLeft].form[i] == oldDiag[0])
		{
			m_pElems[iLeft].form[i] = newDiag[DIM-1];
			m_pElems[iLeft].neig[i] = iNeigLeftBeg;
			m_pElems[iLeft].neig[(i+1)%(DIM+1)] = iRight;
			m_pElems[iLeft].neig[(i+DIM)%(DIM+1)] = iNeigRightBeg;
			break;
		}
	}
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iRight].form[i] == oldDiag[DIM-1])
		{
			m_pElems[iRight].form[i] = newDiag[0];
			m_pElems[iRight].neig[i] = iNeigRightEnd;
			m_pElems[iRight].neig[(i+1)%(DIM+1)] = iLeft;
			m_pElems[iRight].neig[(i+DIM)%(DIM+1)] = iNeigLeftEnd;
		}
	}
	/* The previous edge and the next edge maybe need updating*/
	prevEdge = prevSwapEdge();
	if (prevEdge)
	{
		for (j=0; j<=DIM; j++)
		{
			if (m_pElems[prevEdge->iElemLeft].form[j] == prevEdge->iDiag[0])
			{
				prevEdge->iElemRight = m_pElems[prevEdge->iElemLeft].neig[(j+DIM)%(DIM+1)];
				break;
			}
		}
	}
	nextEdge = nextSwapEdge();
	if (nextEdge)		
	{
		for (j=0; j<=DIM; j++)
		{
			if (m_pElems[nextEdge->iElemRight].form[j] == nextEdge->iDiag[DIM-1])
			{
				nextEdge->iElemLeft = m_pElems[nextEdge->iElemRight].neig[(j+DIM)%(DIM+1)];
				break;
			}
		}
	}
	
	/* Check whether new edge is still to be swapped or not*/
	if (IntersectProp(m_pNodes[newDiag[0]].pt,m_pNodes[newDiag[DIM-1]].pt,
		m_pNodes[m_pBnds[iBnd].beg+INIT_NOD_NUM-1].pt,
		m_pNodes[m_pBnds[iBnd].end+INIT_NOD_NUM-1].pt))
	{
		for (j=0; j<DIM; j++)
		{
			curEdge->iDiag[j] = newDiag[j];
		}
		moveAlong();
	}
	else
	{
		/*Check whether new edge is just the boundary edge that needs recovering*/
		if (newDiag[0] == m_pBnds[iBnd].beg+INIT_NOD_NUM-1 &&
			newDiag[DIM-1] == m_pBnds[iBnd].end+INIT_NOD_NUM-1)
		{
			m_pElems[iLeft].iReserved = 1;
			m_pElems[iRight].iReserved = -1;
		}
		else if (newDiag[0] == m_pBnds[iBnd].end+INIT_NOD_NUM-1 &&
			newDiag[DIM-1] == m_pBnds[iBnd].beg+INIT_NOD_NUM-1)
		{
			m_pElems[iRight].iReserved = 1;
			m_pElems[iLeft].iReserved = -1;
		}	
		delSwapEdge();
	}
	
	return true;
}

/*
 * add the edge need swapping
 */
bool DTIso2D::addSwapEdge(SwapEdge *edge)
{
	m_lstSwapEdge.push_back(edge);
	return true;
}

void DTIso2D::addSwapEdge(INTEGER iLeft, INTEGER iRight, INTEGER beg, INTEGER end)
{
	SwapEdge* pEdge = new SwapEdge();
	pEdge->iElemLeft = iLeft;
	pEdge->iElemRight = iRight;
	pEdge->iDiag[0] = beg;
	pEdge->iDiag[1] = end;
	
	addSwapEdge(pEdge);
}

/*
 * current swapping edge
 */
SwapEdge * DTIso2D::curSwapEdge()
{
	INTEGER i;
	std::list<SwapEdge *>::iterator pos;
	for (pos=m_lstSwapEdge.begin(),i=0; i<m_nCurInd; i++,pos++);
	return *(pos);
}

/*
 * previous swapping edge
 */
SwapEdge * DTIso2D::prevSwapEdge()
{
	INTEGER i;
	std::list<SwapEdge *>::iterator pos;
	if (m_nCurInd == 0)
	{
		return NULL;
	}
	for (pos=m_lstSwapEdge.begin(),i=0; i<m_nCurInd-1; i++,pos++);
	return *(pos);
}

/*
 * next swapping edge
 */
SwapEdge * DTIso2D::nextSwapEdge()
{
	INTEGER i;
	std::list<SwapEdge *>::iterator pos;
	if (m_nCurInd == m_lstSwapEdge.size()-1)
	{
		return NULL;
	}
	for (pos=m_lstSwapEdge.begin(),i=0; i<m_nCurInd+1; i++,pos++);
	return *(pos);
}

/*
 * delete current edge
 */
bool DTIso2D::delSwapEdge()
{
	INTEGER i;
	std::list<SwapEdge *>::iterator pos;
	for (pos=m_lstSwapEdge.begin(),i=0; i<m_nCurInd; i++,pos++);
	m_lstSwapEdge.erase(pos);
	if (m_nCurInd == m_lstSwapEdge.size())
	{
		m_nCurInd = 0;
	}
	if(m_nCurInd>=m_lstSwapEdge.size())
		m_nCurInd = 0;
	return true;
}

/*
 * move to next edge
 */
bool DTIso2D::moveAlong()
{
	m_nCurInd++;
	if (m_nCurInd == m_lstSwapEdge.size())
	{
		m_nCurInd = 0;
	}
	return true;
}

/*
 * clear all outer elements
 */
bool DTIso2D::clrOuterEles()
{
	INTEGER iEle,iNeig,iBnd,iElem,iSrch;
	int i,j,k,loop=0;

	//First, set the flags, 1 for inner , -1 for outer
	for (iBnd=0; iBnd<m_nBnds; iBnd++)
	{
		bool bFind = false;
		for (iEle=0; iEle<m_nElems; iEle++)
		{
			if (isDelEle(iEle))
			{
				continue;
			}
			for (k=0; k<=DIM; k++)
			{
				if (m_pElems[iEle].form[k] == m_pBnds[iBnd].beg+INIT_NOD_NUM-1)
				{
					iElem = iEle;
					bFind = true;
					break;
				}
			}
			if (bFind)
			{
				break;
			}
		}
		if (!bFind)
		{
			cout<<"***Error: Cannot find element for boundary "<<iBnd<<endl;
			continue;
		}
		// then, search the edge
		bFind = false;		
		iSrch = iElem;
		do {
			int index;
			for (j=0; j<=DIM; j++)
			{
				if (m_pElems[iSrch].form[j] == m_pBnds[iBnd].beg+INIT_NOD_NUM-1)
				{
					index = j;
					break;
				}
			}
			if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[iBnd].end+INIT_NOD_NUM-1)
			{
				//Find this bound edge
				m_pBnds[iBnd].ele = iSrch;
				m_pElems[iSrch].iReserved = 1;
				iNeig =	m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
				m_pElems[iNeig].iReserved = -1;
				bFind = true;
				break;
			}
			iSrch = m_pElems[iSrch].neig[(index+1)%(DIM+1)];
			if (iSrch == NULL_NEIG)
			{
				break;
			}
		} while(iSrch != iElem);
		if (!bFind)
		{
			//Try the other direction
			iSrch = iElem;
			do {
				int index;
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iSrch].form[j] == m_pBnds[iBnd].beg+INIT_NOD_NUM-1)
					{
						index = j;
						break;
					}
				}
				if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[iBnd].end+INIT_NOD_NUM-1)
				{
					//Find this bound edge
					m_pBnds[iBnd].ele = iSrch;
					m_pElems[iSrch].iReserved = 1;
					iNeig =	m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
					m_pElems[iNeig].iReserved = -1;
					bFind = true;
					break;
				}
				iSrch = m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
				if (iSrch == NULL_NEIG)
				{
					break;
				}
			} while(iSrch != iElem);
		}
		if (!bFind)
		{
			cout<<"***Error: Boundary edge "<<iBnd<<" is missing!"<<endl;
		}
	}
//	for (iElem=0; iElem<m_nElems; iElem++)
//	{
//		if (m_pElems[iElem].iReserved == -1)
//		{
//			addTreeSearch(iElem);
//		}
//	}

	for (iElem=0; iElem<m_nElems; iElem++)
	{
		if (!isDelEle(iElem) && m_pElems[iElem].iReserved == -1)
		{
			addTreeSearch(iElem);
			addDeleted(iElem);
			while (!isTreeSearchEmpty())
			{
				iEle = pickTreeSearch();				
				for (i=0; i<=DIM; i++)
				{
					iNeig = m_pElems[iEle].neig[i];
					if (iNeig != NULL_NEIG && !isDelEle(iNeig))
					{
						if (m_pElems[iNeig].iReserved == 1)
						{
							for (j=0; j<=DIM; j++)
							{
								if (m_pElems[iNeig].neig[j] == iEle)
								{
									m_pElems[iNeig].neig[j] = NULL_NEIG;
									break;
								}
							}			
						}
						else /*if (m_pElems[iNeig].iReserved == 0)*/
						{
							if (m_pElems[iNeig].iReserved == 0)
							{
								m_pElems[iNeig].iReserved = -1;
							}							
							addTreeSearch(iNeig);
							addDeleted(iNeig);
						}
					}
				}
			}
		}
	}
	updateNewEleLoc();
	clrDelEles();
	clrTstEles();
	return true;
}

/*
 * add inner point
 */
int DTIso2D::addInnerPnt(INTEGER iNod, INTEGER iEle)
{
	/* make sure this point is in the element*/
	int i,k,m;
	INTEGER iCurNod, iNextNod;
	INTEGER iElem, iSrch;
	Elem *pElem = NULL, *pSrch = NULL;
	INTEGER nnew[DIM+1];
	POINT pnew[DIM+1]; 
	POINT cen;
	REAL d_a, rad;
	INTEGER iEmp; /*empty position for new created element*/
	INTEGER iTakEle;
	int nfc = 0, npc = 0;
	int nBegin = 0;
	int nCase = 1;
	float dist;
	bool bValid = true;
	bool bAccept = true;
	if (isDelEle(iEle))
	{
		//rejected, the element has been deleted
		return 3;
	}
	for (i=0; i<=DIM; i++)
	{
		iCurNod = m_pElems[iEle].form[i];
		iNextNod = m_pElems[iEle].form[(i+1)%(DIM+1)];
		if (!Left(m_pNodes[iCurNod].pt, m_pNodes[iNextNod].pt, m_pNodes[iNod].pt))
		{
			bValid = false;
			break;
		}
	}
	if (!bValid)
	{
		//rejected, maybe the element has been deleted and reused
		return 3;
	}
//	addTstEle(iEle);
	addDeleted(iEle);
	addTreeSearch(iEle);
	int ii=0;
//	char str[20];
//	sprintf(str, "%d.dt2", ii++);
//	writeDt2(str);
	while (!isTreeSearchEmpty())
	{
		iElem = pickTreeSearch();
		pElem = &(m_pElems)[iElem];
		for (i = 0; i <= DIM; i++)
		{
			iSrch = (pElem->neig)[i];
			pSrch = &(m_pElems)[iSrch];
			if (iSrch == NULL_NEIG)
			{
				nBegin = 1;
			}
			else if (/*!isTstEle(iSrch) && */!isDelEle(iSrch))
			{
//				addTstEle(iSrch);
				REAL dt = distSqua(pSrch->cen, m_pNodes[iNod].pt); /*distance square between two points*/
				REAL cri = dt - pSrch->rad;
				if (cri < -EPS_ZERO_SQ)
				{/*incircle cirteria is broken*/
					addDeleted(iSrch);
					addTreeSearch(iSrch);
//					sprintf(str, "%d.dt2", ii++);
//					writeDt2(str);
				}		
				else if (cri > EPS_ZERO_SQ)/*incircle criteria is kept, perform tree search*/
				{
					nBegin = 1;
				}//else if (cri > EPS_ZERO_SQ)
				else
				{/* four points in the same circle */
					nCase = 0;
					goto INNER_RECOVERED;
				}
			}/*if (!isTested(pSrch) && isValid(pSrch))*/
			if (nBegin)
			{
				/*find the nonshared node index*/
//				for (j = 0; j <= DIM; j++)
//				{
//					no = (pElem->form)[j];
//					if (no != (pSrch->form)[0] && no != (pSrch->form)[1] && no != (pSrch->form)[2])
//						break;
//				}
//#ifdef _ERROR_CHK
//				if (j > DIM)
//				{
//					printf("can not find uncommon point!!\n");
//					//go to error handle
//				}
//					
//#endif //_ERROR_CHK				
				for (k = 0; k < DIM; k++)
				{
					nnew[k] = (pElem->form)[(i+k+1)%(DIM+1)];//[(j+k+1)%(DIM+1)];
					for (m = 0; m < DIM; m++)
						pnew[k][m] = (m_pNodes[nnew[k]].pt)[m];
				}
				nnew[DIM] = iNod; /*add a node*/
				for (m = 0; m < DIM; m++)
					pnew[DIM][m] = m_pNodes[iNod].pt[m];
				//different from boundary node inserting, must check the density
				for (k=0; k<DIM; k++)
				{
					dist = sqrt(distSqua(m_pNodes[nnew[k]].pt,m_pNodes[iNod].pt));
					if (/*dist<g_alpha*m_pNodes[iNod].density*/ true /*||*/&& dist<g_alpha*m_pNodes[nnew[k]].density)
					{
						bAccept = false;
						break;
					}
				}
				if (!bAccept)
				{
					//rejected
					nCase = 3;
					goto INNER_RECOVERED;
				}
				calcElePar(pnew, &d_a, &cen, &rad);
//#ifdef _ERROR_CHK
				if (d_a < 0.0)
				{
					printf("negative volume\n!!");
					nCase = 3;
					goto INNER_RECOVERED;
					//go to error handle
				}
				if (d_a < EPS_ZERO_SQ)
				{
					printf("volume is near zero!!\n");
					nCase = 3;
					goto INNER_RECOVERED;
					//go to error handle
				}
//#endif //_ERROR_CHK
				nfc = nfc + 1;
				
				iEmp = crtEle(nnew, cen, rad); /*create a new element & return its location*/
				((m_pElems)[iEmp].neig)[DIM]/*[0]*/ = iSrch;
//				sprintf(str, "%d.dt2", ii++);
//				writeDt2("Temp.dt2");
        		
				//maintain a pointer to parent element for inner point.
				//these code can be improved in future
				// TEMP [2/13/2006]
//				for (k=0; k<=DIM; k++)
//					m_pNodes[nnew[k]].iReserved = iEmp; //This causes problem and also just "m_pNodes[nnew[k]].iReserved = iEmp;"

   				/*
   				 * check two nodes of the created element. If they are taken by a new created element,
   				 * then update the neighboring information correspondingly
   				 */
   				for (k = 0; k < DIM; k++)
   				{
   					iTakEle = getNodTakEle(nnew[k]);
   					if (iTakEle >= 0)
   					{//taken
   						((m_pElems)[iEmp].neig)[(k+1)%DIM] = iTakEle;
   						if (((m_pElems)[iTakEle].form)[0] == nnew[k])
   			 				((m_pElems)[iTakEle].neig)[1] = iEmp;
   						else
   			 				((m_pElems)[iTakEle].neig)[0] = iEmp;
   					}
   					else
   					{
   						addNodTakEle(nnew[k], iEmp);
						npc = npc + 1;
   					}
				}
				nBegin = 0;
			}
		}/*for (i = 0; i < DIM; i++)*/
	}/*while (!isTreeSearchEmpty())*/
	
#ifdef _ERROR_CHK
	if (nfc != (DIM - 1) * npc - 4 * (DIM - 2))
	{
		printf("Node %d: the number of faces & the number of nodes is imbanlant!!\n",iNod);
		//go to error handle
	}
#endif //_ERROR_CHK

	/*update neighboring info. of elements near the cavity*/
	updateNeigInfo();
	updateNewEleLoc();
	goto INNER_FINSHED;
INNER_RECOVERED:
	/*error happens, try to recover*/
	undoInsertion();
	clrTreeSearch();
INNER_FINSHED:
	clrAddEles();
	clrDelEles();
	//clear test flags
//	clrTstEles();
	//clear taken flags
	clrTakNods();
	if (nCase == 1) //successful
	{
	}
	return nCase;
}

/*
 * create new inner points, using gravity center
 */
int DTIso2D::crtInnerPnts()
{
	INTEGER iEle,iNod,iForm;
	int j,k;
	REAL dist, spac;
	bool bAccept;
	iNod = m_nNodes;
	m_lstInstInnerNods.clear();
	for (iEle=0; iEle<m_nElems; iEle++)
	{
		if (!isDelEle(iEle)/* && !isTstEle(iEle)*/)//comment this check //  [11/7/2005]
		{
			for (k=0; k<DIM; k++)
			{
				m_pNodes[iNod].pt[k] = 0;
			}
			m_pNodes[iNod].density = 0;
			for (j=0; j<=DIM; j++)
			{
				iForm = m_pElems[iEle].form[j];
				for (k=0; k<DIM; k++)
				{
					m_pNodes[iNod].pt[k] += m_pNodes[iForm].pt[k];
				}
				m_pNodes[iNod].density += m_pNodes[iForm].density;
			}
			for (k=0; k<DIM; k++)
			{
				m_pNodes[iNod].pt[k] /= (DIM+1);
			}
			m_pNodes[iNod].density /= (DIM+1);

			g_nCreatePnts++;		
#if 0
			spac = spacFrmSrc(m_pNodes[iNod].pt);
			m_pNodes[iNod].density = Min(spac, m_pNodes[iNod].density);
#endif
			m_pNodes[iNod].iReserved = 0;
			bAccept = true;
			for (j=0; j<=DIM; j++)
			{
				iForm = m_pElems[iEle].form[j];
				dist = sqrt(distSqua(m_pNodes[iForm].pt,m_pNodes[iNod].pt));
				if (/*dist < g_alpha*m_pNodes[iNod].density
					&&*/ dist < g_alpha*m_pNodes[iForm].density) // add this check //  [11/7/2005]
				{
					bAccept = false;
					//enabEleTst(iEle, true); //comment this //  [11/7/2005]
					break;
				}
			}
			if (bAccept) 
			{
				InnerPnt *ppnt = (InnerPnt *)malloc(sizeof(InnerPnt));
				ppnt->iEle = iEle;
				ppnt->iNod = iNod;
				
//				m_pElems[iEle].iReserved = iNod; //  [2/13/2006]
				
				m_lstInstInnerNods.push_back(ppnt);
				
				iNod++;
				g_nAcceptPnts++;
			}
		}		
	}
	m_nInners = iNod - m_nNodes;
	return m_nInners;
}

int DTIso2D::getMaxiNod()
{
		assert( m_lstInstInnerNods.size() > 0 );
		std::list<InnerPnt *>::iterator it = m_lstInstInnerNods.begin();
		int iMaxNod = (*it)->iNod;
		for(;it!=m_lstInstInnerNods.end();it++)
		{
				if( (*it)->iNod > iMaxNod )
						iMaxNod = (*it)->iNod;
		}
		return iMaxNod;
}

/*
 * inner point insertion
 */
int DTIso2D::innerPntInst()
{
	int iSucc = 0, i, iFail = 0, j = 0;
	INTEGER iInner;
	INTEGER iNod,iEle;
	POINT pnt;
	int nResult;
	int i_NotFind = 0, iMaxNod = 0;
	std::list<InnerPnt *>::iterator it_fir, it, it_old, it_last, it_e; 
	InnerPnt *pInnerPnt,*pLastPnt;

	if (crtInnerPnts() == 0)
		return 0;
	
	iInner = m_nNodes;
	m_nNodes += m_nInners;
 	while (iSucc != m_nInners)
	{
		assert(!m_lstInstInnerNods.empty());
		it_fir = m_lstInstInnerNods.begin();
		pInnerPnt = *(it_fir);
		iNod = pInnerPnt->iNod;
		iEle = pInnerPnt->iEle;
		if (isDelEle(iEle))
		{
			//discard current node and replace with the last node
			if (iNod == m_nNodes-1)
			{
				delete (*it_fir);
				m_lstInstInnerNods.erase(it_fir);	
				m_nInners--;
				m_nNodes--;
			}
			else
			{
				it_last = m_lstInstInnerNods.end();
				it_last--;
				pLastPnt = *it_last;
				assert(it_last != m_lstInstInnerNods.begin());
				assert(pLastPnt->iNod == m_nNodes - 1);
				m_pNodes[iNod] = m_pNodes[pLastPnt->iNod];
				pInnerPnt->iEle = pLastPnt->iEle;
				delete pLastPnt;
				m_lstInstInnerNods.erase(it_last);
				m_nInners--;
				m_nNodes--;
			}
			continue;
		}
		nResult = addInnerPnt(iNod, iEle);

// 		nErrCode = checkGlobalData();
// 		if (nErrCode != 0)
// 		{
// 			printError(nErrCode);
// 			exit(1);
// 		}

		if (nResult == 1)
		{
			delete (*it_fir);
			m_lstInstInnerNods.erase(it_fir);	
			iSucc++;
			iFail = 0;
			recDistPnts();
			m_pNodes[iNod].density = Min(m_pNodes[iNod].density, spacFrmSrc(m_pNodes[iNod].pt));
		}
		else /*if (nResult == 3 )*/
		{
			it_last = m_lstInstInnerNods.end();
			it_last--;
			pLastPnt = *it_last;
			assert(pLastPnt->iNod == m_nNodes - 1);
			if (m_lstInstInnerNods.size() > 1)
			{
				assert(it_last != m_lstInstInnerNods.begin());
				m_pNodes[iNod] = m_pNodes[pLastPnt->iNod];
				pInnerPnt->iEle = pLastPnt->iEle;
				delete pLastPnt;
				m_lstInstInnerNods.erase(it_last);
				m_nInners--;
				m_nNodes--;
			}
			else
			{
				m_pNodes[pInnerPnt->iNod].iReserved = -1;
				delete pInnerPnt;
				m_lstInstInnerNods.erase(it_fir);
				m_nInners--;				
			}
		}
// 		else
// 		{
// 			iFail++;
// 			if (iFail >= m_lstInstInnerNods.size())
// 			{
// 				if (!isDisturbed(iNod))
// 				{
// 					for (i=0; i<DIM; i++)
// 					{
// 						pnt[i] = m_pNodes[iNod].pt[i];
// 					}
// 					addDistInfo(iNod, pnt);
// 				}
// 				//disturb the point
// 				for (i = 0; i < DIM; i++)
// 				{
// 					(m_pNodes[iNod].pt)[i] = (m_pNodes[iNod].pt)[i] * (1.0 + EPS_DISTURB);
// 					(m_pNodes[iNod].pt)[i] = (m_pNodes[iNod].pt)[i] * (1.0 + EPS_DISTURB);
// 				}
// 			}
// 			else //delay the insertion of iNod
// 			{
// 				m_lstInstInnerNods.erase(it_fir);
// 				m_lstInstInnerNods.push_back(pInnerPnt);
// 			}
// 		}
 	}
	clrDelEles();
	clrTstEles();
	return iSucc;
}

REAL DTIso2D::spacFrmSrc(POINT pnt)
{
	INTEGER j;
	REAL st = 1.0e12, spac;

	for (j=0; j<m_nPntSrcNum; j++)
	{
		spac = spacFrmPnt(m_pntSources[j], pnt);
		st = Min(st, spac);
	}
	/* line source control */
	for (j=0; j<m_nLneSrcNum; j++)
	{
		spac = spacFrmLne(m_lineSources[j].points[0], m_lineSources[j].points[1], pnt);
		st = Min(st, spac);
	}
#ifdef _3D_DECOM_
	/* triangle source control */
	for (j=0; j<m_nTriSrcNum; j++)
	{
		spac = spacFrmTri(m_triSources[j].points[0], m_triSources[j].points[1],
			m_triSources[j].points[2], pnt);
		st = Min(st, spac);
	}
#endif

	return st;
}
/*
 * computes the spacing at a point from a source point
 */ 
REAL DTIso2D::spacFrmPnt(PointSource src, POINT pt)
{
	REAL diff, ae;
	REAL dist ;
	const REAL BIG = 50;
	REAL CLG2 = log(2.0);
#ifdef _3D_DECOM_
	POINT3D pt2;
	pt2[0] = pt[0];
	pt2[1] = pt[1];
	pt2[2] = 0;
	
	dist = sqrt(distSqua_3D(src.pt, pt2));
#else
	dist = sqrt(distSqua(src.pt, pt));
#endif
	if (dist <= src.rInnerRad)
	{
		return src.rIntensity;
	}
	else
	{
		dist = dist-src.rInnerRad;
		diff = src.rOuterRad - src.rInnerRad;
		if (diff <= 0)
		{
			cout<<"***Error: Error point source"<<endl;
			return -1;
		}
		else
		{
			diff = CLG2/diff;
			ae = Min(dist*diff, BIG);
			return src.rIntensity*exp(ae);
		}
	}
	return 0;
}

/*
 * computes the spacing at a point from a source line
 */ 
REAL DTIso2D::spacFrmLne(PointSource src1, PointSource src2, POINT pt)
{
	REAL sca = 0;
	REAL w1, w2;
	REAL tolg = 1.0e-5;
	PointSource src3;
	INTEGER i;
	REAL dist;
#ifdef _3D_DECOM_
	dist = sqrt(distSqua_3D(src1.pt, src2.pt));
#else
	dist = sqrt(distSqua(src1.pt, src2.pt));
#endif
	if (dist < tolg)
	{
		cout<<"***Error: Error distance in line source."<<endl;
		return -1;
	}
#ifdef _3D_DECOM_
	for (i=0; i<3; i++)
	{
		src3.pt[i] = (src2.pt[i] - src1.pt[i])/dist;		
	}
	POINT3D pt2;
	pt2[0] = pt[0];
	pt2[1] = pt[1];
	pt2[2] = 0;
	for (i=0; i<3; i++)
	{
		sca += (pt2[i] - src1.pt[i])*src3.pt[i];
	}
#else
	for (i=0; i<2; i++)
	{
		src3.pt[i] = (src2.pt[i] - src1.pt[i])/dist;		
	}
	for (i=0; i<2; i++)
	{
		sca += (pt[i] - src1.pt[i])*src3.pt[i];
	}
#endif
	if (sca <= 0)
	{
		return spacFrmPnt(src1,pt);
	}
	else if (sca >= dist)
	{
		return spacFrmPnt(src2,pt);
	}
	else
	{
		w2 = sca/dist;
		w1 = 1 - w2;
#ifdef _3D_DECOM_
		for (i=0; i<3; i++)
		{
			src3.pt[i] = w1*src1.pt[i] + w2*src2.pt[i];
		}
#else
		for (i=0; i<2; i++)
		{
			src3.pt[i] = w1*src1.pt[i] + w2*src2.pt[i];
		}
#endif
		src3.rInnerRad = w1*src1.rInnerRad + w2*src2.rInnerRad;
		src3.rOuterRad = w1*src1.rOuterRad + w2*src2.rOuterRad;
		src3.rIntensity = w1*src1.rIntensity + w2*src2.rIntensity;
		return spacFrmPnt(src3,pt);
	}
	return 0;
}
#ifdef _3D_DECOM_
/*
 * computes the spacing at a point from a source triangle
 */ 
REAL DTIso2D::spacFrmTri(PointSource src1, PointSource src2, PointSource src3, POINT pt)
{
	REAL sca = 0;
	REAL w1, w2, w3, sh1=0.0, sh3=0.0, h1=0.0, h3=0.0;
	REAL tolg = 1e-5;
	PointSource src4;
	POINT3D zero,x1,x2,x3,x4,x5,x6;
	int i;
	REAL dist;

	for (i=0; i<3; i++)
	{
		zero[i] = 0.0;
	}

	for (i=0; i<3; i++)
	{
		x1[i] = src2.pt[i] - src1.pt[i];
	}
	dist = sqrt(distSqua_3D(zero, x1));
	if (dist < tolg)
	{
		cout<<"***Error: Error distance in triangle source."<<endl;
		return -1;
	}

	for (i=0; i<3; i++)
	{
		x2[i] = src3.pt[i] - src2.pt[i];
	}
	dist = sqrt(distSqua_3D(zero, x2));
	if (dist < tolg)
	{
		cout<<"***Error: Error distance in triangle source."<<endl;
		return -1;
	}

	for (i=0; i<3; i++)
	{
		x3[i] = src1.pt[i] - src3.pt[i];
	}
	dist = sqrt(distSqua_3D(zero, x3));
	if (dist < tolg)
	{
		cout<<"***Error: Error distance in triangle source."<<endl;
		return -1;
	}

	//normal to the triangle xn=x1^x2
	x6[0] = x1[1]*x2[2]-x1[2]*x2[1];
    x6[1] = x1[2]*x2[0]-x1[0]*x2[2];
    x6[2] = x1[0]*x2[1]-x1[1]*x2[0];
	
	//normal to side 1-2 x4=xn^x1
    x4[0] = x6[1]*x1[2]-x6[2]*x1[1];
    x4[1] = x6[2]*x1[0]-x6[0]*x1[2];
    x4[2] = x6[0]*x1[1]-x6[1]*x1[0];

	//normal to side 2-3 x5=xn^x2
    x5[0] = x6[1]*x2[2]-x6[2]*x2[1];
    x5[1] = x6[2]*x2[0]-x6[0]*x2[2];
    x5[2] = x6[0]*x2[1]-x6[1]*x2[0];

	for (i=0; i<3; i++)
	{
		x6[i] = pt[i] - src2.pt[i];
	}
	
	//compute area weights, using the point projected to the triangle 
	for (i=0; i<3; i++)
	{
		h1 += x6[i]*x5[i];
		sh1 += -x1[i]*x5[i];
		h3 += x6[i]*x4[i];
		sh3 += x2[i]*x4[i];
	} 
	w1 = h1/sh1;
	w3 = h3/sh3;
	w2 = 1.0-w1-w3;
	
	if (w1 <= 0)
	{
		return spacFrmLne(src2, src3, pt);
	}
	else if (w2 <= 0)
	{
		return spacFrmLne(src1, src3, pt);
	}
	else if (w3 <= 0)
	{
		return spacFrmLne(src1, src2, pt);
	}
	else
	{
		for (i=0; i<3; i++)
		{
			src4.pt[i] = w1*src1.pt[i]+w2*src2.pt[i]+w3*src3.pt[i];
		}
		src4.rInnerRad = w1*src1.rInnerRad+w2*src2.rInnerRad+w3*src3.rInnerRad;
		src4.rOuterRad = w1*src1.rOuterRad+w2*src2.rOuterRad+w3*src3.rOuterRad;
		src4.rIntensity = w1*src1.rIntensity+w2*src2.rIntensity+w3*src3.rIntensity;
		return spacFrmPnt(src4,pt);
	}
	return 0;
}
#endif

/*
 * computes the spacing at a point from background mesh
 */ 
REAL DTIso2D::spacFrmBkGrnd(POINT pt)
{
	INTEGER iEle;
	INTEGER iCurNod, iPrevNod, iNextNod;
	int i;
	Elem* pElem;
	int nFind = 0;
	REAL area, subArea[DIM+1];
	REAL spac = 0;
	REAL BIG = 10000;
	for (iEle=0; iEle<m_nBkGrndElem; iEle++)
	{
		nFind = 1;
		pElem = &(m_pBkGrndElem[iEle]);
		for (i=0; i<=DIM; i++)
		{
			iCurNod = pElem->form[i] - 1;
			iNextNod = pElem->form[(i+1)%(DIM+1)] - 1;
			if (!Left(m_pBkGrndNode[iCurNod].pt,m_pBkGrndNode[iNextNod].pt,pt))
			{
				nFind = 0;
				break;
			}
		}
		if (nFind)
		{
			break;
		}
	}
	if (nFind)
	{
		area = Area2(m_pBkGrndNode[m_pBkGrndElem[iEle].form[0]-1].pt, 
			m_pBkGrndNode[m_pBkGrndElem[iEle].form[1]-1].pt,
			m_pBkGrndNode[m_pBkGrndElem[iEle].form[2]-1].pt);
		for (i=0; i<=DIM; i++)
		{
			iCurNod = m_pBkGrndElem[iEle].form[i] - 1;
			iPrevNod = m_pBkGrndElem[iEle].form[(i+2)%(DIM+1)] - 1;
			iNextNod = m_pBkGrndElem[iEle].form[(i+1)%(DIM+1)] - 1;
			subArea[i] = Area2(m_pBkGrndNode[iNextNod].pt, m_pBkGrndNode[iPrevNod].pt,pt);
			spac += subArea[i]/area*m_pBkGrndNode[iCurNod].density;
		}
		return spac;
	}
	return BIG;
}

/* recover disturbed points */
bool DTIso2D::recDistPnts()
{
	std::list<DistInfo *>::iterator it;
	DistInfo *pDI;
	INTEGER iNod;
	int i;
	while (!m_lstDistInfo.empty())
	{
		it = m_lstDistInfo.begin();
		pDI = *it;
		iNod = pDI->iNod;
		for (i=0; i<DIM; i++)
		{
			m_pNodes[iNod].pt[i] = pDI->old_pt[i];
		}
		delete pDI;
		m_lstDistInfo.erase(it);
	}
	return true;
}

/*
 * remove empty nodes
 */
int DTIso2D::rmvEmpNods()
{
	INTEGER i,j,iRmvCnt;
	int m;
	int nRmv;

	for (i=INIT_NOD_NUM,j=INIT_NOD_NUM; i<m_nNodes; i++)
	{
		nRmv = 1;
		if (m_pNodes[i].iReserved != -1)
		{
			nRmv = 0;		
		}
		if (nRmv == 1)
		{
			j++;
		}
		else
		{
			m_pNodes[i].iReserved = i-j;
		}
	}
	iRmvCnt = j;

	for (i=INIT_NOD_NUM,j=0; i<m_nNodes; i++)
	{
		if (m_pNodes[i].iReserved != -1)
		{
			m_pNodes[j].density = m_pNodes[i].density;
			for (m=0; m<DIM; m++)
			{
				m_pNodes[j].pt[m] = m_pNodes[i].pt[m];
			}
			j++;
		}
	}
	m_nNodes -= iRmvCnt;

	//update boundary, .beg and .end are correct index of nodes array now
	for (i=0; i<m_nBnds; i++)
	{
		j = m_pBnds[i].beg;
		m_pBnds[i].beg = m_pNodes[j+INIT_NOD_NUM-1].iReserved;
		j = m_pBnds[i].end;
		m_pBnds[i].end = m_pNodes[j+INIT_NOD_NUM-1].iReserved;
	}
	return 1;
}

/*
 * remove empty elements
 */
int DTIso2D::rmvEmpEles()
{
	INTEGER i,j,iNeig,iForm,iBnd,iElem,iRmvCnt;
	int m;
	int nRmv;

	for (i=0,j=0; i<m_nElems; i++)
	{
		nRmv = 1;
		if (!isDelEle(i))
		{
			nRmv = 0;
			for (m=0; m<=DIM; m++)
			{
				if (m_pElems[i].form[m] < INIT_NOD_NUM)
				{
					nRmv = 1;
					break;
				}
			}
		}
		if (nRmv == 1)
		{
			j++;
			m_pElems[i].iReserved = NULL_ELEM;
		}
		else
		{
			m_pElems[i].iReserved = i-j;
		}
	}
	iRmvCnt = j;

	for (i=0,j=0; i<m_nElems; i++)
	{
		if (m_pElems[i].iReserved != NULL_ELEM)
		{
			m_pElems[j].rad = m_pElems[i].rad;
			for (m=0; m<DIM; m++)
			{
				m_pElems[j].cen[m] = m_pElems[i].cen[m];
			}

			for (m=0; m<=DIM; m++)
			{
				iForm = m_pElems[i].form[m];
				m_pElems[j].form[m] = m_pNodes[iForm].iReserved;
				iNeig = m_pElems[i].neig[m];
				if (iNeig != NULL_NEIG)
				{
					m_pElems[j].neig[m] = m_pElems[iNeig].iReserved;
				}
				else
				{
					m_pElems[j].neig[m] = NULL_NEIG;
				}
			}
			j++;
		}
	}
	m_nElems -= iRmvCnt;

	//update boundary's parent element info
	for (iBnd=0; iBnd<m_nBnds; iBnd++)
	{
		iElem = m_pBnds[iBnd].ele;
		if (iElem == NULL_ELEM)
		{
			cout<<"***Error: Boundary edge "<<iBnd<<"'s parent element doesn't exist!"<<endl;
			continue;
		}
		m_pBnds[iBnd].ele = m_pElems[iElem].iReserved;
	}

	return 1;
}

/*
 * smoothing
 */
int DTIso2D::smooth()
{
	int TOTALTIMES = 3;//10;
	int i,m,nEleCnt,index;
	INTEGER iNod,iEle,iSrch,iForm;

	for(iEle=0; iEle<m_nElems; iEle++)
	{
		if (!isDelEle(iEle))
		{
			for (m=0; m<=DIM; m++)
			{
				m_pNodes[m_pElems[iEle].form[m]].iReserved = iEle;
			}
		}
	}//  [2/13/2006]
		
//	for (iNod=0; iNod<m_nNodes; iNod++)
//	{
//		m_pNodes[iNod].iReserved = 0;
//	} //  [2/13/2006]

//	for (iBnd=0; iBnd < m_nBnds; iBnd++)
//	{
//		iNod = m_pBnds[iBnd].beg+INIT_NOD_NUM-1;
//		m_pNodes[iNod].iReserved = NULL_ELEM;
//	} //  [2/13/2006]
	
//	for (iNod=0; iNod<m_nNodes; iNod++)
//	{
//		cout<<iNod<<":\t"<<m_nNodes<<endl;
//		if (m_pNodes[iNod].iReserved != NULL_ELEM)
//		{
//			nFind = 0;
//			for (iEle=0; iEle<m_nElems; iEle++)
//			{
//				for (m=0; m<=DIM; m++)
//				{
//					if (m_pElems[iEle].form[m] == iNod)
//					{
//						nFind = 1;
//						break;
//					}
//				}
//				if (nFind)
//				{
//					m_pNodes[iNod].iReserved = iEle;
//					break;
//				}
//			}
//		}
//	} //these information has been generated during inner point insertion[2/13/2006]

	for (i=0; i<TOTALTIMES; i++)
	{
		for (iNod=m_nBnds+INIT_NOD_NUM/*0*/; iNod<m_nNodes; iNod++) //  [2/13/2006]
		{
			if (m_pNodes[iNod].iReserved != NULL_ELEM)
			{
				nEleCnt = 0;
				for (m=0; m<DIM; m++)
				{
					m_pNodes[iNod].pt[m] = 0;
				}
				iEle = m_pNodes[iNod].iReserved;
#ifdef _ERROR_CHK
				for (m=0; m<=DIM; m++)
				{
					if (m_pElems[iEle].form[m] == iNod)
					{
						break;
					}
				}
				assert(m<=DIM);
#endif
				iSrch = iEle;
				do {
					for (m=0; m<=DIM; m++)
					{
						if (m_pElems[iSrch].form[m] == iNod)
						{
							break;
						}
					}
					assert(m != DIM+1);
					index = (m+1)%(DIM+1);
					iForm = m_pElems[iSrch].form[index];
					for (m=0; m<DIM; m++)
					{
						m_pNodes[iNod].pt[m] += m_pNodes[iForm].pt[m];
					}
					iSrch = m_pElems[iSrch].neig[index];
					assert(iSrch != NULL_ELEM);
					nEleCnt++;
				} while(iSrch != iEle);
				for (m=0; m<DIM; m++)
				{
					m_pNodes[iNod].pt[m] /= nEleCnt;
				}
			}
		}
	}
	return 1;
}

/*
 * output information
 */
int DTIso2D::output()
{
	cout<<"Total Elements: "<<m_nElems<<endl;
	cout<<"Total Nodes: "<<m_nNodes<<endl;
	return 0;
}

int DTIso2D::updateBndParent()
{
	INTEGER iEle,iNeig,iBnd,iElem,iSrch;
	int j,k,loop=0;

	for (iBnd=0; iBnd<m_nBnds; iBnd++)
	{
		bool bFind = false;
		for (iEle=0; iEle<m_nElems; iEle++)
		{
			assert(!isDelEle(iEle));
			for (k=0; k<=DIM; k++)
			{
				if (m_pElems[iEle].form[k] == m_pBnds[iBnd].beg)//Notice that .beg is the correct index of m_pNode array now 
				{
					iElem = iEle;
					bFind = true;
					break;
				}
			}
			if (bFind)
			{
				break;
			}
		}
		if (!bFind)
		{
			cout<<"***Error: Cannot find element for boundary "<<iBnd<<endl;
			continue;
		}
		// then, search the edge
		bFind = false;		
		iSrch = iElem;
		do {
			int index;
			for (j=0; j<=DIM; j++)
			{
				if (m_pElems[iSrch].form[j] == m_pBnds[iBnd].beg)
				{
					index = j;
					break;
				}
			}
			if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[iBnd].end)
			{
				//Find this bound edge
				m_pBnds[iBnd].ele = iSrch;
				m_pElems[iSrch].iReserved = 1;
				iNeig =	m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
				assert(iNeig == NULL_NEIG);
				bFind = true;
				break;
			}
			iSrch = m_pElems[iSrch].neig[(index+1)%(DIM+1)];
			if (iSrch == NULL_NEIG)
			{
				break;
			}
		} while(iSrch != iElem);
		if (!bFind)
		{
			//Try the other direction
			iSrch = iElem;
			do {
				int index;
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iSrch].form[j] == m_pBnds[iBnd].beg)
					{
						index = j;
						break;
					}
				}
				if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[iBnd].end)
				{
					//Find this bound edge
					m_pBnds[iBnd].ele = iSrch;
					m_pElems[iSrch].iReserved = 1;
					iNeig =	m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
					assert(iNeig == NULL_NEIG);
					bFind = true;
					break;
				}
				iSrch = m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
				if (iSrch == NULL_NEIG)
				{
					break;
				}
			} while(iSrch != iElem);
		}
		if (!bFind)
		{
			cout<<"***Error: Boundary edge "<<iBnd<<" is missing!"<<endl;
		}
	}
	return 0;
}

/************************************************************************/
/* Extention interface for PMG_SDG                                      */
/************************************************************************/
int DTIso2D::meshGen()
{
	POINT minW,maxW,minN,maxN;
	minN[0] = minN[1] = -1;
	maxN[0] = maxN[1] = 1;

	calcBox(&minW, &maxW);
	scaleFactor(minW, maxW, minN, maxN);
	scaGeom();
	scaBkGrnd();
//  [1/26/2006]
	calcDens();
//	generator.writeDt2("CalcDens.dt2");

	bndPntInst();
//	generator.writeDt2("bndInsert.dt2");
	recoverBnds();
//	generator.writeDt2("bndRecover.dt2");
	clrOuterEles();
//	generator.writeDt2("outerEleClr.dt2");
	
	int i=0;
	while (innerPntInst()) ;

//	generator.writeDt2("innerPntInst.dt2");

	smooth();
//	generator.writeDt2("smooth.dt2");
	
	rmvEmpNods();
	rmvEmpEles();		

	updateBndParent();

	rescaGeom();//  [1/26/2006]

	output();

	return 0;
}

bool DTIso2D::getBndInfo(double pdBNX[],		/* x coord. of boundary nodes */			
	double pdBNY[],		/* y coord. of boundary nodes */			
	int nBN,	  		/* number of boundary nodes */
	int pnBeg[],		/* start index of boundary edges */
	int pnEnd[] 		/* end index of boundary edges */
)
{
	m_nBnds = nBN;
	INTEGER iNod, iBnd;
	if (m_nBnds==0)
	{
		return false;
	}
	if (m_nBnds > INIT_ALLOC_BND_NUM)
	{
		printf("Error: Not enough memory for input boundaries!\n");
		return false;
	}
	if (m_nBnds + INIT_NOD_NUM > INIT_ALLOC_NOD_NUM)
	{
		printf("Error: Not enough memory for input nodes!\n");
		return false;
	}
//	printf( "bdy_Nods_NUM:%d\n", m_nBnds );
	for (iNod = 0; iNod < m_nBnds; iNod++)
	{
//#ifdef _3D_DECOM_
//	#ifdef _DOUBLE
//			fscanf(fp, "%d%lf%lf%lf", &iTok, &((m_pNodes)[i + INIT_NOD_NUM].pt)[0],
//				&((m_pNodes)[i + INIT_NOD_NUM].pt)[1], &(m_pNodes[i + INIT_NOD_NUM].density));
//	#else
//			fscanf(fp, "%d%f%f", &iTok, &((m_pNodes)[i + INIT_NOD_NUM].pt)[0],
//				&((m_pNodes)[i + INIT_NOD_NUM].pt)[1], &(m_pNodes[i + INIT_NOD_NUM].density));
//	#endif
//#else

		m_pNodes[iNod + INIT_NOD_NUM].pt[0] = pdBNX[iNod];
		m_pNodes[iNod + INIT_NOD_NUM].pt[1] = pdBNY[iNod];
//#endif
	} 
	m_nNodes = m_nBnds + INIT_NOD_NUM;

	for (iBnd=0; iBnd<m_nBnds; iBnd++)
	{
		m_pBnds[iBnd].beg = pnBeg[iBnd]+1;
		m_pBnds[iBnd].end = pnEnd[iBnd]+1;
		m_pBnds[iBnd].curve = 0;
		m_pBnds[iBnd].loop = 1;
		m_pBnds[iBnd].ele = NULL_ELEM;
	}
//	for (iBnd = 0; iBnd < m_nBnds; iBnd++)
//	{
//		m_pBnds[iBnd].beg = 
//		fscanf(fp, "%d%d%d%d%d", &iTok, &(m_pBnds)[i].beg, &(m_pBnds)[i].end,
//				&(m_pBnds)[i].curve, &(m_pBnds)[i].loop);
//		m_pBnds[i].ele = NULL_ELEM;
//	}

	return true;
}

#ifdef _3D_DECOM_
bool DTIso2D::getBackgroundMesh(double pdBMNX[],		/* x coord. of background mesh nodes */		
	double pdBMNY[],		/* y coord. of background mesh nodes */
	double pdBMNZ[],		/* Z coord. of background mesh nodes */	
	double pdBMNSpc[],		/* space values of background mesh */
	int    nBMN,			/* number of background mesh nodes */			
	int    pnBMEFm[],		/* forming points of background mesh elements */
	int    pnBMENg[],		/* neighboring eles. of background mesh elements */	
	int    nBME				/* number of background mesh elements */
)
{
	INTEGER i;
	int j;
	m_nBkGrndNode = nBMN;
	m_nBkGrndElem = nBME;

	if (m_nBkGrndNode > 0)
		m_pBkGrndNode = (Node *)malloc(sizeof(Node)*m_nBkGrndNode);
	else
		m_pBkGrndNode = NULL;

	for (i=0; i<m_nBkGrndNode; i++)
	{
		m_pBkGrndNode[i].pt[0] = pdBMNX[i];
		m_pBkGrndNode[i].pt[1] = pdBMNY[i];
		m_pBkGrndNode[i].pt[2] = pdBMNZ[i];
		m_pBkGrndNode[i].density = pdBMNSpc[i]; 
	}
	if (m_nBkGrndElem > 0)
		m_pBkGrndElem = (Elem *)malloc(sizeof(Elem)*m_nBkGrndElem);
	else
		m_pBkGrndElem = NULL;

	for (i=0; i<m_nBkGrndElem; i++)
	{
		for (j=0; j<=3; j++)
		{
			m_pBkGrndElem[i].form[j] = pnBMEFm[i*(3+1)+j];
			m_pBkGrndElem[i].neig[j] = pnBMENg[i*(3+1)+j];
		}
//			area = ::Area2(m_pBkGrndNode[m_pBkGrndElem[i].form[0]-1].pt,
//				m_pBkGrndNode[m_pBkGrndElem[i].form[1]-1].pt,
//				m_pBkGrndNode[m_pBkGrndElem[i].form[2]-1].pt);
//			if (area < 0)
//			{
//				k = m_pBkGrndElem[i].form[0];
//				m_pBkGrndElem[i].form[0] = m_pBkGrndElem[i].form[DIM];
//				m_pBkGrndElem[i].form[DIM] = k;
//			} // skip the element information
	}
	return true;
}

bool DTIso2D::getSources(double pdCX[],		/* coord. x of center points of sources */
	double pdCY[],		/* coord. y of center points of sources */
	double pdCZ[],		/* coord. z of center points of sources */
	double pdOu[],		/* out radius of sources */
	double pdIn[],      /* inner radius of sources */
	double pdSp[],      /* space values of source */
	int    nPS,			/* number of point sources */
	int    nLS,			/* number of line sources */
	int    nTS			/* number of tri. sources */
)
{
	INTEGER i;
	int j;

	m_nPntSrcNum = nPS;
	m_nLneSrcNum = nLS;
	m_nTriSrcNum = nTS;

	if (m_nPntSrcNum > 0)
		m_pntSources = (PointSource *)malloc(sizeof(PointSource)*m_nPntSrcNum);
	for (i=0; i<m_nPntSrcNum; i++)
	{
		m_pntSources[i].pt[0] = pdCX[i];
		m_pntSources[i].pt[1] = pdCY[i];
		m_pntSources[i].pt[2] = pdCZ[i];
		m_pntSources[i].rOuterRad = pdOu[i];
		m_pntSources[i].rInnerRad = pdIn[i];
		m_pntSources[i].rIntensity = pdSp[i];
		assert(m_pntSources[i].rOuterRad-m_pntSources[i].rInnerRad > EPS_ZERO_SQ);
	}
	if (m_nLneSrcNum > 0)
		m_lineSources = (LineSource *)malloc(sizeof(LineSource)*m_nLneSrcNum);
	for (i=0; i<m_nLneSrcNum; i++)
	{
		for (j=0; j<2; j++)
		{
			m_lineSources[i].points[j].pt[0] = pdCX[m_nPntSrcNum+i*2+j];
			m_lineSources[i].points[j].pt[1] = pdCY[m_nPntSrcNum+i*2+j];
			m_lineSources[i].points[j].pt[2] = pdCZ[m_nPntSrcNum+i*2+j];
			m_lineSources[i].points[j].rOuterRad = pdOu[m_nPntSrcNum+i*2+j];
			m_lineSources[i].points[j].rInnerRad = pdIn[m_nPntSrcNum+i*2+j];
			m_lineSources[i].points[j].rIntensity = pdSp[m_nPntSrcNum+i*2+j];
			assert(m_lineSources[i].points[j].rOuterRad-m_lineSources[i].points[j].rInnerRad > EPS_ZERO_SQ);
		}
	}
	if (m_nTriSrcNum > 0)
		m_triSources = (TriangleSource *)malloc(sizeof(TriangleSource)*m_nTriSrcNum);
	for (i=0; i<m_nTriSrcNum; i++)
	{
		for (j=0; j<3; j++)
		{
			m_triSources[i].points[j].pt[0] = pdCX[m_nPntSrcNum+m_nLneSrcNum+i*3+j];
			m_triSources[i].points[j].pt[1] = pdCY[m_nPntSrcNum+m_nLneSrcNum+i*3+j];
			m_triSources[i].points[j].pt[2] = pdCZ[m_nPntSrcNum+m_nLneSrcNum+i*3+j];
			m_triSources[i].points[j].rOuterRad = pdOu[m_nPntSrcNum+m_nLneSrcNum+i*3+j];
			m_triSources[i].points[j].rInnerRad = pdIn[m_nPntSrcNum+m_nLneSrcNum+i*3+j];
			m_triSources[i].points[j].rIntensity = pdSp[m_nPntSrcNum+m_nLneSrcNum+i*3+j];
			assert(m_triSources[i].points[j].rOuterRad-m_triSources[i].points[j].rInnerRad > EPS_ZERO_SQ);
		}
	}

	return true;
}
#else
bool DTIso2D::getBackgroundMesh(double pdBMNX[],		/* x coord. of background mesh nodes */		
	double pdBMNY[],		/* y coord. of background mesh nodes */
	double pdBMNSpc[],		/* space values of background mesh */
	int    nBMN,			/* number of background mesh nodes */			
	int    pnBMEFm[],		/* forming points of background mesh elements */
	int    pnBMENg[],		/* neighboring eles. of background mesh elements */	
	int    nBME				/* number of background mesh elements */
)
{
	INTEGER i, k;
	int j;
	m_nBkGrndNode = nBMN;
	m_nBkGrndElem = nBME;
	
	if (m_nBkGrndNode > 0)
		m_pBkGrndNode = (Node *)malloc(sizeof(Node)*m_nBkGrndNode);
	else
		m_pBkGrndNode = NULL;

	for (i=0; i<m_nBkGrndNode; i++)
	{
		m_pBkGrndNode[i].pt[0] = pdBMNX[i];
		m_pBkGrndNode[i].pt[1] = pdBMNY[i];
	}
	if (m_nBkGrndElem > 0)
		m_pBkGrndElem = (Elem *)malloc(sizeof(Elem)*m_nBkGrndElem);
	else
		m_pBkGrndElem = NULL;

	for (i=0; i<m_nBkGrndElem; i++)
	{
		for (j=0; j<=DIM; j++)
		{
			m_pBkGrndElem[i].form[j] = pnBMEFm[i*(DIM+1)+j];
			m_pBkGrndElem[i].neig[j] = pnBMENg[i*(DIM+1)+j];
		}
		double	area = Area2(m_pBkGrndNode[m_pBkGrndElem[i].form[0]-1].pt,
				m_pBkGrndNode[m_pBkGrndElem[i].form[1]-1].pt,
				m_pBkGrndNode[m_pBkGrndElem[i].form[2]-1].pt);
		if (area < 0)
		{
			k = m_pBkGrndElem[i].form[0];
			m_pBkGrndElem[i].form[0] = m_pBkGrndElem[i].form[DIM];
			m_pBkGrndElem[i].form[DIM] = k;
		} 
	}
	return true;
}

bool DTIso2D::getSources(double pdCX[],		/* coord. x of center points of sources */
	double pdCY[],		/* coord. y of center points of sources */
	double pdOu[],		/* out radius of sources */
	double pdIn[],      /* inner radius of sources */
	double pdSp[],      /* space values of source */
	int    nPS,			/* number of point sources */
	int    nLS			/* number of line sources */
)
{
	INTEGER i;
	int j;

	m_nPntSrcNum = nPS;
	m_nLneSrcNum = nLS;

	if (m_nPntSrcNum > 0)
		m_pntSources = (PointSource *)malloc(sizeof(PointSource)*m_nPntSrcNum);
	for (i=0; i<m_nPntSrcNum; i++)
	{
		m_pntSources[i].pt[0] = pdCX[i];
		m_pntSources[i].pt[1] = pdCY[i];
		m_pntSources[i].rOuterRad = pdOu[i];
		m_pntSources[i].rInnerRad = pdIn[i];
		m_pntSources[i].rIntensity = pdSp[i];
		assert(m_pntSources[i].rOuterRad-m_pntSources[i].rInnerRad > EPS_ZERO_SQ);
	}
	if (m_nLneSrcNum > 0)
		m_lineSources = (LineSource *)malloc(sizeof(LineSource)*m_nLneSrcNum);
	for (i=0; i<m_nLneSrcNum; i++)
	{
		for (j=0; j<2; j++)
		{
			m_lineSources[i].points[j].pt[0] = pdCX[m_nPntSrcNum+i*2+j];
			m_lineSources[i].points[j].pt[1] = pdCY[m_nPntSrcNum+i*2+j];
			m_lineSources[i].points[j].rOuterRad = pdOu[m_nPntSrcNum+i*2+j];
			m_lineSources[i].points[j].rInnerRad = pdIn[m_nPntSrcNum+i*2+j];
			m_lineSources[i].points[j].rIntensity = pdSp[m_nPntSrcNum+i*2+j];
			assert(m_lineSources[i].points[j].rOuterRad-m_lineSources[i].points[j].rInnerRad > EPS_ZERO_SQ);
		}
	}

	return true;
}
#endif

bool DTIso2D::outputMesh(double **ppMNX,      /* x coord. of mesh nodes */
	double **ppMNY,      /* y coord. of mesh nodes */
	double **ppMNSpc,	 /* space values of mesh nodes */
	int    *pnMN,		 /* number of mesh nodes */
	int    **ppnMEFm,    /* forming points of mesh elements */
	int    **ppnMENg,    /* neighboring eles. of mesh elements */
	int    *pnME,        /* number of mesh elements */
	int    **ppnPrt		 /* parents of boundary segments */
)
{
	INTEGER i;
	int j;

	*pnME = m_nElems;
	*pnMN = m_nNodes;
	*ppMNX = (double *)malloc(sizeof(double)*m_nNodes);
	*ppMNY = (double *)malloc(sizeof(double)*m_nNodes);
	*ppMNSpc = (double *)malloc(sizeof(double)*m_nNodes);

	for (i=0; i<m_nNodes; i++)
	{
		(*ppMNX)[i] = m_pNodes[i].pt[0];
		(*ppMNY)[i] = m_pNodes[i].pt[1];
		(*ppMNSpc)[i] = m_pNodes[i].density;
	}

	*ppnMEFm = (int *)malloc(sizeof(int)*m_nElems*(DIM+1));
	*ppnMENg = (int *)malloc(sizeof(int)*m_nElems*(DIM+1));

	for (i=0; i<m_nElems; i++)
	{
		for (j=0; j<=DIM; j++)
		{
			(*ppnMEFm)[i*(DIM+1)+j] = m_pElems[i].form[j];
			(*ppnMENg)[i*(DIM+1)+j] = m_pElems[i].neig[j];
		}
	}

	*ppnPrt = (int *)malloc(sizeof(int)*m_nBnds);

	for (i=0; i<m_nBnds; i++)
	{
		(*ppnPrt)[i] = m_pBnds[i].ele;
	}

	return true;
}

/*
* Check the validity of the data after ia MYPOINT insertion
* Parameters:
*	 void
* Result:
*   and error code
*/
int DTIso2D::checkGlobalData(bool bInvHole)
{
	INTEGER i;
	int m, l, k, n, iNeig, iNgNg;
	Elem *pElem = NULL, *pNeig = NULL;
	Node *pNode = NULL;
	int nErr = 0;
	
	/*Check if the neighboring relations are symmetry & correct*/
	for (i = 0; i < m_nElems; i++)
	{
		if (!isDelEle(i))
		{
			pElem = &(m_pElems[i]);
			for (m = 0; m <= DIM; m++)
			{
				iNeig = (pElem->neig)[m];
				if (iNeig != NULL_NEIG)
				{
					pNeig = &(m_pElems[iNeig]);
					if (isDelEle(iNeig))
					{//Error
						nErr = ERR_DEL_NEIG;
						return nErr;
					}
					else
					{
						//ensure there are ia neighbor for pNeig MYPOINTing to pElem
						for (l = 0; l <= DIM; l++)
						{
							iNgNg = (pNeig->neig)[l];
							if (iNgNg == i)
								break;
						}
						
						if (l > DIM)
						{//ERROR
							nErr = ERR_NEIG_NOT_SYM;
							return nErr;
						}
						else
						{
							//Check if ia common face exists between two neighbors
							for (n = 0; n <= DIM; n++)
							{
								if (n != l)
								{
									for (k = 0; k <= DIM; k++)
									{
										if (k != m)
										{
											if ((pElem->form)[k] == (pNeig->form)[n])
												break;
										}
									}
									
									if (k > DIM)
									{//ERROR
										nErr = ERR_FACE_NOT_COMM;
										return nErr;
									}
								}
							}//for (n = 0; n <= DIM; n++)
						}
						
					}
				}
			}
		}//if (isDelEle(i))
	}//for (i = 0; i < m_nElems; i++)
	
	return 0;/*checkNeigs(!bInvHole);*/
}

}
#pragma optimize("",on)