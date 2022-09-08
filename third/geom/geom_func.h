/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 面向网格生成的几何函数库 (版本号：0.1)
 * Geometry Function Library for Mesh Generation (Version 0.1)
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2012年11月28日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 28/11/2012
 * 
 * 联系方式
 *   电话：+86-571-87951883
 *   传真：+86-571-87953167
 *   邮箱：chenjj@zju.edu.cn
 * For further information, please conctact
 *  Tel: +86-571-87951883
 *  Fax: +86-571-87953167
 * Mail: chenjj@zju.edu.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
#ifndef __geom_func_h__
#define __geom_func_h__

#include <stdio.h>
//#include<math.h>
//#include<stdlib.h>

namespace GEOM_FUNC
{

// extern void exactinit();

int cmpTwoTriDir(double *,double *,double *,double *,double *,double *);
/* 线面相交判断 */
int lnFacInt(double ln1[], double ln2[], 
			 double fac1[], double fac2[], double fac3[],  
			 double pnt[], int *val,
			 int edgFac1 = -1, int edgFac2 = -1, int edgFac3 = -1, 
			 double *pt1 = NULL, double *pt2 = NULL, double *pt3 = NULL);

double exactinit_threadSafe();
#ifdef __cplusplus
extern "C" {
#endif
/**
 * @brief o3dstaticfilter, ispstaticfilter
 * Static filters for orient3d() and insphere().
 * @author Wang, Junji
 * @date May, 31st, 2018
 */
extern double o3dstaticfilter;
extern double ispstaticfilter;
/* function prototype */

/* extern the predicates */
extern double macheps;
double exactinit();
double fixedSplitPoint(double s1, double s2, double pnt1, double pnt2);
double orient2d(double *pa, double *pb, double *pc);
double orient3d(double *pa, double *pb, double *pc, double *pd);
double orient3dexact(double *pa, double *pb, double *pc, double *pd);
double incircle(double*, double*, double*, double*);
double insphere(double *pa, double *pb, double *pc, double *pd, double *pe);
double insphereexact(double *pa, double *pb, double *pc, double *pd, double *pe);

void norm_3p(double *p1 , double *p2 ,double *p3, double *normal);
extern int tri_tri_overlap_test_3d(double p1[3], double q1[3], double r1[3], 
			    double p2[3], double q2[3], double r2[3]);
extern int one_node_same_tri_tri_overlap_3d(double p1[3], double q1[3], double r1[3], 								 
			    double p2[3], double q2[3], double r2[3]);  //zhaodawei add 2010-07-16

/* -----------------------------------------------------------------------------------
 * 这是一个判定1条线段和三角形是否相交的代码，并返回交点位置 
 * -----------------------------------------------------------------------------------*/
enum LIN_TRI_INT_TYPE 
{
	LTI_INTERSECT_NUL = 0,
	LTI_INTERSECT_NOD,
	LTI_INTERSECT_EDG,
	LTI_INTERSECT_FAC,
	LTI_INTERSECT_INS,			/* 直线完全在平面里 */
	LTI_INTERSECT_DEG_FACE		/* 这是临时设置的枚举变量，今后可能会被替换。
								   表示输入的面片为退化面片，目前暂不判定一个
								   退化的面片是否和线相交（没有这个需求）*/
};

/* -----------------------------------------------------------------------------------
 * 辅助函数，假设p1/p2/p3/p4共面，将p1/p2/p3/p4向使得p1p2p3投影面积最大的坐标平面投影
 * 如果p1p2p3是退化面片，返回0，否则，返回1
 * -----------------------------------------------------------------------------------*/
int proj_four_coplanr_points(double p1[3], double p2[3], double p3[3], double p4[3],
		double proj1[2], double proj2[2], double proj3[2], double proj4[2]);

/* -----------------------------------------------------------------------------------
 * 这是一个判定1条线段和三角形是否相交的代码(2D)，并返回交点位置
 * 相交类型：
 * PNT 相交于1个点（intCod返回点的编号0~2，i代表面的第i个顶点）：
 * EDG 相交于1条边（intCod返回边的编号0~2，i代表(i,(i+1)%3形成的边)
 * FAC 相交于1个面
 * intPnt返回交点的值
 * 注意：在共面情形下，一条线可能会和一个面的两条边都相交，此时，
 * 根据调用该算法的边界恢复过程需求，我们只返回离得最近的那个交点
 * 当将这个代码用于表面可能相交的情形时，我们需要根据具体情况进行更新
 * -----------------------------------------------------------------------------------*/
int lin_tri_intersect2d(double linep[2][2], double facep[3][2], int *intTyp, int *intCod, double intPnt[2]);

/* -----------------------------------------------------------------------------------
 * 这是一个判定1个位于面上的3维点和面的关系
 * 相交类型：
 * PNT 相交于1个点（intCod返回点的编号0~2，i代表面的第i个顶点）：
 * EDG 相交于1条边（intCod返回边的编号0~2，i代表(i,(i+1)%3形成的边)
 * FAC 相交于1个面
 * intPnt返回交点的值
 * -----------------------------------------------------------------------------------*/
extern int pnt_tri_intersect3d(double poinp[3], double facep[3][3], int *intTyp, int *intCod, double intPnt[3]);

extern int lin_tri_intersect3d(double linep[2][3], double facep[3][3], int *intTyp, int *intCod, double intPnt[3], bool bEpsilon = true);

extern int lin_tri_intersect3d_idx(int iline[2], int iface[3], double linep[2][3], double facep[3][3], int *intTyp, int *intCod, double intPnt[3], bool bEpsilon = true);

extern int isintersect_oneSharePoint(int iface[], int iline[], double facep[][3], double linep[][3], int *intTyp, int *intCod, double intPnt[3]);

/* ------------------------------------------------------------------------------------
 * 线面相交判断
 * 这是一个过渡函数，我们尝试用lin_tri_intersect3d去实现它，并逐步替代它
 * ---------------------------------------------------------------------------------- */
int lnFacInt2(double ln1[], double ln2[], 
			 double fac1[], double fac2[], double fac3[],  
			 double pnt[], int *val,
			 int edgFac1 = -1, int edgFac2 = -1, int edgFac3 = -1, 
			 double *pt1 = NULL, double *pt2 = NULL, double *pt3 = NULL, bool bEpsilon = true);

/* -----------------------------------------------------------------------------------
 * 这是一个判定2个三角面片是否相交的代码，它将在small polyhedron reconnection算法中调用
 * -----------------------------------------------------------------------------------*/
extern int tri_tri_intersect3d(int facei1[3], int facei2[3], double facep1[3][3], double facep2[3][3]);
extern int tri_tri_intersect3d_fast(int facei1[3], int facei2[3], double *facep1[3], double *facep2[3]);


extern int lin_tri_intersect3d_check(double linep[2][3], double facep[3][3]);

/* 计算一个四面体单元的形状质量值 */
extern double tetrahedron_gamma(double p1[3], double p2[3], double p3[3], double p4[3]);

#ifdef __cplusplus
}
#endif

}
#endif