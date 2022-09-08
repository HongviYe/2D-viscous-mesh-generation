/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 三维各向同性Delaunay网格生成器 (版本号：0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2005年9月15日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * 联系方式
 *   电话：+86-571-87953165
 *   传真：+86-571-87953167
 *   邮箱：zdchenjj@yahoo.com.cn
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
#ifndef __iso2d_error_h__
#define __iso2d_error_h__

/* 成功 successful */
#define ERR_SUCCESS        0

/* 相邻单元被删除 a neighbor is deleted */
#define ERR_DEL_NEIG       1

/* 相邻关系不对称 neighbor relation is not symmetrical */
#define ERR_NEIG_NOT_SYM   2

/* 相邻单元共享面设置不正确 uncorrect shared face between two neighbors */
#define ERR_FACE_NOT_COMM  3

/* 三角化中存在无效空洞 invalid hole in the triangulation */
#define ERR_INVALID_HOLE   4

/* 最终结果中仍然有未删除的单元 undelete elements exist in the resulting mesh */
#define ERR_DEL_ELEM       5


/* 节点的母单元记录错误 invalid record for the parent element of a node */
#define ERR_INVALID_NOD_PRT 6

/* 某个面片的OuIn设置不正确 invalid out/inner settings for some elements */
#define ERR_INVALID_OU_IN  7

/* 非流型网格拓扑结构 non-manifold mesh */
#define ERR_NON_MANIFOLD_MESH 8

/* 
 * 打印错误信息
 * print error info.
 */
int printError(int nErr);

#endif /*  __iso2d_error_h__ */