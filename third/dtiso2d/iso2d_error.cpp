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
#include "iso2d_error.h"
#include <stdio.h>

/* 
 * 打印错误信息
 * print error info.
 */
int printError(int nErr)
{
	int nRet = 1;
	switch (nErr)
	{
	case ERR_SUCCESS:
		printf("Successful!\n");
		break;
	case ERR_DEL_NEIG:
		printf("A neighbor is deleted!\n");
		break;
	case ERR_NEIG_NOT_SYM:
		printf("Neighboring relations are not symmetry!\n");
		break;
	case ERR_FACE_NOT_COMM:
		printf("No common face exists between two neighboring elements!\n");
		break;
	case ERR_INVALID_HOLE:
		printf("An invalid holes exist in the triangulation!\n");
		break;
	case ERR_DEL_ELEM:
		printf("undelete elements exist in the resulting mesh\n");
		break;
	case ERR_INVALID_NOD_PRT:
		printf("Invalid record for the parent element of a node\n");
		break;
	case ERR_INVALID_OU_IN:
		printf("Invalid out/inner settings for some elements\n");
		break;
	case ERR_NON_MANIFOLD_MESH:
		printf("Non-manifold mesh\n");
		break;
	default:
		nRet = 0;
		break;
	}
	return nRet;
}