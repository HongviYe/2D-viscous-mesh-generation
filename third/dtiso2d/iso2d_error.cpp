/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ά����ͬ��Delaunay���������� (�汾�ţ�0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2005��9��15��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * ��ϵ��ʽ
 *   �绰��+86-571-87953165
 *   ���棺+86-571-87953167
 *   ���䣺zdchenjj@yahoo.com.cn
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
 * ��ӡ������Ϣ
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