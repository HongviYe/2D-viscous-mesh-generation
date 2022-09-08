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
#ifndef __iso2d_error_h__
#define __iso2d_error_h__

/* �ɹ� successful */
#define ERR_SUCCESS        0

/* ���ڵ�Ԫ��ɾ�� a neighbor is deleted */
#define ERR_DEL_NEIG       1

/* ���ڹ�ϵ���Գ� neighbor relation is not symmetrical */
#define ERR_NEIG_NOT_SYM   2

/* ���ڵ�Ԫ���������ò���ȷ uncorrect shared face between two neighbors */
#define ERR_FACE_NOT_COMM  3

/* ���ǻ��д�����Ч�ն� invalid hole in the triangulation */
#define ERR_INVALID_HOLE   4

/* ���ս������Ȼ��δɾ���ĵ�Ԫ undelete elements exist in the resulting mesh */
#define ERR_DEL_ELEM       5


/* �ڵ��ĸ��Ԫ��¼���� invalid record for the parent element of a node */
#define ERR_INVALID_NOD_PRT 6

/* ĳ����Ƭ��OuIn���ò���ȷ invalid out/inner settings for some elements */
#define ERR_INVALID_OU_IN  7

/* �������������˽ṹ non-manifold mesh */
#define ERR_NON_MANIFOLD_MESH 8

/* 
 * ��ӡ������Ϣ
 * print error info.
 */
int printError(int nErr);

#endif /*  __iso2d_error_h__ */