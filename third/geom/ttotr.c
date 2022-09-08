/*
*		  		
*  Triangle-Triangle Overlap Test Routines				
*  July, 2002                                                          
*  Updated December 2003                                                
*                                                                       
*  This file contains C implementation of algorithms for                
*  performing two and three-dimensional triangle-triangle intersection test 
*  The algorithms and underlying theory are described in                    
*                                                                           
* "Fast and Robust Triangle-Triangle Overlap Test 
*  Using Orientation Predicates"  P. Guigue - O. Devillers
*                                                 
*  Journal of Graphics Tools, 8(1), 2003                                    
*                                                                           
*  Several geometric predicates are defined.  Their parameters are all      
*  points.  Each point is an array of two or three double precision         
*  floating point numbers. The geometric predicates implemented in          
*  this file are:                                                            
*                                                                           
*    int tri_tri_overlap_test_3d(p1,q1,r1,p2,q2,r2)                         
*    int tri_tri_overlap_test_2d(p1,q1,r1,p2,q2,r2)                         
*                                                                           
*    int tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2,
*                                     coplanar,source,target)               
*                                                                           
*       is a version that computes the segment of intersection when            
*       the triangles overlap (and are not coplanar)                        
*                                                                           
*    each function returns 1 if the triangles (including their              
*    boundary) intersect, otherwise 0                                       
*                                                                           
*                                                                           
*  Other information are available from the Web page                        
*  http://www.acm.org/jgt/papers/GuigueDevillers03/                         
*                                                                           
*/
#include <math.h>

/* function prototype */

int tri_tri_overlap_test_3d(double p1[3], double q1[3], double r1[3], 
			    double p2[3], double q2[3], double r2[3]);


int coplanar_tri_tri3d(double  p1[3], double  q1[3], double  r1[3],
		       double  p2[3], double  q2[3], double  r2[3],
		       double  N1[3], double  N2[3]);


int tri_tri_overlap_test_2d(double p1[2], double q1[2], double r1[2], 
			    double p2[2], double q2[2], double r2[2]);


int tri_tri_intersection_test_3d(double p1[3], double q1[3], double r1[3], 
				 double p2[3], double q2[3], double r2[3],
				 int * coplanar, 
				 double source[3],double target[3]);

//*****************************************************zhaodawei add 2010-07-16

int one_node_same_tri_tri_overlap_3d(double p1[3], double q1[3], double r1[3], 
									 
									 double p2[3], double q2[3], double r2[3]);  //zhaodawei

int is_one_node_same_2tri_inter( double p1[3], double q1[3], double r1[3], 
								
								double p2[3], double q2[3], double r2[3] ); //zzhaodawei

int is_tri_tri_inter_or_in( double s1[2], double s2[2], double s3[2],
						   double t1[2], double t2[2], double t3[2] ); //zhaodawei

int is1pInTri3( double p[2], double t1[2], double t2[2], double t3[2] ); //zhaodawei

int getMainIndexOfPlaneFrom3p(double p0[3], double p1[3], double p2[3]); //zhaodawei

void prjPtTo2D(double p[3],double pj[2],int index);  //zhaodawei

int is_seg_tri_inter_or_in( double p[2], double t1[2], double t2[2], double t3[2] ); //zhaodawei

double one_node_tri_orient3D( double p[3], double p2[3], double q2[3], double r2[3] ); //zhaodawei

//********************************************************************

/* coplanar returns whether the triangles are coplanar  
*  source and target are the endpoints of the segment of 
*  intersection if it exists) 
*/

//******************************************zhaodawei add 2010-07-16

void vec_2p(double *a , double *b , double *ab)
{
	int i;
	for(i=0; i<3; i++)
		ab[i]=b[i]-a[i];
}
void vec_crop(double *vector1, double *vector2, double *vector3)
{
	vector3[0] = vector1[1]*vector2[2]-vector2[1]*vector1[2];
	vector3[1] = vector2[0]*vector1[2]-vector1[0]*vector2[2];
	vector3[2] = vector1[0]*vector2[1]-vector2[0]*vector1[1];
}

double vec_val(double *vec )
{
	return( sqrt(vec[0] * vec[0] +vec[1] * vec[1] + vec[2] *vec[2] )) ;
}

void norm_3p(double *p1 , double *p2 ,double *p3, double *normal )
{
	double v12[3] , v13[3] ;
	vec_2p(p1, p2, v12 ) ;
	vec_2p(p1, p3, v13 ) ; // the direction of crop has changed ?!
	vec_crop(v12, v13, normal ) ;
}

//****************************************************************



/* some 3D macros */

#define CROSS(dest,v1,v2)                       \
               dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
               dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
               dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
 


#define SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; \
                        dest[1]=v1[1]-v2[1]; \
                        dest[2]=v1[2]-v2[2]; 


#define SCALAR(dest,alpha,v) dest[0] = alpha * v[0]; \
                             dest[1] = alpha * v[1]; \
                             dest[2] = alpha * v[2];



#define CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) {\
  SUB(v1,p2,q1)\
  SUB(v2,p1,q1)\
  CROSS(N1,v1,v2)\
  SUB(v1,q2,q1)\
  if (DOT(v1,N1) > 0.0f) return 0;\
  SUB(v1,p2,p1)\
  SUB(v2,r1,p1)\
  CROSS(N1,v1,v2)\
  SUB(v1,r2,p1) \
  if (DOT(v1,N1) > 0.0f) return 0;\
  else return 1; }

/* Permutation in a canonical form of T2's vertices */
#define TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
  if (dp2 > 0.0f) { \
     if (dq2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2) \
     else if (dr2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
     else CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) }\
  else if (dp2 < 0.0f) { \
    if (dq2 < 0.0f) CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0.0f) CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
    else CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
  } else { \
    if (dq2 < 0.0f) { \
      if (dr2 >= 0.0f)  CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
      else CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2)\
    } \
    else if (dq2 > 0.0f) { \
      if (dr2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
      else  CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
    } \
    else  { \
      if (dr2 > 0.0f) CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
      else if (dr2 < 0.0f) CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2)\
      else return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);\
     }}}
  
/*
*
*  Three-dimensional Triangle-Triangle Overlap Test
*
*/
int tri_tri_overlap_test_3d(double p1[3], double q1[3], double r1[3], 

			    double p2[3], double q2[3], double r2[3])
{
  double dp1, dq1, dr1, dp2, dq2, dr2;
  double v1[3], v2[3];
  double N1[3], N2[3]; 
  
  /* Compute distance signs  of p1, q1 and r1 to the plane of
     triangle(p2,q2,r2) */
  SUB(v1,p2,r2)
  SUB(v2,q2,r2)
  CROSS(N2,v1,v2)

  SUB(v1,p1,r2)
  dp1 = DOT(v1,N2);
  SUB(v1,q1,r2)
  dq1 = DOT(v1,N2);
  SUB(v1,r1,r2)
  dr1 = DOT(v1,N2);
  
  if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return 0; 

  /* Compute distance signs  of p2, q2 and r2 to the plane of
     triangle(p1,q1,r1) */
  SUB(v1,q1,p1)
  SUB(v2,r1,p1)
  CROSS(N1,v1,v2)

  SUB(v1,p2,r1)
  dp2 = DOT(v1,N1);
  SUB(v1,q2,r1)
  dq2 = DOT(v1,N1);
  SUB(v1,r2,r1)
  dr2 = DOT(v1,N1);
  
  if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return 0;

  /* Permutation in a canonical form of T1's vertices */
  if (dp1 > 0.0f) 
  {
    if (dq1 > 0.0f) TRI_TRI_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
    else if (dr1 > 0.0f) TRI_TRI_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)	
    else TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
  } 
  else if (dp1 < 0.0f) 
  {
    if (dq1 < 0.0f) TRI_TRI_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
    else if (dr1 < 0.0f) TRI_TRI_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    else TRI_TRI_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
  } 
  else 
  {
    if (dq1 < 0.0f) 
	{
      if (dr1 >= 0.0f) TRI_TRI_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
      else TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
    }
    else if (dq1 > 0.0f) 
	{
      if (dr1 > 0.0f) TRI_TRI_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
      else TRI_TRI_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    }
    else  
	{
      if (dr1 > 0.0f) TRI_TRI_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
      else if (dr1 < 0.0f) TRI_TRI_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
      else return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);
    }
  }
  return 0;
};



int coplanar_tri_tri3d(double p1[3], double q1[3], double r1[3],
		       double p2[3], double q2[3], double r2[3],
		       double normal_1[3], double normal_2[3]){
  
  double P1[2],Q1[2],R1[2];
  double P2[2],Q2[2],R2[2];

  double n_x, n_y, n_z;

  n_x = ((normal_1[0]<0)?-normal_1[0]:normal_1[0]);
  n_y = ((normal_1[1]<0)?-normal_1[1]:normal_1[1]);
  n_z = ((normal_1[2]<0)?-normal_1[2]:normal_1[2]);


  /* Projection of the triangles in 3D onto 2D such that the area of
     the projection is maximized. */


  if (( n_x > n_z ) && ( n_x >= n_y )) {
    // Project onto plane YZ

      P1[0] = q1[2]; P1[1] = q1[1];
      Q1[0] = p1[2]; Q1[1] = p1[1];
      R1[0] = r1[2]; R1[1] = r1[1]; 
    
      P2[0] = q2[2]; P2[1] = q2[1];
      Q2[0] = p2[2]; Q2[1] = p2[1];
      R2[0] = r2[2]; R2[1] = r2[1]; 

  } else if (( n_y > n_z ) && ( n_y >= n_x )) {
    // Project onto plane XZ

    P1[0] = q1[0]; P1[1] = q1[2];
    Q1[0] = p1[0]; Q1[1] = p1[2];
    R1[0] = r1[0]; R1[1] = r1[2]; 
 
    P2[0] = q2[0]; P2[1] = q2[2];
    Q2[0] = p2[0]; Q2[1] = p2[2];
    R2[0] = r2[0]; R2[1] = r2[2]; 
    
  } else {
    // Project onto plane XY

    P1[0] = p1[0]; P1[1] = p1[1]; 
    Q1[0] = q1[0]; Q1[1] = q1[1]; 
    R1[0] = r1[0]; R1[1] = r1[1]; 
    
    P2[0] = p2[0]; P2[1] = p2[1]; 
    Q2[0] = q2[0]; Q2[1] = q2[1]; 
    R2[0] = r2[0]; R2[1] = r2[1]; 
  }

  return tri_tri_overlap_test_2d(P1,Q1,R1,P2,Q2,R2);
    
};



/*
*                                                                
*  Three-dimensional Triangle-Triangle Intersection              
*
*/

/*
   This macro is called when the triangles surely intersect
   It constructs the segment of intersection of the two triangles
   if they are not coplanar.
*/

#define CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) { \
  SUB(v1,q1,p1) \
  SUB(v2,r2,p1) \
  CROSS(N,v1,v2) \
  SUB(v,p2,p1) \
  if (DOT(v,N) > 0.0f) {\
    SUB(v1,r1,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) <= 0.0f) { \
      SUB(v2,q2,p1) \
      CROSS(N,v1,v2) \
      if (DOT(v,N) > 0.0f) { \
	SUB(v1,p1,p2) \
	SUB(v2,p1,r1) \
	alpha = DOT(v1,N2) / DOT(v2,N2); \
	SCALAR(v1,alpha,v2) \
	SUB(source,p1,v1) \
	SUB(v1,p2,p1) \
	SUB(v2,p2,r2) \
	alpha = DOT(v1,N1) / DOT(v2,N1); \
	SCALAR(v1,alpha,v2) \
	SUB(target,p2,v1) \
	return 1; \
      } else { \
	SUB(v1,p2,p1) \
	SUB(v2,p2,q2) \
	alpha = DOT(v1,N1) / DOT(v2,N1); \
	SCALAR(v1,alpha,v2) \
	SUB(source,p2,v1) \
	SUB(v1,p2,p1) \
	SUB(v2,p2,r2) \
	alpha = DOT(v1,N1) / DOT(v2,N1); \
	SCALAR(v1,alpha,v2) \
	SUB(target,p2,v1) \
	return 1; \
      } \
    } else { \
      return 0; \
    } \
  } else { \
    SUB(v2,q2,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) < 0.0f) { \
      return 0; \
    } else { \
      SUB(v1,r1,p1) \
      CROSS(N,v1,v2) \
      if (DOT(v,N) >= 0.0f) { \
	SUB(v1,p1,p2) \
	SUB(v2,p1,r1) \
	alpha = DOT(v1,N2) / DOT(v2,N2); \
	SCALAR(v1,alpha,v2) \
	SUB(source,p1,v1) \
	SUB(v1,p1,p2) \
	SUB(v2,p1,q1) \
	alpha = DOT(v1,N2) / DOT(v2,N2); \
	SCALAR(v1,alpha,v2) \
	SUB(target,p1,v1) \
	return 1; \
      } else { \
	SUB(v1,p2,p1) \
	SUB(v2,p2,q2) \
	alpha = DOT(v1,N1) / DOT(v2,N1); \
	SCALAR(v1,alpha,v2) \
	SUB(source,p2,v1) \
	SUB(v1,p1,p2) \
	SUB(v2,p1,q1) \
	alpha = DOT(v1,N2) / DOT(v2,N2); \
	SCALAR(v1,alpha,v2) \
	SUB(target,p1,v1) \
	return 1; \
      }}}} 

								

#define TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
  if (dp2 > 0.0f) { \
     if (dq2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2) \
     else if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
     else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) }\
  else if (dp2 < 0.0f) { \
    if (dq2 < 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    else CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
  } else { \
    if (dq2 < 0.0f) { \
      if (dr2 >= 0.0f)  CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
      else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2)\
    } \
    else if (dq2 > 0.0f) { \
      if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
      else  CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    } \
    else  { \
      if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
      else if (dr2 < 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2)\
      else { \
       	*coplanar = 1; \
	return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);\
     } \
  }} }
  

/*
   The following version computes the segment of intersection of the
   two triangles if it exists. 
   coplanar returns whether the triangles are coplanar
   source and target are the endpoints of the line segment of intersection 
*/

int tri_tri_intersection_test_3d(double p1[3], double q1[3], double r1[3], 
				 double p2[3], double q2[3], double r2[3],
				 int * coplanar, 
				 double source[3], double target[3] )
				 
{
  double dp1, dq1, dr1, dp2, dq2, dr2;
  double v1[3], v2[3], v[3];
  double N1[3], N2[3], N[3];
  double alpha;

  // Compute distance signs  of p1, q1 and r1 
  // to the plane of triangle(p2,q2,r2)


  SUB(v1,p2,r2)
  SUB(v2,q2,r2)
  CROSS(N2,v1,v2)

  SUB(v1,p1,r2)
  dp1 = DOT(v1,N2);
  SUB(v1,q1,r2)
  dq1 = DOT(v1,N2);
  SUB(v1,r1,r2)
  dr1 = DOT(v1,N2);
  
  if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return 0; 

  // Compute distance signs  of p2, q2 and r2 
  // to the plane of triangle(p1,q1,r1)

  
  SUB(v1,q1,p1)
  SUB(v2,r1,p1)
  CROSS(N1,v1,v2)

  SUB(v1,p2,r1)
  dp2 = DOT(v1,N1);
  SUB(v1,q2,r1)
  dq2 = DOT(v1,N1);
  SUB(v1,r2,r1)
  dr2 = DOT(v1,N1);
  
  if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return 0;

  // Permutation in a canonical form of T1's vertices


  if (dp1 > 0.0f) {
    if (dq1 > 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
    else if (dr1 > 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
	
    else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
  } else if (dp1 < 0.0f) {
    if (dq1 < 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
    else if (dr1 < 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    else TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
  } else {
    if (dq1 < 0.0f) {
      if (dr1 >= 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
      else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
    }
    else if (dq1 > 0.0f) {
      if (dr1 > 0.0f) TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
      else TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    }
    else  {
      if (dr1 > 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
      else if (dr1 < 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
      else {
	// triangles are co-planar

	*coplanar = 1;
	return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);
      }
    }
  }
};





/*
*
*  Two dimensional Triangle-Triangle Overlap Test    
*
*/


/* some 2D macros */

#define ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))


#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) {\
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f)\
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0f)\
      if (ORIENT_2D(P1,P2,Q1) > 0.0f) {\
	if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) return 1; \
	else return 0;} else {\
	if (ORIENT_2D(P1,P2,R1) >= 0.0f)\
	  if (ORIENT_2D(Q1,R1,P2) >= 0.0f) return 1; \
	  else return 0;\
	else return 0;}\
    else \
      if (ORIENT_2D(P1,Q2,Q1) <= 0.0f)\
	if (ORIENT_2D(R2,Q2,R1) <= 0.0f)\
	  if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) return 1; \
	  else return 0;\
	else return 0;\
      else return 0;\
  else\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) \
      if (ORIENT_2D(Q1,R1,R2) >= 0.0f)\
	if (ORIENT_2D(P1,P2,R1) >= 0.0f) return 1;\
	else return 0;\
      else \
	if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) {\
	  if (ORIENT_2D(R2,R1,Q2) >= 0.0f) return 1; \
	  else return 0; }\
	else return 0; \
    else  return 0; \
 };



#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
        if (ORIENT_2D(P1,Q1,R2) >= 0.0f) return 1; \
        else return 0;} else { \
      if (ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
	if (ORIENT_2D(R1,P1,P2) >= 0.0f) return 1; else return 0;} \
      else return 0; } \
  } else {\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) {\
      if (ORIENT_2D(P1,P2,R1) >= 0.0f) {\
	if (ORIENT_2D(P1,R1,R2) >= 0.0f) return 1;  \
	else {\
	  if (ORIENT_2D(Q1,R1,R2) >= 0.0f) return 1; else return 0;}}\
      else  return 0; }\
    else return 0; }}



int ccw_tri_tri_intersection_2d(double p1[2], double q1[2], double r1[2], 
				double p2[2], double q2[2], double r2[2]) {
  if ( ORIENT_2D(p2,q2,p1) >= 0.0f ) {
    if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) return 1;
      else INTERSECTION_TEST_EDGE(p1,q1,r1,p2,q2,r2)
    } else {  
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
	INTERSECTION_TEST_EDGE(p1,q1,r1,r2,p2,q2)
      else INTERSECTION_TEST_VERTEX(p1,q1,r1,p2,q2,r2)}}
  else {
    if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
	INTERSECTION_TEST_EDGE(p1,q1,r1,q2,r2,p2)
      else  INTERSECTION_TEST_VERTEX(p1,q1,r1,q2,r2,p2)}
    else INTERSECTION_TEST_VERTEX(p1,q1,r1,r2,p2,q2)}
};


int tri_tri_overlap_test_2d(double p1[2], double q1[2], double r1[2], 
			    double p2[2], double q2[2], double r2[2]) {
  if ( ORIENT_2D(p1,q1,r1) < 0.0f )
    if ( ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,q2,r2);
  else
    if ( ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,q2,r2);

};


double one_node_tri_orient3D( double p[3], double p2[3], double q2[3], double r2[3] )
{
	double dp1;
	double v1[3], v2[3];
	double N1[3], N2[3];
	SUB(v1,p2,r2)
	SUB(v2,q2,r2)
	CROSS(N2,v1,v2)
		
	SUB(v1,p,r2)
    dp1 = DOT(v1,N2);  //在这里等于0，注意浮点误差
	return dp1;
}

/*
* 当两个三角形中，有一节点一样的情形下，两三角形是否相交
* 相交返回1，不相交返回0
* p1和p2代表相同的那个点
*/
int one_node_same_tri_tri_overlap_3d(double p1[3], double q1[3], double r1[3], 
									 
			    double p2[3], double q2[3], double r2[3])
{
  double dp1, dq1, dr1, dp2, dq2, dr2, dt;
  double v1[3], v2[3];
  double N1[3], N2[3]; 
  double d11,d12,d21,d22;
  int n_zero1 = 1, n_zero2 = 1;

  double d31,d32,d41,d42;
  
  /* Compute distance signs  of p1, q1 and r1 to the plane of
     triangle(p2,q2,r2) */


  SUB(v1,p2,r2)
  SUB(v2,q2,r2)
  CROSS(N2,v1,v2)

  SUB(v1,p1,r2)
  dp1 = DOT(v1,N2);  //在这里等于0，注意浮点误差
  SUB(v1,q1,r2)
  dq1 = DOT(v1,N2);
  SUB(v1,r1,r2)
  dr1 = DOT(v1,N2);
  d11 = dq1;
  d12 = dr1;
  if ( ( d11 * d12 ) > 0.0f )  return 0;  //不相交

  
  SUB(v1,q1,p1)
  SUB(v2,r1,p1)
  CROSS(N1,v1,v2)

  SUB(v1,p2,r1)
  dp2 = DOT(v1,N1);  //在这里等于0，注意浮点误差
  SUB(v1,q2,r1)
  dq2 = DOT(v1,N1);
  SUB(v1,r2,r1)
  dr2 = DOT(v1,N1);
	
  d21 = dq2;
  d22 = dr2;

  if ( ( d21 * d22 ) > 0.0f ) return 0; //不相交

  //三角形各有两个点在两三角形的交线上
  if( dq1 == 0.0f ) n_zero1++;
  if( dr1 == 0.0f ) n_zero1++;
  if( dq2 == 0.0f ) n_zero2++;
  if( dr2 == 0.0f ) n_zero2++;
  if( n_zero1 == 2 && n_zero2 == 2 ) return 0;

  

  if( ( ( d11 * d12 ) == 0.0f ) && ( ( d21 * d22 ) < 0.0f ) )
  {
	  if( d11 == 0.0f )
	  {
		  if( is1Seg1TriInter( p2, q2, r2, q1 ) ) return 1;
		  else return 0;
	  }
	  if( d12 == 0.0f )
	  {
		  if( is1Seg1TriInter( p2, q2, r2, r1 ) ) return 1;
		  else return 0;
	  }
  }
  else if( ( ( d11 * d12 ) < 0.0f ) && ( ( d21 * d22 ) == 0.0f ) )
  {
	  if( d21 == 0.0f )
	  {
		  if( is1Seg1TriInter( p1, q1, r1, q2 ) ) return 1;
		  else return 0;
	  }
	  if( d22 == 0.0f )
	  {
		  if( is1Seg1TriInter( p1, q1, r1, r2 ) ) return 1;
		  else return 0;
	  }
  }

  if( ( ( d11 * d12 ) < 0.0f ) && ( ( d21 * d22 ) < 0.0f ) )
  {
	  if( d11 > 0.0f )
	  {
		  dt = one_node_tri_orient3D( q1, r1, q2, p2 );
		  if( d11 * dt > 0 ) return 0;
		  else return 1;
	  }
	  else if( d12 > 0.0f )
	  {
		  dt = one_node_tri_orient3D( r1, q1, q2, p2 );
		  if( d12 * dt > 0 ) return 0;
		  else return 1;
	  }
  }


// 下面是共面的情况
	
	if( is_one_node_same_2tri_inter( p1,q1,r1,p2,q2,r2 ) ) return 1;
	else return 0;
}

//p and p is the same node.
//传进来参数的时候，必须保证两个三角形的第一个点是相同的那个点
int is_one_node_same_2tri_inter( double p1[3], double q1[3], double r1[3], 
								 
			    double p2[3], double q2[3], double r2[3] )
{
//判断三角形1是否与三角形2相交或三角形1是否在三角形2内部
	int flag1 = 1, flag2 = 1;
	double p12[2], q12[2], r12[2], p22[2], q22[2], r22[2];
	double tmp1[2], tmp2[2], tmp3[2];

	int index = getMainIndexOfPlaneFrom3p(p2, q2, r2);
	prjPtTo2D(p1, tmp1, index);
	prjPtTo2D(q1, tmp2, index);
	prjPtTo2D(p2, p22, index);
	prjPtTo2D(q2, q22, index);
	prjPtTo2D(r2, r22, index);
	prjPtTo2D(r1, tmp3, index);

	flag1 = is_tri_tri_inter_or_in( tmp1, tmp2, tmp3, p22, q22, r22 );
	
	index = getMainIndexOfPlaneFrom3p(p1, q1, r1);
	prjPtTo2D(p2, tmp1, index);
	prjPtTo2D(q2, tmp2, index);
	prjPtTo2D(p1, p12, index);
	prjPtTo2D(q1, q12, index);
	prjPtTo2D(r1, r12, index);
	prjPtTo2D(r2, tmp3, index);

	flag2 = is_tri_tri_inter_or_in( tmp1, tmp2, tmp3, p12, q12, r12 );

	return ( flag1 || flag2 );
}

//s1 is same to t1.
//判断两个共面的三角形s1s2s3是否与面t1t2t3相交或在面t1t2t3里面
//返回1表示函数成立，返回0表示函数不成立
int is_tri_tri_inter_or_in( double s1[2], double s2[2], double s3[2],
							double t1[2], double t2[2], double t3[2] )
{
	int i_s2, i_s3, ii_s2, ii_s3;
	i_s2 = is1pInTri3( s2, t1, t2, t3 );
	i_s3 = is1pInTri3( s3, t1, t2, t3 );
	if( i_s2 == 2 && i_s3 == 2 ) return 1;
	if( i_s2 == 1 || i_s3 == 1 ) return 1;

	ii_s2 = ( s2, t1, t2, t3 );
	ii_s3 = ( s3, t1, t2, t3 );

	if( ii_s2 || ii_s3 ) return 1;
	else return 0;
}

//判断线段是否与三角形相交
int is_seg_tri_inter_or_in( double p[2], double t1[2], double t2[2], double t3[2] )
{
	double t1t2t3, t2t3p, t1t2p, t1t3p;
	t1t2t3 = orient2d( t1, t2, t3 );
	
	t2t3p = orient2d( t2, t3, p );
	t1t2p = orient2d( t1, t2, p );
	t1t3p = orient2d( t1, t3, p );

	return ( ( t1t2t3 * t2t3p <= 0 ) && ( t1t2p * t1t3p < 0 ) );

}

//test p is or isn't in tri p1 q1 r1
//判断点p是否在三角形t1t2t3里面还是在t1t2t3的边的延长线上
//在延长线上返回2，在里面返回1，否则返回0
int is1pInTri3( double p[2], double t1[2], double t2[2], double t3[2] )
{
	double t1t2p, t2t3p, t3t1p;
	t1t2p = orient2d( t1, t3, p );
	t2t3p = orient2d( t2, t3, p );
	t3t1p = orient2d( t3, t1, p );
	if( ( t1t2p < 0 && t2t3p < 0 && t3t1p < 0 )
		|| ( t1t2p > 0 && t2t3p > 0 && t3t1p > 0 ) )
		return 1;
	else if( t1t2p == 0 || t3t1p == 0 ) return 2;
	else return 0;
}

//判断一条线段与1个三角形共面的情况下是否相交
int is1Seg1TriInter( double t1[3], double t2[3], double t3[3], double p[3] )
{
	double tmp[2], tmp1[2], tmp2[2], tmp3[2];
	int i_p;
	int index = getMainIndexOfPlaneFrom3p(t1, t2, t3);
	prjPtTo2D(t1, tmp1, index);
	prjPtTo2D(t2, tmp2, index);
	prjPtTo2D(t3, tmp3, index);
	prjPtTo2D(p, tmp, index);

	i_p = is1pInTri3( tmp, tmp1, tmp2, tmp3 );
	if( i_p ) return 1;

	return is_seg_tri_inter_or_in( tmp, tmp1, tmp2, tmp3 );
}

int getMainIndexOfPlaneFrom3p(double p0[3], double p1[3], double p2[3])
{
	//貌似是求法向量三个分量，哪个最大，返回它的索引
	//？？？在这里到底是什么意思呢
	double norm[3];
	int index=-1;
	double max=0;
	int i = 0;
	norm_3p(p0, p1, p2, norm);
	
	for(i=0; i<3; i++)
		if(fabs(norm[i]) > max)
		{
			max = fabs(norm[i]);
			index = i;
		}
		return index;
}

//得到除索引index外的另外两个坐标值
void prjPtTo2D(double p[3],double pj[2],int index)
{
	int j=0;
	int i = 0;
	for(i=0; i<3; i++)
		if(i != index)
			pj[j++]=p[i];
}






