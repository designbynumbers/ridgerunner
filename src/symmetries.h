/* 

   symmetries.h : Includes basic matrix types and so forth for 3-space vector operations. 
   We use the cblas standard to store our matrices in row-major order.

*/

typedef double matrix3_t[3][3];

/* We will need a few basic matrix operations: matrix-vector product and matrix-matrix product. 
   Of course, we can get this from cblas, but given the cblas interface, it's convenient to have
   wrapper macros: */

/* void cblas_dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY); */

#define rr_Axy(A,x,y)  /* Computes the matrix-plc_vector product Ax and stores it in y. */    
  cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1.0,A,3,x.c,1,0,y,1);

/*void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc); */

#define rr_AB(A,B,C)     /* Computes the matrix-matrix product AB and stores it in C. */
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, A, 3, B, 3, 0.0, C, 3);

/* We now introduce a data type encoding a symmetry of a link. Such a symmetry has to include
   both a geometric transformation of space AND a corresponding map from each vertex of the 
   plcurve to a target vertex. */

typedef struct plc_vertex_loc {

  int cp;
  int vt;

}

typedef struct plc_symmetry_type {

  matrix3_t transform;
  plCurve   *curve;
  struct plc_vertex_loc **target; /* Array of curve->nc arrays of curve->cp[cp].nv arrays of plc_vertex_loc */ 

} plc_symmetry;

typedef struct plc_symmetry_group_type {

  int n;
  plc_symmetry **sym;

} plc_symmetry_group;

void plc_identity_matrix(matrix3_t *A);
void plc_rotation_matrix(plc_vector axis, plc_vector angle,matrix3_t *A);
void plc_reflection_matrix(plc_vector axis,matrix3_t *A);

/* We now define a high level interface for dealing with symmetries. */

plc_symmetry *plc_symmetry_new(plCurve *model);
void plc_symmetry_free(plc_symmetry *A);

/* This creates a plc_symmetry from a transform by searching to try to figure
   out the "intended" target of each vertex under the transform A. */
plc_symmetry *plc_build_symmetry(matrix3_t A,plCurve *L);

/* This is a combination of matrix multiplication and applying the permutation
   of vertices in the symmetries to build a new symmetry (matrix product BA). 
   Returns NULL on fail. */
plc_symmetry *plc_compose_symmetries(plc_symmetry *A,plc_symmetry *B);

/* We define a couple of standard groups as well. Return NULL if the build fails. */
/* Remember that the curves have to basically have the desired symmetry to start. */
plc_symmetry_group *plc_rotation_group(plCurve *L,plc_vector axis, int n);
plc_symmetry_group *plc_reflection_group(plCurve *L,plc_vector axis);

/* This symmetrizes a plCurve over a group. */
void plc_symmetrize(plCurve *L,plc_symmetry_group *G);




