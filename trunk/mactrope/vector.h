/*
 * 
 * Prototypes for routines in vector.c.
 *
 * $Id: vector.h,v 1.1 2004-08-16 22:35:30 michaelp Exp $
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "octrope.h"

/*
 * Prototypes for vector routines. 
 *
 */

octrope_vector octrope_vplus(octrope_vector A,octrope_vector B);
octrope_vector octrope_vminus(octrope_vector A,octrope_vector B);
octrope_vector octrope_cross(octrope_vector A,octrope_vector B);
octrope_vector octrope_scalarmult(double x,octrope_vector A);

/*
 * Macro replacements (requires some reprogramming)
 *
 */

#define octrope_dot(A,B)    (A.c[0]*B.c[0] + A.c[1]*B.c[1] + A.c[2]*B.c[2])
#define octrope_norm(A)     sqrt(octrope_dot(A,A))
#define octrope_vadd(A,B)  \
                    A.c[0] += B.c[0]; A.c[1] += B.c[1]; A.c[2] += B.c[2];
#define octrope_vsub(A,B)  \
                    A.c[0] -= B.c[0]; A.c[1] -= B.c[1]; A.c[2] -= B.c[2];
#define octrope_vsmult(s,V) V.c[0] *= s; V.c[1] *= s; V.c[2] *= s;
  /* Add a multiple of B to A */
#define octrope_vmadd(A,s,B)  (A).c[0] += (s)*(B).c[0]; \
                              (A).c[1] += (s)*(B).c[1]; \
                              (A).c[2] += (s)*(B).c[2];
  /* A = B + s(C-B)                               *
   * equivalent to                                *
   *   A = C; vsub(A,B); vsmult(s,A); vsadd(A,B); */
#define octrope_vweighted(A,s,B,C)  \
    (A).c[0] = (B).c[0] + s*((C).c[0] - (B).c[0]); \
    (A).c[1] = (B).c[1] + s*((C).c[1] - (B).c[1]); \
    (A).c[2] = (B).c[2] + s*((C).c[2] - (B).c[2]); 
