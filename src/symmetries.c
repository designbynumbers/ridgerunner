/* 

   symmetries.c : This file is part of the ridgerunner distribution. 
   
   The purpose is to allow us to symmetrize plCurves and variations over
   subgroups of SO(3) as a way to enforce symmetries during RR runs.

*/

void plc_rotation_matrix(plc_vector axis, plc_vector theta,matrix3_t *A)

/* Generates the matrix in SO(3) corresponding to rotation around axis with angle theta. */

{
  plc_vector u; /* This will be the normalized axis. */
  u = plc_normalize_vect(axis,NULL);
  
  (*A)[0][0] = cos(theta) + u.c[0]*u.c[0]*(1 - cos(theta));
  (*A)[0][1] = u.c[0]*u.c[1]*(1 - cos(theta)) - u.c[2]*sin(theta);
  (*A)[0][2] = u.c[0]*u.c[2]*(1 - cos(theta)) + u.c[1]*sin(theta);

  (*A)[1][0] = u.c[1]*u.c[0]*(1 - cos(theta)) + u.c[2]*sin(theta);
  (*A)[1][1] = cos(theta) + u.c[1]*u.c[1]*(1 - cos(theta));
  (*A)[1][2] = u.c[1]*u.c[2]*(1 - cos(theta)) - u.c[2]*sin(theta);

  (*A)[2][0] = u.c[2]*u.c[0]*(1 - cos(theta)) - u.c[1]*sin(theta);
  (*A)[2][1] = u.c[2]*u.c[1]*(1 - cos(theta)) + u.c[0]*sin(theta);
  (*A)[2][2] = cos(theta) + u.c[2]*u.c[2]*(1 - cos(theta));
}
  

void plc_reflection_matrix(plc_vector axis, matrix3_t *A)

/* Generates the matrix in O(3) corresponding to reflection across the plane normal to axis. */
/* This turns out to be 

   A[i][j] = delta(i,j) - 2 axis[i]*axis[j]/||axis||^2

*/

{
  plc_vector u;
  int i,j;

  u = plc_normalize_vect(axis,NULL);
  
  for(i=0;i<3;i++) {
   
    for(j=0;j<3;j++) {

      (*A)[i][j] = (i == j ? 1:0) - 2*u.c[i]*u.c[j];

    }
    
  }

}

plc_symmetry *plc_symmetry_new(plCurve *model)

/* Generate an (empty) symmetry which will be applied to the curve in model. 
   We are basically using the model to pick off number-of-components and 
   vertices-per-component data. */

{
  plc_symmetry *sym;
  int cp;

  sym = calloc_or_die(1,sizeof(plc_symmetry));
  sym->curve = model;

  sym->target = calloc(model->nc,sizeof(struct plc_vertex_loc **));

  for(cp=0;cp<model->nc;cp++) {

    sym->target[cp] = calloc(model->cp[cp].nv,sizeof(struct plc_vertex_loc *));

    for(vt=0;vt<model->cp[cp].nv;vt++) {
      
      sym->target[cp][vt] = calloc_or_die(1,sizeof(struct plc_vertex_loc));

    }

  }

}
  
void plc_symmetry_free(plc_symmetry **Sym) 

/* Kills a plc_symmetry. Safe for multiple calls on 
   the same pointer, which is set to NULL after the call. */

{
  int cp, vt;
  plc_symmetry *sym;
  sym = *Sym;

  if (sym != NULL) {

    if (sym->target != NULL) {
      
      for(cp=0;cp<sym->curve->nc;cp++) {
	
	if (sym->target[cp] != NULL) { 

	  for(vt=0;vt<sym->curve->cp[cp].nv;vt++) {

	    if (sym->target[cp][vt] != NULL) {

	      free(sym->target[cp][vt]); 

	    }

	  }

	}
	
      }
    
      free(sym->target);
    
    }

    free(sym);
    
  }

  (*Sym) = NULL;

}   

plc_symmetry *plc_build_symmetry(matrix3_t A,plCurve *L)

/* Uses plc_nearest_vertex to try to figure out the "intended" target of each vertex under the symmetry.
   This can fail if the curve is not symmetric to begin with to within < (1/2) an edgelength. In this case, 
   you should build the symmetry map yourself and then symmetrize the results. */

{
  int cp,afterAcp;
  int vt,afterAvt;
  plc_vector afterA;

  struct plc_nearest_vertex_pc_data *pc_data = NULL;

  /* We are going to mark off "used" vertex targets as they come in. */
  /* If a vertex target is hit twice, we fail. */

  used = calloc(plc_num_verts(L),sizeof(bool));

  /* Now we figure out the (putative) target of each vertex under A */

  plc_symmetry *build;
  build = plc_symmetry_new(L);
  build->transform = A;  
  
  for(cp=0;cp<L->nc;cp++) {

    for(vt=0;vt<L->cp[cp].nv;vt++) {

      rr_Axy(A,L->cp[cp].vt[vt],afterA);
      plc_nearest_vertex(afterA,L,&afterAcp,&afterAvt,&pc_data,NULL);
      
      build->target[cp][vt].cp = afterAcp;
      build->target[cp][vt].vt = afterAvt;

      if (used[plc_vertex_num(L,afterAcp,afterAvt)]) { /* We fail. Cleanup and quit. */
	
	free(used);
	free(build);
	return NULL;

      } else { /* Mark this off for the future */

	used[plc_vertex_num(L,afterAcp,afterAvt)] = true;

      }
     
    }
      
  }

  /* We survived this far, so the build must have worked. */

  free(used);
  return build;

}

plc_symmetry *plc_compose_symmetries(plc_symmetry *A,plc_symmetry *B)

/* Combine two symmetries in the order A first, then B (or matrix multiplication BA). */

{
  plc_symmetry *build;
  plCurve *L;

  L = A->curve;
  build = plc_symmetry_new(L);

  int cp,vt;

  for(cp=0;cp<L->nc;cp++) {

    for(vt=0;vt<L->cp[cp].nv;vt++) {

      build->target[cp][vt] = B->target[A->target[cp][vt].cp][A->target[cp][vt].vt];

    }
    
  }

  rr_AB(B->transform,A->transform,build->transform);

  return build;

}

plc_symmetry_group *plc_rotation_group(plCurve *L,plc_vector axis, int n)

/* Creates the symmetry group Z/nZ of rotations around axis. */

{
  double TWO_PI = 6.2831853071795864769;
  plc_symmetry_group *build;
  build = calloc(1,sizeof(plc_symmetry_group));
  
  build->n = n;
  build->sym = calloc(n,sizeof(plc_symmetry *));

  /* We now try to set up the group. */

  matrix3_t A;
  plc_identity_matrix(&A);
  build->sym[0] = plc_build_symmetry(A,L);

  int i;
  for(i=1;i<n;i++) {
    plc_rotation_matrix(axis,i*(TWO_PI/(n-1)),&A);
    build->sym[i] = plc_build_symmetry(A,L);
  }

  /* We need to scan and make sure the builds succeeded. */

  for(i=0;i<n;i++) {

    if (build->sym[i] == NULL) {

      int j;
      for(j=0;j<n;j++) { plc_symmetry_free(&(build->sym[j])); } /* Uses the fact that plc_symmetry_free is NULL-safe. */
      free(build->sym);
      free(build);

      return NULL;

    }

  }

}

plc_symmetry_group *plc_reflection_group(plCurve *L,plc_vector axis)

/* Creates the Z/2Z group of reflections over the plane normal to axis. */

{
  plc_symmetry_group *build;
  build = calloc(1,sizeof(plc_symmetry_group));
  
  build->n = 2;
  build->sym = calloc(2,sizeof(plc_symmetry *));

  /* We now try to set up the group. */

  matrix3_t A;
  plc_identity_matrix(&A);
  build->sym[0] = plc_build_symmetry(A,L);

  plc_reflection_matrix(axis,&A);
  build->sym[1] = plc_build_symmetry(A,L);

  /* We need to scan and make sure the builds succeeded. */

  for(i=0;i<2;i++) {

    if (build->sym[i] == NULL) {

      int j;
      for(j=0;j<n;j++) { plc_symmetry_free(&(build->sym[j])); } /* Uses the fact that plc_symmetry_free is NULL-safe. */
      free(build->sym);
      free(build);

      return NULL;

    }

  }

}

/* We are at last ready to symmetrize a curve over a group! */
