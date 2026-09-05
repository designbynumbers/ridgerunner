/*
 * Main routines for liboctrope.a
 *
 * $Id: octrope.c,v 1.51 2007-11-28 16:16:41 ashted Exp $
 *
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of liboctrope.
   
liboctrope is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

liboctrope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with liboctrope; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <float.h>
#include <ctype.h>
#include <stdarg.h>

#include"plCurve.h"
#include "octrope.h"
#if SIZEOF_INT < 4
#include "octcnv2.h"
#else
#include "octcnv4.h"
#endif
#include "octmem.h"
#include "octrope_tls.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef min
  #define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef max
  #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

/* Parallelize the edge loop when the curve is large enough to amortize OMP. */
#ifndef OCTROPE_OMP_MIN_EDGES
#define OCTROPE_OMP_MIN_EDGES 64
#endif

typedef struct vertarray_type { /* used for sorting vertices or edges */
  plc_vector v; /* edge midpoint */
  double length;
  int component;
  int vertex;
  int oct[2];
} varr;

struct box_type { /* A 3-dimensional box -- a part of the octree */
  double min_x;
  double min_y;
  double min_z;
  double max_x;
  double max_y;
  double max_z;
  plc_vector corner[8];
  plc_vector center;
  double half_diagonal;
  struct box_type *sub_boxes[8];
  int num_edges;
  varr *edges;
};
typedef struct box_type box;

#ifdef DEBUG
static OCTROPE_TLS int strut_check_count;
static OCTROPE_TLS int rc_count[5];
static OCTROPE_TLS int wc_count[5];
static OCTROPE_TLS int bi_count[2][9];
#endif

/* Per-call / per-thread tree state (reentrant across concurrent octrope calls). */
static OCTROPE_TLS int points_per_box = 1;
static OCTROPE_TLS int octree_levels = 0;
static int requested_levels = 0; /* process-wide preference from octrope_set_levels */
static OCTROPE_TLS varr *by_oct;
static OCTROPE_TLS int num_edges; /* Number of edges in the link */

static int debug_level = 0;

int octrope_error_num;      /* These have to live in a specific SOURCE file, don't they? */
char octrope_error_str[1024];

void octrope_set_debug(int level) 
{
  octrope_error_num = octrope_error_str[0] = 0;
  debug_level = abs(level);
}

inline int octrope_debug_level()
{
  octrope_error_num = octrope_error_str[0] = 0;
  return debug_level;
}

void octrope_set_levels(const int levels) 
{
  octrope_error_num = octrope_error_str[0] = 0;
  requested_levels = abs(levels);
}

void octrope_set_threads(int nthreads)
{
  octrope_error_num = octrope_error_str[0] = 0;
#ifdef _OPENMP
  if (nthreads > 0) {
    omp_set_num_threads(nthreads);
  }
#else
  (void)nthreads;
#endif
}

int octrope_get_threads(void)
{
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

/* Estimate the memory needed to work with an num_edges-edge link */
int octrope_est_mem(const int num_edges)
{
  int my_octree_levels;
  int last_power_of_8;
  int total_boxes;
  int temp_int;
  int log_2;

  octrope_error_num = octrope_error_str[0] = 0;
  my_octree_levels = octree_levels;
  temp_int = num_edges;
  log_2 = 0;
  while (temp_int) {
    log_2++;
    temp_int >>= 1;
  }
  if (my_octree_levels == 0) { /* Caller did not specify number of levels */
    my_octree_levels = log_2 * 3;
    my_octree_levels += 3; /* So that the integer division is a "ceil" */
    my_octree_levels /= 4; /* rather than a "floor".                   */
  } else {
    /* No need to go below 1 point per box */
    my_octree_levels = min(my_octree_levels,log_2);
  }
  /* Make sure octal tag stays in 2 ints */
  my_octree_levels = min(my_octree_levels,2*octree_max_bits+1);

  last_power_of_8 = (log_2+2)/3;
  if (last_power_of_8 <= my_octree_levels) {
    total_boxes = ceil(((1 << 3*last_power_of_8) - 1)/7) +
                  num_edges*(my_octree_levels - last_power_of_8);
  } else {
    total_boxes = ceil(((1 << 3*my_octree_levels) -1)/7);
  }

  return total_boxes*sizeof(box)+my_octree_levels*sizeof(box *) +
         3*num_edges*sizeof(varr *) + num_edges*sizeof(varr);
}

void vsort_test(varr *A,int dim,int l,int h) 
{
  int s;

  for (s = l; s < h-1; s++) {
    printf(" %%%3.3f%%\n",A[s].v.c[dim]); 
    if (A[s].v.c[dim] > A[s+1].v.c[dim]) {
      printf("---- Sort error ----\n");
    }
  }
  printf(" %%%3.3f%%\n",A[s].v.c[dim]);
}

/* 
 * Sort a vertex array 
 *
 */
void sort_varr(varr **A,int dim,int l,int h) 
{ 
  int s,t;
  varr *temp;
  double pivot_val;

  if (h - l <= 10) { /* for sorting the last few, use an insertion sort */
    for (s = l+1; s <= h; s++) {
      temp = A[s];
      t = s;
      while ((t > l) && (A[t-1]->v.c[dim] > temp->v.c[dim])) {
        A[t] = A[t-1];
        t--;
      }
      A[t] = temp;
    }
  } else { /* quicksort using median-of-three to find pivot */
    /* Pick the pivot and put it in the h slot */
    s = (h+l)/2;
    /* Make sure that the "middlest" of the three is in A[h] */
    if (A[h]->v.c[dim] < A[s]->v.c[dim]) {
      if (A[s]->v.c[dim] < A[l]->v.c[dim]) { /* h < s < l */
        temp = A[h]; A[h] = A[s]; A[s] = temp;
      } else if (A[h]->v.c[dim] < A[l]->v.c[dim]) { /* h < l <= s */
        temp = A[h]; A[h] = A[l]; A[l] = temp;
      } 
    } else if (A[l]->v.c[dim] < A[s]->v.c[dim]) { /* l < s <= h */
      temp = A[h]; A[h] = A[s]; A[s] = temp;
    } else if (A[h]->v.c[dim] > A[l]->v.c[dim]) { /* s <= l < h */ 
      temp = A[h]; A[h] = A[l]; A[l] = temp;
    }
    pivot_val = A[h]->v.c[dim]; 
    s = l-1; t = h; 
    for (;;) {
      s++;
      while (A[s]->v.c[dim] < pivot_val && s < h) { s++; }
      t--;
      while (A[t]->v.c[dim] >= pivot_val && t > l) { t--; }
      if (s >= t) break;
      temp = A[s]; A[s] = A[t]; A[t] = temp;
    }
    temp = A[s]; A[s] = A[h]; A[h] = temp;
    t = s;
    while (s > l && A[--s]->v.c[dim] == pivot_val);
    if (l < s) {
      sort_varr(A,dim,l,s);
    }
    s = t;
    while (s < h && A[++s]->v.c[dim] == pivot_val);
    if (s < h) {
      sort_varr(A,dim,s,h);
    }
  }
} /* sort_varr */

/* 
 * Sort the pieces of the octree
 *
 */
#define o_g(A,B) (A[0] > B[0] || (A[0] == B[0] && A[1] > B[1]))
#define o_l(A,B) (A[0] < B[0] || (A[0] == B[0] && A[1] < B[1]))
#define o_ge(A,B) (A[0] > B[0] || (A[0] == B[0] && A[1] >= B[1]))
#define o_le(A,B) (A[0] < B[0] || (A[0] == B[0] && A[1] <= B[1]))
#define o_e(A,B) (A[0] == B[0] && A[1] == B[1])
void sort_oct(varr *A,int l,int h) 
{ 
  int s,t;
  varr temp;
  int pivot_val[2];

  if (h - l <= 10) { /* for sorting the last few, use an insertion sort */
    for (s = l+1; s <= h; s++) {
      temp = A[s];
      t = s;
      while ((t > l) && o_g(A[t-1].oct,temp.oct)) {
        A[t] = A[t-1];
        t--;
      }
      A[t] = temp;
    }
  } else { /* quicksort using median-of-three to find pivot */
    /* Pick the pivot and put it in the h slot */
    s = (h+l)/2;
    if (o_l(A[h].oct,A[s].oct)) {
      if (o_l(A[s].oct,A[l].oct)) { /* h < s < l */
        temp = A[h]; A[h] = A[s]; A[s] = temp;
      } else if (o_l(A[h].oct,A[l].oct)) { /* h < l <= s */
        temp = A[h]; A[h] = A[l]; A[l] = temp;
      } 
    } else if (o_l(A[l].oct,A[s].oct)) { /* l < s <= h */
      temp = A[h]; A[h] = A[s]; A[s] = temp;
    } else if (o_g(A[h].oct,A[l].oct)) { /* s <= l < h */ 
      temp = A[h]; A[h] = A[l]; A[l] = temp;
    }
    s = l-1; t = h; pivot_val[0] = A[t].oct[0]; pivot_val[1] = A[t].oct[1];
    for (;;) {
      s++;
      while (o_l(A[s].oct,pivot_val) && s < h) { s++; }
      t--;
      while (o_ge(A[t].oct,pivot_val) && t > l) { t--; }
      if (s >= t) break;
      temp = A[s]; A[s] = A[t]; A[t] = temp;
    }
    temp = A[s]; A[s] = A[h]; A[h] = temp;
    t = s; 
    s--;
    while (o_e(A[s].oct,pivot_val) && s > l) { s--; }
    if (l < s) {
      sort_oct(A,l,s);
    }
    s = t;
    s++;
    while (o_e(A[s].oct,pivot_val) && s < h) { s++; }
    if (s < h) {
      sort_oct(A,s,h);
    }
  }
} /* sort_oct */

#define clear_subboxes(B) \
  B->sub_boxes[0] = B->sub_boxes[1] = B->sub_boxes[2] = B->sub_boxes[3] = \
  B->sub_boxes[4] = B->sub_boxes[5] = B->sub_boxes[6] = B->sub_boxes[7] = NULL

void initialize_box(box *B, varr *varr_ptr, const int leaf_box, 
                    const plc_vector *vtx0,
                    const plc_vector *vtx1) 
{
  B->min_x = min(vtx0->c[0], vtx1->c[0]);
  B->min_y = min(vtx0->c[1], vtx1->c[1]);
  B->min_z = min(vtx0->c[2], vtx1->c[2]);
  B->max_x = max(vtx0->c[0], vtx1->c[0]);
  B->max_y = max(vtx0->c[1], vtx1->c[1]);
  B->max_z = max(vtx0->c[2], vtx1->c[2]);
  clear_subboxes(B);
  B->num_edges = (leaf_box) ? 1 : 0;
  B->edges = varr_ptr;
}

#define octree_leaf(B) (B->sub_boxes[0] == NULL && \
                        B->sub_boxes[1] == NULL && \
                        B->sub_boxes[2] == NULL && \
                        B->sub_boxes[3] == NULL && \
                        B->sub_boxes[4] == NULL && \
                        B->sub_boxes[5] == NULL && \
                        B->sub_boxes[6] == NULL && \
                        B->sub_boxes[7] == NULL)
void print_octree(box *B,const char *tag,int *edges_found)
{
  int cnt;
  char new_tag[80];
  
  printf("%s %7.3lf %7.3lf, %7.3lf %7.3lf, %7.3lf %7.3lf -- %d\n",
    tag,B->min_x,B->max_x,B->min_y,B->max_y,B->min_z,B->max_z,
    B->num_edges);
  if (octree_leaf(B)) {
    *edges_found += B->num_edges;
    for (cnt = 0; cnt < B->num_edges; cnt++) {
      printf("  (%3.3lf,%3.3lf,%3.3lf) %d:%d (%o|%o)\n",
        B->edges[cnt].v.c[0], B->edges[cnt].v.c[1], B->edges[cnt].v.c[2],
        B->edges[cnt].component, B->edges[cnt].vertex, B->edges[cnt].oct[0],
        B->edges[cnt].oct[1]);
    }  
  } else {
    strcpy(new_tag,tag);
    strcat(new_tag,"-+");
    for (cnt = 0; cnt < 8; cnt++) {
      if (B->sub_boxes[cnt] != NULL) {
        print_octree(B->sub_boxes[cnt],new_tag,edges_found);
      }
    }
  }
} /* print_octree */

void show_varr(varr **A) {
  int vcnt;

  for (vcnt = 0; vcnt < num_edges; vcnt++) {
    printf("(%3.3f,%3.3f,%3.3f) cmp=%d  vert=%d  idx=%o,%o\n",
      A[vcnt]->v.c[0], A[vcnt]->v.c[1], A[vcnt]->v.c[2],
      A[vcnt]->component, A[vcnt]->vertex, A[vcnt]->oct[0], A[vcnt]->oct[1]);
  }
}

void show_octal() {
  int vcnt;

  printf("Edges orderd by Octal Tag:\n");
  for (vcnt = 0; vcnt < num_edges; vcnt++) {
    printf("(%3.3f,%3.3f,%3.3f) cmp=%d  vert=%d  idx=%o,%o\n",
      by_oct[vcnt].v.c[0], by_oct[vcnt].v.c[1], by_oct[vcnt].v.c[2],
      by_oct[vcnt].component, by_oct[vcnt].vertex, by_oct[vcnt].oct[0],
      by_oct[vcnt].oct[1]);
  }
}

void po(box *B) 
{
  int edges;

  edges = 0;
  printf("The Octree:\n");
  print_octree(B,"",&edges);
  printf("Found %d edges.\n",edges);
}

/* sets global num_edges */
box *build_octree(const plCurve *L,const int sl_size,
                  void *mem,const int memsize) {
  box **tree_swath; /* vertical slice of tree */
  int total_boxes;
  varr *varr_ptr;
  varr **by_x,**by_y,**by_z;
  int box_num;
  int change_level;
  int temp_int,log_2;
  int cmp,vert;
  int vcnt,cnt,cnt2;
  int boxes;
  varr temp_varr;
  plc_vector temp_v;
  int oct1,oct2; 
  plc_vector *vtx0,*vtx1;

  num_edges = plc_num_edges(L);
  /* allocate all the space we're going to need (and more) . . . */
  temp_int = num_edges;
  log_2 = 0;
  while (temp_int) {
    log_2++;
    temp_int >>= 1;
  }
  if (requested_levels == 0) { /* Caller did not specify number of levels */
    octree_levels = log_2 * 3;
    octree_levels += 3; /* So that the integer division is a "ceil" */
    octree_levels /= 4; /* rather than a "floor".                   */
  } else {
    /* No need to go below 1 point per box */
    octree_levels = min(requested_levels,log_2);
  }
  /* Make sure octal tag stays in 2 ints */
  octree_levels = min(octree_levels,2*octree_max_bits+1);
#ifdef DEBUG
  if (octrope_debug_level() > 3) {
    printf(
      "Bytes: %zu for tree_swath, %zu for by_x/by_y/by_z, %zu for by_oct.\n",
      octree_levels*sizeof(box *), 3*num_edges*sizeof(varr *),
      num_edges*sizeof(varr));
  }
#endif

  reserve_octmem(octree_levels*sizeof(box *) +
                 3*num_edges*sizeof(varr *) +
                 num_edges*sizeof(varr),mem,memsize);
  if ((tree_swath = (box **)oct_alloc(octree_levels,sizeof(box*))) == NULL) {
    sprintf(octrope_error_str,
      "Unable to allocate %d tree hooks in build_octree.\n", octree_levels);
    octrope_error_num = 1;
    return NULL;
  }
  if ((by_x = (varr **)oct_alloc(num_edges,sizeof(varr *))) == NULL) {
    sprintf(octrope_error_str,
      "Unable to reserve space for %d vertices in by_x.\n", num_edges);
    octrope_error_num = 2;
    return NULL;
  }
  if ((by_y = (varr **)oct_alloc(num_edges,sizeof(varr *))) == NULL) {
    sprintf(octrope_error_str,
      "Unable to reserve space for %d vertices in by_y.\n", num_edges);
    octrope_error_num = 3;
    return NULL;
  }
  if ((by_z = (varr **)oct_alloc(num_edges,sizeof(varr *))) == NULL) {
    sprintf(octrope_error_str,
      "Unable to reserve space for %d vertices in by_z.\n", num_edges);
    octrope_error_num = 4;
    return NULL;
  }
  if ((by_oct = (varr *)oct_alloc(num_edges,sizeof(varr))) == NULL) {
    sprintf(octrope_error_str,
      "Unable to reserve space for %d vertices in by_oct.\n", num_edges);
    octrope_error_num = 5;
    return NULL;
  }
  vcnt = 0;
  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert=0; 
         vert < ((L->cp[cmp].open) ? L->cp[cmp].nv-1 : L->cp[cmp].nv); 
         vert++) {
      temp_v = L->cp[cmp].vt[vert+1];
      /* Store center of edge */
      temp_varr.v.c[0] = (L->cp[cmp].vt[vert].c[0] + temp_v.c[0])/2;
      temp_varr.v.c[1] = (L->cp[cmp].vt[vert].c[1] + temp_v.c[1])/2;
      temp_varr.v.c[2] = (L->cp[cmp].vt[vert].c[2] + temp_v.c[2])/2;
      plc_M_sub_vect(temp_v,L->cp[cmp].vt[vert]);
      temp_varr.length = plc_M_norm(temp_v);
      temp_varr.component = cmp;
      temp_varr.vertex = vert;
      temp_varr.oct[0] = temp_varr.oct[1] = 0; 
      by_oct[vcnt] = temp_varr;
      by_x[vcnt] = by_y[vcnt] = by_z[vcnt] = &by_oct[vcnt];
      vcnt++;
    }
  }
  sort_varr(by_x,0,0,num_edges-1);
#ifdef DEBUG
  if (octrope_debug_level() > 7) {
    printf("Edges ordered by X:\n");
    show_varr(by_x);
  }
#endif
  sort_varr(by_y,1,0,num_edges-1);
#ifdef DEBUG
  if (octrope_debug_level() > 7) {
    printf("Edges ordered by Y:\n");
    show_varr(by_y);
  }
#endif
  sort_varr(by_z,2,0,num_edges-1);
#ifdef DEBUG
  if (octrope_debug_level() > 7) {
    printf("Edges ordered by Z:\n");
    show_varr(by_z);
  }
#endif
  /* Now set the octal tags for the final sort */
  /*
   * Redefine boxes to be the number of leaf boxes in a given direction and
   * find the number of points per "row".
   */
  boxes = 1 << (octree_levels-1); 
  points_per_box = max(ceil(1.0*num_edges/boxes),1);
#ifdef DEBUG
  if (octrope_debug_level() > 0) {
    printf("points_per_box is %d.\n",points_per_box);
  }
#endif
  for (cnt = 1; cnt < boxes; cnt++) { /* Start at 1, since |=0 does nothing */
    for (cnt2 = 0; cnt2 < points_per_box && 
                   cnt * points_per_box + cnt2 < num_edges; cnt2++) {
      vcnt = cnt * points_per_box + cnt2;
      oct1 = bin_to_oct[cnt / octree_max_bin];
      oct2 = bin_to_oct[cnt & (octree_max_bin-1)];
      by_x[vcnt]->oct[0] |= oct1;
      by_x[vcnt]->oct[1] |= oct2;
      by_y[vcnt]->oct[0] |= oct1 << 1;
      by_y[vcnt]->oct[1] |= oct2 << 1;
      by_z[vcnt]->oct[0] |= oct1 << 2;
      by_z[vcnt]->oct[1] |= oct2 << 2;
    }
  }
  sort_oct(by_oct,0,num_edges-1);
  /* Now we reserve enough space for the boxes we need */
  total_boxes = octree_levels; 
  for (cnt = 1; cnt < num_edges; cnt++) {
    if (by_oct[cnt].oct[0] != by_oct[cnt-1].oct[0]) {
      temp_int = by_oct[cnt].oct[0] ^ by_oct[cnt-1].oct[0];
      change_level = octree_max_bits;
    } else {
      temp_int = by_oct[cnt].oct[1] ^ by_oct[cnt-1].oct[1];
      change_level = 0;
    }
    while (temp_int) {
      change_level++;
      temp_int >>= 3;
    }
    total_boxes += change_level;
  }
#ifdef DEBUG
  if (octrope_debug_level() > 0) {
    printf("Reserving %zu bytes for %d boxes on %d levels to hold %d edges.\n",
	   total_boxes*sizeof(box),total_boxes,octree_levels,num_edges);
    // the z modifier means the argument is a size_t (C99)
  }
#endif
  reserve_aux_mem(total_boxes*sizeof(box));
  for (cnt = 0; cnt < octree_levels; cnt++) {
    if ((tree_swath[cnt] = (box *)aux_alloc(1,sizeof(box))) == NULL) {
      sprintf(octrope_error_str,"Unable to allocate box %d in build_octree.\n",cnt);
      octrope_error_num = 6;
      return NULL;
    }
  }
#ifdef DEBUG
  if (octrope_debug_level() > 7) {
    show_octal();
  }

  if (octrope_debug_level() > 1) {
    printf("We have %d edges.\n",num_edges);
  }
#endif

  /* Now build that ol' octree :-) */
  varr_ptr = by_oct;
  /* Preload the boxes in the swath */
  initialize_box(tree_swath[0],varr_ptr,true,
                 &L->cp[varr_ptr->component].vt[varr_ptr->vertex],
                 &L->cp[varr_ptr->component].vt[varr_ptr->vertex+1]);
  for (cnt = 1; cnt < octree_levels; cnt++) {
    initialize_box(tree_swath[cnt],varr_ptr,false,
                   &L->cp[varr_ptr->component].vt[varr_ptr->vertex],
                   &L->cp[varr_ptr->component].vt[varr_ptr->vertex+1]);
    if (cnt > octree_max_bits) {
      box_num = (varr_ptr->oct[0] >> 3*(cnt - octree_max_bits - 1)) & 7;
    } else {
      box_num = (varr_ptr->oct[1] >> 3*(cnt - 1)) & 7;
    }
    tree_swath[cnt]->sub_boxes[box_num] = tree_swath[cnt-1];
  }
  varr_ptr++;
  while (varr_ptr <= &by_oct[num_edges]) {
#ifdef DEBUG
  if (octrope_debug_level() > 8) {
    po(tree_swath[octree_levels-1]);
  }
#endif
    while (varr_ptr <= &by_oct[num_edges-1] && 
           o_e(varr_ptr->oct,tree_swath[0]->edges->oct)) {
      tree_swath[0]->num_edges++;
      vtx0 = &L->cp[varr_ptr->component].vt[varr_ptr->vertex];
      vtx1 = &L->cp[varr_ptr->component].vt[varr_ptr->vertex+1];
      /* Make the box big enough that the whole edge fits inside. */
      if (vtx0->c[0] < vtx1->c[0]) {
        if (vtx0->c[0] < tree_swath[0]->min_x) { 
          tree_swath[0]->min_x = vtx0->c[0];
        }
        if (vtx1->c[0] > tree_swath[0]->max_x) { 
          tree_swath[0]->max_x = vtx1->c[0];
        }
      } else {
        if (vtx1->c[0] < tree_swath[0]->min_x) { 
          tree_swath[0]->min_x = vtx1->c[0];
        }
        if (vtx0->c[0] > tree_swath[0]->max_x) { 
          tree_swath[0]->max_x = vtx0->c[0];
        }
      }
      if (vtx0->c[1] < vtx1->c[1]) {
        if (vtx0->c[1] < tree_swath[0]->min_y) { 
          tree_swath[0]->min_y = vtx0->c[1];
        }
        if (vtx1->c[1] > tree_swath[0]->max_y) { 
          tree_swath[0]->max_y = vtx1->c[1];
        }
      } else {
        if (vtx1->c[1] < tree_swath[0]->min_y) { 
          tree_swath[0]->min_y = vtx1->c[1];
        }
        if (vtx0->c[1] > tree_swath[0]->max_y) { 
          tree_swath[0]->max_y = vtx0->c[1];
        }
      }
      if (vtx0->c[2] < vtx1->c[2]) {
        if (vtx0->c[2] < tree_swath[0]->min_z) { 
          tree_swath[0]->min_z = vtx0->c[2];
        }
        if (vtx1->c[2] > tree_swath[0]->max_z) { 
          tree_swath[0]->max_z = vtx1->c[2];
        }
      } else {
        if (vtx1->c[2] < tree_swath[0]->min_z) { 
          tree_swath[0]->min_z = vtx1->c[2];
        }
        if (vtx0->c[2] > tree_swath[0]->max_z) { 
          tree_swath[0]->max_z = vtx0->c[2];
        }
      }
      varr_ptr++;
    }
    /* Find level at which this differs from the previous edge */
    if (varr_ptr > &by_oct[num_edges-1]) {
      change_level = octree_levels-1; /* we're done building the tree */
    } else {
      if (varr_ptr->oct[0] != tree_swath[0]->edges->oct[0]) {
        temp_int = varr_ptr->oct[0] ^ tree_swath[0]->edges->oct[0];
        change_level = octree_max_bits;
      } else {
        temp_int = varr_ptr->oct[1] ^ tree_swath[0]->edges->oct[1];
        change_level = 0;
      }
      while (temp_int) {
        change_level++;
        temp_int >>= 3;
      }
    }
    if (change_level > octree_levels) {
      sprintf(octrope_error_str,
        "Programming error.  Change_level > Octree_levels.\n");
      octrope_error_num = 7;
      return NULL;
    }
    for (cnt = 1; cnt <= change_level; cnt++) {
      /* Set the corner values */
      for (cnt2 = 0; cnt2 < 8; cnt2++) {
        tree_swath[cnt-1]->corner[cnt2].c[0] = 
           (cnt2 & 1) ? tree_swath[cnt-1]->max_x : tree_swath[cnt-1]->min_x;
        tree_swath[cnt-1]->corner[cnt2].c[1] = 
           (((cnt2 & 2)/2) ^ (cnt2 & 1)) ? 
             tree_swath[cnt-1]->max_y : tree_swath[cnt-1]->min_y;
        tree_swath[cnt-1]->corner[cnt2].c[2] = 
           (((cnt2 & 4)/4) ^ (cnt2 & 1)) ? 
             tree_swath[cnt-1]->max_z : tree_swath[cnt-1]->min_z;
      }
      /* Here we depend on corners 0 and 1 be diagonal from each other */
      tree_swath[cnt-1]->half_diagonal = 0.5*
        plc_norm(plc_vect_diff(tree_swath[cnt-1]->corner[0],
                                    tree_swath[cnt-1]->corner[1]));
      plc_M_vweighted(tree_swath[cnt-1]->center,0.5,
        tree_swath[cnt-1]->corner[0], tree_swath[cnt-1]->corner[1]);

      tree_swath[cnt]->num_edges += tree_swath[cnt-1]->num_edges;
      if (tree_swath[cnt-1]->min_x < tree_swath[cnt]->min_x) {
        tree_swath[cnt]->min_x = tree_swath[cnt-1]->min_x;
      }
      if (tree_swath[cnt-1]->min_y < tree_swath[cnt]->min_y) {
        tree_swath[cnt]->min_y = tree_swath[cnt-1]->min_y;
      }
      if (tree_swath[cnt-1]->min_z < tree_swath[cnt]->min_z) {
        tree_swath[cnt]->min_z = tree_swath[cnt-1]->min_z;
      }
      if (tree_swath[cnt-1]->max_x > tree_swath[cnt]->max_x) {
        tree_swath[cnt]->max_x = tree_swath[cnt-1]->max_x;
      }
      if (tree_swath[cnt-1]->max_y > tree_swath[cnt]->max_y) {
        tree_swath[cnt]->max_y = tree_swath[cnt-1]->max_y;
      }
      if (tree_swath[cnt-1]->max_z > tree_swath[cnt]->max_z) {
        tree_swath[cnt]->max_z = tree_swath[cnt-1]->max_z;
      }
      if (tree_swath[cnt-1]->num_edges <= points_per_box && cnt > 1) { 
        /* Make this a leaf node */
        clear_subboxes(tree_swath[cnt-1]);
      }
    }
    if (varr_ptr <= &by_oct[num_edges-1]) {
      for (cnt = change_level; cnt > 0; cnt--) {
        if ((tree_swath[cnt-1] = (box *)aux_alloc(1,sizeof(box))) == NULL) {
          sprintf(octrope_error_str,
            "Unable to allocate box %d in build_octree.\n",cnt-1);
          octrope_error_num = 8;
          return NULL;
        }
        initialize_box(tree_swath[cnt-1],varr_ptr,(cnt == 1),
          &L->cp[varr_ptr->component].vt[varr_ptr->vertex],
          &L->cp[varr_ptr->component].vt[varr_ptr->vertex+1]);
        if (cnt > octree_max_bits) {
          box_num = (varr_ptr->oct[0] >> 3*(cnt - octree_max_bits - 1)) & 7;
        } else {
          box_num = (varr_ptr->oct[1] >> 3*(cnt - 1)) & 7;
        }
        if (tree_swath[cnt]->sub_boxes[box_num] != NULL) {
          sprintf(octrope_error_str,
            "Programming error in octrope tree build, box full!\n");
          fprintf(stderr,"Context: cnt=%d box_num=%d octree_max_bits=%d\n",
                  cnt,box_num,octree_max_bits);
          fprintf(stderr,"         varr_ptr->oct = %d %d\n", varr_ptr->oct[0],
                  varr_ptr->oct[1]);
          octrope_error_num = 9;
          return NULL;
        }
        tree_swath[cnt]->sub_boxes[box_num] = tree_swath[cnt-1];
      }
    }
    varr_ptr++;
  }
  for (cnt = 0; cnt < 8; cnt++) {
    tree_swath[octree_levels-1]->corner[cnt].c[0] = (cnt & 1) ? 
      tree_swath[octree_levels-1]->max_x : tree_swath[octree_levels-1]->min_x;
    tree_swath[octree_levels-1]->corner[cnt].c[1] = 
      (((cnt & 2)/2) ^ (cnt & 1)) ?
      tree_swath[octree_levels-1]->max_y : tree_swath[octree_levels-1]->min_y;
    tree_swath[octree_levels-1]->corner[cnt].c[2] = 
      (((cnt & 4)/4) ^ (cnt & 1)) ? 
      tree_swath[octree_levels-1]->max_z : tree_swath[octree_levels-1]->min_z;
  }

  return tree_swath[octree_levels-1]; 
} /* build_octree */

static OCTROPE_TLS varr *wedgecheck_edge[2];
int wedge_check(varr *edge, const plc_vector v[3], 
                const plc_vector far_edge[3], double percent,
                const int which, const int last_edge) 
{
  static OCTROPE_TLS plc_vector snlc[2];
  static OCTROPE_TLS plc_vector rnlc[2];
  static OCTROPE_TLS double snc[2];
  static OCTROPE_TLS double rnc[2];
  plc_vector slab_normal;
  plc_vector ramp_normal;
  double slab_norm;
  double ramp_norm;
  double dot_product;
  plc_vector far_end;

  if (edge == wedgecheck_edge[which]) {
    slab_normal = snlc[which];
    ramp_normal = rnlc[which];
    slab_norm = snc[which];
    ramp_norm = rnc[which];
  } else {
    slab_normal = v[2];
    plc_M_sub_vect(slab_normal,v[1]);
    snlc[which] = slab_normal;
    ramp_normal = v[0];
    plc_M_sub_vect(ramp_normal,v[1]); /* points in opposite direction so we don't 
                                       have to calculate the_point again */
    rnlc[which] = ramp_normal;
    snc[which] = slab_norm = plc_M_norm(slab_normal);
    rnc[which] = ramp_norm = plc_M_norm(ramp_normal);
    wedgecheck_edge[which] = edge;
#ifdef DEBUG
    if (octrope_debug_level() > 2) {
      wc_count[4]++;
    }
#endif
  }
  plc_M_vweighted(far_end,percent,far_edge[1],far_edge[2]);
#ifdef DEBUG
  if (octrope_debug_level() > 6) {
    printf("Wedgecheck: (%3.3f,%3.3f,%3.3f)\n",
           far_end.c[0],far_end.c[1],far_end.c[2]);
  }
#endif

  plc_M_sub_vect(far_end,v[1]);
  dot_product = plc_M_dot(slab_normal,far_end)/slab_norm;
  if (dot_product > 0) { /* Inside slab */
#ifdef DEBUG
    if (octrope_debug_level() > 2) {
      wc_count[1]++;
    }
#endif
    return false;
  }
  dot_product = plc_M_dot(ramp_normal,far_end)/ramp_norm;
  if (dot_product > 0) { /* Behind ramp */
#ifdef DEBUG
    if (octrope_debug_level() > 2) {
      wc_count[3]++;
    }
#endif
    return false;
  }
#ifdef DEBUG
  if (octrope_debug_level() > 2) {
    wc_count[2]++; /* Inside ramp */
  }
#endif
  return true;
}

/*
 * ramp_check takes four vectors, a double and a boolean
 *   . One vector which, with the first vector below, determines the "ramp".
 *   . Two vectors which define a single segement.  This segment defines the 
 *     main body of the ramp, "the slab" as it were.  The segment runs from 
 *     the first vector to the second and the front of the slab is the region
 *     beyond the second vector.
 *   . One vector which specifies the point to be checked.
 *   . The double is half the length of the longest edge in the box
 *   . last_edge is true if this is the final edge in an open component.  In
 *     that case, there is no "front"--points which fall in front are
 *     considered to be within the slab.
 * and returns one of the following:
 *  1 -- point is in front of the slab
 *  2 -- point is within the slab 
 *  4 -- point is within the ramp (but not the slab)
 *  8 -- point is behind the slab and ramp
 *
 */
static OCTROPE_TLS varr *rampcheck_edge;
int ramp_check(varr *edge, const plc_vector v[3], 
               const plc_vector v3, const int last_edge) 
{
  static OCTROPE_TLS plc_vector slab_normal;
  static OCTROPE_TLS plc_vector ramp_normal;
  static OCTROPE_TLS double slab_norm;
  static OCTROPE_TLS double ramp_norm;
  plc_vector the_point;
  double dot_product;

#ifdef DEBUG
  if (octrope_debug_level() > 6) {
    printf("Rampcheck: (%3.3f,%3.3f,%3.3f)\n",v3.c[0],v3.c[1],v3.c[2]);
  }
#endif
  if (edge != rampcheck_edge) {
    slab_normal = v[2];
    plc_M_sub_vect(slab_normal,v[1]);
    ramp_normal = v[0];
    plc_M_sub_vect(ramp_normal,v[1]); /* points in opposite direction so we don't 
                                       have to calculate the_point again */
    slab_norm = plc_M_norm(slab_normal);
    ramp_norm = plc_M_norm(ramp_normal);
    rampcheck_edge = edge;
#ifdef DEBUG
    if (octrope_debug_level() > 2) {
      rc_count[4]++;
    }
#endif
  }
  the_point = v3;
  plc_M_sub_vect(the_point,v[1]);
  dot_product = plc_M_dot(slab_normal,the_point)/slab_norm;
  if (!last_edge && dot_product > slab_norm) { /* In front of slab */
#ifdef DEBUG
    if (octrope_debug_level() > 2) {
      rc_count[0]++;
    }
#endif
    return 1;
  }
  if (dot_product > 0) { /* Inside slab */
#ifdef DEBUG
    if (octrope_debug_level() > 2) {
      rc_count[1]++;
    }
#endif
    return 2;
  }
  dot_product = plc_M_dot(ramp_normal,the_point)/ramp_norm;
  if (dot_product > 0) { /* Behind ramp */
#ifdef DEBUG
    if (octrope_debug_level() > 2) {
      rc_count[3]++;
    }
#endif
    return 8;
  }
#ifdef DEBUG
  if (octrope_debug_level() > 2) {
    rc_count[2]++; /* Inside ramp */
  }
#endif
  return 4;
}

void octree_to_skel(box *B,FILE *outfile) 
{
  int cnt;
  
  if (outfile == NULL) {
    outfile = fopen("boxes.vect","w");
    if (outfile == NULL) {
      sprintf(octrope_error_str,"Can't open boxes.vect for write.\n");
      octrope_error_num = 10;
      return;
    }
    fprintf(outfile,"VECT\n");
  }
  if (octree_leaf(B)) {
    for (cnt = 0; cnt < 8; cnt++) {
      fprintf(outfile,"%3.3f %3.3f %3.3f\n",
        B->corner[cnt].c[0], B->corner[cnt].c[1], B->corner[cnt].c[2]);
    }
    fprintf(outfile,"\n");
  } else {
    /* Duplicated here so that it can be commented out if only */
    /* leaf boxes are desired */
    for (cnt = 0; cnt < 8; cnt++) {
      fprintf(outfile,"%3.3f %3.3f %3.3f\n",
        B->corner[cnt].c[0], B->corner[cnt].c[1], B->corner[cnt].c[2]);
    }
    fprintf(outfile,"\n");
    for (cnt = 0; cnt < 8; cnt++) {
      if (B->sub_boxes[cnt] != NULL) {
        octree_to_skel(B->sub_boxes[cnt],outfile);
      }
    }
  }
}

#define LU_DOUBLE_TYPE long double

bool octrope_solveMatrixLU(double A[2][2],double b[2],double x[2])

 /* Procedure uses Crout's method to perform an LU decomposition
    on A and solve the system Ax = b. If the matrix is singular,
    will return false. Otherwise, returns true. */
  
/* This is essentially copied directly from Numerical Recipes */
{ 
  LU_DOUBLE_TYPE B[2][2];
  int index[2];
  LU_DOUBLE_TYPE D;
 
  int      i,imax=0,j,k;
  LU_DOUBLE_TYPE   big,dum,temp,sum;
  LU_DOUBLE_TYPE   vv[3];
  LU_DOUBLE_TYPE   TINY = {1.0e-20};
    
  for (i=0;i<2;i++) for(j=0;j<2;j++) { B[i][j] = A[i][j]; }
  
  /* Now we get to work. */
    
  D = 1.0;                              /* No changes yet. */
  
  for(i=0;i<2;i++) {                     /* Get scaling for each row. */
    
    big = 0.0;
    
    for(j=0;j<2;j++) {
      
      if ((temp = fabsl(B[i][j])) > big) big = temp;
      
    }
    
    if (big == 0.0) {
      
      /* The matrix is really singular*/
      return false;
      
    }
    
    vv[i] = 1.0/big;
    
  }
  
  for(j=0;j<2;j++) {                       /* Loop over columns */
    
    for(i=0;i<j;i++) {
      
      sum = B[i][j];
      for(k=0;k<i;k++) {  sum -= B[i][k] * B[k][j]; }      
      B[i][j] = sum;
      
    }
    
    big = 0.0;                             /* Now we search for a pivot. */
    
    for(i=j;i<2;i++) {
      
      sum = B[i][j];
      
      for(k=0;k<j;k++) {
	
	sum -= B[i][k] * B[k][j];
	
      }
      
      B[i][j] = sum;
      
      if ((dum = vv[i] * fabsl(sum)) >= big) {   
	
	big = dum;
	imax = i;
	
      }
      
    }
    
    if (j != imax) {                        /* Must we interchange rows? */
      
      for(k=0;k<2;k++) {                    /* Perform the interchange.. */
	
	dum = B[imax][k];
	B[imax][k] = B[j][k];
	B[j][k] = dum;
	
      }
      
      D = -(D);                           /* Change parity of D. */
      
      //dum = vv[imax];                       /* Swap scale factors... */
      vv[imax] = vv[j];                     
      //vv[j] = dum;
      
    }
      
    index[j] = imax;
    
    if (B[j][j] == 0.0) B[j][j] = TINY;
    
    if (j != 1) {
      
      dum = 1.0/(B[j][j]);
      
      for (i=j+1;i<2;i++) {
	
	B[i][j] *= dum;
	
      }
      
    }
    
  }

  /* We are now done. We go ahead and solve the system. */

  int ii=-1,ip;

  /* First, we change x to the right hand side. */
    
  for(i=0;i<=1;i++) { x[i] = b[i]; }
  
  for(i=0;i<=1;i++) {   /* Now, we are ready to work.... */

    ip = index[i];
    sum = x[ip];
    x[ip] = x[i];

    if (ii >= 0) {

      for (j=ii;j<=i-1;j++) { sum -= B[i][j]*x[j]; }

    } else if (sum) { ii = i; }

    x[i] = sum;

  }

  for (i=1;i>=0;i--) { 

    sum = x[i];
    for (j=i+1;j<2;j++) { sum -= B[i][j]*x[j]; }
    x[i] = sum/B[i][i];

  }

  /* Last, we have to pass some backwards error analysis */

  LU_DOUBLE_TYPE bcheck[2];

  bcheck[0] = A[0][0]*x[0] + A[0][1]*x[1];
  bcheck[1] = A[1][0]*x[0] + A[1][1]*x[1];

  if (fabsl(bcheck[0] - b[0]) < 1e-8 && fabsl(bcheck[1] - b[1]) < 1e-8) { return true; }

  LU_DOUBLE_TYPE delta[2];

  delta[0] = b[0] - bcheck[0];
  delta[1] = b[1] - bcheck[1];

  /* First, we change x to the right hand side. */

  ii = -1;
  
  for(i=0;i<=1;i++) {   /* Now, we are ready to work.... */

    ip = index[i];
    sum = delta[ip];
    delta[ip] = delta[i];

    if (ii >= 0) {

      for (j=ii;j<=i-1;j++) { sum -= B[i][j]*delta[j]; }

    } else if (sum) { ii = i; }

    delta[i] = sum;

  }

  for (i=1;i>=0;i--) { 

    sum = delta[i];
    for (j=i+1;j<2;j++) { sum -= B[i][j]*delta[j]; }
    delta[i] = sum/B[i][i];

  }

  /* Now we add delta to x */

  x[0] += delta[0]; x[1] += delta[1];

  /* And recheck our solution */

  bcheck[0] = A[0][0]*x[0] + A[0][1]*x[1];
  bcheck[1] = A[1][0]*x[0] + A[1][1]*x[1];

  return (fabsl(bcheck[0] - b[0]) < 1e-8 && fabsl(bcheck[1] - b[1]) < 1e-8);

}

void find_poca(plc_vector edgeA[2], plc_vector edgeB[2],
	       double *t, double *u, double *dist)

/* This is a reference implementation of an algorithm of LUMELSKY (see the /doc directory) */
/* It outperforms the previous Ashton/Cantarella implementation of find_poca by about 20%, */
/* and is more accurate on all but a tiny fraction of the numerically difficult cases. */

{
  plc_vector d1,d2,d12;
  double T,U;

  //A = edgeA[0]; B = edgeA[1]; C = edgeB[0]; D = edgeB[1];

  plc_M_vect_diff(d1,edgeA[1],edgeA[0]);
  plc_M_vect_diff(d2,edgeB[1],edgeB[0]);
  plc_M_vect_diff(d12,edgeB[0],edgeA[0]);

  double D1,D2,R,S1,S2;

  D1 = plc_M_dot(d1,d1); 
  D2 = plc_M_dot(d2,d2);
  R  = plc_M_dot(d1,d2); 
  S1 = plc_M_dot(d1,d12); 
  S2 = plc_M_dot(d2,d12);

  double denom;

  denom = D1*D2 - R*R;

  //assert(fabs(D1) > 1e-10 && fabs(D2) > 1e-10);

  if (denom < 1e-10) {

    /* The situation is somewhat numerically unstable. The right thing to do is to repeat
       the calculation in long double to see if we can salvage things. */

    long double LD1,LD2,LR,LS1,LS2;
    long double LT,LU;

    LD1 = plc_M_dot(d1,d1); 
    LD2 = plc_M_dot(d2,d2);
    LR  = plc_M_dot(d1,d2); 
    LS1 = plc_M_dot(d1,d12); 
    LS2 = plc_M_dot(d2,d12);
    
    long double Ldenom;
    
    Ldenom = LD1*LD2 - LR*LR;
    
    if (Ldenom < 1e-25) {
      
      LT = 0; // We assume that the line segments are parallel. 
	
    } else {
      
      LT= (LS1*LD2 - LS2*LR)/Ldenom; // When it is accurate, this is a lot faster than LU decomposition.
      
      if (LT > 1) { LT = 1; }
      else if (LT < 0) { LT = 0; }
      
    }
    
    LU = ((LT)*LR - LS2)/LD2;
    
    if (LU > 1) { 
      
      LU = 1; 
      LT = ((LU)*LR + LS1)/LD1;
      if (LT > 1) { LT = 1; }
      else if (LT < 0) { LT = 0; }
      
    } else {
      
      if (LU < 0 ) { 
	
	LU = 0; 
	LT = ((LU)*LR + LS1)/LD1;
	if (LT > 1) { LT = 1; }
	else if (LT < 0) { LT = 0; }
	
      }
    }
    
    *dist = sqrt(
		 (LT*d1.c[0] - LU*d2.c[0] - d12.c[0]) * (LT*d1.c[0] - LU*d2.c[0] - d12.c[0]) + 
		 (LT*d1.c[1] - LU*d2.c[1] - d12.c[1]) * (LT*d1.c[1] - LU*d2.c[1] - d12.c[1]) +  
		 (LT*d1.c[2] - LU*d2.c[2] - d12.c[2]) * (LT*d1.c[2] - LU*d2.c[2] - d12.c[2]));
    
    *t = LT;
    *u = LU;

    return;
    
  } else {  /* In the branch where the numerics were difficult, we
	       solved the problem in higher precision and returned.
	       So if we've gotten to this point, the numerics are ok
	       in double precision and we can continue knowing that
	       the denominator is not too small to compute T and U
	       accurately. */
  
    T= (S1*D2 - S2*R)/denom; // When it is accurate, this is a lot faster than LU decomposition.
    
    if (T > 1) { T = 1; }
    else if (T < 0) { T = 0; }
    
  }
  
  U = ((T)*R - S2)/D2;
  
  if (U > 1) { 
    
    U = 1; 
    T = ((U)*R + S1)/D1;
    if (T > 1) { T = 1; }
    else if (T < 0) { T = 0; }
    
  } else {

    if (U < 0 ) { 

      U = 0; 
      T = ((U)*R + S1)/D1;
      if (T > 1) { T = 1; }
      else if (T < 0) { T = 0; }
    
    }
  }

  *dist = sqrt(
	       (T*d1.c[0] - U*d2.c[0] - d12.c[0]) * (T*d1.c[0] - U*d2.c[0] - d12.c[0]) + 
	       (T*d1.c[1] - U*d2.c[1] - d12.c[1]) * (T*d1.c[1] - U*d2.c[1] - d12.c[1]) +  
	       (T*d1.c[2] - U*d2.c[2] - d12.c[2]) * (T*d1.c[2] - U*d2.c[2] - d12.c[2]));
  
  *t = T;
  *u = U;
  
}  

/* 
 * Discover whether this box intersects the edge's "ramp".
 *
 */
int box_intersects(box *B,varr *edge,plc_vector verts[3],int last_edge)
{
  int cnt;
  int ramp_status = 0;

#ifdef DEBUG
  if (octrope_debug_level() > 7) {
    if (B->num_edges == 1) {
      printf("Box Intersects? (%d:%d) %3.3f x %3.3f x %3.3f\n",
        edge->component, edge->vertex,
        B->corner[0].c[0], B->corner[0].c[1], B->corner[0].c[2]);
    } else {
      printf(
        "Box Intersects? (%d:%d) [%3.3f,%3.3f] x [%3.3f,%3.3f] x [%3.3f,%3.3f] (%d)\n",
        edge->component, edge->vertex,
        B->corner[0].c[0], B->corner[7].c[0],
        B->corner[0].c[1], B->corner[7].c[1],
        B->corner[0].c[2], B->corner[7].c[2],B->num_edges);
    }
  }
#endif
  ramp_status = ramp_check(edge,verts,B->corner[0],last_edge);
  if (ramp_status & 6) { 
#ifdef DEBUG
    if (octrope_debug_level() > 2) {
      bi_count[0][0]++; 
    }
#endif
    return true; 
  }

  for (cnt = 1; cnt < 8; cnt++) {
    ramp_status |= ramp_check(edge,verts,B->corner[cnt],last_edge);
    if ((ramp_status & 6) || ramp_status == 9) { 
#ifdef DEBUG
      if (octrope_debug_level() > 2) {
        bi_count[ramp_status == 9][cnt]++;
      }
#endif
      return true; 
    }
  }

#ifdef DEBUG
  if (octrope_debug_level() > 2) {
    bi_count[1][8]--;
  }
#endif
  return false;
} /* box_intersects */

/*
 * Remove struts which are longer than <strutmax> from the list.
 *
 */
void cleanup_strutlist(octrope_strut *st, int *strut_count, 
                       const double strutmax)
{
  int st_idx;

  /* First, take any too-long struts off the end of the list, so that the last
   * strut is a valid one. */
  while ((*strut_count > 0) && (st[*strut_count - 1].length > strutmax)) { 
    (*strut_count)--;
  }
  /* Then back through the list.  Whenever a too-long strut is found, remove it
   * and move the one at the end into its slot */
  for (st_idx = *strut_count - 1; st_idx >= 0; st_idx--) {
    if (st[st_idx].length > strutmax) {
      st[st_idx] = st[--*strut_count];
    }
  }
} /* cleanup_strutlist */

#define v_neq(A,B) (A.c[0] != B.c[0] || A.c[1] != B.c[1] || A.c[2] != B.c[2])
#define final_pv(Pl,vertex) (Pl.open && vertex == Pl.nv - 1)
#define dist_up_bound ((cutoff > 0) ? cutoff : *shortest+epsilon)
/*
 * Take a box and either check for struts with those edges contained in it
 * (if it is a leaf node) or send the sub_boxes off to be checked.  Called
 * only when box is known to intersect ramp.
 *
 */
void check_boxes(const plCurve *L, box *B, plc_vector verts[3],
                 int main_edge, int last_edge, double *shortest, 
                 const double cutoff, const double epsilon, 
                 octrope_strut *st,int sl_size,int *strut_count,
                 const int pi_count) 
{
  int box_num;
  plc_vector check_edge[3];
  int edge;
  double s,t,dist;
  octrope_strut *st_ptr;
  int t_last; /* the check_edge is the final edge of an open component */
#ifdef DEBUG
  plc_vector s_end,t_end; /* the two ends of the poca */
  plc_vector tv[2];
#endif
  varr* B_ee; /* B->edges[edge] */
  varr* bo_me; /* by_oct[main_edge] */
  plc_vector cc_vect;
  double cc_sq_dist;
  plc_vector scale_factors;

  bo_me = &by_oct[main_edge];
  if (octree_leaf(B)) { /* Check edges */
    for (edge = 0; edge < B->num_edges; edge++) {
      B_ee = &B->edges[edge];
      check_edge[1] = L->cp[B_ee->component].vt[B_ee->vertex];
      check_edge[2] = L->cp[B_ee->component].vt[B_ee->vertex+1];
      if (B_ee < bo_me &&  /* Only if in by_oct order... */
          ((B_ee->component != bo_me->component) ||  /* ...and not neighbors */
           (abs(B_ee->vertex - bo_me->vertex) >= pi_count && 
           (L->cp[B_ee->component].open || 
            abs(B_ee->vertex - bo_me->vertex) <=  /* around end, if closed */
              L->cp[B_ee->component].nv-pi_count)))) {
        check_edge[0] = L->cp[B_ee->component].vt[B_ee->vertex-1];
#ifdef DEBUG
        if (octrope_debug_level() > 4) {
          if ((bo_me->component > B_ee->component) ||
              ((bo_me->component == B_ee->component) && 
               (bo_me->vertex > B_ee->vertex))) {
            printf("Poca check between %d:%d and %d:%d.\n",
              B_ee->component,B_ee->vertex,bo_me->component,bo_me->vertex);
          } else {
            printf("Poca check between %d:%d and %d:%d.\n",
              bo_me->component,bo_me->vertex,B_ee->component,B_ee->vertex);
          }
        }
#endif
        find_poca(&verts[1],&check_edge[1],&s,&t,&dist);
#ifdef DEBUG
        if (octrope_debug_level() > 2) {
          strut_check_count++;
        }
#endif

#ifdef DEBUG
        if (octrope_debug_level() > 4) {
          plc_M_vweighted(s_end,s,verts[1],verts[2]);
          plc_M_vweighted(t_end,t,check_edge[1],check_edge[2]);
          printf(
            "  from (%3.3f,%3.3f,%3.3f) to (%3.3f,%3.3f,%3.3f) |%3.3f| %d:%d\n",
            s_end.c[0],s_end.c[1],s_end.c[2],t_end.c[0],t_end.c[1],t_end.c[2],
            dist,bo_me->component,bo_me->vertex);
        }
#endif
        t_last = final_pv(L->cp[B_ee->component],B_ee->vertex+1);
        if (dist <= dist_up_bound && 
            (s < 1 || last_edge) && (t < 1 || t_last) &&
            (s > 0 || wedge_check(bo_me,verts,check_edge,t,0,last_edge)) &&
            (t > 0 || wedge_check(B_ee,check_edge,verts,s,1,t_last))) {
#ifdef DEBUG
          if (octrope_debug_level() > 4) {
            if (*shortest < DBL_MAX) {
              printf("-- STRUT FOUND dist:%3.3f cutoff:%3.3f shortest:%3.3f "
                "epsilon:%3.3f --\n",dist,cutoff,*shortest,epsilon);
            } else {
              printf("-- STRUT FOUND dist:%3.3f cutoff:%3.3f shortest: None "
                "epsilon:%3.3f --\n",dist,cutoff,epsilon);
            }
            printf(
              "   : (%3.3f,%3.3f,%3.3f) -- (%3.3f,%3.3f,%3.3f) : %3.3f-%3.3f\n",
              s_end.c[0],s_end.c[1],s_end.c[2],
              t_end.c[0],t_end.c[1],t_end.c[2],s,t);
          }
#endif
          if (sl_size > 0) {
            /* This strut is so short that all of the struts found so far *
             * will be thrown out in transfer_struts.  Throw them out now */
            if ((cutoff == 0) && (epsilon < *shortest - dist)) { 
              *strut_count = 0;
            }
            if (*strut_count < sl_size) {  /* Use next available slot */
              st_ptr = &st[(*strut_count)++];
            } else { /* compact strut array and try again */
              if (dist < *shortest) {
                *shortest = dist;
              }
              cleanup_strutlist(st,strut_count,dist_up_bound);
              if (*strut_count < sl_size) { 
                st_ptr = &st[(*strut_count)++];
              } else { /* No more space, arbitrarily throw out strut 0 */
#ifdef DEBUG
                if (octrope_debug_level() > 3) {
                  printf("Having to throw out strut due to lack of space.\n");
                }
#endif
                st_ptr = st;
              }
            }
            st_ptr->component[0] = bo_me->component;
            st_ptr->lead_vert[0] = bo_me->vertex;
            st_ptr->position[0] = s;
            st_ptr->component[1] = B_ee->component;
            st_ptr->lead_vert[1] = B_ee->vertex;
            st_ptr->position[1] = t;
            st_ptr->length = dist;
#ifdef DEBUG
            if (octrope_debug_level() > 6) {
              octrope_strut_ends(L,st_ptr,tv);
              printf("Poca between (%3.3f,%3.3f,%3.3f) and (%3.3f,%3.3f,%3.3f)",
                tv[0].c[0],tv[0].c[1],tv[0].c[2],
                tv[1].c[0],tv[1].c[1],tv[1].c[2]);
              printf(" %d;%d %d;%d  (%3.3f,%3.3f).\n",
                st_ptr->component[0],st_ptr->lead_vert[0],
                st_ptr->component[1],st_ptr->lead_vert[1],s,t);
            }
#endif
          }
          if (dist < *shortest) {
            *shortest = dist;
          }
        }
      }
    }
  } else { 
    for (box_num = 0; 
         box_num < 8 && (B->sub_boxes[box_num] == NULL ||
                         B->sub_boxes[box_num]->edges < bo_me);
         box_num++) {
      if (B->sub_boxes[box_num] != NULL) {
        /* See if the box is close enough for a strut to be possible */
        if ((cutoff > 0) || (*shortest < DBL_MAX)) {
          scale_factors = cc_vect = B->sub_boxes[box_num]->center;
          plc_M_sub_vect(cc_vect,bo_me->v);
          /* cc_vect points from the center of the main edge to the center of
           * the box to be checked. */
          plc_M_sub_vect(scale_factors,B->sub_boxes[box_num]->corner[0]);
          scale_factors.c[0] += (dist_up_bound) + (bo_me->length/2);
          scale_factors.c[1] += (dist_up_bound) + (bo_me->length/2);
          scale_factors.c[2] += (dist_up_bound) + (bo_me->length/2);
          /* scale_factors now has length equal to the distance from the center
           * of the box to one corner plus sqrt(3)*(either the cutoff or the
           * length of the shortest strut) plus sqrt(3)*(half the length of the
           * edge in question). */
          cc_sq_dist = 
            (cc_vect.c[0]*cc_vect.c[0])/(scale_factors.c[0]*scale_factors.c[0])
          + (cc_vect.c[1]*cc_vect.c[1])/(scale_factors.c[1]*scale_factors.c[1])
          + (cc_vect.c[2]*cc_vect.c[2])/(scale_factors.c[2]*scale_factors.c[2]);
        } else {
          cc_sq_dist = 0; /* Guaranteed less than 3 ;-) */
        }
        if (cc_sq_dist <= 3 &&
            box_intersects(B->sub_boxes[box_num],bo_me,verts,last_edge)) {
          check_boxes(L,B->sub_boxes[box_num],verts,main_edge,last_edge,
                      shortest,cutoff,epsilon,st,sl_size,strut_count,pi_count); 
        }
      }
    }
  }
} /* check_boxes */

/* Return the maximum turning angle between any two edges of a link */
double max_turning_angle(const plCurve *L)
{
  int vert;
  int cmp;
  plc_strand *comp;
  plc_vector e1,e2;
  double max_angle = 0.0;
  double angle;
  
  for (cmp = 0; cmp < L->nc; cmp++) {
    comp = &L->cp[cmp];
    if (!comp->open) { 
      /* Closed loop, check first and last edges against "virtual edge" */
      e2 = comp->vt[comp->nv-1];   
      e1 = comp->vt[0];
      plc_M_sub_vect(e2,e1);
      plc_M_sub_vect(e1,comp->vt[1]);
      angle = 
        acos(plc_M_dot(e1,e2)/sqrt(plc_M_dot(e1,e1)*plc_M_dot(e2,e2)));
      max_angle = max(angle,max_angle);
      e1 = comp->vt[comp->nv-2];
      plc_M_sub_vect(e1,comp->vt[comp->nv-1]);
      angle = 
        acos(plc_M_dot(e1,e2)/sqrt(plc_M_dot(e1,e1)*plc_M_dot(e2,e2)));
      max_angle = max(angle,max_angle);
    }
    for (vert = 1; vert < comp->nv-1; vert++) {
      e1 = comp->vt[vert];
      e2 = comp->vt[vert+1];
      plc_M_sub_vect(e2,e1);
      plc_M_sub_vect(e1,comp->vt[vert-1]);
      angle = 
        acos(plc_M_dot(e1,e2)/sqrt(plc_M_dot(e1,e1)*plc_M_dot(e2,e2)));
      max_angle = max(angle,max_angle);
    }
  }
  return max_angle;
} /* max_turning_angle */

int octrope_struts(plCurve *L, 
                   const double cutoff, const double epsilon, 
                   octrope_strut *strutlist, int sl_size,
                   double *shortest, void *mem, const int memsize) 
{
  int cmp,vert;
  plc_vector verts[3];
  int edge;
  box *B;
  double short_strut;
  int strut_count;
  double max_turn;
  int pi_count;

#ifdef DEBUG
  int edges_found;
#endif

  octrope_error_num = octrope_error_str[0] = 0;
#ifdef DEBUG
  if (octrope_debug_level() > 8) {
    printf("octrope_struts:called.\n");
  }
  if (octrope_debug_level() > 2) {
    strut_check_count = 0;
    rc_count[0] = rc_count[1] = rc_count[2] = rc_count[3] = rc_count[4] = 0;
    wc_count[0] = wc_count[1] = wc_count[2] = wc_count[3] = wc_count[4] = 0;
    bi_count[0][0] = bi_count[0][1] = bi_count[0][2] = bi_count[0][3] = 0;
    bi_count[0][4] = bi_count[0][5] = bi_count[0][6] = bi_count[0][7] = 0;
    bi_count[1][0] = bi_count[1][1] = bi_count[1][2] = bi_count[1][3] = 0;
    bi_count[1][4] = bi_count[1][5] = bi_count[1][6] = bi_count[1][7] = 0;
    bi_count[0][8] = bi_count[1][8] = 0;
  }
#endif
  rampcheck_edge = NULL;
  wedgecheck_edge[0] = wedgecheck_edge[1] = 0;
  /* Sanity check */
  if (L == NULL) {
    sprintf(octrope_error_str,"Link pointer is NULL in octrope_struts().\n");
    octrope_error_num = 13;
    return -1;
  }
  if (epsilon < 0) {
    sprintf(octrope_error_str,"Epsilon less than zero in octrope_struts().\n");
    octrope_error_num = 14;
    return -1;
  }
  if (cutoff < 0) {
    sprintf(octrope_error_str,"Cutoff less than zero in octrope_struts().\n");
    octrope_error_num = 16;
    return -1;
  }
  if (strutlist == NULL && sl_size != 0) {
    sprintf(octrope_error_str,
      "Strut list pointer is NULL with sl_size nonzero in octrope_struts().\n");
    octrope_error_num = 15;
    return -1;
  }

  /* Build the Octree */
  plc_fix_wrap(L);
  B = build_octree(L,sl_size,mem,memsize);
  max_turn = max_turning_angle(L);

  /* We now set pi_count, which is a lower bound on the number of edges 
     required to turn by pi radians. In general, we floor for safety when
     pi_count is large. However, for curves with very large turning angle,
     we can get pi_count = 1 from this computation. Hilarity ensues. So 
     we make sure that pi_count >= 2 afterwards. */

  pi_count = floor(acos(-1)/max_turn); 
  if (pi_count < 2) { pi_count = 2; }  

#ifdef DEBUG
  if (octrope_debug_level() > 3) {
    printf("Maximum turning angle: %13.13f\n",max_turn);
  }
  if (octrope_debug_level() > 6) {
    edges_found = 0;
    print_octree(B,"",&edges_found);
    printf("In traversing the tree we found %d edges.\n",edges_found);
  }
#endif

  /* Set shortest large so as to catch the first one */
  short_strut = DBL_MAX;

  /* Mark the strut list as empty */
  strut_count = 0;

  /* Don't check first edge as struts can only be *
   * from an edge to one before it.               */
  {
    varr *shared_by_oct = by_oct;
    int shared_num_edges = num_edges;
    int use_omp = 0;
#ifdef _OPENMP
    use_omp = (omp_get_max_threads() > 1 &&
               shared_num_edges >= OCTROPE_OMP_MIN_EDGES);
#endif

    if (!use_omp) {
      for (edge = shared_num_edges - 1; edge > 0; edge--) {
        cmp = shared_by_oct[edge].component;
        vert = shared_by_oct[edge].vertex;
        verts[0] = L->cp[cmp].vt[vert-1];
        verts[1] = L->cp[cmp].vt[vert];
        verts[2] = L->cp[cmp].vt[vert+1];

        check_boxes(L,B,verts,edge,final_pv(L->cp[cmp],vert+1),
          &short_strut,cutoff,epsilon,strutlist,sl_size,&strut_count,pi_count);
      }
    } else {
#ifdef _OPENMP
      int nthreads = omp_get_max_threads();
      int tid;
      octrope_strut *local_lists = NULL;
      int *local_counts = NULL;
      double *local_short = NULL;
      int local_sl = (sl_size > 0) ? sl_size : 0;

      local_short = (double *)malloc((size_t)nthreads * sizeof(double));
      local_counts = (int *)calloc((size_t)nthreads, sizeof(int));
      if (local_short == NULL || local_counts == NULL) {
        free(local_short);
        free(local_counts);
        sprintf(octrope_error_str,
          "Unable to allocate OpenMP scratch in octrope_struts().\n");
        octrope_error_num = 17;
        free_octmem();
        return -1;
      }
      if (local_sl > 0) {
        local_lists = (octrope_strut *)malloc(
          (size_t)nthreads * (size_t)local_sl * sizeof(octrope_strut));
        if (local_lists == NULL) {
          free(local_short);
          free(local_counts);
          sprintf(octrope_error_str,
            "Unable to allocate OpenMP strut buffers in octrope_struts().\n");
          octrope_error_num = 17;
          free_octmem();
          return -1;
        }
      }
      for (tid = 0; tid < nthreads; tid++) {
        local_short[tid] = DBL_MAX;
      }

#pragma omp parallel default(none) \
  shared(L,B,shared_by_oct,shared_num_edges,cutoff,epsilon,pi_count, \
         local_lists,local_counts,local_short,local_sl,short_strut) \
  private(edge,cmp,vert,verts)
      {
        int mytid = omp_get_thread_num();
        double my_short = DBL_MAX;
        int my_count = 0;
        octrope_strut *my_list = (local_sl > 0)
          ? &local_lists[mytid * local_sl] : NULL;

        /* Workers start with empty TLS; inherit the master's tree. */
        by_oct = shared_by_oct;
        num_edges = shared_num_edges;
        rampcheck_edge = NULL;
        wedgecheck_edge[0] = wedgecheck_edge[1] = NULL;

#pragma omp for schedule(dynamic, 8) nowait
        for (edge = shared_num_edges - 1; edge > 0; edge--) {
          double bound_short;
          double shared_short;

          cmp = shared_by_oct[edge].component;
          vert = shared_by_oct[edge].vertex;
          verts[0] = L->cp[cmp].vt[vert-1];
          verts[1] = L->cp[cmp].vt[vert];
          verts[2] = L->cp[cmp].vt[vert+1];

          /* Share improved bounds across threads for better pruning. */
#pragma omp atomic read
          shared_short = short_strut;
          bound_short = my_short;
          if (shared_short < bound_short) {
            bound_short = shared_short;
          }

          check_boxes(L,B,verts,edge,final_pv(L->cp[cmp],vert+1),
            &bound_short,cutoff,epsilon,my_list,local_sl,&my_count,pi_count);

          if (bound_short < my_short) {
            my_short = bound_short;
#pragma omp critical (octrope_shortest)
            {
              if (my_short < short_strut) {
                short_strut = my_short;
              }
            }
          }
        }
        local_counts[mytid] = my_count;
        local_short[mytid] = my_short;
      }

      for (tid = 0; tid < nthreads; tid++) {
        if (local_short[tid] < short_strut) {
          short_strut = local_short[tid];
        }
      }

      if (local_sl > 0 && strutlist != NULL) {
        double bound = (cutoff > 0) ? cutoff : short_strut + epsilon;
        for (tid = 0; tid < nthreads; tid++) {
          int i;
          for (i = 0; i < local_counts[tid]; i++) {
            octrope_strut *s = &local_lists[tid * local_sl + i];
            if (s->length > bound) {
              continue;
            }
            if (strut_count < sl_size) {
              strutlist[strut_count++] = *s;
            } else {
              cleanup_strutlist(strutlist,&strut_count,bound);
              if (strut_count < sl_size) {
                strutlist[strut_count++] = *s;
              }
            }
          }
        }
      }

      free(local_lists);
      free(local_counts);
      free(local_short);

      /* Restore master TLS tree pointers (freed via free_octmem below). */
      by_oct = shared_by_oct;
      num_edges = shared_num_edges;
#endif /* _OPENMP */
    }
  }

  if (sl_size > 0 && strut_count > 0) {
    cleanup_strutlist(strutlist,&strut_count,
      (cutoff > 0) ? cutoff : short_strut + epsilon);
  }

#ifdef DEBUG
  if (octrope_debug_level() > 2) {
    printf("Did %d strut checks.\n",strut_check_count);
    printf("Rampcheck counts: %d %d %d %d %d\n",rc_count[0],rc_count[1],
      rc_count[2],rc_count[3],rc_count[4]);
    printf("Wedgecheck counts: %d %d %d %d %d\n",wc_count[0],wc_count[1],
      wc_count[2],wc_count[3],wc_count[4]);
    printf("BoxIntersects counts: %d %d %d %d %d %d %d %d %d\n",
      bi_count[0][0], bi_count[0][1], bi_count[0][2], bi_count[0][3],
      bi_count[0][4], bi_count[0][5], bi_count[0][6], bi_count[0][7],
      bi_count[0][8]);
    printf("                      %d %d %d %d %d %d %d %d %d\n",
      bi_count[1][0], bi_count[1][1], bi_count[1][2], bi_count[1][3],
      bi_count[1][4], bi_count[1][5], bi_count[1][6], bi_count[1][7],
      bi_count[1][8]);
  }
#endif

  free_octmem();
  if (shortest != NULL) { *shortest = short_strut; }
  return strut_count;
} /* octrope_struts */

inline double octrope_curvelength(const plCurve *L)
{
  double tot_length;
  int cmp,nv,vert;
  plc_strand *cp;
  plc_vector temp_vect;

  octrope_error_num = octrope_error_str[0] = 0;
  if (L == NULL) {
    sprintf(octrope_error_str,"Link pointer NULL in octrope_curvelength().\n");
    octrope_error_num = 28;
    return -1;
  }

  tot_length = 0;
  for (cmp = 0; cmp < L->nc; cmp++) {
    cp = &L->cp[cmp];
    nv = (cp->open) ? cp->nv-1 : cp->nv;
    for (vert = 0; vert < nv; vert++) {
      temp_vect = cp->vt[vert+1];
      plc_M_sub_vect(temp_vect,cp->vt[vert]);
      tot_length += plc_M_norm(temp_vect);
    }
  }
 
  return tot_length;
}

/*
 * Find the length of the shortest strut.
 *
 */
inline double octrope_poca(plCurve *L, void *mem, const int memsize) 
{
  double min_strut_length;

  /* Don't set octrope_error_num, trust octrope_minrad to take care of that */
  octrope_struts(L,0,0,NULL,0,&min_strut_length,mem,memsize);
  return min_strut_length;
}

/* 
 * Find the value of MinRad.
 *
 */
inline double octrope_minradval(plCurve *L)
{
  /* Don't set octrope_error_num, trust octrope_minrad to take care of that */
  return octrope_minrad(L,0,0,NULL,0,NULL);
}

/* 
 * Find the thickness, where thickness is given by
 * 
 * thickness = min(min_strutlength/2,MinRad/lambda)
 *
 * If lambda is 0, the denominator is just 
 *   min_strutlength/2
 *
 */
inline double octrope_thickness(plCurve *L, void *mem, 
                                const int memsize, const double lambda) 
{
  double thickness;

  /* Don't set octrope_error_num, trust octrope_minrad to take care of that */
  octrope(L,                        NULL    /* ropelength       */,
    &thickness,                     NULL    /* curve_len_ret    */,
    NULL    /* min_rad_ret      */, NULL    /* min_strut_ret    */,
    0       /* mr_cutoff        */, 0       /* mr_epsilon       */, 
    NULL    /* min_rad_locs     */, 0       /* mr_size          */, 
    NULL    /* num_min_rad_locs */, 0       /* strut_cutoff     */,
    0       /* strut_epsilon    */, NULL    /* strutlist        */,
    0       /* sl_size          */, NULL    /* num_strut_ret    */,
    mem, memsize, lambda);
  return thickness;
}

/* 
 * Find the ropelength, where ropelength is given by
 * 
 *                        total length
 * ropelength = ------------------------------------
 *              min(min_strutlength/2,MinRad/lambda)
 *
 * If lambda is 0, the denominator is just 
 *   min_strutlength/2
 *
 */
inline double octrope_ropelength(plCurve *L, void *mem, 
                                 const int memsize, const double lambda) 
{
  double ropelength;

  /* Don't set octrope_error_num, trust octrope_minrad to take care of that */
  octrope(L,                        &ropelength,
    NULL    /* thickness        */, NULL    /* curve_len_ret    */,
    NULL    /* min_rad_ret      */, NULL    /* min_strut_ret    */,
    0       /* mr_cutoff        */, 0       /* mr_epsilon       */, 
    NULL    /* min_rad_locs     */, 0       /* mr_size          */, 
    NULL    /* num_min_rad_locs */, 0       /* strut_cutoff     */,
    0       /* strut_epsilon    */, NULL    /* strutlist        */,
    0       /* sl_size          */, NULL    /* num_strut_ret    */,
    mem, memsize, lambda);
  return ropelength;
}

void octrope(plCurve *L,

             double *ropelength,
             double *thickness_ret,

             double *curve_len_ret,

             double *min_rad_ret, 
             double *min_strut_ret, 
     
             const double mr_cutoff,
             const double mr_epsilon, 
             octrope_mrloc *min_rad_locs, 
             const int mr_size,
             int *num_min_rad_locs, 

             const double strut_cutoff,
             const double strut_epsilon, 
             octrope_strut *strutlist, 
             const int sl_size,
             int *num_strut_ret,
     
             void *mem, const int memsize,

             const double lambda)
{
  double min_strutlength;
  double min_rad = 0.0;
  double tot_length = 0.0;
  double thickness = 0.0;
  int num_struts;

  /* Sanity check */

  octrope_error_num = octrope_error_str[0] = 0;
  if (L == NULL) {
    sprintf(octrope_error_str,"Link pointer is NULL in octrope().\n");
    octrope_error_num = 11;
    return;
  }
  if (lambda < 0) { /* lambda==0 means minrad doesn't figure into thickness */
    sprintf(octrope_error_str,"Lambda is negative in octrope().\n");
    octrope_error_num = 12;
    return;
  }

  /* 
   * We now gather the min_rad location information. 
   */

  if ((ropelength  != NULL || thickness_ret    != NULL || 
       min_rad_ret != NULL || num_min_rad_locs != NULL)
      && lambda > 0) 

    /* If lambda == 0, we never record any minradlocs. */

    {
      min_rad = octrope_minrad(L,mr_cutoff,mr_epsilon,min_rad_locs,mr_size,
			       num_min_rad_locs);
#ifdef DEBUG
      if (octrope_debug_level() > 0) {
	printf("Minrad: %3.3f\n",min_rad);
      }
#endif
      if (min_rad_ret != NULL) { *min_rad_ret = min_rad; }
      
    } 

  /* We set special values for min_rad_ret and num_min_rad_locs if we are in the
     zero stiffness case */

  if (lambda == 0) { 

    if (min_rad_ret != NULL) { *min_rad_ret = DBL_MAX; }
    if (num_min_rad_locs != NULL) { *num_min_rad_locs = 0; }

  }
  
  /* Now we gather the strut data. */

  if (ropelength    != NULL || thickness_ret != NULL || 
      min_strut_ret != NULL || num_strut_ret != NULL) { 
    num_struts = octrope_struts(L,strut_cutoff,strut_epsilon,strutlist,
                                sl_size,&min_strutlength,mem,memsize);
 
    if (min_strut_ret != NULL) { *min_strut_ret = min_strutlength; }
    if (num_strut_ret != NULL) { *num_strut_ret = num_struts; }
#ifdef DEBUG
    if (octrope_debug_level() > 0) {
      printf("Mininum strutlength: %3.3f\n",min_strutlength);
    }
#endif

    if (num_struts == 0) { /* The cutoff and epsilon mean that we didn't record min_strutlength. */
      /* But we'll use it to calculate ropelength below, which is BAD. We now recompute with no  */
      /* cutoff requirements. */

      octrope_struts(L,DBL_MAX,0,NULL,0,&min_strutlength,NULL,0);

    }
      
  }
  
  /* Now we need the thickness, so we find the min of twice minrad/lambda 
     and min_strutlength */
  if (ropelength != NULL || thickness_ret != NULL) {
    if (lambda == 0) {
      thickness = min_strutlength/2.0;
    } else {
      thickness = min(min_strutlength/2.0,min_rad/lambda);
    } 
    if (thickness_ret != NULL) { *thickness_ret = thickness; }
  }

  /* We want to know the actual length of the curve */
  if (ropelength != NULL || curve_len_ret != NULL) {
    tot_length = octrope_curvelength(L);
    if (curve_len_ret != NULL) { *curve_len_ret = tot_length; }
  }
  
  /* And with everything in hand we can do the division and find ropelength. */
  if (ropelength != NULL) {
#ifdef DEBUG
    if (octrope_debug_level() > 0) {
      printf("Total length: %3.3f,  Ropelength: %3.3f\n",tot_length,
        tot_length/thickness);
    }
#endif
    if (thickness == 0) {
      (*ropelength) = DBL_MAX;
    } else {
      (*ropelength) = tot_length/thickness;
    }
  }

  return;
}

/*
 * Translate the strut-format information into two vectors which are the
 * ends of the strut in question.  S points to the strut, se should
 * be the address of a vector[2].
 *
 */
void octrope_strut_ends(const plCurve *L, const octrope_strut *S, 
                              plc_vector se[2]) 
{  
  int i;
  plc_strand *cp;
  plc_vector *vert1,*vert2;

  octrope_error_num = octrope_error_str[0] = 0;
  if (S == NULL) { 
    sprintf(octrope_error_str,"Strut pointer NULL in octrope_strut_ends().\n");
    octrope_error_num = 17;
    return;
  } else if (se == NULL) { 
    sprintf(octrope_error_str,"Vector pointer NULL in octrope_strut_ends().\n");
    octrope_error_num = 18;
    return;
  }

  for (i=0; i < 2; i++) {
    cp = &L->cp[S->component[i]];
    vert1 = &cp->vt[S->lead_vert[i]];
    vert2 = &cp->vt[S->lead_vert[i]+1];
    plc_M_vweighted(se[i],S->position[i],*vert1,*vert2);
  }
} /* octrope_strut_ends */

/* Procedure positions the file pointer on next non-whitespace character,   *
 * returning false if EOF happens first. We skip anything between a # and a *
 * newline.                                                                 */
static int skip_whitespace_and_comments(FILE *infile)
{
  int thischar,commentflag = {false};

  /* First, we check to make sure that infile looks legit. */
  if (infile == NULL) {
    octrope_error_num = 61;
    sprintf(octrope_error_str,"skip_whitespace_and_comments: infile is a null pointer.\n");
    return -1;
  }
  
  /* Now we start to work. */
  for(;;) {
    thischar = fgetc(infile);

    if (thischar == EOF) {  /* Reached end of file before a non-space, non-comment */
      return 0;
    } else if (thischar == '#') { /* Started a comment. */
      commentflag = true;
    } else if (thischar == '\n' && commentflag) { /* End a comment. */
      commentflag = false;
    } else if (!isspace(thischar) && !commentflag) { /* Found a hit! */
      ungetc(thischar,infile);
      return 1;
    } /* It must have been a space or a non-space in a comment. */
  }
}

/* Procedure scans for nfloats floating point (or double) numbers, ignoring  *
 * whitespace and comments between them. We expect the variable length       *
 * arguments to contain a collection of pointers to doubles. If not, there's *
 * trouble.                                                                  */
static int scandoubles(FILE *infile,int ndoubles, ...)
{
  int nconverted = 0,i;
  va_list ap;
  double *thisdouble;

  /* First, we check for overall sanity. */

  if (infile == NULL) {
    octrope_error_num = 71;
    sprintf(octrope_error_str,"scandoubles: infile is a null pointer.\n");
    return -1;
  }

  if (ndoubles < 1) {
    octrope_error_num = 72;
    sprintf(octrope_error_str,"scandoubles: ndoubles (%d) is less than one.\n",ndoubles);
    return -1;
  }

  va_start(ap,ndoubles);

  /* Now we're ready to work. */

  for (i=0;i<ndoubles;i++) {    /* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */
      return nconverted;
    }

    thisdouble = va_arg(ap,double *);
    if (fscanf(infile,"%lf",thisdouble) != 1) { /* We couldn't scan. */
      return nconverted;        /* So give up here */
    } else {                    /* Else record our victory */
      nconverted++;
    }
  }
  va_end(ap);

  return nconverted;
}

/* Procedure scans for nints integers, ignoring whitespace and     *
 * comments between them. We expect the variable length arguments  *
 * to contain a collection of pointers to ints. If not,            *
 * there's trouble.                                                */
static int scanints(FILE *infile,int nints, ...)
{
  int nconverted = 0,i;
  va_list ap;
  int *thisint;

  /* First, we check for overall sanity. */

  if (infile == NULL) {
    octrope_error_num = 73;
    sprintf(octrope_error_str,"scanints: infile is a null pointer.\n");
    return -1;
  }

  if (nints < 1) {
    octrope_error_num = 74;
    sprintf(octrope_error_str,"scanints: nints (%d) is less than one.\n",nints);
    return -1;
  }

  va_start(ap,nints);

  /* Now we're ready to work. */
  for (i=0;i<nints;i++) {       /* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */
      return nconverted;
    }
    thisint = va_arg(ap,int *);

    if (fscanf(infile,"%d",thisint) != 1) {     /* We couldn't scan. */
      return nconverted;        /* So give up here */
    } else {                    /* Else record our victory */
      nconverted++;
    }
  }
  va_end(ap);
  return nconverted;
}
/* Reads a disk file of stored strut and mrloc records for a given link L */
/* During this process, we go through a fair bit of error checking to make */
/* sure that the struts make sense for the link L. */

void octrope_strutfile_read(plCurve *L,
                            int *n_struts,octrope_strut **strutlist,
                            int *n_mrloc,octrope_mrloc **mrlist,
                            FILE *file)

{
  int i;
  octrope_strut ts;
  octrope_mrloc tmr;
  
  /* First, we check for the 'STRUTS' keyword. */
  
  octrope_error_num = octrope_error_str[0] = 0;
  if (fscanf(file," STRUTS ") == EOF) {
    
    octrope_error_num = 20;
    sprintf(octrope_error_str,
      "octrope_strutfile_read: Couldn't find STRUTS keyword.");    
    return;
  }

  /* Now we read two integers giving the number of struts and mrlocs. */

  if (scanints(file,2,n_struts,n_mrloc) != 2) {
  
    octrope_error_num = 30;
    sprintf(octrope_error_str,
            "octrope_strutfile_read: Couldn't parse <struts> <mrlocs> line");
    return;
  }

  if (*n_struts < 0 || *n_mrloc < 0 || *n_struts > 1E6 || *n_mrloc > 1E6) {

    octrope_error_num = 40;
    sprintf(octrope_error_str,
      "octrope_strutfile_read: number of struts (%d) or minrad locations (%d)"
      " seems unlikely.\n"
      "                        Only 1,000,000 struts are allowed.\n",
      *n_struts,*n_mrloc);
    return;
  }

  /* Now we allocate space for the struts and mrlocs. */

  (*strutlist) = (octrope_strut *)calloc(*n_struts,sizeof(octrope_strut));
  (*mrlist)    = (octrope_mrloc *)calloc(*n_mrloc,sizeof(octrope_mrloc));

  if (strutlist == NULL || mrlist == NULL) {

    octrope_error_num = 50;
    sprintf(octrope_error_str,
      "octrope_strutfile_read: Couldn't allocate space for "
      "%d struts and %d mrloc records.\n",*n_struts,*n_mrloc);
    return;

  }

  /* Now we start attempting to read lines of struts. */

  for(i=0;i<*n_struts;i++) {

    if (scanints(file, 4, 
                 &((*strutlist)[i].component[0]), 
                 &((*strutlist)[i].component[1]), 
                 &((*strutlist)[i].lead_vert[0]),
                 &((*strutlist)[i].lead_vert[1])) != 4) {

      octrope_error_num = 64;
      sprintf(octrope_error_str,
        "octrope_strutfile_read: Couldn't parse strut %d "
        " component and vertex information.\n",i);
      return;

    }

    ts = (*strutlist)[i];

    if (ts.component[0] >= L->nc || ts.component[1] >= L->nc) {

      octrope_error_num = 66;
      sprintf(octrope_error_str,
        "octrope_strutfile_read: Component information (%d,%d) "
        " for strut %d doesn't match L, which has %d components.\n",
        ts.component[0],ts.component[1],i,L->nc);
      return;

    }

    if (ts.lead_vert[0] >= L->cp[ts.component[0]].nv ||
        ts.lead_vert[1] >= L->cp[ts.component[1]].nv) {

      octrope_error_num = 67;
      sprintf(octrope_error_str,
        "octrope_strutfile_read: Vertices (%d,%d) in strut %d "
        "don't match component (%d,%d) of L, which have (%d,%d) verts.\n",
        ts.lead_vert[0],ts.lead_vert[1], i,
        ts.component[0],ts.component[1],
        L->cp[ts.component[0]].nv, L->cp[ts.component[1]].nv);
      return;
      
    }
    
    if (scandoubles(file,4,
      &((*strutlist)[i].position[0]),&((*strutlist)[i].position[1]),
      &((*strutlist)[i].length),&((*strutlist)[i].compression)) != 4) {

      octrope_error_num = 65;
      sprintf(octrope_error_str,
        "octrope_strutfile_read: Couldn't parse strut %d "
        " position, length, and compression information.\n",i);
      return;

    }
    
  }

  /* We now try to read the minrad information */

  for(i=0;i<*n_mrloc;i++) {

    if (scanints(file,2,&((*mrlist)[i].component),&((*mrlist)[i].vert)) != 2) {

      octrope_error_num = 92;
      sprintf(octrope_error_str,
        "octrope_strutfile_read: Couldn't parse mrloc %d "
        " component and vertex information.\n",i);
      return;

    } 
    
    tmr = (*mrlist)[i];

    if (tmr.component >= L->nc || tmr.component < 0) {

      octrope_error_num = 93;
      sprintf(octrope_error_str,
        "octrope_strutfile_read: Mrloc %d "
        " component information %d doesn't match L, which has %d components.\n"
              ,i,tmr.component,L->nc);
      return;

    }

    if (tmr.vert >= L->cp[tmr.component].nv || tmr.vert < 0) {

      octrope_error_num = 94;
      sprintf(octrope_error_str,"octrope_strutfile_read: Mrloc %d "
              " vertex information %d doesn't match component %d of L, which "
              " has %d vertices.\n",i,tmr.vert,tmr.component,
              L->cp[tmr.component].nv);
      return;

    }

    if (scandoubles(file,2,&((*mrlist)[i].mr),&((*mrlist)[i].compression)) != 2) {

      octrope_error_num = 95;
      sprintf(octrope_error_str,
        "octrope_strutfile_read: Couldn't parse mrloc %d "
        " mr and compression information.\n",i);
      return;

    } 
                
  }

}


/*
 *  Writes a list of struts to the disk file in file. Aside from basic
 *  NULL-pointer checking, there's not much error checking or correction 
 *  that's possible here. We won't sanity-check the strutlist or mrlist here,
 *  and so we do not need the link that all this (hopefully) refers to.
 */

void octrope_strutfile_write(int n_struts,octrope_strut *strutlist,
                             int n_mrloc,octrope_mrloc *mrlist,
                             FILE *file)

{
  int i;

  octrope_error_num = octrope_error_str[0] = 0;
  if (file == NULL) {
    
    octrope_error_num = 36;
    sprintf(octrope_error_str,"octrope_strutfile_write: Passed NULL pointer "
            "for output file.\n");
    return;
    
  }
  
  if (n_struts < 0 || (n_struts != 0 && strutlist==NULL)) {
    
    octrope_error_num = 37;
    sprintf(octrope_error_str,"octrope_strutfile_write: n_struts (%d) "
            "or strutlist ptr (%p) appear corrupt.\n", n_struts,strutlist);
    return;
    
  }
  
  if (n_mrloc < 0 || (n_mrloc != 0 && mrlist == NULL)) {
    
    octrope_error_num = 38;
    sprintf(octrope_error_str,"octrope_strutfile_write: n_mrloc (%d) "
            " mrlist ptr (%p) appear corrupt.\n",n_mrloc,mrlist);
    return;
    
  }

  /* We have now checked everything as best we can. It's time to write 
     things to disk. */

  fprintf(file,"STRUTS \n");
  fprintf(file,"%d %d \n",n_struts,n_mrloc);
  
  for(i=0;i<n_struts;i++) {

    fprintf(file,"%d %d %d %d %16g %16g %16g %16g\n",
            strutlist[i].component[0],
            strutlist[i].component[1],
            strutlist[i].lead_vert[0],
            strutlist[i].lead_vert[1],
            strutlist[i].position[0],
            strutlist[i].position[1],
            strutlist[i].length,
            strutlist[i].compression);

  }

  for(i=0;i<n_mrloc;i++) {

    fprintf(file,"%d %d %16g %16g \n",
            mrlist[i].component,
            mrlist[i].vert,
            mrlist[i].mr,
            mrlist[i].compression);

  }
            
}

/* Print or store library version number. */

void octrope_version( char *version, size_t strlen) {
  if (version == NULL) {
    printf("Octrope Version: %s\n",PACKAGE_VERSION);
  } else {
    (void)snprintf(version,strlen,PACKAGE_VERSION);
  }
}
