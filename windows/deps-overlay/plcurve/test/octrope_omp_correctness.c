/*
 * Compare octrope shortest-strut / thickness with 1 vs N OpenMP threads.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "plCurve.h"
#include "octrope.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static plCurve *make_trefoil(int nv)
{
  int i, cc = 0;
  bool open = false;
  double theta, t_step;
  double pi = 3.14159265358979;
  plCurve *L = plc_new(1, &nv, &open, &cc);
  if (L == NULL) {
    return NULL;
  }
  t_step = 4.0 * pi / (double)nv;
  for (i = 0, theta = t_step / 2.0; i < nv; i++, theta += t_step) {
    L->cp[0].vt[i].c[1] = (1 + 0.66 * cos(1.5 * theta)) * cos(theta);
    L->cp[0].vt[i].c[0] = (1 + 0.66 * cos(1.5 * theta)) * sin(theta);
    L->cp[0].vt[i].c[2] = 0.66 * sin(1.5 * theta);
  }
  return L;
}

static plCurve *make_hopf(int total_verts)
{
  int i, nv[2], ccarray[2] = {0, 0};
  bool open[2] = {false, false};
  double theta, t_step;
  double pi = 3.14159265358979;
  int vpc = total_verts / 2;
  plCurve *L;

  nv[0] = nv[1] = vpc;
  L = plc_new(2, nv, open, ccarray);
  if (L == NULL) {
    return NULL;
  }
  t_step = 2.0 * pi / (double)vpc;
  for (i = 0, theta = t_step / 2.0; i < vpc; i++, theta += t_step) {
    L->cp[0].vt[i].c[1] = (1 / cos(t_step / 2.0)) * cos(theta);
    L->cp[0].vt[i].c[0] = (1 / cos(t_step / 2.0)) * sin(theta);
    L->cp[0].vt[i].c[2] = 0;
    L->cp[1].vt[i].c[0] = 0;
    L->cp[1].vt[i].c[1] = -(1 / cos(t_step / 2.0)) * cos(theta) + 1.0;
    L->cp[1].vt[i].c[2] = (1 / cos(t_step / 2.0)) * sin(theta);
  }
  return L;
}

static int nearly_equal(double a, double b, double rel)
{
  double scale = fabs(a) + fabs(b) + 1.0;
  return fabs(a - b) <= rel * scale;
}

static int check_curve(const char *name, plCurve *L, int nthreads)
{
  double short1, shortN, thi1, thiN;
  int memsize = octrope_est_mem(plc_num_edges(L));
  void *mem = malloc((size_t)memsize);
  int nstrut1, nstrutN;
  int sl = 6 * plc_num_edges(L);
  octrope_strut *slist = (octrope_strut *)malloc((size_t)sl * sizeof(octrope_strut));
  int ok = 1;

  if (mem == NULL || slist == NULL) {
    fprintf(stderr, "%s: malloc failed\n", name);
    free(mem);
    free(slist);
    return 0;
  }

  octrope_set_threads(1);
  nstrut1 = octrope_struts(L, 0, 1e-6, slist, sl, &short1, mem, memsize);
  thi1 = octrope_thickness(L, mem, memsize, 1.0);
  if (octrope_error_num != 0 || nstrut1 < 0) {
    fprintf(stderr, "%s: serial octrope failed: %s\n", name, octrope_error_str);
    ok = 0;
    goto done;
  }

  octrope_set_threads(nthreads);
  nstrutN = octrope_struts(L, 0, 1e-6, slist, sl, &shortN, mem, memsize);
  thiN = octrope_thickness(L, mem, memsize, 1.0);
  if (octrope_error_num != 0 || nstrutN < 0) {
    fprintf(stderr, "%s: omp octrope failed: %s\n", name, octrope_error_str);
    ok = 0;
    goto done;
  }

  if (!nearly_equal(short1, shortN, 1e-12)) {
    fprintf(stderr, "%s: shortest mismatch serial=%g omp=%g (threads=%d)\n",
            name, short1, shortN, nthreads);
    ok = 0;
  }
  if (!nearly_equal(thi1, thiN, 1e-12)) {
    fprintf(stderr, "%s: thickness mismatch serial=%g omp=%g (threads=%d)\n",
            name, thi1, thiN, nthreads);
    ok = 0;
  }

  if (ok) {
    printf("PASS %s edges=%d short=%g thi=%g threads=%d (struts %d vs %d)\n",
           name, plc_num_edges(L), short1, thi1, nthreads, nstrut1, nstrutN);
  }

done:
  free(mem);
  free(slist);
  return ok;
}

int main(void)
{
  int nthreads = 4;
  int fails = 0;
  plCurve *L;
#ifdef _OPENMP
  if (omp_get_max_threads() < 2) {
    nthreads = 2;
  }
  printf("OpenMP available, testing with %d threads (max=%d)\n",
         nthreads, omp_get_max_threads());
#else
  printf("OpenMP not compiled in; verifying serial path only\n");
  nthreads = 1;
#endif

  L = make_trefoil(128);
  if (L == NULL || !check_curve("trefoil128", L, nthreads)) {
    fails++;
  }
  plc_free(L);

  L = make_trefoil(300);
  if (L == NULL || !check_curve("trefoil300", L, nthreads)) {
    fails++;
  }
  plc_free(L);

  L = make_hopf(200);
  if (L == NULL || !check_curve("hopf200", L, nthreads)) {
    fails++;
  }
  plc_free(L);

  if (fails) {
    fprintf(stderr, "octrope_omp_correctness: %d failure(s)\n", fails);
    return 1;
  }
  printf("octrope_omp_correctness: all passed\n");
  return 0;
}
