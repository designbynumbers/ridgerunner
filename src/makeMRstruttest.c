
/* makeMRstruttest.c: Generates a file to be used to test the MR strut architecture.
   Distributed, but not installed. */

#include <portability.h>
#include <plCurve.h>

int main() {

  double PI = 3.14159265359;
  double s1,s2,c1,c2;
  bool TRUE = (1 == 1);

  c1 = cos(2*PI/5.0);
  c2 = cos(PI/5.0);
  s1 = sin(2*PI/5.0);
  s2 = sin(4*PI/5.0);

  plCurve *thiscurve;
  int  nv = 5;
  bool open = TRUE;
  int  cc = 1;
  
  thiscurve = plc_new(1,&nv,&open,&cc);

  thiscurve->cp[0].vt[0] = plc_build_vect(0,4,0);
  thiscurve->cp[0].vt[1] = plc_build_vect(4*s1,4*c1,0);
  thiscurve->cp[0].vt[2] = plc_build_vect(4*s2,-4*c2,0);
  thiscurve->cp[0].vt[3] = plc_build_vect(-4*s2,-4*c2,0);
  thiscurve->cp[0].vt[4] = plc_build_vect(-4*s1,4*c1,0);

  /* Fix both endpoints and constrain to xy plane. */

  plc_set_fixed(thiscurve,0,0,thiscurve->cp[0].vt[0]);
  plc_set_fixed(thiscurve,0,4,thiscurve->cp[0].vt[4]);
  plc_constrain_to_plane(thiscurve,0,1,3,plc_build_vect(0,0,1),0);
  
  FILE *outfile;

  outfile = fopen("MRtest_pentagon.vect","w");
  if (outfile == NULL) { printf("Could not open MRtest_pentagon.vect.\n"); exit(1); }

  plc_write(outfile,thiscurve);
  
  printf("makeMRstruttest.c\n"
	 "Generated 5 vertex test polygon.\n"
	 "Run the test with runMRtest.\n");

  fclose(outfile);
  plc_free(thiscurve);
  exit(0);

}
  
