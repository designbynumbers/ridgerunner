/* This is a short self-test program which tests placeMinRad2Sparse
   functionality.  At the moment, it's undeveloped and just contains
   some raw code from an old version. This won't work as written. */

static void
glom( plCurve* inLink, search_state* inState )
{
  octrope_mrloc strut;
  
  strut.component = 0;
  strut.vert = 1;
  strut.mr = octrope_minradval(inLink);
  
  int Acols = (3*3)*(1);
  double* A = (double*)calloc(Acols, sizeof(double)); // calloc zeros A
  
  double* blah = (double*)calloc(Acols, sizeof(double));
  
  double dpsum = 0, dx;
  double diff;
  int i, mItr;
  
  for( i=0; i<1000; i++ )
    {
      plc_vector motion[3];
      plc_vector foo;
      double initialMR;
      
      dpsum = 0;
      
      for( mItr=0; mItr<3; mItr++ )
	{
	  motion[mItr].c[0] = (drand48()-0.5)*2;
	  motion[mItr].c[1] = (drand48()-0.5)*2;
	  motion[mItr].c[2] = (drand48()-0.5)*2;
	}
      
      initialMR = octrope_minradval(inLink);
      
      placeMinradStruts2(A, inLink, &strut, 1, inState, 0);
      
      foo.c[0] = A[0];
      foo.c[1] = A[1];
      foo.c[2] = A[2];
      dpsum += plc_M_dot(foo,motion[0]);
      foo.c[0] = A[3];
      foo.c[1] = A[4];
      foo.c[2] = A[5];
      dpsum += plc_M_dot(foo,motion[1]);
      foo.c[0] = A[6];
      foo.c[1] = A[7];
      foo.c[2] = A[8];
      dpsum += plc_M_dot(foo,motion[2]);
      
      for( dx=1e-1; dx>1e-10; dx *= 0.1 )
	{
	  plCurve* workerLink = plc_copy(inLink);
	  plc_fix_wrap(workerLink);
	  step( workerLink, dx, motion, inState );
	  
	  // if everything's right: (minrad of workerLink) - (minrad of inLink) <--- ~= dpsum*dx. we'll see!
	  diff = octrope_minradval(workerLink) - octrope_minradval(inLink);
	  printf( "diff: %e / dpsum*dx: %e / (diff - dpsum*dx): %e / still!!: %e\n", 
		  diff, dpsum*dx, diff-(dpsum*dx), (diff-(dpsum*dx))/dx );
	  
	  plc_free(workerLink);
	}
    }
  
  free(A);
}

