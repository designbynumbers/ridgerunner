

#include "gradient/eqedge.h"
#include "errors.h"

void init_display()
{
	gclpipe = fopen("/tmp/OOGL","w");
	fatalifnull_(gclpipe);
	
	system("rm -f rmov.*");

	system("rm -f /tmp/struts.vect");
	system("touch /tmp/struts.vect");
	
	system("rm -f /tmp/dl.vect");
	system("touch /tmp/dl.vect");
	
	system("rm -f /tmp/dVdt.vect");
	system("touch /tmp/dVdt.vect");
	
	fprintf(gclpipe,"(normalization g0 none)\n");
	fprintf(gclpipe,"(bbox-draw g0 no)\n");
}

void
shutdown_display()
{
	fclose(gclpipe);
}

void refresh_display(octrope_link *L) 

     /* Procedure displays the curve on the linked copy of Geomview. */
     /* We note that OOGLPIPE must be running in Geomview, and we must */
     /* be able to write to /tmp. */

{

  FILE *curvepipe;

  /* Open the curve file and update it. */

   curvepipe = fopen("/tmp/curvepipe.oogl","w");
  octrope_link_draw(curvepipe,L);
  fclose(curvepipe);

  /* Now tell Geomview about the situation. */

  fprintf(gclpipe,"(geometry curve { < /tmp/curvepipe.oogl})\n");
  fprintf(gclpipe,"(geometry struts { < /tmp/struts.vect})\n");
  fprintf(gclpipe,"(geometry dl { < /tmp/dl.vect})\n");
  fprintf(gclpipe,"(geometry dVdt { < /tmp/dVdt.vect})\n");
    
  fprintf(gclpipe,"(look g0)\n");

  fflush(gclpipe);

}

void 
export_struts(octrope_link* inLink, octrope_strut* strutSet, int inSize, double* compressions, double time)
{
	if( inSize == 0 )
		return;

	FILE* fp = fopen("/tmp/struts.vect", "w");
	int i=0;
	double maxCompression;
	
	FILE* kp = fopen("/tmp/struts.ascii", "w");
	
	if( inSize == 0 )
	{
		fclose(fp);
		return;
	}
	
	fprintf(fp, "VECT\n");
	fprintf(fp, "%d %d %d\n", inSize, 2*inSize, inSize);
	fprintf(fp, "2");
	for( i=1; i<inSize; i++ )
		fprintf(fp, " 2");
	fprintf(fp, "\n");
	
	fprintf(fp, "1");
	for( i=1; i<inSize; i++ )
		fprintf(fp, " 1");
	fprintf(fp, "\n");
	
	for( i=0; i<inSize; i++ )
	{
		octrope_vector  points[2];
		
		octrope_strut_ends( inLink, &strutSet[i], points );
			
		fprintf(fp, "%lf %lf %lf\n",
			points[0].c[0], 
			points[0].c[1], 
			points[0].c[2] );
		
		fprintf(kp, "%lf %lf %lf\n",
			points[0].c[0], 
			points[0].c[1], 
			points[0].c[2] );
		
		fprintf(fp, "%lf %lf %lf\n",
			points[1].c[0], 
			points[1].c[1], 
			points[1].c[2] );
			
		fprintf(kp, "%lf %lf %lf\n\n",
			points[1].c[0], 
			points[1].c[1], 
			points[1].c[2] );
	}
	
	maxCompression = 0;
	for( i=0; i<inSize; i++ )
	{
		if( compressions[i] > maxCompression )
			maxCompression = compressions[i];
	}
	
	if( maxCompression == 0 )
		maxCompression = 1e-6;
	
	for( i=0; i<inSize; i++ )
	{
		fprintf(fp, "%f,%f,%f,0\n", compressions[i]/maxCompression, 0.0, 0.0);
	}
	
	fclose(kp);
	fclose(fp);
}

void
exportVect( const octrope_vector* dl, octrope_link* link, const char* fname )
{
	int totalVerts = 0, cItr, vItr, i;
	FILE* fp = fopen( fname, "w" );
	
	for( cItr=0; cItr<link->nc; cItr++ )
		totalVerts += link->cp[cItr].nv;
	
	fprintf(fp, "VECT\n");
	fprintf(fp, "%d %d 1\n", totalVerts, 2*totalVerts );
	
	fprintf(fp, "2");
	for( i=1; i<totalVerts; i++ )
		fprintf( fp, " 2");
	fprintf(fp, "\n");
	
	fprintf(fp, "1");
	for( i=1; i<totalVerts; i++ )
		fprintf( fp, " 0");
	fprintf(fp, "\n");
	
	for( cItr=0, i=0; cItr<link->nc; cItr++ )
	{
		for( vItr=0; vItr<link->cp[cItr].nv; vItr++, i++ )
		{
			fprintf( fp, "%lf %lf %lf\n", 
				link->cp[cItr].vt[vItr].c[0],
				link->cp[cItr].vt[vItr].c[1],
				link->cp[cItr].vt[vItr].c[2] );
			fprintf( fp, "%lf %lf %lf\n", 
				link->cp[cItr].vt[vItr].c[0] + dl[i].c[0],
				link->cp[cItr].vt[vItr].c[1] + dl[i].c[1],
				link->cp[cItr].vt[vItr].c[2] + dl[i].c[2] );
		}
	}
	
	// Deep sky blue <- according to Rawdon 8)
	if( strcmp( fname, "/tmp/dVdt.vect") == 0 )
		fprintf( fp, "0,0,0,0\n");
	else
		fprintf( fp, "0,0.7461,0.9661,.5\n" );
	
	fclose(fp);
}
