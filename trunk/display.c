

#include "gradient/eqedge.h"
#include "errors.h"

void init_display()
{
	
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

extern int gSuppressOutput;

void
export_ted(octrope_link* inLink, octrope_strut* strutSet, 
			int inSize, octrope_mrloc* minradSet, int minradLocs, 
			double* compressions, search_state* inState)
{
	char	fname[1024];

	if( gSuppressOutput != 0 )
		return;
	
	preptmpname(fname, "struts_ted.txt", inState);
	FILE* ted = fopen(fname,"w");
	
	fprintf(ted, "STRUTS\n");
	fprintf(ted, "%d %d\n", inSize, minradLocs);
	
	int i;
	for( i=0; i<inSize; i++ )
	{
		fprintf(ted, "%d %d %d %d %3.10lf %3.10lf %3.10lf %e\n", 
			strutSet[i].component[0], strutSet[i].component[1],
			strutSet[i].lead_vert[0], strutSet[i].lead_vert[1],
			strutSet[i].position[0], strutSet[i].position[1],
			strutSet[i].length, compressions[i] );
	}
	for( i=0; i<minradLocs; i++ )
	{
		fprintf(ted, "%d %d %3.10lf %e\n", minradSet[i].component,
					minradSet[i].vert, minradSet[i].mr,
					compressions[i+inSize]);
	}
	
	fclose(ted);
}

extern int gSurfaceBuilding;

short gSurfaceIndex = 0; // the next filename for the surface building
short gSurfaceItr = 0; // how long till next surface output

void 
export_struts(octrope_link* inLink, octrope_strut* strutSet, int inSize, double* compressions, search_state* inState)
{
	if( inSize == 0 || compressions == NULL )
		return;
	
	if( gSuppressOutput == 1 )
		return;
		
	double time = inState->time;

	char	fname[1024];

	preptmpname(fname, "struts.vect", inState);
	FILE* fp = fopen(fname, "w");
	int i=0;
	double maxCompression;
		
	preptmpname(fname, "struts_meta.txt", inState);
	FILE* meta = fopen(fname,"w");
	
	preptmpname(fname, "struts_ted.txt", inState);
	FILE* ted = fopen(fname, "w");
	
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
		octrope_vector  points[2];
		
		octrope_strut_ends( inLink, &strutSet[i], points );
		
		fprintf( meta, "%d %d %d %d %f\n", strutSet[i].component[0], 
										strutSet[i].lead_vert[0],
										strutSet[i].component[1],
										strutSet[i].lead_vert[1],
										compressions[i]/maxCompression );
					
		fprintf(fp, "%lf %lf %lf\n",
			points[0].c[0], 
			points[0].c[1], 
			points[0].c[2] );
				
		fprintf(fp, "%lf %lf %lf\n",
			points[1].c[0], 
			points[1].c[1], 
			points[1].c[2] );
			
	}
		
	for( i=0; i<inSize; i++ )
	{
		fprintf(fp, "%f,%f,%f,0\n", compressions[i]/maxCompression, 0.0, 0.0);
	}
	
	fclose(fp);
	fclose(meta);
	fclose(ted);
	
	if( gSurfaceBuilding != 0 )
	{
		if( gSurfaceItr == 0 )
		{
			char cmd[1024];
			printf("updating strut surface %d\n", gSurfaceIndex);
			sprintf(cmd, "cp /tmp/struts.vect /tmp/struts%d.vect", gSurfaceIndex);
			system(cmd);
			sprintf(cmd, "cp /tmp/struts_meta.txt /tmp/struts_meta%d.txt", gSurfaceIndex);
			system(cmd);
			if( gSurfaceIndex == 5 )
				gSurfaceIndex = 0;
			else
				gSurfaceIndex++;
			gSurfaceItr = 0; //(rand() % 1);
		}
		else 
		{
			gSurfaceItr--;
		}
	}
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
				link->cp[cItr].vt[vItr].c[0] + .25*dl[i].c[0],
				link->cp[cItr].vt[vItr].c[1] + .25*dl[i].c[1],
				link->cp[cItr].vt[vItr].c[2] + .25*dl[i].c[2] );
		}
	}
	
	// Deep sky blue <- according to Rawdon 8)
/*	if( strcmp( fname, "/tmp/dVdt.vect") == 0 )
		fprintf( fp, "0,0,0,0\n");
	else*/
		fprintf( fp, "0,0.7461,0.9661,.5\n" );
	
	fclose(fp);
}
