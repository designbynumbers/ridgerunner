#include <sys/types.h>
#include <sys/stat.h>
#include <sys/dir.h>

#include <iostream>
#include <fstream>

using namespace std;

void asciiIt( char* path );

int main (int argc, char * const argv[]) 
{
    DIR* workingDir = opendir(".");
    dirent* dp;
    
    system("mkdir -p asciis");
    
    if( workingDir != NULL ) 
    {
        while( (dp = readdir(workingDir)) != NULL )
        {
            char	path[512];
            sprintf( path, "%s", dp->d_name );
      //      cout << path << endl;

            struct stat finfo;
            if( stat( path, &finfo ) == 0 )
            {
                if( !(finfo.st_mode & S_IFDIR) )	// if this isn't a directory...
                {
                    if( strstr(path, "rmov.") != NULL && 
                    		strstr(path, ".ascii") == NULL )
                    {
                    		asciiIt(path);
                    
					cout << "load asciis/" << path << ".ascii" << endl;
					cout << "psmode=2" << endl;
					cout << "cyl-rad=0.5" << endl;
					cout << "smart=false" << endl;
				//	cout << "rotate x .1" << endl;
				//	cout << "rotate y .5" << endl;
				//	cout << "rotate z .5" << endl;
					cout << "display true" << endl;
					cout << "ppmout " << path << " 1" << endl;
					//cout << "psout " << path << endl;
                    }
                }
            }
        }
    }
    
    return 0;
}

void
asciiIt( char* path )
{
	ifstream f(path);
	char		line[1024];
	int		verts = 0;
	int		dummy, comps;
	int		compsizes[20];
	char		asciiname[512];
	
	// eat vect line
	f.getline(line,1023);
	
	
	f.getline(line, 1023);
	sscanf( line, "%d %d %d", &comps, &verts, &dummy );
	
	// eat next 2 useless lines
	for( int i=0; i<comps; i++ )
	{
		f >> compsizes[i];
		compsizes[i] *= -1; // for some reason this is negative
	}
	// eat last line
	f.getline(line, 1023);
	f.getline(line, 1023);
	
	// strip extensions
	path[strlen(path)-5] = '\0';
	
	sprintf(asciiname, "asciis/%s.ascii", path);
	
	ofstream o(asciiname);
	
	int toCheck = 0;
	int compTracker = 0, good;
	for( int i=0; i<verts; i++ )
	{
		good = 0;
		f.getline(line,1023);
		o << line << endl;
		for( int j=0; j<comps; j++ )
		{
			if( compTracker == compsizes[toCheck]-1 )
			{
				o << endl;
				compTracker = 0;
				good = 1;
				toCheck++;
				compsizes[j] = 650000; // in case we'd hit this again
			}
		}
		if( !good )
			compTracker++;	
	}
}
