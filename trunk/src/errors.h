/*
 *  errors.h
 *  ridgerunner3
 *
 *  Created by Michael Piatek on Fri Jan 16 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _H_errors
#define _H_errors

#include "octrope.h"
#include "plCurve.h"

#ifndef true
#define true (1==1)
#endif
#ifndef false
#define false (1==0)
#endif

enum 
{
    kNoErr              = 0,
    kNULLPointer        = 100,
    kFileIOErr          = 200,
    kUnsupportFileType  = 201,
    kGenericNumerics    = 300,
    kSplineProblem      = 400

};

void DebugThrow( int inErr, const char* inFile, long inLine );
void DebugWarning( int inErr, const char* inFile, long inLine );

void error_write( plCurve* inLink );

#define fatal_(err) DebugThrow( (err), __FILE__, __LINE__ )
#define warning_(err) DebugWarning( (err), __FILE__, __LINE__ )

// we need to wrap these extended defines in do-while statements so the 
// semicolons make sense and code blocks are self-contained
#define fatalifnull_(ptr)                   \
            do {                            \
                if( (ptr) == NULL )         \
                    fatal_(kNULLPointer);   \
            } while( false ) 
            
#define warningifnull_(ptr)                 \
            do {                            \
                if( (ptr) == NULL )         \
                    warning_(kNULLPointer); \
            } while( false ) 

#define fataliferr_(err)                    \
            do {                            \
                int __theErr = err;         \
                if( __theErr != kNoErr )    \
                {                           \
                    fatal_(__theErr);       \
                }                           \
            } while( false )  

#define warningiferr_(err)                  \
            do {                            \
                int __theErr = err;         \
                if( __theErr != kNoErr )     \
                {                           \
                    warning_(__theErr);     \
                }                           \
            } while( false )  



#endif // _H_errors
