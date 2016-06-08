#ifndef MULTIGRID_C
#define MULTIGRID_C

#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "../ops/gecrsmv.h"
#include "multigrid.h"

void
mgcrs(	pcrs *A,                   /* in    */
    	prealvector x,             /* in/out */
      	pcrealvector b,            /* in     */
		const int Level,			/* in     */
		const int nLevel,			/* in     */
		const int nPreSmoothing,	/* in     */
		const int nPostSmoothing,	/* in     */
		const int gamma, 			/* in     */
		pindexmatrix *f2s,	 		/* in     */
		real tol)               	/* in     */
{
	index k;
	
	if( (Level+1) == nLevel-1){ /*If coarsest grid is achieved -> solve exactly */
		(void) cgcrs(A[Level], x, b, 1e-6, x->length);
	}else{
    	
		prealvector r;
		prealvector r_restricted;
		prealvector v;
		transpose t = notrans;
		

		/*Perform some Pre-Smoothing steps */
		gscrs(A[Level], x, b, nPreSmoothing);
		
		/* Compute residual, r = b-A*x */
		r = new_realvector(A[Level]->numc);
    	copy_realvector(r, b);
    	gecrsmv(t, -1., A[Level], x, 1., r);
		
		/* Restrict the residual -> coarser grid on Level+1 */
		r_restricted = new_realvector(A[Level+1]->numc);
		v = new_realvector(A[Level+1]->numc);
		fill_realvector(v,0.0);
		restriction(r,r_restricted,f2s[Level]);
		
		/* Compute the correction on the coarser grid (Level+1) */
		for(k=0;k<(index) gamma;k++){ 
			mgcrs(A, v, r_restricted, Level+1, nLevel, nPreSmoothing, nPostSmoothing, gamma, f2s, tol); 
		}
		
		/* Prolongate the correction on the finer grid (Level) */
		prolongation(v,r,f2s[Level]); /* We can overwrite r, we don't need it anymore...and it has already the right size */
		
		/* Add now the prolongated correction */
		axpy_realvector(1.0, r, x);
		
		/*Perform some Post-Smoothing steps */
		gscrs(A[Level], x, b, nPostSmoothing);
	}
	
}


void
mgcrs_constrains(pcrs *A,                    /* in    */
                 prealvector x,              /* in/out */
                 pcrealvector b,             /* in     */
                 pindexvector *fixedNodes,   /* in     */
				 const int Level,			 /* in	  */
				 const int nLevel,  		 /* in	  */
				 const int nPreSmoothing,	 /* in	  */
				 const int nPostSmoothing,	 /* in	  */
				 const int gamma, 			 /* in     */
				 pindexmatrix *f2s,		 	 /* in     */
				 real tol)                   /* in     */
{
    
	index k, i;
	
	if( (Level+1) == nLevel-1){ /*If coarsest grid is achieved -> solve exactly */
		(void) cgcrs_constrains(A[Level], x, b, fixedNodes[Level], 1e-6, x->length);
	}else{
    	
		prealvector r;
		prealvector r_restricted;
		prealvector v;
		transpose t = notrans;
		pindexvector fixedMask = new_indexvector(x->length);
		
		/* Build fixedMask for gs */
		for (k=0; k<fixedMask->length; ++k) {
	        fixedMask->vals[k] = 0;
	        for (i=0; i<fixedNodes[Level]->length; ++i) {
	            if (k==fixedNodes[Level]->vals[i]-INDEX_BASE) ++fixedMask->vals[k];
	        }
	    }
		

		/*Perform some Pre-Smoothing steps */
		gscrs_constrains(A[Level], x, b, fixedMask, nPreSmoothing);
		
		/* Compute residual, r = b-A*x */
		r = new_realvector(A[Level]->numc);
    	copy_realvector(r, b);
    	gecrsmv(t, -1., A[Level], x, 1., r);
		
    	/* Incorporate constrains */
    	for ( k =  0 ; k < fixedNodes[Level]->length; ++k){
        	r->vals[fixedNodes[Level]->vals[k]-INDEX_BASE] = 0.;
    	}
		
		/* Restrict the residual -> coarser grid on Level+1 */
		r_restricted = new_realvector(A[Level+1]->numc);
		v = new_realvector(A[Level+1]->numc);
		fill_realvector(v,0.0);
		restriction(r,r_restricted,f2s[Level]);
		
		/* Compute the correction on the coarser grid (Level+1) */
		for(k=0;k<(index) gamma;k++){ 
			mgcrs_constrains(A, v, r_restricted, fixedNodes, Level+1, nLevel, nPreSmoothing, nPostSmoothing, gamma, f2s, tol); 
		}
		
		/* Prolongate the correction on the finer grid (Level) */
		prolongation(v,r,f2s[Level]); /* We can overwrite r, we don't need it anymore...and it has already the right size */
		
		/* Again we have to pay attention to the fixedNodes -> set them to 0 */
		for ( k =  0 ; k < fixedNodes[Level]->length; ++k){
        	r->vals[fixedNodes[Level]->vals[k]-INDEX_BASE] = 0.0;
    	}
		
		/* Add now the prolongated correction */
		axpy_realvector(1.0, r, x);
		
		/*Perform some Post-Smoothing steps */
		gscrs_constrains(A[Level], x, b, fixedMask, nPostSmoothing);
	}
}



#endif
