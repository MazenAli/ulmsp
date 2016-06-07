/*! \file multigrid.h
 *  The following functions perform a Multigrid-solver for a CRS matrix.
 */

#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/indexvector.h"
#include "../Lib/mesh/mesh.h"
#include "../storage/crs.h"
#include "gscrs.h"
#include "cgcrs.h"


/*!
 * Solve \f$ Ax=b \f$ with \f$ A \f$ in CRS format.
 * Returns number of iterations.
 */
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
		real tol);               	/* in     */


/*!
 * Solve \f$ Ax=b \f$ with \f$ A \f$ in CRS format.
 * Iterations are not performed on fixedNodes corresponding to
 * dirichlet boundary conditions.
 * Returns number of iterations.
 */
void
mgcrs_constrains(pcrs *A,                   /* in    */
                 prealvector x,             /* in/out */
                 pcrealvector b,            /* in     */
                 pindexvector *fixedNodes,  /* in     */
				 const int Level,			/* in     */
				 const int nLevel,			/* in     */
				 const int nPreSmoothing,	/* in     */
				 const int nPostSmoothing,	/* in     */
				 const int gamma, 			/* in     */
				 pindexmatrix *f2s,	 		/* in     */
                 real tol);               	/* in     */


#endif
