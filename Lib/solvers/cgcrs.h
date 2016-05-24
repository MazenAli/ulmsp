/*! \file cgcrs.h
 *  Conjugate Gradient Method (CG) for a CRS matrix.
 */

#ifndef CGCRS_H
#define CGCRS_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/crs.h"

/*!
 * Solve \f$ Ax=b \f$ with \f$ A \f$ in CRS format.
 * Returns number of iterations.
 */
index
cgcrs(pccrs A,
      prealvector x,
      pcrealvector b,
      real tol,
      index maxit);

/*!
 * Solve \f$ Ax=b \f$ with \f$ A \f$ in CRS format.
 * Iterations are not performed on fixedNodes corresponding to
 * dirichlet boundary conditions.
 * Returns number of iterations.
 */
index
cgcrs_constrains(pccrs A,                   /* in     */
                 prealvector x,             /* in/out */
                 pcrealvector b,            /* in     */
                 pcindexvector fixedNodes,  /* in     */
                 real tol,                  /* in     */
                 index maxit);              /* in     */

/*!
 * Solve \f$ Ax=b \f$ with \f$ A \f$ in CRS format.
 * Iterations are not performed on fixedNodes corresponding to
 * dirichlet boundary conditions.
 * The Gauss-Seidel method is used as a preconditioner.
 * Returns number of iterations.
 */
index
pcggscrs_constrains(pccrs A,                   /* in     */
                    prealvector x,             /* in/out */
                    pcrealvector b,            /* in     */
                    pcindexvector fixedNodes,  /* in     */
                    real tol,                  /* in     */
                    index itgs,                /* in     */
                    index maxit);              /* in     */

#endif
