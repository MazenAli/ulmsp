/*! \file cgcrs.h
 *  Conjugate Gradient Method (CG) for a CRS matrix.
 */

#ifndef CGCRS_H
#define CGCRS_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/indexvector.h"
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
 * Uses the Jacobi preconditioner.
 * Returns number of iterations.
 */
index
pcgdiagcrs(pccrs A,                   /* in     */
           prealvector x,             /* in/out */
           pcrealvector b,            /* in     */
           real tol,                  /* in     */
           index maxit);              /* in     */

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
 * Uses the Jacobi preconditioner.
 * Returns number of iterations.
 */
index
pcgdiagcrs_constrains(pccrs A,                   /* in     */
                      prealvector x,             /* in/out */
                      pcrealvector b,            /* in     */
                      pcindexvector fixedNodes,  /* in     */
                      real tol,                  /* in     */
                      index maxit);              /* in     */


/*!
 * Solve \f$ Ax=b \f$ with \f$ A \f$ in CRS format.
 * Iterations are not performed on fixedNodes corresponding to
 * dirichlet boundary conditions.
 * Uses the symmetric GS preconditioner.
 * Returns number of iterations.
 */
index
pcgsymgscrs_constrains(pccrs A,                   /* in     */
                       prealvector x,             /* in/out */
                       pcrealvector b,            /* in     */
                       pcindexvector fixedNodes,  /* in     */
                       real tol,                  /* in     */
                       index maxit);              /* in     */

#endif
