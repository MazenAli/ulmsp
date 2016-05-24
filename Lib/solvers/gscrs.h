/*! \file gscrs.h
 *  Gauss-Seidel iterative method for a CRS matrix. Performs fixed number of
 *  iterations. Assumes diagonal entries are sorted first.
 */

#ifndef GSCRS_H
#define GSCRS_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/indexvector.h"
#include "../storage/crs.h"

/*!
 * Performs fixed number of iterations on \f$ Ax=b \f$
 * with \f$ A \f$ in CRS format.
 */
void
gscrs(pccrs A,
      prealvector x,
      pcrealvector b,
      index it);

/*!
 * Performs fixed number of iterations on \f$ Ax=b \f$
 * with \f$ A \f$ in CRS format.
 * Iterations are not performed on fixed nodes corresponding to
 * dirichlet boundary conditions. fixedMask defines fixed nodes
 * by 1 for fixed.
 */
void
gscrs_constrains(pccrs A,
                 prealvector x,
                 pcrealvector b,
                 pcindexvector fixedMask,
                 index it);

#endif
