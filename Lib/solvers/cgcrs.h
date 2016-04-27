/**
 *  Conjugate Gradient Method (CG) for a CRS matrix.
 */

#ifndef CGCRS_H
#define CGCRS_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/crs.h"

index
cgcrs(pccrs A,
      prealvector x,
      pcrealvector b,
      real tol,
      index maxit);

#endif
