/**
 *  Conjugate Gradient Method (CG) for a CRS matrix.
 */

#ifndef CGCOO_H
#define CGCOO_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/coo.h"

index
cgcoo(pccoo A,
      prealvector x,
      pcrealvector b,
      real tol,
      index maxit);

#endif
