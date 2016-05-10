/*! \file gausscrs.h
 *  Gauss method for a CRS matrix.
 */

#ifndef GAUSSCRS_H
#define GAUSSCRS_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/indexvector.h"
#include "../storage/crs.h"

/*!
 * Solve \f$ Ax=b \f$ with \f$ A \f$ in CRS format.
 * Store result in b.
 */
void
gausscrs(pccrs A,
         prealvector b,
         pcindexvector fixed);

#endif
