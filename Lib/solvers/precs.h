/*! \file precs.h
 *  Preconditioners for iterative methods.
 *  Assume diagonal entries are sorted first.
 */

#ifndef PRECS_H
#define PRECS_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/indexvector.h"
#include "../storage/crs.h"

/*!
 * Diagonal (Jacobi).
 */
void
diagcrs(pccrs A,
        prealvector x,
        pcrealvector b);

/*!
 * Diagonal (Jacobi). Performs solve on nodes not
 * marked in fixedMask by 1.
 */
void
diagcrs_constrains(pccrs A,
                   prealvector x,
                   pcrealvector b,
                   pcindexvector fixedMask);

/*!
 * Forward GS. Performs solve on nodes not
 * marked in fixedMask by 1.
 */
void
frwgscrs_constrains(pccrs A,
                    prealvector x,
                    pcrealvector b,
                    pcindexvector fixedMask);

/*!
 * Backward GS. Performs solve on nodes not
 * marked in fixedMask by 1.
 */
void
bkwgscrs_constrains(pccrs A,
                   prealvector x,
                   pcrealvector b,
                   pcindexvector fixedMask);

/*!
 * Symmetric GS. Performs solve on nodes not
 * marked in fixedMask by 1.
 */
void
symgscrs_constrains(pccrs A,
                   prealvector x,
                   pcrealvector b,
                   pcindexvector fixedMask);

#endif
