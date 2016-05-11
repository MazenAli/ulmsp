/*! \file settings.h
 *  Some typedefs and macros.
 */

#ifndef SETTINGS_H
#define SETTINGS_H

#include <stddef.h>

typedef size_t index; /*!< Data type for indexing. */
typedef double real;  /*!< Data type for reals. */

/*! Used in matrix vector operations to perform \f$ A \f$ or \f$ A^T \f$. */
typedef enum
{
    notrans = 0,
    trans   = 1
} transpose;

#ifndef INDEX_BASE
    #define INDEX_BASE 0 /*!< Base for indexing vectors and matrices.
                              Set in options.inc. Default is 0. */
#endif

/*! Check if i is valid against INDEX_BASE. */
void
check_base(index i);

#endif
