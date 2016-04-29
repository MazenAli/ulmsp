/**
 *  Typedefs and macros
 */

#ifndef SETTINGS_H
#define SETTINGS_H

#include <stddef.h>

typedef size_t index;
typedef double real;

typedef enum
{
    notrans = 0,
    trans   = 1
} transpose;

#ifndef INDEX_BASE
    #define INDEX_BASE 0
#endif

#endif
