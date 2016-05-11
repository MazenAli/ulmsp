#ifndef SETTINGS_C
#define SETTINGS_C

#include <assert.h>

#include "settings.h"

void
check_base(index i)
{
    index base = INDEX_BASE;
    assert(i>=base);
    (void) base;
    (void) i;
}

#endif
