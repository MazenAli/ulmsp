#ifndef INDEXVECTOR_C
#define INDEXVECTOR_C

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* Conflicting typedefs for "index" */
#define index string_index
    #include <string.h>
#undef index

#include <math.h>

#include "indexvector.h"

pindexvector
new_indexvector(index length)
{
    pindexvector x;
    x = (pindexvector) malloc(sizeof(indexvector));
    init_indexvector(x, length);

    return x;
}


void
init_indexvector(pindexvector x, index length)
{
    assert(x);

    if (length>0) {
        x->vals = (index*) malloc(length*sizeof(index));
    } else {
        x->vals = NULL;
    }
    x->length = length;
}


void
resize_indexvector(pindexvector x, index length)
{
    assert(x);

    if (length>0) {
        x->vals = (index*) realloc(x->vals,length*sizeof(index));
    } else {
        if (x->vals) free(x->vals);
        x->vals = NULL;
    }
    x->length = length;
}


void
copy_indexvector(pindexvector dest, pcindexvector src)
{
    assert(dest);
    assert(src);

    if (dest->length!=src->length) {
        del_indexvector(dest);
        init_indexvector(dest, src->length);
    }

    memcpy(dest->vals, src->vals, dest->length*sizeof(index));
}

void
swap_indexvector(pindexvector x, pindexvector y)
{
    index tmp_length, *tmp_vals;

    assert(x);
    assert(y);

    tmp_length = x->length;
    x->length  = y->length;
    y->length  = tmp_length;

    tmp_vals = x->vals;
    x->vals  = y->vals;
    y->vals  = tmp_vals;
}


void
del_indexvector(pindexvector x)
{
    assert(x);

    if (x->vals) free(x->vals);
    x->length = 0;
}


index
getentry_indexvector(pcindexvector x, index i)
{
    assert(x);
    check_base(i);
    assert(i<x->length+INDEX_BASE);

    return x->vals[i-INDEX_BASE];
}


void
setentry_indexvector(pindexvector x, index i, index entry)
{
    assert(x);
    check_base(i);
    assert(i<x->length+INDEX_BASE);

    x->vals[i-INDEX_BASE] = entry;
}


void
addentry_indexvector(pindexvector x, index i, index entry)
{
    assert(x);
    check_base(i);
    assert(i<x->length+INDEX_BASE);

    x->vals[i-INDEX_BASE] += entry;
}


void
fill_indexvector(pindexvector x, index val)
{
    index i;

    assert(x);

    for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
        setentry_indexvector(x, i, val);
    }
}


void
print_indexvector(pcindexvector x)
{
    index i;

    assert(x);

    for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
        (void) printf("%lu\n", getentry_indexvector(x, i));
    }
}


pindexvector
load_indexvector(char *fname)
{
  FILE *file;
  index cnt, a;
  pindexvector x;

  file = fopen(fname,"r");

  if (file == NULL){
    printf("\n fopen() Error!!!\n\n");
    return NULL;
  }

  cnt = 0;
  while (fscanf(file,"%lu",&a) != EOF) cnt++;

  x = (pindexvector) malloc(sizeof(indexvector));
  init_indexvector(x, cnt);

  fseek(file,0L,SEEK_SET);
  cnt = 0;
  while (fscanf(file,"%lu",&(x->vals[cnt])) != EOF) cnt++;
  fclose(file);
  return x;
}


int
write_indexvector(char *fname, pindexvector v)
{
  FILE *file;
  index i, *vx;

  file = fopen(fname,"w");

  if (file == NULL){
    printf("\n write_indexvector() Error!!!\n\n");
    return 0;
  }

  vx = v->vals;
  for ( i = 0 ; i < v->length ; i++ ) {
    fprintf(file,"%lu\n", vx[i]);
  }
  fclose(file);
  return 1;
}

#endif
