/**
 *  Refine the triangle (0,0), (1,0), (0,1).
 */

#include <stdio.h>
#include<time.h>

#include "settings.h"
#include "mesh.h"
#include "gecrsmv.h"
#include "gecoomv.h"
#include "gausscrs.h"

#define ELAPSED(t0,t1) ((int) ((t1 - t0) / (double) CLOCKS_PER_SEC * 1000))


real VolForce(real *x){
  return x[0];  /* f(x,y) = x */
}

real uD(real *x){
  (void) x;
  return 0.0;   /* u(x,y) = 0, value on Dirichlet boundary */
}


int
main()
{
    int i;
    index j;
    index nBdry = 1;
    long TIME[9];
    pcoo S = new_coo(1,1,1);
    pcrs A = new_crs(1,1,1);
    prealvector rhs = new_realvector(1);
    pindexmatrix elements;
    pindexvector material;
    pindexmatrix elements2edges;
    pindexmatrix edgeno;
    pindexvector bdrytyp = new_indexvector(1);
    pindexvector fixedNodes = new_indexvector(0);
    pindexmatrix bdrylist[1];
    pindexvector bdry2edgeslist[1];
    transpose t   = notrans;
    prealvector y;
    prealvector x;


    /* coordinates */
    prealmatrix coordinates = new_realmatrix(2,3);

    coordinates->vals[0] = 0.0;  coordinates->vals[1] = 0.0;
    coordinates->vals[2] = 1.0;  coordinates->vals[3] = 0.0;
    coordinates->vals[4] = 0.0;  coordinates->vals[5] = 1.0;

    printf("coordinates:\n");
    print_realmatrix(coordinates);

    /* elements */
    elements = new_indexmatrix(3,1);

    elements->vals[0] = INDEX_BASE+0;  elements->vals[1] = INDEX_BASE+1;  elements->vals[2] = INDEX_BASE+2;

    printf("elements:\n");
    print_indexmatrix(elements);

    /* material */
    material = new_indexvector(1);

    material->vals[0] = INDEX_BASE+13;

    printf("material:\n");
    print_indexvector(material);

    /* elements2edges */
    elements2edges = new_indexmatrix(3,1);

    elements2edges->vals[0] = INDEX_BASE+2;  elements2edges->vals[1] = INDEX_BASE+0;  elements2edges->vals[2] = INDEX_BASE+1;

    printf("elements2edges:\n");
    print_indexmatrix(elements2edges);

    /* edgeno */
    edgeno = new_indexmatrix(0,0);

    /* bdrylist */
    bdrytyp->vals[0] = 0;
    bdrylist[0] = new_indexmatrix(2,3);
    bdrylist[0]->vals[0] = INDEX_BASE+0; bdrylist[0]->vals[1] = INDEX_BASE+1;
    bdrylist[0]->vals[2] = INDEX_BASE+1; bdrylist[0]->vals[3] = INDEX_BASE+2;
    bdrylist[0]->vals[4] = INDEX_BASE+2; bdrylist[0]->vals[5] = INDEX_BASE+0;
    bdry2edgeslist[0] = new_indexvector(3);
    bdry2edgeslist[0]->vals[0] = INDEX_BASE+2; 
    bdry2edgeslist[0]->vals[1] = INDEX_BASE+0; 
    bdry2edgeslist[0]->vals[2] = INDEX_BASE+1; 

    printf("bdrylist:\n");
    print_indexmatrix(bdrylist[0]);

    for (i=0; i<6; i++){  /* 11 scales fine, 14 too large for 8 GByte PC */
      TIME[0] = clock();
      refine_uniform(coordinates,       /* in / out */
                     elements,          /* in / out */
                     material,          /* in / out */
                     elements2edges,    /* in / out */
                     edgeno,            /* out      */
                     bdrylist,          /* in / out */
                     bdry2edgeslist,    /* in / out */
                     nBdry);            /* in       */

      TIME[1] = clock();
      buildStiffness(coordinates,       /* in   */
                     elements,          /* in   */
                     S);                /* out  */

      TIME[2] = clock();
      init_coo2crs(A, S);


      TIME[3] = clock();

      getFixed(coordinates->cols, bdrytyp, bdrylist, fixedNodes);

      buildRhs(coordinates, elements, VolForce, rhs); 
      setDirichletData2Rhs(coordinates, fixedNodes, uD, rhs);

      TIME[4] = clock();

      y = new_realvector(A->numr);
      x = new_realvector(A->numc);
      for (j=0 ; j<x->length; j++) x->vals[j] = 1.0;

      TIME[5] = clock();

      gecrsmv(t,(real) 1.0, A,x,(real) 1.0,y);

      TIME[6] = clock();

      for (j=0 ; j<x->length; j++) x->vals[j] = 1.0;

      TIME[7] = clock();

      gausscrs(A, rhs, fixedNodes);

      TIME[8] = clock();


      printf("-- %2i ------------------\n",i+1);
      printf("Time for Refinement = %i ms\n", ELAPSED(TIME[0],TIME[1]));
      printf("Time for Assembling = %i ms\n", ELAPSED(TIME[1],TIME[2]));
      printf("Time for COO -> CRS = %i ms\n", ELAPSED(TIME[2],TIME[3]));
      printf("Time for RHS        = %i ms\n", ELAPSED(TIME[3],TIME[4]));
      printf("Time for alloc x,y  = %i ms\n", ELAPSED(TIME[4],TIME[5]));
      printf("Time for A * x      = %i ms\n", ELAPSED(TIME[5],TIME[6]));
      printf("Time for (re)init x = %i ms\n", ELAPSED(TIME[6],TIME[7]));
      printf("Time for solving    = %i ms\n", ELAPSED(TIME[7],TIME[8]));
      printf("Degrees of elements / freedom %lu / %lu\n", elements->cols, coordinates->cols);
      printf("Degrees on boundary %lu \n", fixedNodes->length);

    }

    printf("---------------------\n");

    del_realmatrix(coordinates);
    del_indexmatrix(elements);
    del_indexvector(material);
    del_indexmatrix(elements2edges);
 
    return 0;
}
