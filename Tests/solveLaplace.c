/**
 *  Refine the triangle (0,0), (1,0), (0,1).
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "settings.h"
#include "realvector.h"
#include "indexvector.h"
#include "realmatrix.h"
#include "indexmatrix.h"
#include "mesh.h"
#include "gecrsmv.h"
#include "gecoomv.h"
#include "cgcrs.h"
#include "gscrs.h"

#define ELAPSED(t0,t1) ((int) ((t1 - t0) / (double) CLOCKS_PER_SEC * 1000))


real
f1(real x[2], real m)
{
    (void) m;
    (void) x;
    return 1.;  /* f(x,y) =  */
}


real*
f2(real x[2], real m)
{
    (void) x;
    (void) m;
    static real ret[2];
    ret[0] = (real) 0;
    ret[1] = (real) 0;

    return ret;
}


real
g(real *x, real m)
{
    (void) m;
    (void) x;
    return -0.;   /* value on Neumann boundary */
}


real
uD(real *x)
{
    (void) x;
    return 0.;   /* value on Dirichlet boundary */
}


int
main(int argc, char **argv)
{
    int i, nItCG, nRef = 2;
    long TIME[10];

    pcoo S = new_coo(1,1,1);
    pcrs A = new_crs(1,1,1);
    pindexmatrix edgeno = new_indexmatrix(0,0);
    pindexvector fixedNodes = new_indexvector(0);
    prealvector material
            = load_realvector("./Tests/Example1/material.dat");
    prealmatrix coordinates
            = load_realmatrix("./Tests/Example1/coordinates.dat", 2, 1);
    pindexmatrix elements
            = load_indexmatrix("./Tests/Example1/elements.dat",3, 1);
    pindexmatrix elements2edges
            = load_indexmatrix("./Tests/Example1/elements2edges.dat",3, 1);
    int nBdry = 2;
    pindexvector bdrytyp = new_indexvector(2);
    pindexmatrix bdrylist[2];
    pindexvector bdry2edgeslist[2];
    prealvector rhs;
    prealvector solcg;
    prealvector solpcg;

    /* Number of refinements */
    if (argc>1){
      if (atoi(argv[1]) < 12) nRef = atoi(argv[1]);
    }
    /* Load geometry */
    bdrytyp->vals[0] = 0;
    bdrytyp->vals[1] = 1;
    bdrylist[0] = load_indexmatrix("./Tests/Example1/Dirichlet.dat",2,1);
    bdrylist[1] = load_indexmatrix("./Tests/Example1/Neumann.dat",2,1);
    bdry2edgeslist[0] = load_indexvector("./Tests/Example1/Dirichlet2edges.dat");
    bdry2edgeslist[1] = load_indexvector("./Tests/Example1/Neumann2edges.dat");


    /* Show geometry */
    printf("====================\n");
    printf("coordinates:\n");
    printf("====================\n");
    print_realmatrix(coordinates);
    printf("====================\n");
    printf("elements:\n");
    printf("====================\n");
    print_indexmatrix(elements);
    printf("====================\n");
    printf("elements2edges:\n");
    printf("====================\n");
    print_indexmatrix(elements2edges);
    printf("====================\n");
    printf("material:\n");
    printf("====================\n");
    print_realvector(material);
    printf("====================\n");
    printf("bdrylist:\n");
    printf("====================\n");
    print_indexmatrix(bdrylist[0]);
    printf("--------------------\n");
    print_indexmatrix(bdrylist[1]);
    printf("====================\n");
    printf("bdry2edgeslist:\n");
    printf("====================\n");
    print_indexvector(bdry2edgeslist[0]);
    printf("--------------------\n");
    print_indexvector(bdry2edgeslist[1]);
    printf("====================\n");

    TIME[0] = clock();
    /* Refine mesh uniformly */
    for (i=0; i<nRef; i++){  /* 11 scales fine, 14 too large for 8 GByte PC */
      refine_uniform(coordinates,       /* in / out */
                     elements,          /* in / out */
                     material,          /* in / out */
                     elements2edges,    /* in / out */
                     edgeno,            /* out      */
                     bdrylist,          /* in / out */
                     bdry2edgeslist,    /* in / out */
                     nBdry);
    }
    TIME[1] = clock();

    write_realmatrix("./Tests/Example1/coordinates_fine.dat",coordinates,1);
    write_indexmatrix("./Tests/Example1/elements_fine.dat",elements,1);
    write_indexmatrix("./Tests/Example1/elements2edges_fine.dat",elements,1);
    write_realvector("./Tests/Example1/material_fine.dat",material);
    write_indexmatrix("./Tests/Example1/Dirichlet_fine.dat",bdrylist[0],1);
    write_indexmatrix("./Tests/Example1/Neumann_fine.dat",bdrylist[1],1);
    write_indexvector("./Tests/Example1/Dirichlet2edges_fine.dat",bdry2edgeslist[0]);
    write_indexvector("./Tests/Example1/Neumann2edges_fine.dat",bdry2edgeslist[1]);

    rhs    = new_realvector(coordinates->cols);
    solcg  = new_realvector(coordinates->cols);
    solpcg = new_realvector(coordinates->cols);
    fill_realvector(solcg, 0.);
    fill_realvector(solpcg, 0.);
    fill_realvector(rhs, 0.);

    TIME[2] = clock();

    buildStiffness(coordinates, elements, S);

    TIME[3] = clock();

    init_coo2crs(A, S);

    TIME[4] = clock();

    getFixed(coordinates->cols, bdrytyp, bdrylist, fixedNodes);

    TIME[5] = clock();

    buildRhs(coordinates, elements, bdrylist, material, f1, f2, g, rhs);
    setDirichletData2Rhs(coordinates, fixedNodes, uD, solcg);
    setDirichletData2Rhs(coordinates, fixedNodes, uD, solpcg);

    TIME[6] = clock();

    nItCG = cgcrs_constrains(A, solcg, rhs, fixedNodes, 1e-6, coordinates->cols);

    printf("No. iterations CG %i\n", nItCG);

    TIME[7] = clock();

    nItCG = pcggscrs_constrains(A, solpcg, rhs, fixedNodes, 1e-6, 5,
                                coordinates->cols);

    printf("No. iterations PCG %i\n", nItCG);

    TIME[8] = clock();

    write_realvector("./Tests/Example1/sol_fine.dat",solpcg);

    TIME[9] = clock();

    printf("---------------------\n");
    printf("Time for refinement       = %i ms\n", ELAPSED(TIME[0],TIME[1]));
    printf("Time for saving ref. mesh = %i ms\n", ELAPSED(TIME[1],TIME[2]));
    printf("Time for assembling       = %i ms\n", ELAPSED(TIME[2],TIME[3]));
    printf("Time for COO -> CRS       = %i ms\n", ELAPSED(TIME[3],TIME[4]));
    printf("Time for fixedNodes       = %i ms\n", ELAPSED(TIME[4],TIME[5]));
    printf("Time for RHS              = %i ms\n", ELAPSED(TIME[5],TIME[6]));
    printf("Time for solving CG       = %i ms\n", ELAPSED(TIME[6],TIME[7]));
    printf("Time for solving PCG      = %i ms\n", ELAPSED(TIME[7],TIME[8]));
    printf("Time store solution       = %i ms\n\n", ELAPSED(TIME[8],TIME[9]));
    printf("Degrees of elements / freedom %lu / %lu\n", elements->cols, coordinates->cols);

    printf("---------------------\n");

    del_coo(S);
    del_crs(A);
    del_realvector(rhs);
    del_realvector(solcg);
    del_realvector(solpcg);
    del_indexmatrix(edgeno);
    del_indexvector(fixedNodes);

    del_realmatrix(coordinates);
    del_indexmatrix(elements);
    del_indexmatrix(elements2edges);
    del_realvector(material);

    return 0;
}
