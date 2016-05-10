/*! \file mesh.h
 *  Functions for generating discretized FEM problem.
 */

#ifndef MESH_H
#define MESH_H

#include "../settings.h"
#include "../storage/coo.h"
#include "../storage/crs.h"
#include "../storage/realvector.h"
#include "../storage/realmatrix.h"
#include "../storage/indexvector.h"
#include "../storage/indexmatrix.h"

void refine_uniform(prealmatrix coordinates,        /* in / out */
                    pindexmatrix elements,          /* in / out */
                    pindexvector material,          /* in / out */
                    pindexmatrix elements2edges,    /* in / out */
                    pindexmatrix edgeno,            /* out      */
                    pindexmatrix *bdrylist,         /* in / out */
                    pindexvector *bdry2edgeslist,   /* in / out */
                    const int nBdry);             /* in       */

void prolongation(pcrealvector x,                   /* in  */
                  prealvector y,                    /* out */
                  pcindexmatrix edgeno);            /* in  */

void restriction(pcrealvector x,                    /* in  */
                 prealvector y,                     /* out */
                 pcindexmatrix edgeno);             /* in  */

void stima_laplace(real p1[2],                      /* in  */ 
                   real p2[2],                      /* in  */ 
                   real p3[2],                      /* in  */ 
                   real m[3][3]);                   /* out */ 

void buildStiffness(prealmatrix coordinates,        /* in  */
                    pindexmatrix elements,          /* in  */
                    pcoo S);                        /* out */

void rhs_laplace(real p1[2],                      /* in  */ 
                   real p2[2],                      /* in  */ 
                   real p3[2],                      /* in  */
                   real (*f)(real*),                /* in  */
                   real m[3]);                      /* out */ 

void buildRhs(prealmatrix coordinates,              /* in  */
              pindexmatrix elements,                /* in  */
              real (*f)(real*),                     /* in  */
              prealvector b);                       /* out  */

void getFixed(const index     nCoord,               /* in  */
              pcindexvector   bdrytyp,              /* in  */
              pindexmatrix   *bdrylist,             /* in  */
              pindexvector    fixedNodes);          /* out */

void setDirichletData2Rhs(prealmatrix coordinates,  /* in  */
                          pindexvector fixedNodes,  /* in  */
                          real (*f)(real*),         /* in  */
                          prealvector rhs);         /* in/out */

void buildStiffnessByInterpolation(pcrs    Ain,     /* in  */
                              pcrs         Aout,    /* out */
                              pindexmatrix s2p);    /* in  */

int
create_hierarchy(const int       nLevel,            /* in       */
                 prealmatrix     coordinates,       /* in / out */
                 pindexmatrix    elements,          /* in / out */
                 pindexvector    material,          /* in / out */
                 pindexmatrix    elements2edges,    /* in       */
                 pindexmatrix   *f2s,               /* out      */
                 const int       nBdry,             /* in       */
                 pcindexvector   bdrytyp,           /* in       */
                 pindexmatrix   *bdrylist,          /* in / out */
                 pindexvector   *bdry2edgeslist,    /* in / out */
                 pcrs           *A,                 /* in / out */
                 pindexvector   *fixed);            /* out      */

#endif
