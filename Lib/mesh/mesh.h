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

/*! Uniform (RGB) mesh refinement. */
void refine_uniform(prealmatrix  coordinates,        /* in / out */
                    pindexmatrix elements,          /* in / out */
                    prealvector  material,          /* in / out */
                    pindexmatrix elements2edges,    /* in / out */
                    pindexmatrix edgeno,            /* out      */
                    pindexmatrix *bdrylist,         /* in / out */
                    pindexvector *bdry2edgeslist,   /* in / out */
                    const int nBdry);             /* in       */

/*! Prolongation operator, required for multilevel methods. */
void prolongation(pcrealvector x,                   /* in  */
                  prealvector y,                    /* out */
                  pcindexmatrix edgeno);            /* in  */

/*! Restriction operator, required for multilevel methods. */
void restriction(pcrealvector x,                    /* in  */
                 prealvector y,                     /* out */
                 pcindexmatrix edgeno);             /* in  */

/*! Assemble local stiffness matrix. */
void stima_laplace(real p1[2],                      /* in  */ 
                   real p2[2],                      /* in  */ 
                   real p3[2],                      /* in  */ 
                   real m[3][3]);                   /* out */ 

/*! Assemble global stiffness matrix. */
void buildStiffness(prealmatrix coordinates,        /* in  */
                    pindexmatrix elements,          /* in  */
                    pcoo S);                        /* out */

/*! Assemble local right hand side force contribution. */
void
rhs_laplace(real  p1[2],                /* in  */
            real  p2[2],                /* in  */
            real  p3[2],                /* in  */
            real  mat,
            real  (*f1)(real[2], real), /* in  */
            real* (*f2)(real[2], real), /* in  */
            real  m[3]);                /* out */

/*! Assemble local right hand side neumann contribution. */
void
rhs_neumann(real  p1[2],                /* in  */
            real  p2[2],                /* in  */
            real  mat,                  /* in  */
            real  (*g)(real[2], real),  /* in  */
            real* (*f2)(real[2], real), /* in  */
            real  m[2]);                /* out */

/*! Assemble global right hand side. */
void
buildRhs(pcrealmatrix        coordinates, /* in  */
         pcindexmatrix       elements,    /* in  */
         const pindexmatrix  *bdrylist,   /* in  */
         pcrealvector  material,          /* in  */
         real  (*f1)(real[2], real),      /* in  */
         real* (*f2)(real[2], real),      /* in  */
         real  (*g)(real[2], real),       /* in  */
         prealvector b);                  /* out  */

/*! Get Dirichlet nodes. */
void getFixed(const index     nCoord,               /* in  */
              pcindexvector   bdrytyp,              /* in  */
              pindexmatrix   *bdrylist,             /* in  */
              pindexvector    fixedNodes);          /* out */

/*! Set Dirichlet data for solution. */
void setDirichletData2Rhs(prealmatrix coordinates,  /* in  */
                          pindexvector fixedNodes,  /* in  */
                          real (*f)(real*),         /* in  */
                          prealvector rhs);         /* in/out */

/*! [TODO] Build stiffness matrix by interpolation. */
void buildStiffnessByInterpolation(pcrs Ain,          /* in  */
                                   pcrs Aout,         /* out */
                                   pindexmatrix s2p); /* in  */

/*! Create hierarchy for multilevel methods. */
int
create_hierarchy(const int       nLevel,            /* in       */
                 prealmatrix     coordinates,       /* in / out */
                 pindexmatrix    elements,          /* in / out */
                 prealvector     material,          /* in / out */
                 pindexmatrix    elements2edges,    /* in       */
                 pindexmatrix   *f2s,               /* out      */
                 const int       nBdry,             /* in       */
                 pcindexvector   bdrytyp,           /* in       */
                 pindexmatrix   *bdrylist,          /* in / out */
                 pindexvector   *bdry2edgeslist,    /* in / out */
                 pcrs           *A,                 /* in / out */
                 pindexvector   *fixed);            /* out      */

#endif
