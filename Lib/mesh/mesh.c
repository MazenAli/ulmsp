#ifndef MESH_C
#define MESH_C

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>

#include "mesh.h"

void refine_uniform(prealmatrix   coordinates,      /* in / out */
                    pindexmatrix  elements,         /* in / out */
                    prealvector   material,         /* in / out */
                    pindexmatrix  elements2edges,   /* in / out */
                    pindexmatrix  edgeno,           /* out      */
                    pindexmatrix *bdrylist,         /* in / out */
                    pindexvector *bdry2edgeslist,   /* in / out */
                    const int     nBdry){           /* in       */
    /* Allocate local variables */
    int iBdry;
    index *B, *B2e, *nB, *nB2e;
    index nEdges, nTria, nCoord;
    index *E, *E2e, *nE, *nE2e, *nEno;
    real *M, *nM;
    real *C, *nC;
    index i, i0, i1, k, p;
    int isucc[3] = {1,2,0}, iprae[3] = {2,0,1};
    prealmatrix new_coordinates;
    prealvector new_material;
    pindexmatrix new_elements, new_elements2edges, new_bdrylist, new_edgeno;
    pindexvector new_bdry2edgeslist;

    /* Check proper input */
    assert(coordinates);  assert(elements);
    assert(material);     assert(elements2edges);
    assert(edgeno);       assert(bdrylist);
    assert(nBdry);

    /* Set mesh sizes */
    nCoord = coordinates->cols;    /* Number of coordinates */
    nTria  = elements->cols;       /* Number of triangles */

    E2e = elements2edges->vals;
    nEdges = 0;                    /* Number of edges */
    for (i=0 ; i<3*(elements2edges->cols); ++i)
       if ( E2e[i]> nEdges) nEdges = E2e[i];
    nEdges+=1-INDEX_BASE;

    /* Allocate storage for refined mesh */
    new_coordinates = new_realmatrix(2,nCoord+nEdges);
    new_elements = new_indexmatrix(3,4*nTria);
    new_material = new_realvector(4*nTria);
    new_elements2edges = new_indexmatrix(3,4*nTria);
    new_edgeno = new_indexmatrix(2,nEdges);

    /* Declare local variables for convinience */
    C = coordinates->vals; nC = new_coordinates->vals;
    E = elements->vals; nE = new_elements->vals;
    M = material->vals; nM = new_material->vals;
    nE2e = new_elements2edges->vals;
    nEno = new_edgeno->vals;

    /* Get endpoints for each edge, i.e. compute edgeno */
    for (i=0 ; i<nTria ; ++i){
      for (k=0 ; k<3 ; ++k){
        p=3*i+k;
        nEno[2*(E2e[p]-INDEX_BASE)  ] = E[p];
        nEno[2*(E2e[p]-INDEX_BASE)+1] = E[3*i+isucc[k]];
      }
    }

    /* Copy old coordinates */
    for (i=0 ; i<2*nCoord ; i++) nC[i] = C[i];
    /* Compute new coordinates */
    for (i=0 ; i<nEdges ; i++){
      i0 = nEno[2*i]  -INDEX_BASE;
      i1 = nEno[2*i+1]-INDEX_BASE;
      nC[2*nCoord+2*i  ] = 0.5*(C[2*i0  ] + C[2*i1  ]);
      nC[2*nCoord+2*i+1] = 0.5*(C[2*i0+1] + C[2*i1+1]);
    }


    /* Loop over all elements */
    for (i=0 ; i<nTria ; i++){
      for (k=0 ; k<3 ; k++){
        nE[12*i+3*k+0] = E[3*i+k];
        nE[12*i+3*k+1] = nCoord + E2e[3*i+k];
        nE[12*i+3*k+2] = nCoord + E2e[3*i+iprae[k]];

        nE2e[12*i+3*k+0] = 2*E2e[3*i+k]-INDEX_BASE+ (E[3*i+k]>E[3*i+isucc[k]]);
        nE2e[12*i+3*k+1] = 2*nEdges+3*i+k+INDEX_BASE;
        nE2e[12*i+3*k+2] = 2*E2e[3*i+iprae[k]]-INDEX_BASE+ (E[3*i+k]>E[3*i+iprae[k]]);

        nM[4*i+k] = M[i];
      }
      nE[12*i+ 9] = nCoord + E2e[3*i  ];
      nE[12*i+10] = nCoord + E2e[3*i+1];
      nE[12*i+11] = nCoord + E2e[3*i+2];

      nE2e[12*i+ 9] = 2*nEdges+3*i+1+INDEX_BASE;
      nE2e[12*i+10] = 2*nEdges+3*i+2+INDEX_BASE;
      nE2e[12*i+11] = 2*nEdges+3*i+0+INDEX_BASE;

      nM[4*i+3] = M[i];
    }

    /* Refine boundary */
    new_bdrylist = new_indexmatrix(0,0);
    new_bdry2edgeslist = new_indexvector(0);
    /* Loop over all boundary pieces */
    for (iBdry=0 ; iBdry<nBdry ; iBdry++){
      resize_indexmatrix(new_bdrylist,2,2*bdrylist[iBdry]->cols);
      resize_indexvector(new_bdry2edgeslist,2*bdrylist[iBdry]->cols);
      B = bdrylist[iBdry]->vals; B2e = bdry2edgeslist[iBdry]->vals;
      nB = new_bdrylist->vals; nB2e = new_bdry2edgeslist->vals;
      for (i=0 ; i<bdrylist[iBdry]->cols ; i++){
        /* Refine bdrylist[iBdry]*/
        nB[4*i  ] = B[2*i];
        nB[4*i+1] = nCoord+B2e[i];
        nB[4*i+2] = nCoord+B2e[i];
        nB[4*i+3] = B[2*i+1];

        /* Refine bdry2edgeslist[iBdry]*/
        nB2e[2*i  ] = 2*B2e[i  ] + (B[2*i  ] > B[2*i+1]) - INDEX_BASE;
        nB2e[2*i+1] = 2*B2e[i  ] + (B[2*i+1] > B[2*i  ]) - INDEX_BASE;
      }
      swap_indexmatrix(bdrylist[iBdry], new_bdrylist);
      swap_indexvector(bdry2edgeslist[iBdry], new_bdry2edgeslist);
    }

    swap_realmatrix(coordinates, new_coordinates);
    swap_realvector(material, new_material);
    swap_indexmatrix(edgeno, new_edgeno);
    swap_indexmatrix(elements, new_elements);
    swap_indexmatrix(elements2edges, new_elements2edges);

    del_realmatrix(new_coordinates);
    del_realvector(new_material);
    del_indexmatrix(new_edgeno);
    del_indexmatrix(new_elements);
    del_indexmatrix(new_elements2edges);

    del_indexmatrix(new_bdrylist);
    del_indexvector(new_bdry2edgeslist);
}

void prolongation(pcrealvector x,        /* in  */
                  prealvector y,         /* out */
                  pcindexmatrix edgeno){ /* in  */
  index i, ndim, *ptr;
  real *in, *out;

  ndim = x->length;
  in = x->vals; out = y->vals; ptr = edgeno->vals;

  for (i=0; i<ndim; i++) out[i] = in[i];
  for (i=0; i<edgeno->cols; i++){									//edgeno->cols = "number of sons"
    out[ndim+i] = 0.5 * (in[ ptr[2*i  ] - INDEX_BASE]
                       + in[ ptr[2*i+1] - INDEX_BASE]);
  }
  return;
}

void restriction(pcrealvector x,        /* in  */
                 prealvector y,         /* out */
                 pcindexmatrix edgeno){ /* in  */
  index i, ndim, *ptr;
  real *in, *out;

  ndim = y->length;
  in = x->vals; out = y->vals; ptr = edgeno->vals;

  for (i=0; i<ndim; i++) out[i] = in[i];
  for (i=0; i<edgeno->cols; i++){
    out[ ptr[2*i  ]-INDEX_BASE] += 0.5*in[ndim+i];
    out[ ptr[2*i+1]-INDEX_BASE] += 0.5*in[ndim+i];
  }
  return;
}

void stima_laplace(real p1[2],    /* in  */
                   real p2[2],    /* in  */
                   real p3[2],    /* in  */
                   real m[3][3])  /* out */
{
  int i, j;
  real d[3][2], fac;

  for (i = 0 ; i < 2 ; ++i){
     d[0][i] = p3[i]-p2[i];
     d[1][i] = p1[i]-p3[i];
     d[2][i] = p2[i]-p1[i];
  }
  fac = 1./(2.*(d[1][0]*d[2][1]-d[2][0]*d[1][1]));
  for ( i = 0 ; i < 3 ; ++i){
    for ( j = 0 ; j < i ; ++j){
      m[i][j] = fac * (d[i][0]*d[j][0] + d[i][1]*d[j][1]);
      m[j][i] = m[i][j];
    }
    m[i][i] = fac * (d[i][0]*d[i][0] + d[i][1]*d[i][1]);
  }
}

void buildStiffness(prealmatrix coordinates,        /* in  */
                    pindexmatrix elements,          /* in  */
                    pcoo S)                         /* out  */
{
  index i, j, k, nz, nT, nC, *E, *Si, *Sj;
  real M[3][3], *C, *Sx;

  nT = elements->cols;
  nC = coordinates->cols;

  resize_coo(S, 9*nT, nC, nC);
  Sx = S->vals; Si = S->rows; Sj = S->cols;
  C = coordinates->vals; E = elements->vals;

  nz = 0;
  for (k=0; k < nT; ++k){
    stima_laplace((real*)(C+2*(E[3*k+0]-INDEX_BASE)),
                  (real*)(C+2*(E[3*k+1]-INDEX_BASE)),
                  (real*)(C+2*(E[3*k+2]-INDEX_BASE)),M);
    for (i=0; i<3; i++){
      for (j=0; j<3; j++){
        Sx[nz] = M[i][j];
        Si[nz] = E[3*k+i];
        Sj[nz] = E[3*k+j];
        ++nz;
      }
    }
  }
}


void
rhs_laplace(real  p1[2],                /* in  */
            real  p2[2],                /* in  */
            real  p3[2],                /* in  */
            real  mat,                  /* in */
            real  (*f1)(real[2], real), /* in  */
            real* (*f2)(real[2], real), /* in  */
            real  m[3])                 /* out */
{
    int i;
    real midpoint[2], d[2][2], facf1, facf2, *_f2;

    for (i = 0 ; i < 2 ; ++i){
        d[0][i] = p1[i]-p3[i];
        d[1][i] = p2[i]-p1[i];
    }
    midpoint[0] = (p1[0] + p2[0] + p3[0])/3.;
    midpoint[1] = (p1[1] + p2[1] + p3[1])/3.;
    facf1 = f1(midpoint, mat)/6.*(d[0][0]*d[1][1]-d[1][0]*d[0][1]);
    for ( i = 0 ; i < 3 ; ++i){
        m[i]  = facf1;
    }

    _f2 = f2(midpoint, mat);

    facf2  = (p2[1]-p3[1])*_f2[0];
    facf2 += (p3[0]-p2[0])*_f2[1];
    m[0]  -= 0.5*facf2;

    facf2  = (p3[1]-p1[1])*_f2[0];
    facf2 += (p1[0]-p3[0])*_f2[1];
    m[1]  -= 0.5*facf2;

    facf2  = (p1[1]-p2[1])*_f2[0];
    facf2 += (p2[0]-p1[0])*_f2[1];
    m[2]  -= 0.5*facf2;
}


void
rhs_neumann(real  p1[2],                /* in  */
            real  p2[2],                /* in  */
            real  mat,                  /* in  */
            real  (*g)(real[2], real),  /* in  */
            real* (*f2)(real[2], real), /* in  */
            real  m[2])                 /* out */
{
    real midpoint[2], normal[2], E, *_f2;

    midpoint[0] = (p1[0] + p2[0])/2.;
    midpoint[1] = (p1[1] + p2[1])/2.;

    normal[0] = p2[1]-p1[1];
    normal[1] = p1[0]-p2[0];

    E = normal[0]*normal[0]+normal[1]*normal[1];
    E = sqrt(E);

    normal[0] /= E;
    normal[1] /= E;

    _f2 = f2(midpoint, mat);

    m[0]  = _f2[0]*normal[0]+_f2[1]*normal[1];
    m[0] += g(midpoint, mat);
    m[0] *= 0.5*E;
    m[1]  = m[0];
}


void
buildRhs(pcrealmatrix  coordinates,    /* in  */
         pcindexmatrix elements,       /* in  */
         const pindexmatrix *bdrylist, /* in  */
         pcrealvector  material,       /* in  */
         real  (*f1)(real[2], real),   /* in  */
         real* (*f2)(real[2], real),   /* in  */
         real  (*g)(real[2], real),    /* in  */
         prealvector b)                /* out  */
{
    index k, nT, nC, idx[3], *E;
    real M[3], *C, *bx;
    pcindexmatrix neumann = bdrylist[1];

    nC = coordinates->cols;
    nT = elements->cols;

    resize_realvector(b, nC);

    bx  = b->vals; C = coordinates->vals; E = elements->vals;

    for (k=0; k < nT; ++k){
        idx[0] = E[3*k+0]-INDEX_BASE;
        idx[1] = E[3*k+1]-INDEX_BASE;
        idx[2] = E[3*k+2]-INDEX_BASE;
        rhs_laplace((real*)(C+2*idx[0]),
                    (real*)(C+2*idx[1]),
                    (real*)(C+2*idx[2]),
                    material->vals[k],
                    f1, f2,
                    M);
        bx[idx[0]] += M[0];
        bx[idx[1]] += M[1];
        bx[idx[2]] += M[2];
    }

    for (k=0; k<neumann->cols; ++k) {
        idx[0] = neumann->vals[2*k+0]-INDEX_BASE;
        idx[1] = neumann->vals[2*k+1]-INDEX_BASE;

        rhs_neumann((real*)(C+2*idx[0]),
                    (real*)(C+2*idx[1]),
                    material->vals[k],
                    g, f2,
                    M);

        bx[idx[0]] += M[0];
        bx[idx[1]] += M[1];
    }
}


/*
 * Collect all nodes of typ 0 in fixedNodes
 */
void getFixed(const index     nCoord,     /* in       */
              pcindexvector   bdrytyp,    /* in       */
              pindexmatrix   *bdrylist,   /* in       */
              pindexvector    fixedNodes) /* out      */
{
  int i, nBdry;
  index j, nz, cnt;
  bool *flag;

  nBdry = (int) bdrytyp->length;

  flag = (bool*) calloc(nCoord,sizeof(bool));
  cnt = 0;
  for (i=0 ; i<nBdry ; i++){
    if (!bdrytyp->vals[i]){
      cnt += (bdrylist[i]->cols)+1;
      for (j=0 ; j<2*bdrylist[i]->cols ; j++){
        flag[bdrylist[i]->vals[j]-INDEX_BASE] = 1;
      }
    }
  }

  resize_indexvector(fixedNodes,cnt);
  nz = 0;
  for (j=0; j<nCoord; j++) {
    if (flag[j]) 
      fixedNodes->vals[nz++] = j+INDEX_BASE;
  }


  resize_indexvector(fixedNodes,nz);
  if (!flag){ 
	free(flag);
  	return;
	}

	free(flag);
}

void setDirichletData2Rhs(prealmatrix coordinates,
                          pindexvector fixedNodes,
                          real (*f)(real*),
                          prealvector rhs)
{
  index i, ptr;
  for ( i=0 ; i < fixedNodes->length ; i++ ) {
    ptr = fixedNodes->vals[i]-INDEX_BASE;
    rhs->vals[ptr] = f(&(coordinates->vals[2*ptr]));
  }
}


void buildStiffnessByInterpolation(pcrs Ain, pcrs Aout, pindexmatrix s2p)
{
  index i, j, nS, nH, nE, nz;
  pcoo AxP_coo;
  pcoo S;
  pcrs AxP_crs;

  nH = Ain->numc;    /* Huge dimension */
  nE = (s2p->cols);  /* Number of son -> parents relations */
  nS = nH - nE;      /* Small dimension */

  /* Let interpolation operator be P = [P_c, P_f] 
   * and A = [A_cc, A_cf / A_fc, A_ff]
   *
   * Notice, A is expected to be symmetric and P_c is the identity matrix.
   */

  /* Compute [A_cf / A_ff] * P_f */
  AxP_coo = new_coo( 2 * ( Ain->rowptr[nH] - Ain->rowptr[nS] ), nH, nS);
  nz = 0;
  for (i=0; i< nE; i++){
    for(j = Ain->rowptr[nS+i]; j < Ain->rowptr[nS+i+1]; j++){
      AxP_coo->vals[nz] = 0.5 * Ain->vals[j];
      AxP_coo->cols[nz] = s2p->vals[2*i  ];
      AxP_coo->rows[nz] = Ain->colind[j];
      nz++;
      AxP_coo->vals[nz] = 0.5 * Ain->vals[j];
      AxP_coo->cols[nz] = s2p->vals[2*i+1  ];
      AxP_coo->rows[nz] = Ain->colind[j];
      nz++;
    }
  }

  /* Convert to crs format */
  AxP_crs= new_crs(1,1,1);
  init_coo2crs(AxP_crs, AxP_coo);
  del_coo(AxP_coo);

  /* Allocate memory to store P^T * A * P in coo format */
  nz = Ain->rowptr[nS] + 2 * AxP_crs->rowptr[nH];
  S = new_coo(Ain->rowptr[nS],nS,nS);
  /* Copy entries of A_cc to S */
  nz = 0;
  for (i=0; i< nS; i++){
    for(j = Ain->rowptr[i]; j < Ain->rowptr[i+1]; ++j){
      if (Ain->colind[j]-INDEX_BASE<nS){
        S->vals[nz] = Ain->vals[j];
        S->rows[nz] = i+INDEX_BASE;
        S->cols[nz] = Ain->colind[j];
        nz++;
      }
    }
  }
  resize_coo(S,nz+ 2 * AxP_crs->rowptr[nH],nS,nS);

  /* Copy entries of A_cf * P_f and (A_cf * P_f)^T to S */
  for (i=0; i< nS; i++){
    for(j = AxP_crs->rowptr[i] ; j < AxP_crs->rowptr[i+1]; j++){
      S->vals[nz] = AxP_crs->vals[j];
      S->cols[nz] = i+INDEX_BASE;
      S->rows[nz] = AxP_crs->colind[j];
      nz++;
      S->vals[nz] = AxP_crs->vals[j];
      S->cols[nz] = AxP_crs->colind[j];
      S->rows[nz] = i+INDEX_BASE;
      nz++;
    }
  }

  /* Compute P_f^T * A_cf * P_f and store to S */
  for (i=0; i< nE; i++){
    for( j = AxP_crs->rowptr[nS+i] ; j < AxP_crs->rowptr[nS+i+1] ; j++){
      S->vals[nz] = 0.5 * AxP_crs->vals[j];
      S->rows[nz] = s2p->vals[2*i  ];
      S->cols[nz] = AxP_crs->colind[j];
      nz++;
      S->vals[nz] = 0.5 * AxP_crs->vals[j];
      S->rows[nz] = s2p->vals[2*i+1];
      S->cols[nz] = AxP_crs->colind[j];
      nz++;
    }
  }

  /* Convert to crs format */
  del_crs(AxP_crs);
  init_coo2crs(Aout, S);
  del_coo(S);

  return;
}


int create_hierarchy(const int       nLevel,          /* in       */
                     prealmatrix     coordinates,     /* in / out */
                     pindexmatrix    elements,        /* in / out */
                     prealvector     material,        /* in / out */
                     pindexmatrix    elements2edges,  /* in       */
                     pindexmatrix   *f2s,             /* out      */
                     const int       nBdry,           /* in       */
                     pcindexvector   bdrytyp,         /* in       */
                     pindexmatrix   *bdrylist,        /* in / out */
                     pindexvector   *bdry2edgeslist,  /* in / out */
                     pcrs           *A,               /* in / out */
                     pindexvector   *fixed)           /* out      */
{
  int i;
  pcoo S = new_coo(1,1,1);
  pindexmatrix edgeno;

  edgeno = new_indexmatrix(1,1);

  /* Store fixed for coarsest mesh*/
  fixed[nLevel-1] = new_indexvector(0);
  getFixed(coordinates->cols, bdrytyp, bdrylist, fixed[nLevel-1]);
  
  /* Refine mesh (nLevel-1) times */
  for (i=nLevel-2; i>=0; i--){  
    refine_uniform(coordinates, elements, material,
                   elements2edges, edgeno, bdrylist,
                   bdry2edgeslist, nBdry);
 	fixed[i] = new_indexvector(1);
    getFixed(coordinates->cols, bdrytyp, bdrylist, fixed[i]);
    f2s[i] = new_indexmatrix(0,0);
    copy_indexmatrix(f2s[i],edgeno);
  }


    /* Build stiffness matrix for finest mesh */
    buildStiffness(coordinates, elements, S);
    A[0] = new_crs(1,1,1);
    init_coo2crs(A[0], S); 
	del_coo(S);
	
    /* Build stiffness matrix for coarser meshes */
    for (i=1; i<nLevel; i++){
	  A[i] = new_crs(1,1,1);
      buildStiffnessByInterpolation(A[i-1],A[i],f2s[i-1]);
    }
  
	
  del_indexmatrix(edgeno);
  return 1;
}

#endif
