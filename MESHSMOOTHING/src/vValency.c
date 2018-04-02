/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             EGADS Tessellation using wv with Quad Tessellation
 *
 *      Copyright 2011-2018, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include <string.h>
#include <unistd.h>		// usleep
#include <nlopt.h>
#include <time.h>
#ifdef  WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <winsock2.h>
#endif
#include "egads.h"
#include "wsserver.h"
int TOTSWAPS     = 0;
int OVERALLSWAPS = 0;
#define EPS  1.E-03
#define DEPS 1.E-03
//#define DEBUG_AREA
int NOUT = 10000;
#define NLOPTMAXEVAL 100000

#define DTOL 1.e-9
#define TRIPLEOPTI 1
#define PI              3.1415926535897931159979635
#define DEBUG
int GLO = 0 ;
int QUADTESS  = 0;
int OPTCOUNT  = 0;
int TESSCOUNT = 0;
#define ANGPASS 0.3
#define ERRCTT  3.554147
#define  FRONTFACE 0
int SCOUNT = 0;
#define FORCETRIANGULATION
#define MAX_LINKS 100
#define FORCEQUADS
/* structure to hold on to the EGADS triangulation per Body */
typedef struct {
  ego    *faces;
  ego    *edges;
  ego    body;
  ego    tess;
  int    mtype;
  int    nfaces;
  int    nedges;
  int    plen;
  double *uvs;      /* face #1 updated uvs */
} bodyData;


/* globals used in these functions */
static wvContext *cntxt;
static bodyData  *bodydata;
static int       nbody;
static float     focus[4];

typedef struct{
  int *verts;
  int  nV; // nV = origin + peaks
  int *quads;
  int  nQ;
} vStar;

typedef struct {
  int     *quadIdx, *quadAdj, **valence, *vType;
  int      nQ, nS, totVerts;
  ego      face;
  double  *xyzs, *uv;
  vStar   **star;
} quadMap;

void NLOPTtermination(int n) {
  switch(n) { 
  case 1: printf(" NLOPT SUCCEESS!\n");
  break;
  case 2: printf(" NLOPT STOPVALUE REACHED!\n");
  break;
  case 3: printf(" NLOPT ABS OR RELATIVE FUNCTION TOL REACHED!\n");
  break;
  case 4: printf(" NLOPT ABS OR RELATIVE XTOL REACHED!\n");
  break;
  case 5: printf(" NLOPT MAX EVALUATIONS REACHED!\n");
  break;
  default: 
    printf(" NLOPT FAILED!\n");
    break;
  }
}


#ifdef DEBUG

int FCOUNT = 1;
static void
print_mesh(quadMap *qm) {
  int i, vID, k;
  FILE *fout;
  char buffer[33];
  snprintf(buffer, sizeof(char) * 32, "MESH%i.txt", FCOUNT);
  printf(" WRITING ON FILE %s\n",buffer);
  ++FCOUNT;
  fout = fopen(buffer,"w");
  if (fout == NULL ) return;
  for ( i = 0 ; i < qm ->nQ; ++i) {
    if ( qm ->quadIdx[4*i] != -1 ) {
      //printf("====== QUAD %d =======\n",i + 1);
      for ( k = 0; k < 4; ++k) {
        vID = qm->quadIdx[4*i + k] - 1;
        //printf(" V ID %d WITH COORDS %lf %lf %lf\n",vID + 1, qm->xyzs[3*vID], qm->xyzs[3*vID +1],qm->xyzs[3*vID + 2]);
        fprintf(fout, "%lf  %lf  %lf \n",qm->xyzs[3*(qm->quadIdx[4*i + k] - 1)], qm->xyzs[3*(qm->quadIdx[4*i + k] - 1) +1], qm->xyzs[3*(qm->quadIdx[4*i + k] - 1) + 2]);
      }
      fprintf(fout, "%lf  %lf  %lf \n",qm->xyzs[3*(qm->quadIdx[4*i] - 1)], qm->xyzs[3*(qm->quadIdx[4*i] - 1) +1], qm->xyzs[3*(qm->quadIdx[4*i] - 1) + 2]);
      fprintf(fout,"\n");
    }
  }
  fclose(fout);
}
static void
printMeshStats(quadMap *qm, int sweep) {
  int fix, move, i,len ;
  int intVal[MAX_LINKS], boundVal[MAX_LINKS];
  FILE *fout;
  char buffer[33];
  snprintf(buffer, sizeof(char) * 32, "MESH_STATS_%d.txt", sweep);
  printf(" WRITING ON FILE %s\n",buffer);
  fout = fopen(buffer, "w");
  len = qm ->totVerts;
  for ( i = 0; i < MAX_LINKS; ++i) {
    intVal[i]   = 0;
    boundVal[i] = 0;
  }
  for ( i = 0; i < len; ++i) {
    if ( qm -> vType[i] == -1){
      ++intVal[qm ->valence[i][0]];
    } else {
      ++boundVal[qm ->valence[i][0]];
    }
  }
  fprintf(fout,"--------------------- MESH STATS AFTER %d SWEEP --------------\n",sweep);
  fprintf(fout," INTERIOR VERTICES\n");
  for ( i = 0 ; i < MAX_LINKS; ++i) {
    if ( intVal[i]  > 0 ) fprintf(fout," VALENCY %d = %d VERTICES\n", i, intVal[i]);
  }
  fprintf(fout," BOUNDARY VERTICES\n");
  for ( i = 0 ; i < MAX_LINKS; ++i) {
    if ( boundVal[i]  > 0 ) fprintf(fout," VALENCY %d = %d VERTICES\n", i, boundVal[i]);
  }
  fclose(fout);
  return ;
}
static void
print_quad_specs(quadMap *qm, int id) {
  --id;
  if (id < 0 ) return ;
  int i;
  printf(" QUAD %d HAS VERTICES ",id +1);
  for ( i = 0 ; i < 4; ++i ) printf(" %d ",qm -> quadIdx[4*id + i]);
  printf("\t AND ADJACENT QUADS ");
  for ( i = 0 ; i < 4; ++i ) printf(" %d ",qm -> quadAdj[4*id + i]);
  printf("\n");
}
static void
print_star(vStar *star) {

  int i, j;
  printf(" ===== STAR CENTRED AT VERTEX %d HAS  %d VALENCIES\n", star->verts[0], star->verts[1]);
  printf(" ROSA  DE %d PICOS \n",star->nV);
  for ( i = 1 ; i < star->nV; ++ i) {
    if ( i%2 == 0)
      printf(" OPP  PICO %d CON VALENCIA %d\n",star->verts[2*i], star->verts[2*i+1]);
    else
      printf(" LINK  PICO %d CON VALENCIA %d\n",star->verts[2*i], star->verts[2*i+1]);
  }
  for ( i = 0 ; i < star ->nQ; ++i) printf(" QUAD %d = %d\n",i,star->quads[i]);
}
static void
print_starFile(vStar *star, quadMap *qm) {
  int i, vID, k;
  FILE *fout;
  char buffer[33];
  snprintf(buffer, sizeof(char) * 32, "STAR%i_SCOUNT_%i_OPTIN_%i.txt", FCOUNT, star->verts[0], OPTCOUNT);
  printf(" ------WRITING ON FILE %s\n",buffer);
  ++SCOUNT;
  fout = fopen(buffer,"w");
  if (fout == NULL ) return;
  for ( i = 0 ; i < star->nQ; ++i) {
    fprintf(fout, "%lf  %lf  %lf \n",qm->xyzs[3*(star->verts[0] - 1)],qm->xyzs[3*(star->verts[0] - 1) + 1],qm->xyzs[3*(star->verts[0] - 1) + 2]);
    for ( k = 1 ; k <= 3; ++k){
      vID =  (2*i + k);
      if ( vID >= star->nV) vID = 1;
      vID = star->verts[2*vID] - 1;
      fprintf(fout, "%lf  %lf  %lf %d  %d\n", qm->xyzs[3*vID],qm->xyzs[3*vID + 1],qm->xyzs[3*vID + 2],(k+2*i)%star->nV,vID);
    }
    fprintf(fout, "%lf  %lf  %lf \n",qm->xyzs[3*(star->verts[0] - 1)],qm->xyzs[3*(star->verts[0] - 1) + 1],qm->xyzs[3*(star->verts[0] - 1) + 2]);
    fprintf(fout,"\n");
  }
  for ( i = 0 ; i < star->nV; ++i) {
    vID = star->verts[2*i];
    fprintf(fout, "# VERT %d = %lf  %lf  %lf\n",star->verts[2*i], qm ->xyzs[3*(vID - 1)], qm ->xyzs[3*(vID - 1) +1], qm ->xyzs[3*(vID - 1) + 2]);
  }
  fclose(fout);
}

#endif
static int EG_buildStar(quadMap *qM, vStar **star, int vID );
static double star_algebraic_area(quadMap *qm, vStar *star) ;
double NLOPT_objEqualAnglesTess(unsigned n, const double *uvs, double *grad, void *inQM);
static int optimize_angles(quadMap *qm, int nV, int *vID);
static int quad_algebraic_area(quadMap *qm, int qID, int vID,  int *validArea, double *intAngle ) ;
static int checkInvalidElement(quadMap *qm, int qID, int v0ID,int pullVertex) {
  int i, j, k, v[4], stat, aux, v1ID = -1, adjQ = -1, area, it = 0;
  double dir[3], dirOpp[3], uvo[2], uvOri[2], vXYZ[18], penalty, angle;
  int optiAng[1];
  printf(" CHECK INVALID QUAD FOR VERTEX %d IN QUAD %d\n", v0ID, qID);
  if ( v0ID == -1  || qID == -1) return EGADS_SUCCESS;
  print_quad_specs(qm, qID);
  vStar *star = NULL;
  stat = EG_buildStar(qm, &star, v0ID -1);
  if ( stat != EGADS_SUCCESS ) {
    EG_free(star);
    return EGADS_GEOMERR;
  }
  if ( qm -> vType[v0ID -1 ] >= 0 ) {
    stat = quad_algebraic_area(qm, qID, v0ID, &area, &angle);
    EG_free(star);
    if ( area == 1) return EGADS_SUCCESS;
    else return EGADS_GEOMERR;
  }
  printf(" BUILD STAR FOR VERTEX %d   STAT %d\n",v0ID, stat);
  penalty = star_algebraic_area(qm , star);
  printf(" FOR VERTEX %d  PENALTY %lf\n",v0ID, penalty);
  if (pullVertex == 0 || qm -> vType[v0ID] != -1) {
    EG_free(star);
    if ( penalty >= 1 ) return EGADS_GEOMERR;
    else             return EGADS_SUCCESS;
  }
  // INTERIOR VERTEX: MOVE IT AROUND
  // GET ADJACENT QUAD AND APPROPRIATE EDGE TO PULL FROM
  for ( k = 0 ; k < 4; ++k)
    if (qm -> quadIdx[4*(qID -1) + k] == v0ID) break;
  if ( qm -> quadAdj[4*(qID -1) + k] != -1)
    adjQ = qm ->quadAdj[4*(qID -1) +  k];
  else {
    printf(" I CAN't FIND ANY NON VOID ADJACENT QUAD FOR VERTEX %d\n",v0ID + 1);
    print_quad_specs(qm, qID);
    EG_free(star);
    return EGADS_GEOMERR;
  }
  print_quad_specs(qm, adjQ);
  for ( k = 0 ; k < 4; ++k) {
    if ( qm->quadIdx[4*(adjQ - 1) + k  ] == v0ID ) {
      v1ID = qm->quadIdx[4*(adjQ - 1) + (int)(k +2)%4  ];
      break;
    }
  }
  printf(" USING FOR %d AS CENTRE %d  %d\n",v0ID, k + 2, v1ID);
  dir  [0] = 1.e-04*(qm -> uv[2*(v1ID - 1)    ] - qm -> uv[2*(v0ID - 1)    ]);
  dir  [1] = 1.e-04*(qm -> uv[2*(v1ID - 1) + 1] - qm -> uv[2*(v0ID - 1) + 1]);
  uvo  [0] = qm -> uv[2*(v0ID - 1)    ];
  uvo  [1] = qm -> uv[2*(v0ID - 1) + 1];
  uvOri[0] = qm -> uv[2*(v0ID - 1)    ];
  uvOri[1] = qm -> uv[2*(v0ID - 1) + 1];
  optiAng[0] = v0ID - 1;
  printf(" ORI VERTEX %lf  %lf  PUSHED TOWARDS %lf  %lf \n",qm -> xyzs[3*(v0ID - 1)], qm -> xyzs[3*(v0ID - 1) + 1], 
      qm -> xyzs[3*(v1ID - 1)],qm -> xyzs[3*(v1ID - 1) +1]);
  printf(" SPEED %lf  %lf\n", dir[0], dir[1]);

  while ( penalty >= 1.0 )
  {
    stat    = optimize_angles(qm, 1, optiAng);
    uvo[0] += dir[0];
    uvo[1] += dir[1];
    printf(" Vo %lf %lf  -> %lf %lf\t",uvOri[0], uvOri[1], uvo[0], uvo[1]);
    qm -> uv[2*(v0ID - 1)    ] = uvo[0];
    qm -> uv[2*(v0ID - 1) + 1] = uvo[1];
    stat                       = EG_evaluate(qm -> face, &qm -> uv[2*(v0ID - 1)], vXYZ);
    printf(" POSITION %lf %lf\n",qm -> xyzs[3*(v0ID - 1)], qm -> xyzs[3*(v0ID - 1) +1]);
    if ( stat != EGADS_SUCCESS) {
      printf(" IN FUNCTION CHECKINVALIDELEMENT  EG_evaluate %d\n",stat);
      EG_free(star);
      return stat;
    }
    qm -> xyzs[3*(v0ID - 1)    ] = vXYZ[0];
    qm -> xyzs[3*(v0ID - 1) + 1] = vXYZ[1];
    qm -> xyzs[3*(v0ID - 1) + 2] = vXYZ[2];
    for ( i = 0 ; i < 2; ++ i)
      dirOpp[i] = qm -> uv[2*(v0ID - 1) + i] - qm -> uv[2*(v0ID - 1) + i];
    if ( (dirOpp[0]*dir[0] + dirOpp[1]*dir[1]) < 0) goto reset;
    penalty = star_algebraic_area(qm , star);
    printf(" STAT %d  PENALTY %lf  \n",stat, penalty);
    if ( it %100 == 0 ) print_mesh(qm);
    if ( it > 10000 || stat != EGADS_SUCCESS) exit (1);//goto reset;
    ++it;
  }
  printf(" GOOD! LEAVING WITH PENALTY %lf \n",penalty);
  EG_free(star);
  return EGADS_SUCCESS;
  reset:
  qm -> uv[2*(v0ID - 1)    ] = uvOri[0];
  qm -> uv[2*(v0ID - 1) + 1] = uvOri[1];
  stat                       = EG_evaluate(qm -> face, &qm -> uv[2*(v0ID - 1)], vXYZ);
  if ( stat != EGADS_SUCCESS) {
    printf(" IN FUNCTION CHECKINVALIDELEMENT  EG_evaluate %d\n",stat);
    return stat;
  }
  qm -> xyzs[3*(v0ID - 1)    ] = vXYZ[0];
  qm -> xyzs[3*(v0ID - 1) + 1] = vXYZ[1];
  qm -> xyzs[3*(v0ID - 1) + 2] = vXYZ[2];
  printf(" BAD! LEAVING WITH PENALTY %lf \n",penalty);
  EG_free(star);
  return EGADS_GEOMERR;

}


/*


static int
quadInteriorAngle(quadMap *qm, int qID, int vID, double *theta) {
  double alpha[2], dot[2], norm[2], vA[3], vB[3], vC[3];
  int vPos = -1, i, qV[4], k;
  for ( i = 0 ; i < 4; ++i) {
    if (qm ->quadIdx[4*(qID - 1) + i] == vID ) {
      vPos = i;
      break;
    }
  }
  if ( vPos == -1 ) {
    printf(" VERTEX %d is STAR CENTRE AND SHOULD BELONG TO QUAD %d\n",vID, qID);
    return EGADS_INDEXERR;
  }
  for ( i = 0 ; i < 4; ++i) {
    qV[i] = qm ->quadIdx[4*(qID - 1) + (i+vPos)%4] - 1;
  }
  //print_quad_specs(qm, qID);

  // GET INTERNAL ANGLE AS A SUM OF TWO ANGLES
  for ( k = 1 ; k <= 2; ++k) {
    dot[k - 1] = 0.0;
    norm[0]    = 0.0; norm[1] = 0.0;
    for ( i = 0 ; i < 3; ++i){
      vA[i]     = qm -> xyzs[3*qV[k    ] + i ];
      vB[i]     = qm -> xyzs[3*qV[0    ] + i ];
      vC[i]     = qm -> xyzs[3*qV[k + 1] + i ];
      dot[k-1] += (vA[i] - vB[i])*(vC[i] - vB[i]);
      norm[0]  += (vA[i] - vB[i])*(vA[i] - vB[i]);
      norm[1]  += (vC[i] - vB[i])*(vC[i] - vB[i]);
    }
    norm[0] = sqrt(norm[0]); norm[1] = sqrt(norm[1]);
    if ( norm[0] < EPS || norm[1] < EPS) {
 *theta = 2.0*M_PI;
      return EGADS_SUCCESS;

    }
    else dot[k-1]  /= (norm[0]*norm[1]);
    if      ( fabs(dot[k-1] - 1.0) < DTOL ) {
      printf(" DEGENERATE : IMPOSE 2pi\n");
 *theta = 2*M_PI;
      return EGADS_SUCCESS;
    }
    else if ( fabs(dot[k-1] + 1.0) < DTOL ) alpha[k-1] = M_PI;
    else alpha[k-1] = acos(dot[k-1]);
  }
 *theta = alpha[0] + alpha[1];
  return EGADS_SUCCESS;
}
 */

static double
triang_algebraic_area(double A[], double B[], double C[]){
  int i;
  double AB[3], BC[3], CA[3];
  for ( i = 0 ; i < 3; ++i){
    //  printf(" A = %lf   B = %lf  C = %lf\n",A[i], B[i], C[i]);
    AB[i] = B[i] - A[i];
    BC[i] = C[i] - B[i];
    CA[i] = A[i] - A[i];
  }
  return 0.5*(AB[0]*BC[1] - BC[0]*AB[1] );
}

static int 
computePenalty(quadMap *qm, int qID, int vID, double *penalty, int *validArea ) {
  double alpha, dot, z, zerr;
  int stat;
  stat = quad_algebraic_area(qm, qID, vID, &(*validArea), &alpha);
  if ( stat != EGADS_SUCCESS) {
    printf(" ERROR COMPUTING PENALTY FUNCTION  %d",stat);
    return stat;
  }
  dot = cos(alpha - 0.5*M_PI);
  if ( alpha > 2.0*M_PI - acos(ANGPASS)) dot = -1.0;
  z         = ERRCTT*(dot - ANGPASS*0.5)/ANGPASS;
  zerr      = erfc(z);
  *penalty  = zerr;//exp(zerr);
  if ( validArea == 0)
    *penalty = zerr + exp(ERRCTT*(1 + dot));
#ifdef DEBUG_AREA
  printf("\n APPLiYING PENALTY TO QUAD %d:: z %lf erfc  %lf  penalty %lf  DOT %lf  %lf  NEG AREA %d\n",qID, z,zerr, *penalty, dot, alpha*180.0/M_PI, *validArea);
#endif
  return EGADS_SUCCESS;

}

static int
quad_algebraic_area(quadMap *qm, int qID, int vID,  int *validArea, double *intAngle ) {
  int i, vPos, v, qV[4], k, tri, piv[4], vL, vR, stat, negArea[2];
  double vA[3], vB[3], vC[3],  area[2*2], z, zerr, pen0, alpha, dot, norm[2], theta[2];
#ifdef DEBUG_AREA
  printf("\n-------------------------- LOOKING AT QUAD %d  AND CENTRE  %d -------------\n", qID, vID);
#endif
  for ( i = 0 ; i < 4; ++i) {
    if (qm ->quadIdx[4*(qID - 1) + i] == vID ) {
      vPos = i;
      break;
    }
  }
  if ( vPos == -1 ) {
    printf(" VERTEX %d is STAR CENTRE AND SHOULD BELONG TO QUAD %d\n",vID, qID);
    return EGADS_INDEXERR;
  }
  for ( k = 0 ; k < 2; ++k) { // CHeck additional area to ensure that the quad is correct
    negArea[k] = 0;
    for ( tri = 0 ; tri < 2; ++ tri) {
      piv[0]  = (k  + vPos)%4;
      piv[1]  = (k  + vPos + 1 + tri)%4;
      piv[2]  = (k  + vPos + 2 + tri)%4;
      for ( i = 0 ; i < 3; ++i){
        vA[i] = qm -> xyzs[3*qV[piv[0]] + i ];
        vB[i] = qm -> xyzs[3*qV[piv[1]] + i ];
        vC[i] = qm -> xyzs[3*qV[piv[2]] + i ];
      }
      area[2*k + tri] = triang_algebraic_area(vA, vB, vC);
#ifdef DEBUG_AREA
      printf("\n K %d  TRI %d\n AREA BETWEEN VERTEX %d %d %d \n", k,tri,qV[piv[0]] + 1,qV[piv[1]] + 1,qV[piv[2]] + 1);
      printf(" TRIANG AREA FOR A= %lf %lf B= %lf %lf C= %lf %lf\n", vA[0], vA[1], vB[0], vB[1], vC[0], vC[1] );
      printf(" AREA TRI= %lf \n",area[2*k+tri]);
#endif
      if ( k == 0 ) { // Compute internal angle
        dot = 0.0; norm[0] = 0.0; norm[1] = 0.0;
        for ( i = 0 ; i < 3; ++i) {
          dot     += (vB[i] - vA[i])*(vC[i] - vA[i]);
          norm[0] += (vB[i] - vA[i])*(vB[i] - vA[i]);
          norm[1] += (vC[i] - vA[i])*(vC[i] - vA[i]);
        }
        norm[0] = sqrt(norm[0]); norm[1] = sqrt(norm[1]);
        if ( norm[0] < EPS || norm[1] < EPS) theta[tri] = 2.*M_PI; 
        else dot /= (norm[0]*norm[1]);
        if ( fabs(dot - 1.0) < DTOL )  theta[tri] = M_PI;
        if ( area[2*k + tri ] > 0 ) theta[tri] = 2.0 * M_PI -theta[tri];
      }
    }
    if     ((area[2*k ] * area[2*k + 1]) < 0.0 ) negArea[k] = -1;
    else if(fabs(area[2*k]) < DEPS || fabs(area[2*k + 1]) < DEPS) negArea[k] = 2;
  }
  *intAngle = theta[0] + theta[1];
  if      ( negArea[0] == -1 || negArea[1] == -1) *validArea = -1;
  else if ( negArea[0] == 2  || negArea[1] ==  2) *validArea = 0;
  else *validArea = 1;
  return EGADS_SUCCESS;
}
static double
star_algebraic_area(quadMap *qm, vStar *star) {
  int i, stat, area;
  double penalty, ret = 0.0;;
#ifdef DEBUG_AREA
  printf(" ======= AREA FOR STAR CENTRE %d TOTAL QUADS %d\n",star->verts[0], star -> nQ);
#endif
  for ( i = 0 ; i < star -> nQ; ++i) {
    stat = computePenalty(qm, star -> quads[i],star->verts[0], &penalty, &area);
#ifdef DEBUG_AREA
    printf(" AREA FOR QUAD %d = %d  PENALTY %lf\n",star->quads[i], area, penalty);
#endif
    if ( ret < penalty)
      ret = penalty;
  }
  return ret;
}



static int
EG_buildStar(quadMap *qM, vStar **star, int vID ){
  int i, j, id0 = -1, k = 0, itQ, vROSE[MAX_LINKS], vQUADS[MAX_LINKS];
  int auxm,auxV, auxQ,  r, quad0, quadID, prevQuad, nPeaks, adjQ;
  int qLoop[5] = {0, 1, 2, 3, 0};
  printf(" BUILDING STAR WITH CENTRE %d AND VALENCY  %d  and quad %d\n",vID + 1, qM->valence[vID][0], qM->valence[vID][1] );
  // if (qM -> vType[vID] != -1) {
  //  printf(" THIS IS A VERTEX BELONGING TO THE DOMAIN BOUNDARY CORNER. WE SHOULDN't BE TOUCHING IT !!!\n");
  //    return EGADS_NOTTOPO;
  //}
  vROSE[0]      = vID + 1;
  vROSE[1]      = qM -> valence[vID][0];
  quad0         = qM -> valence[vID][1] - 1;
  quadID        = quad0;
  prevQuad      = -5;
  itQ           =  0;
  vQUADS[itQ++] = quad0 + 1; // Accumulates the quads that we are going to touch
  for ( j = 0 ; j < 4; ++j) {
    if ((qM->quadIdx[4*quad0 + j ] == (vID + 1) )) {
      id0 = j;
      break;
    }
  }
  if ( id0 == -1) {
    printf("I couldn't find the star centre (Vertex %d) in Quad %d\n. quadMap build UNSUCCESSFUL. Check updateUVs function!\n", vROSE[0], quad0);
    return (EGADS_INDEXERR);
  }
  //print_quad_specs(qM, quad0);
  r        =  2;
  vROSE[r] = -1;
  // for ( k = 0; k < qM ->valence[vID][0]; ++k) printf(" VALENCY %d = %d\n",k + 1, qM -> valence[vID][2+k]);
  do {
    ////  printf(" AT QUAD %d \n",quadID + 1);
    for ( j = 0 ; j <= 2; ++j) {
      auxm = (id0 + j)%4;
      if ( qM -> quadAdj[4*quadID + auxm] != (prevQuad +1) && qM -> quadAdj[4*quadID + auxm] != (quad0+1)) {
        if ( vROSE[2] != qM -> quadIdx[4*quadID + qLoop[auxm + 1] ]) { // At last quad, this link will be repeated; vROSE[2] is the first link
          vROSE[r++] = qM -> quadIdx[4*quadID + qLoop[auxm + 1] ];
          auxm       = qM -> valence[vROSE[r-1] - 1][0];
          vROSE[r++] = auxm;
        }
      }
    }
    id0      = (int)(id0 + 3)%4;
    prevQuad = quadID;
    quadID   = qM -> quadAdj[4*prevQuad + id0] - 1;
    //   printf(" ADJ %d \n",quadID + 1);
    //for ( k = 0 ; k < r/2; ++k) printf(" STAR %d = %d\n",k,vROSE[2*k]);
    if ( quadID < 0) { // find the other boundary vertex. We will loop around the valences
      for ( k = 0; k < qM ->valence[vID][0]; ++k) {
        auxV = qM -> valence[vID ][k + 2] - 1;
        auxQ = qM -> valence[auxV][1    ] - 1;
        //   printf(" LINK %d USING QUAD %d \n",auxV+1, auxQ + 1);
        for ( i = 0 ; i < r/2; ++i) {
          if ( (auxV + 1) == vROSE[2*i] ) {
            auxQ = -2; i = r;
          }
        }
        if ( auxQ >= 0) {
          // Find the quad common to vID andauxV
          //    printf(" ASS QUAD %d\n",auxQ + 1);
          auxm   = -1;
          for ( i = 0 ; i < 4; ++i) {
            if ( qM -> quadIdx[4*auxQ + i ] == vID + 1) auxm = auxQ;
          }
          if ( auxm == -1) { // Wrong quad. Look around its adjacents to find the one that contains vID
            //    printf(" WE NEED TO FIND THE APPROPRIATE QUAD\n");
            for ( i = 0 ; i < 4; ++i) {
              adjQ = qM -> quadAdj[4*auxQ + i ] - 1;
              //   printf(" ADJ %d \n",adjQ + 1);
              if ( adjQ >= 0) {
                for ( j = 0 ; j < 4; ++j) {
                  if ( qM -> quadIdx[4*adjQ + j ] == vID + 1) {
                    auxm = adjQ;
                  }
                }
              }
              //      printf(" AUXM %d \n",auxm + 1);
            }
          }
          if ( auxm != -1 ) {
            for ( i = 0 ; i < 4; ++i) {
              if ( qM ->quadAdj[4*auxm + i] == -1 ) goto getId0;
            }
          }
        }
        //   printf(" k %d\n",k);
      }
      //  printf(" WE HAVE LINKED ALL THE VERTICES. STAR IS DONE\n ");
      auxm = quad0;
    }
    getId0:
    if ( quadID < 0) {
      vROSE[r++]    = -1; vROSE[r++] = 0;
      vQUADS[itQ++] = -1;
      quadID        = auxm;
    }
    id0 = -1;
    //   printf(" LOOKING  FOR ID0 IN QUAD %d\n",quadID + 1);
    for ( j = 0 ; j < 4; ++j) {
      if ( ( qM->quadIdx[4*quadID + j ] == (vID + 1) ) ) {//&&  (qM ->quadAdj[4*quadID + j] == (prevQuad + 1) ) ) {
        id0  = j;
        break;
      }
    }
    if ( id0 == -1) {
      printf(" WE MESSED UP. BUG ON FINDING ADJ QUAD!!\n");
      exit(1);
    }
    if (quadID != quad0 ) vQUADS[itQ++] = quadID+1;

  }
  while(quadID != quad0);
  if ( itQ != vROSE[1] ) {
    printf(" BUG! NUMBER OF QUADS AND VALENCY SHOULD BE THE SAME!  %d != %d\n",itQ, vROSE[1]);
    exit (1);
  }
  nPeaks = r/2;
  *star            = (vStar*)EG_alloc(sizeof(vStar));
  if ((*star) == NULL ) return EGADS_MALLOC;
  (*star) -> verts = (int*)EG_alloc(r*sizeof(int));
  (*star) -> nV    = nPeaks;
  (*star) -> quads = (int*)EG_alloc(itQ*sizeof(int));
  (*star) -> nQ    = itQ;
  if ( (*star) -> verts == NULL || (*star) -> quads == NULL) {
    EG_free((*star));
    return EGADS_MALLOC;
  }
  for ( r = 0 ; r < itQ; ++r)      (*star) -> quads[r] = vQUADS[r];
  for ( r = 0 ; r < 2*nPeaks; ++r) {
    (*star) -> verts[r] = vROSE[r];
  }
  print_star((*star));
  return EGADS_SUCCESS;
}

static int
EG_removeDoublet(quadMap *qm, int vID, int qID) {
  int val, vA, vB, valA, i, j, jM, jP, k, adjQ;
  int locVal[MAX_LINKS];
  static int loop[5] = {0, 1, 2, 3, 0};
  if ( qm -> valence[vID][0] != 2) return EGADS_RANGERR;
  // First: Remove link to other vertices
  for ( k = 0 ; k < 2; ++k) {
    vA   = qm -> valence[vID][2 + k] - 1;
    valA = qm -> valence[vA ][0    ];
    for ( i = 0 ; i < valA + 1; ++ i) locVal[i] = qm -> valence[vA][i + 1];
    --valA;
    qm -> valence[vA]    = EG_reall(qm -> valence[vA], (valA + 2)*sizeof(int));
    qm -> valence[vA][0] = valA;
    qm -> valence[vA][1] = locVal[0];
    for ( j = i = 0 ; i < valA; ++ i) {
      if ( locVal[i + 1] != vID + 1) {
        qm -> valence[vA][2 + j] = locVal[i];
        ++j;
      }
    }
  }
  // Now build new quad
  adjQ = - 1;
  // Get Adjacent quad
  for ( j = 0 ; j < 4; ++j)
    if ( qm -> quadIdx[4*qID + j] == vA +1) {
      adjQ = qm -> quadAdj[4*qID + j] - 1;
      break;
    }
  if ( adjQ == -1) return EGADS_NULLOBJ;
  // GEt new vertex
  for ( i = j = 0 ; j < 4; ++j) {
    if ( qm -> quadIdx[4*qID + j] == qm -> quadIdx[4*adjQ + i]) {
      ++i; j = -1 ;
    }
  }
  // Remove doublet and insert vertex i
  for ( j = 0 ; j < 4; ++j)
    if ( qm -> quadIdx[4*qID + j] == vID + 1 ) {
      qm -> quadIdx[4*qID + j] = qm -> quadIdx[4*adjQ + i];
      break;
    }
  // Reallocate neighbors for new quad
  jM = (j - 1);   jP = (j + 1);
  if ( jM == -1 ) jM = 3;
  if ( jP ==  4 ) jP = 0;
  for ( k = 0 ; k < 4; ++k) {
    vA = qm -> quadIdx[4*adjQ  + loop[k]   ];
    vB = qm -> quadIdx[4*adjQ  + loop[k +1]];
    if ( ((vA == qm -> quadIdx[4*qID + jM]) && (vB == qm -> quadIdx[4*qID + j ] )) ||
        ((vA  == qm -> quadIdx[4*qID + j ]) && (vB == qm -> quadIdx[4*qID + jM] ))) {
      qm -> quadAdj[4*qID + jM] = qm -> quadAdj[4*adjQ + loop[k]];
    }
    if ( ((vA == qm -> quadIdx[4*qID + jP]) && (vB == qm -> quadIdx[4*qID + j ] )) ||
        ((vA  == qm -> quadIdx[4*qID + j])  && (vB == qm -> quadIdx[4*qID + jP] ))){
      qm -> quadAdj[4*qID + j] = qm -> quadAdj[4*adjQ + loop[k]];
    }
  }
  // Delete quad
  for ( j = 0 ; j < 4; ++j) {
    qm -> quadIdx[4*adjQ + j ] = -1; // set its vertices to 1
    qm -> quadAdj[4*adjQ + j ] = -2;
  }
  // Delete Vertex in map
  qm -> valence[vID]    = EG_reall(qm -> valence[vID], sizeof (int));
  qm -> valence[vID][0] = -1;
  qm -> valence[vID][1] = -2;
  // Reasign correct neighbors
  for ( k = 0 ; k < 4; ++k) {
    vA = qm -> quadAdj[4*qID + k ] - 1;
    for ( j = 0; j < 4; ++ j)
      if ( qm -> quadAdj[4*vA + j ] == (adjQ + 1) ) qm -> quadAdj[4*vA + j ] = qID + 1;
  }
  //print_quad_specs(qm, qID +1 );
  return EGADS_SUCCESS;
}

static int
EG_vertexCollapse(quadMap **qmIN, vStar **starIN,int id0) {
  int        i, j, k,  cID, vC, vC2, vM, vP, rV, nPeaks;
  int        aux, auxID, quadID, stat, adjQ,  val, vCval,  qCount = 0;
  int        poly[2*4], locVal[MAX_LINKS], adjQuads[7]; // local information about vertices and valence
  static int loop[4] = {3, 2, 1, 0}; // for collapsing edges
  double     uv[2], vPos[18];
  quadMap *qm;
  vStar  *star;
  qm   = *qmIN;
  star = *starIN;
  if ( id0%2 == 0) {
    printf(" For COLLAPSING, you should start on a link. I will use the next link: id0 + 1\n");
    id0++;
  }
  poly[0] = star -> verts[0];
  poly[1] = star -> verts[1];
  nPeaks  = qm -> valence[poly[0]][0];
  for ( i = 0; i < 3; ++i){
    aux = id0 + i;
    if ( aux == star ->nV) aux = 1;
    poly[2*(i+1)   ]       = star -> verts[2*aux   ];
    poly[2*(i+1) +1]       = star -> verts[2*aux+ 1];
  }
  for ( i = 0 ; i < 4; ++ i) printf(" COLL POLY %d = %d WIWTH VAL %d\n",i, poly[2*i],poly[2*i+1]);
  aux     = (id0 - 1)/2;
  quadID  = star -> quads[aux];
  printf( " AT QUAD %d\n",quadID + 1);
  // Check that there are at least 3 bad valencies
  aux = 0;
  for ( i = 0 ; i < 4; ++ i) {
    if ( poly[2*i +1] == 4) {
      rV = i; // keep the regular vertex to decide which sides should collapse
      ++aux;
    }
    if ( aux == 2 ) return EGADS_SUCCESS;
  }
  // Look for 3,3 pairs
  if      ( poly[2*0 + 1] == 3 && poly[2*2 + 1] == 3) cID = 0;
  else if ( poly[2*0 + 1] == 3 && poly[2*2 + 1] == 3) cID = 1;
  else cID = (rV + 1)%4;
  printf(" COLLAPSE VERTEX %d = %d\n",cID, poly[2*cID]);
  // Get Quad Centre to place collapsed vertex
  uv[0] = 0.0; uv[1] = 0.0;
  for ( i = 0 ; i < 4; ++i) {
    uv[0] += qm->uv[ 2*(qm->quadIdx[4*quadID + i] - 1)    ];
    uv[1] += qm->uv[ 2*(qm->quadIdx[4*quadID + i] - 1) + 1];
  }
  uv[0] *= 0.25;
  uv[1] *= 0.25;
  stat   = EG_evaluate(qm->face, uv, vPos);
  printf(" VERTEX AT %lf %lf \n",vPos[0], vPos[1]);
  // Collapsing cID with cID +2
  vC              = qm->quadIdx[4*quadID +  cID     ] - 1;
  vC2             = qm->quadIdx[4*quadID + (cID+2)%4] - 1;
  printf(" COLLAPSING VERTEX %d  %d\n",vC+1,vC2 +1);
  qm->uv  [2*vC ] = uv  [0]; qm->uv  [2*vC  + 1] = uv  [1];
  qm->uv  [2*vC2] = uv  [0]; qm->uv  [2*vC2 + 1] = uv  [1];
  qm->xyzs[3*vC ] = vPos[0]; qm->xyzs[3*vC  + 1] = vPos[1]; qm->xyzs[3*vC  + 2] = vPos[2];
  qm->xyzs[3*vC2] = vPos[0]; qm->xyzs[3*vC2 + 1] = vPos[1]; qm->xyzs[3*vC2 + 2] = vPos[2];
  // Modify Quads around vC and point vC2 to vC
#ifdef DEBUG
  adjQuads[qCount] = quadID;
  ++qCount;
  for ( j = 0 ; j < 4; ++j) {
    adjQ = qm ->quadAdj[4*quadID + j] - 1;
    adjQuads[qCount] = adjQ ;
    ++qCount;
#endif
    // point to the new vertex and right neighbor
    //print_quad_specs(qm, adjQ+1);
    for ( i = 0 ; i < 4 ; ++i) {
      if(qm ->quadIdx[4*adjQ + i] == vC2    + 1) qm -> quadIdx[4*adjQ + i] = vC + 1;
      printf(" %d  %d ",qm ->quadAdj[4*adjQ + i],quadID + 1);
      if(qm ->quadAdj[4*adjQ + i] == quadID + 1) qm -> quadAdj[4*adjQ + i] = qm -> quadAdj[4*quadID + loop[j]];
    }
    //print_quad_specs(qm, adjQ+1);
  }
  // Update Valencies
  // Pas valencies of vC2 to vC
  val   = qm -> valence[vC2][0];
  vCval = qm -> valence[vC ][0];
  vM    = cID - 1; if (aux == -1) aux = 3;
  vP    = cID + 1; if (aux ==  4) aux = 0;
  vM    = qm -> quadIdx[4*quadID + vM];
  vP    = qm -> quadIdx[4*quadID + vP];
  for ( i = 0 ; i < vCval + 1; ++i) locVal[i] = qm -> valence[vC][1 + i];
  for ( i = 0 ; i < val; ++i) {
    if ( qm -> valence[vC2][i + 2] != vM && qm -> valence[vC2][i + 2] != vP) {
      locVal[vCval] = qm -> valence[vC2][i + 2];
      ++vCval;
    }
  }
  qm -> valence[vC]    = EG_reall(qm -> valence[vC],(vCval + 2)*sizeof(int));
  qm -> valence[vC][0] = vCval;
  for ( i = 0 ; i < vCval + 1; ++i) qm -> valence[vC][2 + i] = locVal[i];
  // Now remove valence from neighbors
  for ( k = 0 ; k < 2 ; ++ k) {
    if ( k == 0 ) auxID = vM - 1;
    else          auxID = vP - 1;
    val = qm -> valence[auxID][0];
    for ( i = 0 ; i < val + 1; ++i) locVal[i] = qm -> valence[auxID][i+1];
    --val;
    qm -> valence[auxID] = EG_reall(qm -> valence[auxID],(val + 2)*sizeof(int));
    if (qm -> valence[auxID] == NULL) return EGADS_MALLOC;
    for (i = j = 0; i <= val; ++i){
      if( locVal[i + 1] != vC2 + 1) {
        qm -> valence[auxID][2 + j]    = locVal[2*i];
        ++j;
      }
    }
  }
  qm -> valence[auxID][1] =  qm ->valence[ qm -> valence[auxID][2] -1 ][1] - 1;  //COPY A VALID QUAD POSITION
  // Eliminate vertex quadID:
  for ( i = 0 ; i < 4; ++ i) {
    qm -> quadIdx[4*quadID + i ] = -1;
    qm -> quadAdj[4*quadID + i ] = -2;
  }
  // delente vertex vC2
  qm -> valence[vC2]    = EG_reall( qm ->valence[vC2], sizeof(int));
  qm -> valence[vC2][0] = -1;
#ifdef DEBUG
  for ( i = 0 ; i < 4; ++i)
    if ( qm -> quadIdx[4*adjQuads[0] + i ] != -1 || qm -> quadAdj[4*adjQuads[0] + i ] != -2) return EGADS_GEOMERR;
  for ( i = 1; i < 6; ++i){
    for ( j = 0 ; j < 4; ++ j) {
      if ( qm -> quadIdx[4*adjQuads[i] + j ] == vC2 ) {
        printf(" Vertex %d should have disappeared. It is now %d\n",vC2, vC);
        return EGADS_GEOMERR;
      }
      if ( qm -> quadAdj[4*adjQuads[i] + j ] == quadID + 1){
        printf(" Quad %d has disappeared. Nobody should point at it!\n ",quadID + 1);
        return EGADS_GEOMERR;
      }
      if ( qm->quadAdj[4*adjQuads[i] + j] > 0 ) {
        auxID = qm->quadAdj[4*adjQuads[i] + j] - 1; //j-neighbor
        aux = 0;
        for ( k = 0 ; k < 4; ++k) {
          if ( qm -> quadAdj[4*auxID + k ] == adjQuads[i] + 1) {
            aux = 1;
          }
        }
        if ( aux == 0 ) {
          printf(" Quads %d and %d should be neighbors!\n", qm -> quadAdj[4*auxID + k ], adjQuads[i] + 1);
          return EGADS_GEOMERR;
        }
      }
    }
  }

#endif
  *qmIN   = qm;
  EG_free(star);
  stat = EG_buildStar(qm, &star,vC);
  *starIN = star;
  return EGADS_SUCCESS;
}

static int
swappingOperation(quadMap *qm, vStar *star, int poly[], int link, int swapIdx, int *qID, int pullVertex) {
  int stat, i, j, k, aux, auxVAL, auxID, pos, it, adjQ, areaType, area, centreValency;
  int v[4], addV[4], v03[2], OK;
  int perm[6] = {0, 1, 2, 3, 4, 5}, remV[4] = {0, 3, 3, 0}, qL[6]   = { 3, 0, 1, 2, 3, 0};
  int Q1[4], Q2[4], idx1[4], idx2[4], edges[6], q0ID[6], q1ID[6], valence[MAX_LINKS];
  addV[0] = swapIdx; addV[1] = swapIdx + 3;
  addV[2] = addV[1]; addV[3] = addV[0];
  aux     = (link - 1)/2;  // quad location in relation to the vertex
  qID[0]  = star -> quads[aux] - 1;
  ++aux;
  if ( aux == star -> verts[1] ) aux = 0;
  qID[1] = star->quads[aux] - 1;
#ifdef DEBUG
  printf(" ========= EDGE SWAP FUNCTION FOR STAR AROUND VERTEX %d ===========\n",poly[0]);
  printf(" IN QUADS %d AND %d\n",qID[0]+1, qID[1]+1);
  print_quad_specs(qm,qID[0]+1);
  print_quad_specs(qm,qID[1]+1);
  printf(" SWAP INDEX %d\n",swapIdx);
  printf(" WE ARE SWAPPING SIDES %d %d ::: EDGE %d-%dFOR %d-%d\n", swapIdx, swapIdx + 3,poly[0], poly[2*3], poly[2*swapIdx], poly[2*(swapIdx+3)]);
  printf(" VALENCY PAIRS %d %d -- %d %d\n",poly[1], poly[2*3 + 1], poly[2*swapIdx + 1], poly[2*(swapIdx+3) + 1]);
#endif
  /*int vID = poly[0] - 1;
  printf(" ID %d \n",vID +1);
  for ( k = 0; k < qm ->valence[vID][0]; ++k) printf(" VALENCY %d = %d\n",k + 1, qm -> valence[vID][2+k]);
  vID = poly[2*3] - 1;
  printf(" ID %d \n",vID +1);
  for ( k = 0; k < qm ->valence[vID][0]; ++k) printf(" VALENCY %d = %d\n",k + 1, qm -> valence[vID][2+k]);
  vID = poly[2*swapIdx] - 1;
  printf(" ID %d \n",vID +1);
  for ( k = 0; k < qm ->valence[vID][0]; ++k) printf(" VALENCY %d = %d\n",k + 1, qm -> valence[vID][2+k]);
  vID = poly[2*(swapIdx + 3)] - 1;
  printf(" ID %d \n",vID +1);
  for ( k = 0; k < qm ->valence[vID][0]; ++k) printf(" VALENCY %d = %d\n",k + 1, qm -> valence[vID][2+k]);
   */
  // We store the adjacent quad to each edge (ordered). LINKED EDGE 0-3
  for (j = 0 ; j < 4; ++j) if ( qm -> quadIdx[4*qID[0] + j] == poly[2*0] ) break;
  for (k = 0 ; k < 4; ++k) if ( qm -> quadIdx[4*qID[1] + k] == poly[2*3] ) break;
  for ( i = 0 ; i <= 2 ; ++i) {
    edges[i    ] = qm -> quadAdj[4*qID[0] + (j + i)%4];
    edges[i + 3] = qm -> quadAdj[4*qID[1] + (k + i)%4];
  }
  v03[0] = poly[0];
  v03[1] = poly[2*3];
  // Choose which pair to swap
  // MODIFY CORRESPONDING QUADS IN MAP
  for ( i = 0 ; i <= 3; ++ i) {
    Q1[i] = perm[(i + addV[0])%6];
    Q2[i] = perm[(i + addV[1])%6];
  }
  for ( i = 0; i <=2; ++i) {
    q0ID[i    ] = qID[0] + 1;
    q0ID[i + 3] = qID[1] + 1;
    q1ID[Q1[i]] = qID[0] + 1;
    q1ID[Q2[i]] = qID[1] + 1;
  }
  for ( j = 0 ; j < 4; ++j) {
    qm -> quadIdx[4*qID[0] + j]         = poly[2*Q1[j]];
    qm -> quadIdx[4*qID[1] + j]         = poly[2*Q2[j]];
    qm -> valence[poly[2*Q1[j]] -1][1] = qID[0] + 1;
    qm -> valence[poly[2*Q2[j]] -1][1] = qID[1] + 1;
    qm -> quadAdj[4*qID[0] + j]         = edges[Q1[j]];
    qm -> quadAdj[4*qID[1] + j]         = edges[Q2[j]];
  }
  qm -> quadAdj[4*qID[0] + 3] = qID[1] + 1;
  qm -> quadAdj[4*qID[1] + 3] = qID[0] + 1;
  for ( j = 0; j < 6; ++j) {
    for ( i = 0; i < 4; ++i) {
      if ( edges[j] != -1 ) {
        if ( qm -> quadAdj[4*(edges[j] - 1) + i ] == q0ID[j]) qm -> quadAdj[4*(edges[j] - 1) + i ] = q1ID[j];
      }
    }
  }
#ifdef DEBUG
  for ( i = 0; i < 6; ++ i) {
    OK = 0;
    if (edges[i] >= 0) {
      for ( k = 0; k < 4; ++k) {
        if (qm ->quadAdj[4*(edges[i] - 1)  + k ] == q1ID[i] ) {
          OK = 1; k  = 5;
        }
      }
      if (OK == 0 ) {
        printf("Quad at edge %d = %d \n", i, edges[i]);
        print_quad_specs(qm, edges[i]);
        printf(" ONE ADJACENT SHOULD POINT AT %d !!!\n",q1ID[i]);
      }
    }
  }
#endif
  // Now update the valencies : REMOVE PAIR 0-3 AND ADD  PAIR 1-4 ( 2-5)
  for ( j = 0 ; j < 2; ++j) {
    auxID  = poly[2*remV[2*j]   ] - 1; //central vertex
    auxVAL = poly[2*remV[2*j] +1]    ;
    if ( auxID < 0) {
      printf(" Found a problem whilst swapping. We are not going to swap that edge\n");
      return EGADS_INDEXERR;
    }
    for ( i = 0 ; i < auxVAL; ++i) valence[i] = qm -> valence[auxID][2 + i];
    --auxVAL;
    qm -> valence[auxID]    = EG_reall(qm -> valence[auxID], (auxVAL + 2)*sizeof(int));
    qm -> valence[auxID][0] = auxVAL;  // TOTAL VALENCY
    for ( k = i = 0 ; i <= auxVAL; ++i) {
      if ( valence[i] != poly[2*remV[2*j+1]]){
        qm -> valence[auxID][2 + k] = valence[i];
        ++k;
      }
    }
    auxID  = poly[2*addV[2*j]    ] - 1; //central vertex
    auxVAL = poly[2*addV[2*j] + 1];
    for ( i = 0 ; i < auxVAL; ++i) valence[i] = qm -> valence[auxID][2 + i];
    ++auxVAL;
    qm -> valence[auxID]    = EG_reall(qm -> valence[auxID ], (auxVAL + 2)*sizeof(int));
    qm -> valence[auxID][0] = auxVAL;
    for ( i = 0 ; i < auxVAL ; ++i) qm -> valence[auxID][i + 2] = valence[i];
    qm ->valence[auxID ][2 + auxVAL - 1] = poly[2*addV[2*j + 1]];
  }
  // Finally: Check if we have created a doublet
  for ( j = 0 ; j < 2; ++j) {
    auxID = poly[2*remV[2*j] ] - 1; //central vertex
    if (qm -> valence[auxID][0] == 2 && qm -> vType[auxID] != 0) {
      return EGADS_GEOMERR;
      for ( i = 0 ; i < 4 ; ++ i ) {
        if ( qm -> quadIdx[4*qID[1] + i] == auxID + 1) {
          qID[0] = qID[1];
          i    = 4;
        }
      }
      printf(" CALL DOUBLET FOR VERTEX %d IN %d WITH VAL %d\n", auxID + 1, qID[0] + 1, qm -> valence[auxID][0]);
      stat = EG_removeDoublet(qm, auxID,qID[0]);
      if ( stat != EGADS_SUCCESS) return stat;
    }
  }
  // FIND SIDE FROM WICH TO PULL
  adjQ = -2;
  for ( k = 0 ; k < 4; ++k){
    if (qm -> quadIdx[4*qID[0] + k] == v03[0] ) break;
    else if (qm -> quadIdx[4*qID[0] + k] == v03[1] ) {
      aux    = qID[0];
      qID[0] = qID[1];
      qID[1] = aux;
      break;
    }
  }
  stat = 0 ;
  //printf(" NEW QUADS %d AND %d\n",qID[0]+1, qID[1]+1);
  //print_quad_specs(qm,qID[0]+1);
  //print_quad_specs(qm,qID[1]+1);
  for ( j = 0 ; j < 2 ; ++ j) {
    stat += checkInvalidElement(qm, qID[j] + 1, v03[j], pullVertex);
  }
  return stat;
}


static int
EG_edgeSwap(quadMap *qm, vStar *star)
{
  int stat, i, j, k, aux, vL, movingVerts[3], nMove = 0, centreValency;
  int v[4], vIrr[6], qID[2], optiVerts[3], v03[2], link;
  int *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, **locVal = NULL, poly[2*6];
  double fOpti;
  if (star ->nV == 0) return EGADS_EMPTY;
  locIdx = (int*)   EG_alloc(4*qm -> nQ       *sizeof(int ));
  locAdj = (int*)   EG_alloc(4*qm -> nQ       *sizeof(int ));
  locVal = (int **) EG_alloc(  qm -> totVerts *sizeof(int*));
  if ( locIdx == NULL || locAdj == NULL || locVal == NULL) return EGADS_MALLOC;
  for ( j = 0; j< qm -> totVerts; ++j) {
    locVal[j] = (int *) EG_alloc((2 + MAX_LINKS)*sizeof(int));
    if (locVal[j] == NULL ) {
      EG_free(locAdj);
      EG_free(locIdx);
      EG_free(locVal);
      return EGADS_MALLOC;
    }
    for ( i = 0 ; i < qm -> valence[j][0] + 2; ++i) locVal[j][i] = qm -> valence[j][i];
  }
  for ( i = 0 ; i < qm ->nQ; ++i) {
    for ( j = 0 ; j < 4; ++j){
      locIdx[4*i + j] = qm -> quadIdx[4*i + j];
      locAdj[4*i + j] = qm -> quadAdj[4*i + j];
    }
  }
  poly[0] = star->verts[0];
  poly[1] = star->verts[1];
  printf(" 88888888888888      EDGE SWAP FUNCTION     8888888888888888\n");
  // CHECK AROUND ALL THE STAR VERTICES AND TRY TO SWAP ONE THAT DOESN'T PRODUCE INVALID ELEMENTS
  centreValency = qm -> valence[poly[0] - 1][0];
  swapPtr       = (int *) EG_alloc(centreValency *sizeof(int));
  if (swapPtr == NULL) goto cleanup;
  for ( link = 0; link < centreValency; ++link) {
    vL = 2*(2*link + 1);
    //printf(" START AT %d\n",vL);
    for ( i = 0 ; i  < 5; ++i) {
      aux = vL + 2*i;
      //printf(" AUX %d TOT %d  ,   %d\n",aux, 2*star -> nV, star -> nV -1);
      if ( aux >= 2*star ->nV) aux -= 2*(star->nV -1);
      //printf( " AUX NOW %d ---> %d\n",aux, star->verts[aux    ]);
      poly[2*(i+1)   ] = star->verts[aux    ];
      poly[2*(i+1) +1] = star->verts[aux + 1];
      if ( poly[2*(i+1)] < 0) { //This is a virtual quad used when the star centre is at the boundary. Skip
        i  = 5;
        vL = -1;
      }
    }
    swapPtr[link] = -1;
    if ( vL != -1 ) {
      for ( i = 0; i < 6; ++i) {
        // printf(" P(%d) = %d\n",poly[2*i], poly[2*i + 1]);
        if (qm -> vType[ poly[2*i] -1] != -1) vIrr[i] = 100;
        else {
          if      ( poly[2*i+1]  < 4 ) vIrr[i] = -1;
          else if ( poly[2*i+1] == 4 ) vIrr[i] =  0;
          else vIrr[i] =  poly[2*i + 1];
        }
      }
      // CHECK IF WE SHOULD SWAP
      if (qm -> vType[ poly[0] -1] != -1) {
        printf(" SWAPPING FORCED \n");
        if ( (vIrr[1] + vIrr[4]) < (vIrr[2] + vIrr[5]) )  swapPtr[link] = 1;
        else swapPtr[link] = 2;
        printf(" LINK %d\n",swapPtr[link]);
      } else {
        if   ( (vIrr[0] + vIrr[3]) >  0) {
          if      ( (vIrr[1] + vIrr[4]) == -2)  swapPtr[link] = 1;
          else if ( (vIrr[2] + vIrr[5]) == -2)  swapPtr[link] = 2;
          else if (( vIrr[0] + vIrr[3]) >  5) { // two high valencies
            if ((vIrr[1] == -1 && vIrr[4] == 0) ||
                (vIrr[4] == -1 && vIrr[1] == 0)) swapPtr[link] = 1;
            else if ((vIrr[2] == -1 && vIrr[5] == 0) ||
                (vIrr[5] == -1 && vIrr[2] == 0)) swapPtr[link] = 1;
          }
        }
      }
      if ( swapPtr[link] != -1 ) {
        printf(" WE ARE SWAPPING: POSITION %d\n",swapPtr[link]);
        stat = swappingOperation(qm, star, poly, 2*(link + 1), swapPtr[link], qID, 0);
        if ( stat == EGADS_SUCCESS ) {
          link = MAX_LINKS;
          ++TOTSWAPS;
          break;
        } else {
          for ( i = 0 ; i < qm -> nQ; ++i) {
            for ( k = 0 ; k  <4; ++k) {
              qm -> quadIdx[4*i + k] = locIdx[4*i+ k];
              qm -> quadAdj[4*i + k] = locAdj[4*i+ k];
            }
          }
          for ( j = 0 ; j < qm -> totVerts; ++j) {
            qm -> valence[j] = EG_reall(qm -> valence[j], (locVal[j][0] + 2)*sizeof(int));
            for ( i = 0 ; i < locVal[j][0] + 2; ++i)  {
              qm -> valence[j][i] = locVal[j][i];
            }
          }
        }
      }
      else EGADS_INDEXERR;
    }
    else stat = EGADS_INDEXERR;
  }
  if (stat != EGADS_SUCCESS){
    printf(" All changes seem to require moving the vertex. Go ahead and move it\n");
    for ( link = 0 ; link < centreValency; ++link) {
      if (swapPtr[link] != -1) {
        vL = 2*(2*link + 1);
        for ( i = 0 ; i  < 5; ++i) {
          aux = vL + 2*i;
          if ( aux >= 2*star ->nV) aux -= 2*(star->nV -1);
          poly[2*(i+1)   ] = star->verts[aux   ];
          poly[2*(i+1) +1] = star->verts[aux+ 1];
          if ( poly[2*(i+1)] < 0) { //This is a virtual quad used when the star centre is at the boundary. Skip
            i  = 5;
            vL = -1;
          }
        }
        if ( vL != -1) {
          for ( i = 0 ; i  < 6; ++i) printf(" P(%d) = %d\n",poly[2*i], poly[2*i + 1]);
          printf(" SWAPPING %d\n",swapPtr[link]);
          stat = swappingOperation( qm, star, poly, 2*(link + 1), swapPtr[link], qID, 1);
          if ( stat == EGADS_SUCCESS ) {
            ++TOTSWAPS;
            link = MAX_LINKS;
            break;
          }
          else {
            for ( i = 0 ; i < qm -> nQ; ++i) {
              for ( k = 0 ; k < 4; ++k) {
                qm -> quadIdx[4*i + k] = locIdx[4*i+ k];
                qm -> quadAdj[4*i + k] = locAdj[4*i+ k];
              }
            }
            for ( j = 0 ; j < qm -> totVerts; ++j) {
              qm -> valence[j] = EG_reall(qm -> valence[j], (locVal[j][0] + 2)*sizeof(int));
              for ( i = 0 ; i < locVal[j][0] + 2; ++i)  qm -> valence[j][i] = locVal[j][i];
            }
            goto cleanup;
          }
        }
      }
    }
  }
  if (stat != EGADS_SUCCESS) goto cleanup;
  printf(" WE FOUND A VALID SWAP. OPTIMIZE AND LEAVE\n");
  for ( i = 0 ; i < 6; ++i ) printf(" VIRR %d %d\n",poly[2*i],poly[2*i + 1]);
  for ( j = 0 ; j < 2 ; ++ j) {
#if TRIPLEOPTI
    for ( k = 0 ; k < 4; ++k)
      if (qm ->quadIdx[4*qID[j] + k ] == poly[2*3*j]) break;
    i = k - 1; if (i < 0 ) i = 3;
    optiVerts[0] = qm ->quadIdx[4*qID[j] + i ] - 1;
    optiVerts[1] = qm ->quadIdx[4*qID[j] + k ] - 1;
    i = k + 1; if (i >3 ) i = 0;
    optiVerts[2] = qm ->quadIdx[4*qID[j] + i ] - 1;
    nMove = 0 ;
    for ( i = 0 ; i < 3; ++i) {
      if( qm-> vType[ optiVerts[i] ] == -1) movingVerts[nMove++] = optiVerts[i];
    }
    stat = optimize_angles(qm, nMove, movingVerts);

#else
    if (qm -> vType[poly[2*3*j]-1] == -1 ) {
      movingVerts[0] = poly[2*3*j] - 1;
      //stat = optimize_angles(qm, 1, movingVerts);
    }
#endif
  }
  cleanup:
  for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
  EG_free(locVal);
  EG_free(locIdx);
  EG_free(locAdj);
  EG_free(swapPtr);
  printf("EDGE SWAP SUCCESSFUL\n");
  return EGADS_SUCCESS;
}


static int updateUVs(ego tess, int iFace, ego face, double *uvs, double *xyzs)
{
  int          nSweeps = 3, itQ, nPeaks, i, j, k, count, kOK, id0, auxID, auxQ, nseg, newUV, stat, nVars, quadID, prevQuad;
  int          kk,  m, n, iper, len, nTri, optVec[MAX_LINKS], nOpt;
  vStar        *star = NULL;
  const int    *tris, *tric, *ptype, *pindex;
  double       uvbox[4], *lowerBounds = NULL, *upperBounds = NULL, fOpti, *uvOpt=NULL;
  const double *xyzTess, *uvTess;
  static int   qV[6]    = { 0, 1, 2, 5, 0, 1};
  static int   qLoop[5] = { 0, 1, 2, 3, 0};
  quadMap      *qM = NULL;
  stat = EG_getTessFace(tess, iFace, &len,
      &xyzTess, &uvTess, &ptype, &pindex, &nTri, &tris, &tric);
  if (stat != EGADS_SUCCESS) {
    return stat;
  }
  qM = (quadMap*) EG_alloc(sizeof(quadMap));
  if ( qM == NULL ) {
    stat =  EGADS_MALLOC;
    goto cleanup;
  }
  qM -> face     = face;
  qM -> totVerts = len;
  qM -> nQ       = nTri/2;
  qM -> xyzs     = (double*) EG_alloc(3*len    *sizeof(double));
  qM -> uv       = (double*) EG_alloc(2*len    *sizeof(double));
  qM -> vType    = (int*)    EG_alloc(len      *sizeof(   int));
  qM -> quadIdx  = (int*)    EG_alloc(4*qM ->nQ*sizeof(   int));
  qM -> quadAdj  = (int*)    EG_alloc(4*qM ->nQ*sizeof(   int));
  qM -> nS       = 0;
  if ( qM ->quadIdx == NULL || qM ->quadAdj == NULL || qM -> xyzs == NULL || qM -> uv == NULL || qM ->vType == NULL) {
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  for (j = 0; j < len; j++) {
    qM -> xyzs[3*j] = xyzTess[3*j]; qM ->xyzs[3*j + 1] = xyzTess[3*j +1];qM ->xyzs[3*j +2] = xyzTess[3*j +2];
    qM -> uv  [2*j] = uvTess [2*j]; qM ->uv  [2*j + 1] = uvTess [2*j +1];
    qM -> vType [j] = ptype[j];
  }
  for (j = 0; j < qM -> nQ; j++)
    for ( k = 0; k < 4; ++k)  qM -> quadIdx[4*j + k ] = tris[6*j + qV[k+1]];
  qM -> valence = (int **) EG_alloc(len*sizeof(int*));
  if ( qM -> valence == NULL ) {
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  for ( j = 0; j < len; ++j) {
    qM -> valence[j] = (int *) EG_alloc((2 + MAX_LINKS)*sizeof(int));
    if (qM -> valence[j] == NULL ) {
      stat = EGADS_MALLOC;
      goto cleanup;
    }
    qM -> valence[j][0] = 0;
  }
  for ( j = 0; j < qM -> nQ; j++) {
    kk  = 0;
    kOK = 0;
    for ( i = 0 ; i < qM -> nQ; ++i) {
      if ( i == j || i == -1) ++i;
      for ( k = 0 ; k < 4; ++k ) {
        if ( i == j || i == -1) ++i;
        if ( (qM -> quadIdx[4*j + qLoop[kk    ]] == qM -> quadIdx[4*i + qLoop[k]] || qM -> quadIdx[4*j + qLoop[kk    ]] == qM -> quadIdx[4*i + qLoop[k + 1]] )&&
            ( qM -> quadIdx[4*j + qLoop[kk + 1]] == qM -> quadIdx[4*i + qLoop[k]] || qM -> quadIdx[4*j + qLoop[kk + 1]] == qM -> quadIdx[4*i + qLoop[k + 1]]) ) {
          qM -> quadAdj[4*j + kk] = i+1;
          ++kk;
          i   = -1;
          kOK = 1;
          k   = 4;
          if (kk == 4)  i = qM -> nQ;
        }
      }
      if ( (kOK == 0) && (i == qM -> nQ -1) ){
        qM -> quadAdj[4*j + kk] = -1;
        i                       = -1;
        ++kk;
        if (kk == 4) i = qM -> nQ;
      }
      else  kOK = 0 ;
    }
  }
  for ( j = 0 ; j < qM ->nQ; ++j) {
    for ( k = 0 ; k < 4; ++k) {
      if ( qM -> quadAdj[4*j +k] > j +1  || qM -> quadAdj[4*j +k] == -1) {
        auxID = qM -> quadIdx [4*j + qLoop[k]] -1;
        qM   -> valence[auxID][1] = j + 1;
        qM   -> valence[auxID][2 + qM -> valence[auxID][0] ] = qM -> quadIdx[4*j+ qLoop[k+1]];
        ++qM -> valence[auxID][0];
        auxID = qM -> quadIdx [4*j + qLoop[k+1]] -1;
        qM   -> valence[auxID][1] = j + 1;
        qM   -> valence[auxID][2 + qM -> valence[auxID][0] ] = qM -> quadIdx[4*j+ qLoop[k]];
        ++qM -> valence[auxID][0];
      }
    }
  }
#ifdef DEBUG
  printf(" STARTING QUAD MESH \n");
  for ( j = 0 ; j < qM -> nQ; ++j) print_quad_specs(qM,j+1);
  print_mesh(qM);
  printf(" STARTING VALENCES \n");
  for ( j = 0 ; j < qM -> totVerts; ++j) {
    printf(" VERTEX %d\n",j+1);
    for ( i = 0 ; i < qM -> valence[j][0]; ++i) printf(" VALENCE %d = %d\n",i+1, qM ->valence[j][2+ i]);
  }
#endif
  count = 0;
  id0   = 0;
  printMeshStats(qM,0);
  for ( k = 0 ; k < nSweeps; ++k) {
    // DOMAIN CORNERS 
    printf("------------------------\n  DOMAIN CORNERS  \n ----------------------\n");
    for (i = 0 ; i < len; ++i) {
      if (ptype[i] == 0  && qM -> valence[i][0] > 2) {
        printf(" IN UV SWEEP %d  CALLING FOR STAR AROUND %d\n",k, i+1);
        stat = EG_buildStar(qM, &star,i);
        if ( stat != EGADS_SUCCESS) goto cleanup;
        count = TOTSWAPS;
        stat = EG_edgeSwap(qM, star);
        if ( TOTSWAPS > count) print_mesh(qM);
        if ( stat != EGADS_SUCCESS) goto cleanup;
        //stat = EG_vertexCollapse(&qM, &star, i + 1);
      }
    }
    // DOMAIN EDGES
    printf("------------------------\n  DOMAIN EDGES  \n ----------------------\n");
    for (i = 0 ; i < len; ++i) {
      if ( (ptype[i] > 0)  && (qM -> valence[i][0] > 3) ) {
        printf(" IN UV SWEEP %d  CALLING FOR STAR AROUND %d\n",k, i+1);
        stat = EG_buildStar(qM, &star,i);
        if ( stat != EGADS_SUCCESS) goto cleanup;
        count = TOTSWAPS;
        stat = EG_edgeSwap(qM, star);
        if ( TOTSWAPS > count) print_mesh(qM);
        if ( stat != EGADS_SUCCESS) goto cleanup;
        //stat = EG_vertexCollapse(&qM, &star, i + 1);
      }
    }
    printf("------------------------\n  INTERIOR VERTICES  \n ----------------------\n");
    for (i = 0 ; i < len; ++i) {
      if (ptype[i] == -1) {
        printf(" IN UV SWEEP %d  CALLING FOR STAR AROUND %d\n",k, i+1);
        stat = EG_buildStar(qM, &star,i);
        if ( stat != EGADS_SUCCESS) goto cleanup;
        count = TOTSWAPS;
        stat = EG_edgeSwap(qM, star);
        if ( TOTSWAPS > count) print_mesh(qM);
        if ( stat != EGADS_SUCCESS) goto cleanup;
        //stat = EG_vertexCollapse(&qM, &star, i + 1);
      }
    }
    printf(" DONE ONE ROUND. NOW OPTIMIZE FOR THE WHOLE MESH ");
    for ( nOpt = i = 0 ; i < len; ++i)
      if ( qM -> vType[i] == -1)  optVec[nOpt++] = i;
    stat = optimize_angles(qM, nOpt, optVec);
    printf(" TOTAL SWAPS IN THIS ROUND   %d   \n", TOTSWAPS);
    OVERALLSWAPS += TOTSWAPS;
    TOTSWAPS      = 0 ;
    printMeshStats(qM, k + 1);
  }
  for ( nOpt = i = 0 ; i < len; ++i)
    if ( qM -> vType[i] == -1)  optVec[nOpt++] = i;
  stat = optimize_angles(qM, nOpt, optVec);
  print_mesh(qM);
  cleanup:
  EG_free(qM -> quadIdx);
  EG_free(qM -> quadAdj);
  EG_free(qM -> uv);
  EG_free(qM -> xyzs);
  for ( i = 0 ; i < len; ++i) EG_free(qM ->valence[i]);
  EG_free(qM -> valence);
  EG_free(qM -> vType);
  EG_free(qM);
  EG_free(star ->verts);
  EG_free(star ->quads);
  EG_free(star);
  EG_free(tess);
  printf(" OVERALL SWAPPING OPERATION:   %d\n",OVERALLSWAPS);
  return stat;
}



double NLOPT_objEqualAnglesTess(unsigned n, const double *uvs, double *grad, void *inQM){
  int    stat, i, j, k, vFix, vA, vB, aux, oID, peaks;
  double thetaOpt, theta;
  double vXYZ[18], oldXYZ[18], fUV = 0.0, fVert = 0.0, area;
  quadMap *qm = (quadMap*)inQM;
  if (OPTCOUNT == 0 ) SCOUNT = 0 ;
  ++OPTCOUNT;
  for ( i = 0 ; i < qm->nS; ++i) {
    oID  = qm -> star[i] -> verts[0] -1 ;
    stat = EG_evaluate( qm->face, &uvs[2*i], vXYZ);
    qm  -> xyzs[3*oID    ] = vXYZ[0];
    qm  -> xyzs[3*oID + 1] = vXYZ[1];
    qm  -> xyzs[3*oID + 2] = vXYZ[2];
#ifdef DEBUG
    if (OPTCOUNT%NOUT == 0 )
      printf(" UPDATE COORDS %lf  %lf  \n", vXYZ[0], vXYZ[1]);
#endif
  }
  for ( i = 0 ; i < qm->nS; ++i) {
    fVert    = 0.0;
    oID      = qm -> star[i] -> verts[0] -1;
    peaks    = qm -> valence[oID][0];
    thetaOpt = 2.0*PI/(double)(peaks);
    stat     = EG_evaluate( qm->face, &uvs[2*i], vXYZ);
    area     = star_algebraic_area(qm, qm -> star[i]);
    for (j = 0 ; j < qm->star[i] ->nQ; ++j) {
      stat = quadInteriorAngle(qm, qm->star[i]->quads[j], qm ->star[i]->verts[0], &theta);
      fVert += (thetaOpt - theta)*(thetaOpt - theta);
    }
    fUV += fVert + area;
#ifdef DEBUG
    if (OPTCOUNT%NOUT == 0 )
      printf(" IN STAR LOOP COORDIS %lf %lf FUNCTION THETA OPT %lf PICOS  %d  THETA %lf %lf TOTAL %lf  PENALTY  %lf\n", uvs[2*i], uvs[2*i +1], thetaOpt, peaks, theta,fVert, fUV, area);
    //   print_starFile(qm->star[i], qm);
#endif

  }
#ifdef DEBUG
  if (OPTCOUNT%NOUT == 0 )
    printf(" TOTAL FUNCTION %d %lf\n",OPTCOUNT, fUV);
#endif
  return fUV;
}


static int optimize_angles(quadMap *qm, int nV, int *vID)
{
  int i, j, k, stat, iper, peaks, kk, aux;
  double  XYZ[18], uvbox[4], *lowerBounds = NULL, *upperBounds = NULL, fOpti, *uvOpt=NULL;
  double  epsREL = 1.e-08, *epsABS = NULL;
  qm ->nS     = nV;
  qm ->star   = (vStar**) EG_alloc(  nV*sizeof(vStar*));
  uvOpt       = (double*) EG_alloc(2*nV*sizeof(double));
  lowerBounds = (double*) EG_alloc(2*nV*sizeof(double));
  upperBounds = (double*) EG_alloc(2*nV*sizeof(double));
  stat = EG_getRange(qm->face, uvbox, &iper);
  if (lowerBounds == NULL || upperBounds == NULL || uvOpt == NULL || qm ->star == NULL)  return EGADS_MALLOC;
  for ( i = 0; i < nV; ++i){
    stat = EG_buildStar(qm, &qm ->star[i], vID[i]);
    if ( stat != EGADS_SUCCESS) goto cleanup;
    uvOpt[2*i    ]     = qm -> uv[2*vID[i]    ];
    uvOpt[2*i + 1]     = qm -> uv[2*vID[i] + 1];
    lowerBounds[2*i  ] = uvbox[0]; lowerBounds[2*i+1] = uvbox[2];
    upperBounds[2*i  ] = uvbox[1]; upperBounds[2*i+1] = uvbox[3];
    peaks              = (qm->star[i]->nV-1)/2;
    for ( j = 0; j < peaks; ++j) {
      k  = 2*j + 1;
      kk = qm->star[i]->verts[2*k] - 1;
    }
  }
  OPTCOUNT = 0;
#ifdef DEBUG
  printf(" NLOPT OPTIMIZER FOR THE VERTICES:\n");
  for ( i = 0; i < nV; ++i) printf(" V = %d UVS %lf %lf WITH COORDS %lf %lf %lf \n", vID[i] + 1, uvOpt[2*i], uvOpt[2*i+1],qm-> xyzs[3*vID[i]],qm-> xyzs[3*vID[i]+1],qm-> xyzs[3*vID[i]+2]);

  for ( i = 0; i < nV; ++i) {
    int area    = star_algebraic_area(qm, qm -> star[i]);
    // print_starFile(qm->star[i], qm);
  }
  printf(" OPTIMIZE VERTEX LOCATION\n");
#endif
  // constraints for optimizer, u,v in face
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getRange in optimize_angles = %d\n",stat);
    goto cleanup;
  }
  nlopt_opt opt;
  opt = nlopt_create(NLOPT_LN_COBYLA, 2*nV);
  epsABS = (double * ) EG_alloc(2*nV*sizeof(double));
  if ( epsABS == NULL ) {
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  for ( i = 0 ; i < 2*nV; ++i) epsABS[i] = 1.e-08;
  nlopt_set_lower_bounds (opt, lowerBounds);
  nlopt_set_upper_bounds (opt, upperBounds);
  nlopt_set_min_objective(opt, NLOPT_objEqualAnglesTess, qm);
  nlopt_set_stopval(opt,  0.1);
  nlopt_set_maxeval(opt, NLOPTMAXEVAL);
  //nlopt_set_ftol_rel(opt, epsREL);
  nlopt_set_xtol_rel(opt, epsREL);
  //nlopt_set_ftol_abs(opt, epsABS);
  nlopt_set_xtol_abs(opt, epsABS);
  nlopt_result res = nlopt_optimize(opt, uvOpt, &fOpti);
  NLOPTtermination(res);
  if (res >=0 )  stat = EGADS_SUCCESS;
  else {
    stat = EGADS_GEOMERR;
    goto cleanup;
  }
  for (j = 0; j < nV; ++j) {
    qm ->uv[2*vID[j]]     = uvOpt[2*j];
    qm ->uv[2*vID[j] + 1] = uvOpt[2*j + 1];
    stat                  = EG_evaluate(qm ->face, &uvOpt[2*j], XYZ);
    qm -> xyzs[3*vID[j]]  = XYZ[0]; qm -> xyzs[3*vID[j]+1] = XYZ[1]; qm -> xyzs[3*vID[j]+2] = XYZ[2];
  }
#ifdef DEBUG
  printf(" NLOPT EXITED WITH %lf  VALUE AFTER %d ITERATIONS RES = %d\n",fOpti, OPTCOUNT, res);
  for (j = 0; j < nV; ++j)
    printf(" NEW COORDINATES %lf  %lf  %lf\n", qm -> xyzs[3*vID[j]], qm -> xyzs[3*vID[j] + 1], qm -> xyzs[3*vID[j] +2] );
#endif
  nlopt_destroy(opt);
  EG_free(lowerBounds);
  EG_free(upperBounds);
  cleanup:
  for ( i = 0 ; i < nV; ++i){
    EG_free(qm->star[i]->verts);
    EG_free(qm->star[i]->quads);
    EG_free(qm->star[i]);
  }
  EG_free(qm->star);
  EG_free(uvOpt);
  return  stat;
}









int main(int argc, char *argv[])
{
  int          i, j, k, ibody, stat, oclass, mtype, len, ntri, sum;
  int          nseg, ngp, atype, alen, quad, *segs, *senses;
  const int    *tris, *tric, *ptype, *pindex, *ints;
  float        arg, color[3];
  double       box[6], size, params[3];
  const double *xyzs, *uvs, *ts, *reals;
  char         gpname[33], *startapp;
  const char   *OCCrev, *string;
  ego          context, tess, model, geom, *bodies, *dum;
  wvData       items[5];
  float        eye[3]      = {0.0, 0.0, 7.0};
  float        center[3]   = {0.0, 0.0, 0.0};
  float        up[3]       = {0.0, 1.0, 0.0};
  static int   sides[3][2] = {{1,2}, {2,0}, {0,1}};
  static int   sideq[4][2] = {{1,2}, {2,5}, {5,0}, {0,1}};
  static int   neigq[4]    = {  0,     3,     4,     2};
  /* get our starting application line
   *
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");

  if ((argc != 2) && (argc != 5)) {
    printf("\n Usage: vAttr filename [angle maxlen sag]\n\n");
    return 1;
  }

  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

  /* initialize */
  printf(" EG_open           = %d\n", EG_open(&context));
  stat = 0;
  if(argv[1])  stat = EG_loadModel(context, 0, argv[1], &model);
  printf(" EG_loadModel      = %d\n", stat);
  printf(" EG_getBoundingBox = %d\n", EG_getBoundingBox(model, box));
  printf("       BoundingBox = %lf %lf %lf\n", box[0], box[1], box[2]);
  printf("                     %lf %lf %lf\n", box[3], box[4], box[5]);
  printf(" \n");

  size = box[3]-box[0];
  if (size < box[4]-box[1]) size = box[4]-box[1];
  if (size < box[5]-box[2]) size = box[5]-box[2];

  focus[0] = 0.5*(box[0]+box[3]);
  focus[1] = 0.5*(box[1]+box[4]);
  focus[2] = 0.5*(box[2]+box[5]);
  focus[3] = size;

  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
      &bodies, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n", stat);
    return 1;
  }
  printf(" EG_getTopology:     nBodies = %d\n", nbody);
  bodydata = (bodyData *) malloc(nbody*sizeof(bodyData));
  if (bodydata == NULL) {
    printf(" MALLOC Error on Body storage!\n");
    return 1;
  }

  params[0] =  0.1 *size;
  params[1] =  0.01*size;
  params[2] =  20.0;
  if (argc == 5) {
    sscanf(argv[2], "%f", &arg);
    params[2] = arg;
    sscanf(argv[3], "%f", &arg);
    params[0] = arg;
    sscanf(argv[4], "%f", &arg);
    params[1] = arg;
    printf(" Using angle = %lf,  relSide = %lf,  relSag = %lf\n",
        params[2], params[0], params[1]);
    params[0] *= size;
    params[1] *= size;
  }

  /* fill our structure a body at at time */
  for (ibody = 0; ibody < nbody; ibody++) {
#ifdef FORCETRIANGULATION
    stat = EG_attributeAdd(bodies[ibody], ".qParams", ATTRSTRING, 4, NULL, NULL,
        "off");
    if (stat != EGADS_SUCCESS)
      printf(" Body %d: attributeAdd = %d\n", ibody, stat);
#endif


    mtype = 0;
    EG_getTopology(bodies[ibody], &geom, &oclass,
        &mtype, NULL, &j, &dum, &senses);
    bodydata[ibody].body  = bodies[ibody];
    bodydata[ibody].mtype = mtype;
    bodydata[ibody].uvs   = NULL;
    if (mtype == WIREBODY) {
      printf(" Body %d: Type = WireBody\n", ibody+1);
    } else if (mtype == FACEBODY) {
      printf(" Body %d: Type = FaceBody\n", ibody+1);
    } else if (mtype == SHEETBODY) {
      printf(" Body %d: Type = SheetBody\n", ibody+1);
    } else {
      printf(" Body %d: Type = SolidBody\n", ibody+1);
    }

    stat = EG_getBodyTopos(bodies[ibody], NULL, FACE,
        &bodydata[ibody].nfaces, &bodydata[ibody].faces);
    i    = EG_getBodyTopos(bodies[ibody], NULL, EDGE,
        &bodydata[ibody].nedges, &bodydata[ibody].edges);
    if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
      printf(" EG_getBodyTopos Face = %d\n", stat);
      printf(" EG_getBodyTopos Edge = %d\n", i);
      return 1;
    }

    stat = EG_makeTessBody(bodies[ibody], params, &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTessBody %d = %d\n", ibody, stat);
      continue;
    }
    stat = EG_makeTessBody(bodies[ibody], params, &bodydata[ibody].tess);

#ifdef FORCEQUADS
    tess = bodydata[ibody].tess;
    stat = EG_quadTess(tess, &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_quadTess %d = %d  -- reverting...\n", ibody, stat);
      bodydata[ibody].tess = tess;
      continue;
    }
    EG_deleteObject(tess);
    QUADTESS = 1;
#endif
  }
  printf(" \n");

  /* create the WebViewer context */
  cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
  if (cntxt == NULL) {
    printf(" failed to create wvContext!\n");
    for (ibody = 0; ibody < nbody; ibody++) {
      EG_deleteObject(bodydata[ibody].tess);
      EG_free(bodydata[ibody].edges);
      EG_free(bodydata[ibody].faces);
    }
    free(bodydata);

    printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
    printf(" EG_close          = %d\n", EG_close(context));
    return 1;
  }

  /* make the scene */
  for (ngp = sum = stat = ibody = 0; ibody < nbody; ibody++) {
    printf(" NGE P %d\n",ngp);
    quad = 0;
    stat = EG_attributeRet(bodydata[ibody].tess, ".tessType", &atype,
        &alen, &ints, &reals, &string);
    if (stat == EGADS_SUCCESS)
      if (atype == ATTRSTRING)
        if (strcmp(string, "Quad") == 0) quad = 1;

    /* get faces */
    //for (i = 0; i < bodydata[ibody].nfaces; i++) {
    for (i = FRONTFACE; i < FRONTFACE+1; i++) {
      stat = EG_getTessFace(bodydata[ibody].tess, i+1, &len,
          &xyzs, &uvs, &ptype, &pindex, &ntri,
          &tris, &tric);
      if (i == FRONTFACE) {
        bodydata[ibody].plen = len;
        bodydata[ibody].uvs  = (double *) EG_alloc(2*len*sizeof(double));
        if (bodydata[ibody].uvs != NULL)
          for (j = 0; j < 2*len; j++) bodydata[ibody].uvs[j] = uvs[j];
      }
      if (stat != EGADS_SUCCESS) continue;
      sprintf(gpname, "Body %d Face %d", ibody+1, i+1);
      stat = wv_setData(WV_REAL64, len, (void *) xyzs,  WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_INT32, 3*ntri, (void *) tris, WV_INDICES, &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      color[0]  = 1.0;
      color[1]  = ibody;
      color[1] /= nbody;
      color[2]  = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      if (quad == 0) {
        for (nseg = j = 0; j < ntri; j++) {
          for (k = 0; k < 3; k++) {
            if (tric[3*j+k] < j+1) {
              nseg++;
            }
          }
        }
      } else {
        for (nseg = j = 0; j < ntri/2; j++) {
          for (k = 0; k < 4; k++) {
            if (tric[6*j+neigq[k]] < 2*j+1) {
              nseg++;
            }
          }
        }
      }
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d Sides!\n", nseg);
        continue;
      }
      if (quad == 0) {
        for (nseg = j = 0; j < ntri; j++)
          for (k = 0; k < 3; k++)
            if (tric[3*j+k] < j+1) {
              segs[2*nseg  ] = tris[3*j+sides[k][0]];
              segs[2*nseg+1] = tris[3*j+sides[k][1]];
              nseg++;
            }
      } else {
        for (nseg = j = 0; j < ntri/2; j++)
          for (k = 0; k < 4; k++)
            if (tric[6*j+neigq[k]] < 2*j+1) {
              segs[2*nseg  ] = tris[6*j+sideq[k][0]];
              segs[2*nseg+1] = tris[6*j+sideq[k][1]];
              nseg++;
            }
      }
      stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_LINDICES, &items[3]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 3!\n", i, gpname);
      free(segs);
      //    color[0] = color[1] = color[2] = 0.8;
      color[0] = color[1] = color[2] = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_LCOLOR, &items[4]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 4!\n", i, gpname);
      stat = wv_addGPrim(cntxt, gpname, WV_TRIANGLE,
          WV_ON|WV_ORIENTATION, 5, items);
      if (stat < 0)
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
      if (stat > 0) ngp = stat+1;
      sum += ntri;
      printf("BODY DATA VERTEX %d \n",bodydata[0].plen);
    }
    // get edges //
    color[0] = color[1] = 0.0;
    color[2] = 1.0;
    for (i = 0; i < bodydata[ibody].nedges; i++) {
      stat  = EG_getTessEdge(bodydata[ibody].tess, i+1, &len, &xyzs, &ts);
      if (stat != EGADS_SUCCESS) continue;
      if (len == 0) continue;
      nseg = len-1;
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d segments for Body %d Edge %d\n",
            nseg, ibody, i+1);
        continue;
      }
      for (j = 0; j < len-1; j++) {
        segs[2*j  ] = j + 1;
        segs[2*j+1] = j + 2;
      }

      sprintf(gpname, "Body %d Edge %d", ibody+1, i+1);
      stat = wv_setData(WV_REAL64, len, (void *) xyzs, WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS,   &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_INDICES,  &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      free(segs);
      stat = wv_addGPrim(cntxt, gpname, WV_LINE, WV_ON, 3, items);
      if (stat < 0) {
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
      } else {
        if (cntxt != NULL)
          if (cntxt->gPrims != NULL) {
            cntxt->gPrims[stat].lWidth = 1.5;
            if (wv_addArrowHeads(cntxt, stat, 0.05, 1, &nseg) != 0)
              printf(" wv_addArrowHeads Error\n");
          }
        ngp = stat+1;
      }
    }
  }

  printf(" ** %d gPrims with %d triangles **\n", ngp, sum);
  // start the server code //

  stat = 0;
  if (wv_startServer(7681, NULL, NULL, NULL, 0, cntxt) == 0) {

    // we have a single valid server -- stay alive a long as we have a client //
    while (wv_statusServer(0)) {
      usleep(500000);
      if (stat == 0) {
        if (startapp != NULL) system(startapp);
        stat++;
      }
    }
  }


  wv_cleanupServers();
  /* finish up */
  for (ibody = 0; ibody < nbody; ibody++) {
    EG_deleteObject(bodydata[ibody].tess);
    EG_free(bodydata[ibody].edges);
    EG_free(bodydata[ibody].faces);
    EG_free(bodydata[ibody].uvs);
  }
  free(bodydata);

  printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  printf(" EG_close          = %d\n", EG_close(context));


  return 0;
}




/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *wsi, char *text, /*@unused@*/ int lena)
{
  int          i, j, iBody, ient, stat, nattr, atype, alen;
  const int    *pints;
  char         tag[5], gpname[33];
  const char   *name, *pstr;
  double       *xyzs, result[18];
  const double *preals;
  ego          obj;
  wvData       item;
  FILE *gpfile;
  if (strncmp(text,"Picked: ", 8) == 0) {
    iBody = 0;
    sscanf(&text[12], "%d %s %d", &iBody, tag, &ient);
    if (iBody != 0) {
      printf(" Picked: iBody = %d, type = %s, index = %d\n", iBody, tag, ient);
      if (strcmp(tag,"Face") == 0) {
        obj = bodydata[iBody-1].faces[ient-1];
      } else {
        obj = bodydata[iBody-1].edges[ient-1];
      }
      nattr = 0;
      stat  = EG_attributeNum(obj, &nattr);
      if ((stat == EGADS_SUCCESS) && (nattr != 0)) {
        for (i = 1; i <= nattr; i++) {
          stat = EG_attributeGet(obj, i, &name, &atype, &alen,
              &pints, &preals, &pstr);
          if (stat != EGADS_SUCCESS) continue;
          printf("   %s: ", name);
          if ((atype == ATTRREAL) || (atype == ATTRCSYS)) {
            for (j = 0; j < alen; j++) printf("%lf ", preals[j]);
          } else if (atype == ATTRSTRING) {
            printf("%s", pstr);
          } else {
            for (j = 0; j < alen; j++) printf("%d ", pints[j]);
          }
          printf("\n");
        }
      }
    }
  }
  printf("BODY DATA VERTEX %d \n",bodydata[0].plen);
  xyzs = (double *) EG_alloc(3*bodydata[0].plen*sizeof(double));
  for (i = 0; i < bodydata[0].plen; i++) {
    stat = EG_evaluate( bodydata[0].faces[FRONTFACE], &bodydata[0].uvs[2*i], result);
    if (stat != EGADS_SUCCESS)
      printf(" Error: EG_evaluate = %d for %d %d\n", stat, iBody+1, i+1);
    xyzs[3*i  ] = result[0];
    xyzs[3*i+1] = result[1];
    xyzs[3*i+2] = result[2];
    printf(" V %d (%lf,%lf) = %lf  %lf  %lf\n",i+1,bodydata[0].uvs[2*i],bodydata[0].uvs[2*i+1],xyzs[3*i],xyzs[3*i+1],xyzs[3*i+2] );
  }

  /* play around with first Face UVs */
  if ((strcmp(text,"finer") == 0) || (strcmp(text,"coarser")== 0)) {
    for (iBody = 0; iBody < nbody; iBody++) {
      if (bodydata[iBody].uvs == NULL) continue;
      sprintf(gpname, "Body %d Face %d", iBody+1, FRONTFACE +1);
      ient = wv_indexGPrim(cntxt, gpname);
      if (ient < 0) continue;
      updateUVs(bodydata[iBody].tess, FRONTFACE +1, bodydata[iBody].faces[FRONTFACE],
          bodydata[iBody].uvs,xyzs);
      EG_free(xyzs);
      exit(1);
      char buffer[32];
      snprintf(buffer, sizeof(char) * 32, "TESSELATIONS/TESS%i.txt", TESSCOUNT);
      gpfile = fopen(buffer,"w");
      for (i = 0; i < bodydata[iBody].plen; i++) {
        stat = EG_evaluate( bodydata[iBody].faces[FRONTFACE],
            &bodydata[iBody].uvs[2*i], result);
        if (stat != EGADS_SUCCESS)
          printf(" Error: EG_evaluate = %d for %d %d\n", stat, iBody+1, i+1);
        xyzs[3*i  ] = result[0];
        xyzs[3*i+1] = result[1];
        xyzs[3*i+2] = result[2];
        fprintf(gpfile,"%lf %lf %lf %lf %lf\n",bodydata[iBody].uvs[2*i],bodydata[iBody].uvs[2*i+1],xyzs[3*i],xyzs[3*i+1],xyzs[3*i+2]);
      }
      fclose(gpfile);
      ++TESSCOUNT;
      stat = wv_setData(WV_REAL64, bodydata[iBody].plen, (void *) xyzs,
          WV_VERTICES, &item);
      if (stat < 0) printf(" wv_setData = %d for %s!\n", iBody, gpname);
      wv_adjustVerts(&item, focus);
      stat = wv_modGPrim(cntxt, ient, 1, &item);
      if (stat < 0) printf(" wv_modGPrim = %d for %s!\n", iBody, gpname);
    }
  }
  EG_free(xyzs);
}


