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
#include "wsserver.h"
#include"egads.h"
#define PI              3.1415926535897931159979635
#define ERRCTT  3.554147
#define TRYSWAPPING         0
#define TRYCOLLAPSING       5
#define TRYSPLITTING        1
#define TRYSWAPANDSPLIT     4
#define TRYDOUBLESPLIT      2
#define TRYCOLLAPSESPLIT    3
#define OPTYPES 6
#define MAXVALENCE 100
#define ANGLEPASS 0.95*M_PI
int AV = 0 ;
int SO = 0;
int RES = 0 ;
int PLANAR = 0;
int ADDP   = 0;
int DB = 0 ;
int QC = 0 ;
int TOTSWAPS         = 0;
int OVERALLSWAPS     = 0;
int TOTCOLLAPSES     = 0;
int OVERALLCOLLAPSES = 0;
int TOTSPLITS        = 0;
int OVERALLSPLITS    = 0;
int DS = 0;
int PRINTAREA = 0;
int INVALIDCOUNT = 0;
int PLANESURFACE = 0;
int NSWAP = 0 ;
#define EPSK 0.001
double  GEPS, KMIN, KMAX, HMAX, HMIN;
#define MAXRATIO 2.0
int OPTCOUNT  = 0;
int FCOUNT = -1;
int CONVC = 0;
double ANGLE_ERR  = 0.95 * M_PI;
#define EPSAREA 1.E-02
#define DEPS 1.E-08
#define EPSDIR 1.e-06
#define OPTITOL 0.1
#define DTOL 1.e-9

int  IOSTATUS = 0;
//#define DEBUG_AREA
#define DEBUGG
int NOUTMESH = 1000;
int NOUT     = 100;

#define NLOPTMAXEVAL 100000

#define MULTIOPTI 1
int VS = 0;

#define MAX_LINKS 100
#define EPS  1.E-05

typedef struct{
	int    *verts, *quads;
	int  nV, nQ; // nV = origin(1) + peaks (n)
} vStar;

typedef struct {
	int     *quadIdx, *quadAdj, **valence, *vType, *remQuads, *remVerts, curvIdx[2];
	int      nQ, nS, q0, v0, totVerts,setIndex;
	ego      face;
	double  *xyzs, *oriXYZ, *uvs;
	vStar   **star;
} quadMap;


typedef struct {
	ego      *faces;
	ego      *edges;
	ego       body;
	ego       tess;
	int       mtype;
	int       nfaces;
	int       nedges;
	int       plen;
	quadMap   **qm;
} bodyData;


#define DEBUG
int GLO = 0 ;

int QUADTESS  = 0;
int TESSCOUNT = 0;
#define ANGPASS 0.4
#define  FRONTFACE 0
int      SCOUNT = 0;
#define FORCETRIANGULATION
#define FORCEQUADS

/* globals used in these functions */
static wvContext *cntxt;
static bodyData  *bodydata;
static int       nbody, sscount = 0 ;
static float     focus[4];
FILE *filOpti;
char optiName[33], buffer[33];


static void
print_quad_specs(quadMap *qm, int id) ;

static void sampleNormalPlane(double normal[], double point[], int vID, quadMap *qm) ;
static void NLOPTtermination(int n) ;

static void print_starFile(vStar *star, quadMap *qm) ;
static void print_star(vStar *star) ;
static void
printMeshStats(quadMap *qm, int sweep);
static void
print_mesh(quadMap *qm,char[] );
static void
EG_unitVector(double *v, double *norm) {
	double n;
	n     = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	n     = sqrt(n);
	if(n > DEPS) {
		v[0] /=n; v[1] /=n; v[2] /=n;
	}
	else {
		v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
	}
	*norm = n;
	return;
}


static int EG_createQuadMap(bodyData *bodydata);
static int EG_createQuadMap(bodyData *bodydata)
{
	int           nseg, f, stat, j, q, i, auxID, k, kk, kOK,  iper;
	const int    *tris, *tric, *ptype, *pindex;
	int           len, ntri, nquad, vID0, vID1, vS, qS;
	double        dist, surfArea, perimetre, aux, uvbox[2*2], uvLims[2*4],  verts[18*4], daux, dplane = 0.0;
	const double *xyzs, *uvs;
	static int    qV[6]       = { 0, 1, 2, 5, 0, 1};
	static int    qLoop[5]    = { 0, 1, 2, 3, 0};
	bodydata->qm = (quadMap**) EG_alloc(bodydata->nfaces*sizeof(quadMap*));
	if (bodydata->qm == NULL ) {
		return  EGADS_MALLOC;
	}
	for ( f = 0 ; f < bodydata->nfaces; ++f) {
		bodydata->qm[f] = (quadMap*) EG_alloc(sizeof(quadMap));
		if (bodydata->qm[f] == NULL )
			return EGADS_MALLOC;
		stat = EG_getTessFace(bodydata->tess, f+1, &len,
				&xyzs, &uvs, &ptype, &pindex, &ntri,
				&tris, &tric);
		nquad                       = (int)ntri/2;
		vS                          = 2 * len;
		qS                          = 2 * nquad;
		bodydata->qm[f] -> xyzs     = (double*) EG_alloc(3*vS *sizeof(double));
		bodydata->qm[f] -> oriXYZ   = (double*) EG_alloc(3*vS *sizeof(double));
		bodydata->qm[f] -> uvs      = (double*) EG_alloc(2*vS *sizeof(double));
		bodydata->qm[f] -> vType    = (int*)    EG_alloc(  vS *sizeof(   int));
		bodydata->qm[f] -> quadIdx  = (int*)    EG_alloc(4*qS *sizeof(   int));
		bodydata->qm[f] -> quadAdj  = (int*)    EG_alloc(4*qS *sizeof(   int));
		bodydata->qm[f] -> remQuads = (int*)    EG_alloc(  qS *sizeof(   int));
		bodydata->qm[f] -> remVerts = (int*)    EG_alloc(  vS *sizeof(   int));
		bodydata->qm[f] -> star     = (vStar**) EG_alloc(  vS *sizeof(vStar*));
		bodydata->qm[f] -> valence  =  (int **) EG_alloc(  vS *sizeof(  int*));
		if (bodydata->qm[f]-> quadIdx  == NULL  ||bodydata->qm[f] ->  quadAdj == NULL ||bodydata->qm[f] ->  xyzs     == NULL ||
				bodydata->qm[f] -> oriXYZ   == NULL ||bodydata->qm[f] ->  uvs     == NULL ||bodydata->qm[f] ->  vType    == NULL ||
				bodydata->qm[f] ->  remQuads == NULL ||
				bodydata->qm[f] -> remVerts == NULL ||bodydata->qm[f] ->  star    == NULL ||bodydata->qm[f] ->  valence  == NULL )
			return EGADS_MALLOC;
		bodydata->qm[f] -> face        = bodydata->faces[f];
		bodydata->qm[f] -> totVerts    = len;
		bodydata->qm[f] -> v0          = len;
		bodydata->qm[f] -> nQ          = nquad;
		bodydata->qm[f] -> q0          = nquad;
		bodydata->qm[f] -> nS          = 0;
		bodydata->qm[f] -> remQuads[0] = 0;
		bodydata->qm[f] -> remVerts[0] = 0;
		for (j = 0; j < len; j++) {
			bodydata->qm[f] -> valence[j] = (int *) EG_alloc((2 + MAXVALENCE)*sizeof(int));
			if (bodydata->qm[f] -> valence[j] == NULL ) {
				return EGADS_MALLOC;
			}
			bodydata->qm[f] ->  valence[j][0]   = 0;
			bodydata->qm[f] ->  xyzs  [3*j    ] = xyzs  [3*j   ];
			bodydata->qm[f] ->  xyzs  [3*j + 1] = xyzs  [3*j +1];
			bodydata->qm[f] ->  xyzs  [3*j + 2] = xyzs  [3*j +2];
			dplane                             += bodydata->qm[f] ->  xyzs  [3*j + 2];
			bodydata->qm[f] ->  oriXYZ[3*j    ] = bodydata->qm[f] -> xyzs[3*j    ];
			bodydata->qm[f] ->  oriXYZ[3*j + 1] = bodydata->qm[f] -> xyzs[3*j + 1];
			bodydata->qm[f] ->  oriXYZ[3*j + 2] = bodydata->qm[f] -> xyzs[3*j + 2];
			bodydata->qm[f] ->  uvs   [2*j    ] = uvs   [2*j   ];
			bodydata->qm[f] ->  uvs   [2*j + 1] = uvs   [2*j +1];
			bodydata->qm[f] ->  vType [j      ] = ptype [j     ];
		}
		if ( dplane < 1.e-06) PLANAR = 1;
		for (j = 0; j < nquad; j++)
			for ( k = 0; k < 4; ++k)
				bodydata->qm[f] ->  quadIdx[4*j + k ] = tris[6*j + qV[k+1]];
		for ( j = 0; j < nquad; j++) {
			kk  = 0;
			kOK = 0;
			for ( q = 0 ; q < nquad; ++q) {
				if ( q == j || q == -1) continue;
				for ( k = 0 ; k < 4; ++k ) {
					if ( (bodydata->qm[f] -> quadIdx[4*j + qLoop[kk    ]] == bodydata->qm[f] -> quadIdx[4*q + qLoop[k    ]] ||
							bodydata->qm[f] -> quadIdx[4*j + qLoop[kk    ]] == bodydata->qm[f] -> quadIdx[4*q + qLoop[k + 1]]
					) && (bodydata->qm[f] -> quadIdx[4*j + qLoop[kk + 1]] == bodydata->qm[f] -> quadIdx[4*q + qLoop[k    ]] ||
							bodydata->qm[f] -> quadIdx[4*j + qLoop[kk + 1]] == bodydata->qm[f] -> quadIdx[4*q + qLoop[k + 1]]) )
					{
						bodydata->qm[f] ->  quadAdj[4*j + kk] = q + 1;
						++kk;
						q   = -1;
						kOK =  1;
						k   =  4;
						if (kk == 4)  q = nquad;
					}
				}
				if ( (kOK == 0) && (q == nquad -1) ){
					bodydata->qm[f] ->  quadAdj[4*j + kk] = -1;
					q                       = -1;
					++kk;
					if (kk == 4) break;
				}
				else  kOK = 0 ;
			}
		}
		daux  = sqrt (1.0 - 1.2 * EPSK );
		GEPS  = sqrt ( 40.0*( 1.0 - daux) ) ;
		GEPS *= (1.0 - EPSK);
		HMIN  = 100000.0;
		HMAX  = 0.0;
		for ( j = 0 ; j < nquad; ++j) {
			for ( k = 0 ; k < 4; ++k ) {
				if  (bodydata->qm[f] -> quadAdj[4*j + k] == -1 ) {
					vID0 = bodydata -> qm[f] -> quadIdx[4*j +  k       ] - 1;
					vID1 = bodydata -> qm[f] -> quadIdx[4*j + (k + 1)%4] - 1;
					dist = 0.0;
					for ( i = 0 ; i < 3; ++i)
						dist += (bodydata -> qm[f] -> xyzs[3*vID1 + i ] - bodydata -> qm[f] -> xyzs[3*vID0 + i ]) * (bodydata -> qm[f] -> xyzs[3*vID1 + i ] - bodydata -> qm[f] -> xyzs[3*vID0 + i ]);
					dist = sqrt(dist);
					if ( dist < HMIN) HMIN = dist;
					if ( dist > HMAX) HMAX = dist;
				}
			}
		}
		if ( fabs(HMIN - HMAX ) < DEPS) HMIN *=0.5;
		printf(" H MIN MAX %lf %lf\n",HMIN, HMAX);
		KMIN  = GEPS/HMAX;
		KMAX  = GEPS/HMIN;
		printf(" CURVATURE MIN MAX %lf %lf\n", KMIN, KMAX);
		stat = EG_getArea (bodydata->qm[f] -> face, NULL, &surfArea);
		surfArea  /= (double)nquad;
		printf(" ELEMENT AREA %lf \n",surfArea);
		for ( j = 0 ; j < nquad; ++j) {
			for ( k = 0 ; k < 4; ++k) {
				if ( bodydata->qm[f] -> quadAdj[4*j +k] > j +1  || bodydata->qm[f] -> quadAdj[4*j +k] == -1) {
					//print_quad_specs(bodydata->qm[f], j+1);
					auxID = bodydata->qm[f] -> quadIdx [4*j + qLoop[k]] -1;
					bodydata->qm[f]   -> valence[auxID][1] = j + 1;
					bodydata->qm[f]   -> valence[auxID][2 + bodydata->qm[f] -> valence[auxID][0] ] = bodydata->qm[f] -> quadIdx[4*j+ qLoop[k+1]];
					++bodydata->qm[f] -> valence[auxID][0];
					//  printf(" ASSIGNING %d  TO QUAD %d\n",auxID + 1, j+ 1);
					auxID = bodydata->qm[f] -> quadIdx [4*j + qLoop[k+1]] -1;
					bodydata->qm[f]   -> valence[auxID][1] = j + 1;
					bodydata->qm[f]   -> valence[auxID][2 + bodydata->qm[f] -> valence[auxID][0] ] = bodydata->qm[f] -> quadIdx[4*j+ qLoop[k]];
					++bodydata->qm[f] -> valence[auxID][0];
					//    printf(" ASSIGNING %d TO QUAD %d\n",auxID + 1, j+ 1);
				}
			}
		}
	}
	// set vertex movement based on surface area
	return EGADS_SUCCESS;
}

static void EG_destroyQuadMap(quadMap **qm, int nfaces ) {
	int i, j;
	for ( i = 0 ; i < nfaces; ++i) {
		if (qm[i]) {
			if (qm[i] -> valence) {
				for ( j = 0 ; j < qm[i] -> totVerts; ++j)
					EG_free( qm[i] -> valence[j]) ;
				EG_free( qm[i] -> valence) ;
			}
			if (qm[i] -> star[j]){
				EG_free(qm[i] -> star[j] -> verts    );
				EG_free(qm[i] -> star[j] -> quads    );
				EG_free(qm[i] -> star[j]    );
			}
			EG_free(qm[i] -> quadIdx );
			EG_free(qm[i] -> quadAdj );
			EG_free(qm[i] -> uvs     );
			EG_free(qm[i] -> xyzs    );
			EG_free(qm[i] -> oriXYZ  );
			EG_free(qm[i] -> remQuads);
			EG_free(qm[i] -> remVerts);
			EG_free(qm[i] -> vType   );

			EG_free(qm[i]            );
		}
	}
	EG_free(qm);
}



static int EG_averageCoords(quadMap *qm, int vID ) {
	int i, auxID, stat;
	double uvc[2], eval[18];
	uvc[0] = 0.0; uvc[1] = 0.0;
	printf(" RECALCULATING VERTEX %d\n",vID);
	for ( i = 0 ; i < qm ->valence[vID - 1][0]; i++){
		auxID   = qm -> valence[vID - 1][2 + i] - 1;
		uvc[0] += qm -> uvs    [2*auxID    ];
		uvc[1] += qm -> uvs    [2*auxID + 1];
	}
	qm -> uvs[2*(vID -1 )    ] = uvc[0] / (double)qm ->valence[vID - 1][0];
	qm -> uvs[2*(vID -1 ) + 1] = uvc[1] / (double)qm ->valence[vID - 1][0];
	stat                   = EG_evaluate(qm->face, &qm ->uvs[2*(vID -1)], eval);
	if ( stat != EGADS_SUCCESS) {
		printf(" EG_EVALUATE IN AVERAGE COORDS = %d\n", stat);
		return stat;
	}
	qm -> xyzs[3*(vID -1)    ]  = eval[0];
	qm -> xyzs[3*(vID -1) + 1]  = eval[1];
	qm -> xyzs[3*(vID -1) + 2]  = eval[2];
	return EGADS_SUCCESS;
}


static int
EG_splittingOperation(quadMap *qm, vStar *star, int poly[], int id0, int distanceToId0, int *);
static int EG_removeDoublet(quadMap *qm, int vID) ;
static int EG_swapAndSplit(quadMap *qm, int vID);
static int EG_edgeSwap(quadMap *qm, int vID);
static int EG_vertexCollapse(quadMap *qm, int vID);
static int EG_vertexSplit(quadMap *qm, int vID);
static int EG_collapseAndSplit(quadMap *qm, int qID) ;
static int EG_buildStar(quadMap *qm, vStar **star, int vID );
double NLOPT_objEqualAnglesTess(unsigned n, const double *uvs, double *grad, void *inQM);
static int optimize_angles(quadMap *qm, int nV, int *vID);
static int quad_algebraic_area(quadMap *qm, int qID, int vID,  int *validArea, double *angles,double *area) ;

static int
EG_getCommonEdge(quadMap *qm, int q1, int q2, int *edge) {
	int i, j, e = 0, v0, v1;
	edge[0] = -1; edge[1] = -1;
	edge[2] = -1; edge[3] = -1;
	for ( i = 0 ; i < 4; ++i) {
		if ( qm -> quadAdj[4*(q1 - 1) + i ] == q2) {
			edge[e++] = i;
			v1 = qm -> quadIdx[4*(q1 - 1) + i    ];  //quad q2 loops in opposite direction. AB in q1 is BA in q2
			v0 = qm -> quadIdx[4*(q1 - 1) + (i + 1)%4];
			for ( j = 0 ; j < 4; ++j) {
				if ( (qm -> quadAdj[4*(q2 - 1) + j ] == q1) && (qm -> quadIdx[4*(q2 - 1) + j ] == v0) && (qm -> quadIdx[4*(q2 - 1) + (j + 1)%4] == v1) )
					edge[e++] = j;
			}
		}
	}
	if ( edge[0] != -1 && edge[1] != -1)  {
		return EGADS_SUCCESS;
	}
	return EGADS_EMPTY;
}

static int
EG_getAdjacentToVertexPair(quadMap *qm, int qID, int v1, int v2, int *adj) {
	int i, j, auxQ, id[2];
	// print_quad_specs(qm, qID);
	for ( i = 0 ; i < 4; ++i) {
		id[0] = - 1; id[1] = -1;
		if      ( qm -> quadIdx[4*(qID - 1) + i ] == v1 ) {
			id[0] = v1; id[1] = v2;
		}
		else if ( qm -> quadIdx[4*(qID - 1) + i ] == v2 ) {
			id[0] = v2; id[1] = v1;
		}
		if ( id[0] != -1 ) {
			auxQ = qm -> quadAdj[4*(qID - 1) + i ] - 1;
			if ( auxQ >= 0 ) {
				//  print_quad_specs(qm, auxQ + 1);
				for ( j = 0 ; j < 4; ++j)
					if ( qm -> quadIdx[4*auxQ + j ] == id[1]) {
						adj[0] = i; adj[1] = auxQ + 1;
						return EGADS_SUCCESS;
					}
			}
		}
	}
	return EGADS_EMPTY;
}

static int
resizeQm(quadMap *qm);
static int
backupQm(quadMap *qm, int *nV, int *nQ, int **quadIdx, int **quadAdj, int **vtype, int ***valences, int **Qrem, int **Vrem, double **uvs) ;

static int
resetQm(quadMap *qm, int nV, int nQ, int *quadIdx, int *quadAdj, int *vtype, int **valences, int *Qrem, int *Vrem, double *uvs) {
	int i, j, k, stat;
	double eval[18];
	qm -> nQ       = nQ;
	qm -> totVerts = nV;
	for ( j = 0 ; j < qm -> nQ; ++j) {
		for ( k = 0 ; k  <4; ++k) {
			qm -> quadIdx[4*j + k] = quadIdx[4*j+ k];
			qm -> quadAdj[4*j + k] = quadAdj[4*j+ k];
		}
	}
	for ( j = 0 ; j < qm -> totVerts; ++j) {
		stat = EG_evaluate(qm -> face, &uvs[2*j], eval);
		if ( stat != EGADS_SUCCESS) return stat;
		qm -> vType[  j    ] = vtype[  j    ];
		qm -> uvs  [2*j    ] = uvs  [2*j    ];
		qm -> uvs  [2*j + 1] = uvs  [2*j + 1];
		qm -> xyzs [3*j    ] = eval [0];
		qm -> xyzs [3*j + 1] = eval [1];
		qm -> xyzs [3*j + 2] = eval [2];
		if ( qm -> vType[j] == -2) {
			qm -> valence[j][0] = -1;
		} else {
			for ( i = 0 ; i < valences[j][0] + 2; ++i)
				qm -> valence[j][i] = valences[j][i];
		}
	}
	for ( i = 0 ; i <= Qrem[0]; i++ ) {
		qm -> remVerts[i] = Vrem[i];
		qm -> remQuads[i] = Qrem[i];
	}
	snprintf(buffer,33,"RESETP%d",RES);
	print_mesh(qm, buffer);
	++RES;
	return EGADS_SUCCESS;
}

static int
resizeQm(quadMap *qm) {
	int stat, nV, vRem, nQ, qRem,  i, j, k, *vpiv = NULL, *qpiv = NULL,vB, qB, aux;
	int  *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, *locVtype, **locVal = NULL, *locQrem = NULL, *locVrem = NULL;
	double *uv0= NULL, eval[18];
	nV = 0 ; vRem = 0 , nQ = 0 , qRem = 0;
	for ( i = 0 ; i < qm -> totVerts; i++) {
		if(qm -> vType[i] != -2 ) nV++;
		else vRem++;
	}
	for ( i = 0 ; i < qm -> nQ; i++) {
		if(qm -> quadIdx[4*i] != -2 ) nQ++;
		else qRem++;
	}
	if ( vRem != qRem ) {
		if ( qRem > 0 ) {
			printf(" I HAVE %d REMOVED VERTICES BUT ACTUALLY %d REMOVED QUADS!!!!!\n ", vRem, qRem);exit(1);
		}
	}
	if ( vRem == 0 )return EGADS_SUCCESS;
	stat = backupQm ( qm, &vB, &qB, &locIdx, &locAdj, &locVtype, &locVal,&locQrem, &locVrem,  &uv0);
	if ( stat != EGADS_SUCCESS) return stat;
	vpiv = (int*) EG_alloc(qm -> totVerts * sizeof(int));
	qpiv = (int*) EG_alloc(qm -> nQ       * sizeof(int));
	if (vpiv == NULL || qpiv == NULL) return EGADS_MALLOC;
	for ( j = i = 0 ; i < qm -> nQ; i++) {
		if(qm -> quadIdx[4*i] != -2 ) {
			qpiv[i] = j;
			++j;
		} else qpiv[i] = -2;
	}
	for ( j = i = 0 ; i < qm -> totVerts; i++) {
		if(qm -> vType[i] != -2 ) {
			vpiv[i] = j;
			++j;
		} else vpiv[i] = -2;
	}
	nV             = qm -> totVerts - vRem;
	nQ             = qm -> nQ       - qRem;
	qm -> xyzs     = EG_reall(qm -> xyzs,     3*nV*sizeof(double));
	qm -> uvs      = EG_reall(qm -> uvs ,     2*nV*sizeof(double));
	qm -> vType    = EG_reall(qm -> vType,      nV*sizeof(int));
	qm -> quadIdx  = EG_reall(qm -> quadIdx,  4*nQ*sizeof(int));
	qm -> quadAdj  = EG_reall(qm -> quadAdj,  4*nQ*sizeof(int));
	qm -> valence =  EG_reall(qm -> valence ,   nV*sizeof(  int*));
	if (qm->  quadIdx  == NULL ||qm ->  quadAdj == NULL ||qm ->  xyzs     == NULL ||
			qm -> uvs      == NULL ||qm ->  vType   == NULL || qm ->  star    == NULL ||qm ->  valence  == NULL ) {
		for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
		EG_free(locVal);
		EG_free(locIdx);
		EG_free(locAdj);
		EG_free(locVtype);
		EG_free(locQrem);
		EG_free(locVrem);
		EG_free(uv0);
		return EGADS_MALLOC;
	}
	for ( j = i = 0 ; i < qm -> nQ; i++){
		//	printf(" Q %d  -> ACTUALLY %d \n",i, qpiv[i]);
		if ( qpiv[i] == -2) continue;
		for ( k = 0 ; k < 4; ++k) {
			qm -> quadIdx[4*qpiv[i] + k ] = vpiv[ locIdx[4*i + k] - 1] + 1;
			if ( locAdj[4*i + k] == -1 ) qm -> quadAdj[4*qpiv[i] + k ] = - 1;
			else  qm -> quadAdj[4*qpiv[i] + k ] = qpiv[ locAdj[4*i + k] - 1] + 1;
			//	printf(" Q %d  VERT %d = %d   \n", qpiv[i] +1, qm -> quadIdx[4*qpiv[i] + k ] , qm -> quadAdj[4*qpiv[i] + k ] );

		}
	}
	for ( i = 0 ; i < qm -> totVerts; i++) {
		if ( vpiv[i] == -2 ) continue;
		j = vpiv[i];
		printf("i %d  j %d TOT %d TOT OLD %d \n",i, j, nV, qm -> totVerts);
		qm -> vType[j      ] = locVtype[  i    ];
		qm -> uvs  [2*j    ] = uv0     [2*i    ];
		qm -> uvs  [2*j + 1] = uv0     [2*i + 1];
		stat                 = EG_evaluate(qm -> face, &qm -> uvs[2*j], eval);
		qm -> xyzs[3*j]      = eval[0]; qm -> xyzs[3*j + 1] = eval[1]; qm -> xyzs[3*j + 2] = eval[2];
		qm -> valence[j]     = EG_reall (qm->valence[j], (locVal[i][0] + 2) * sizeof(int));
		qm -> valence[j][0]  = locVal[i][0];
		qm -> valence[j][1]  = qpiv[locVal[i][1] - 1] + 1;
		for ( k = 0 ; k < locVal[i][0]; ++k) qm ->valence[j][2 + k] = vpiv[locVal[i][2 + k] -1] + 1;
	}
	//#ifdef DEBUG
	int OK;
	printf(" DEBUG RESIZE\n");
	for ( i = 0 ; i < nQ; i++) {
		print_quad_specs(qm, i+1);
		for ( j = 0 ; j < 4; ++j ) {
			//		printf(" Q = %d  ADJ %d  POINTS BACK ?\t", i+1, qm -> quadAdj[4*i + j]);
			if ( qm -> quadAdj[4*i + j ] == -1 ) continue;
			OK = 0 ;
			for ( k = 0 ; k < 4; ++k ) {
				if ( qm -> quadAdj[ 4* ( qm -> quadAdj[4*i + j] -1) + k] == i + 1) {
					OK = 1;
					break;
				}
			}
			if ( OK == 0) {
				printf(" ========== ERROR =============\n");
				printf(" QUAD %d POINTS AT %d BUT NOT THE OTHER WAY AROUND \n", i + 1, qm -> quadAdj[4*i + j ]);
				print_quad_specs(qm, i+1);
				print_quad_specs(qm, qm -> quadAdj[4*i + j]);
				exit (1);
			}
		}
	}
	for ( i = 0 ; i < nV; i++) {
		OK = 0 ;
		for ( j = 0 ; j < qm -> valence[i][0]; j++) {
			aux = qm -> valence[i][2 + j] - 1;
			for ( k = 0 ; k < qm -> valence[aux][0] ; ++k ) {
				if ( qm -> valence[aux][k + 2] == i + 1) {
					OK = 1;
					break;
				}
			}
			if ( OK == 0) {
				printf(" ========== ERROR =============\n");
				printf(" VERTEX %d has valence %d BUT NOT THE OTHER WAY AROUND \n", i + 1, aux + 1);
				for ( k = 0 ; k < qm -> valence[aux][0] ; ++k )
					printf(" VALENCE %d = %d \n", k, qm -> valence[aux][k + 2]);
				exit (1);
			}
		}
	}
	//#endif
	EG_free(vpiv);
	EG_free(qpiv);
	for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
	EG_free(locVal);
	EG_free(locIdx);
	EG_free(locAdj);
	EG_free(locVtype);
	EG_free(locQrem);
	EG_free(locVrem);
	EG_free(uv0);
	qm -> totVerts    = nV;
	qm -> nQ          = nQ;
	qm -> remQuads[0] = 0 ;
	qm -> remVerts[0] = 0 ;
	print_mesh(qm,"AFTERRED");
	return EGADS_SUCCESS;
}




static int
backupQm(quadMap *qm, int *nV, int *nQ, int **quadIdx, int **quadAdj, int **vtype, int ***valences, int **Qrem, int **Vrem, double **uvs) {
	int i, j, v;
	snprintf(buffer,33,"BACKUP%d",RES);
	print_mesh(qm, buffer);
	++RES;
	if ( qm -> remQuads[0] != qm -> remVerts[0]){
		printf(" MISS MATCH OF REMOVED QUADS AND VERTS. THIS IS A TOTAL DISSASTER!\n");
		return EGADS_MALLOC;
	}
	*nV        = qm -> totVerts;
	*nQ        = qm -> nQ;
	*quadIdx   = (int*)    EG_alloc(4*qm -> nQ              *sizeof(int   ));
	*quadAdj   = (int*)    EG_alloc(4*qm -> nQ              *sizeof(int   ));
	*vtype     = (int*)    EG_alloc(  qm -> totVerts        *sizeof(int   ));
	*valences  = (int **)  EG_alloc(  qm -> totVerts        *sizeof(int*  ));
	*uvs       = (double*) EG_alloc(2*qm -> totVerts        *sizeof(double));
	*Qrem      = (int *)   EG_alloc( (qm -> remQuads[0] + 1)*sizeof(int   ));
	*Vrem      = (int *)   EG_alloc( (qm -> remVerts[0] + 1)*sizeof(int   ));
	if ( (*quadIdx) == NULL || (*quadAdj) == NULL || (*valences) == NULL
			|| (*vtype) == NULL || (*uvs) == NULL || (*Vrem == NULL) || (*Qrem == NULL)) return EGADS_MALLOC;
	for ( j = 0; j< qm -> totVerts; ++j) {
		(*valences)[j] = (int *) EG_alloc((2 + 100)*sizeof(int));
		if ((*valences)[j] == NULL ) {
			EG_free(*quadIdx);
			EG_free(*quadAdj);
			EG_free(*vtype);
			EG_free(*valences);
			EG_free(*Vrem);
			EG_free(*Qrem);
			return EGADS_MALLOC;
		}
		for ( i = 0 ; i < qm -> valence[j][0] + 2; ++i) (*valences)[j][i] = qm -> valence[j][i];
	}
	for ( i = 0 ; i < qm ->nQ; ++i) {
		for ( j = 0 ; j < 4; ++j){
			v = qm -> quadIdx[4*i + j];
			if ( v > qm -> totVerts  ) {
				printf(" TOT V %d BACKUP QM = %d  idx  %d = %d -> ", qm -> totVerts, i + 1, j, v);
				exit(1);
			}
			(*quadIdx)[4*i + j] = qm -> quadIdx[4*i + j];
			(*quadAdj)[4*i + j] = qm -> quadAdj[4*i + j];
		}
	}
	for ( i = 0 ; i < qm ->totVerts; ++i) {
		(*uvs)[2*i] = qm -> uvs[2*i]; (*uvs)[2*i + 1] = qm -> uvs[2*i + 1];
		(*vtype)[i] = qm -> vType[i];
	}
	for ( i = 0 ; i <= qm -> remQuads[0]; ++i){
		(*Qrem)[i] = qm -> remQuads[i];
		(*Vrem)[i] = qm -> remVerts[i];
	}
	return EGADS_SUCCESS;
}









static int checkInvalidElement(quadMap *qm, int qID) {
	int    i, v[4], stat, area;
	double uvxyz[5], penalty, angles[4], quadArea;
	PRINTAREA = 1;

	if ( qID == -1) {
		printf(" ARE YOU SURE YOU WANT TO LOOK HERE???\n");
	}
	if ( qm -> quadIdx[4*(qID-1)] == -2 ) {
		printf(" IT DOESN't EXIST. DON'T CHECK\n");
		return EGADS_SUCCESS; // empty quad. probably has been removed elsewhere
	}
	printf("\n\n INVALID ELEMENT Q   %d   \n", qID);
	print_quad_specs(qm, qID);
	for ( i = 0; i < 4; ++i) v[i] = qm -> quadIdx[4*(qID-1) + i];
	snprintf(buffer,33, "AVCOORD_%d",AV);
	++AV;
	print_mesh(qm, buffer);;

	stat = quad_algebraic_area(qm, qID, v[0], &area, angles, &quadArea);
	//printf(" STAT IN QUAD ALG AREA = %d   AREA %d\n", stat, area);
	if ( stat != EGADS_SUCCESS) {
		printf(" STAT %d \n", stat);
		return stat;
	}
	if ( area == 1) {
		printf(" AREA OK. OUT \n");
		return EGADS_SUCCESS;
	}
	snprintf(buffer,33, "AVCOORD_%d",AV);
	++AV;
	print_mesh(qm, buffer);;
	for  ( i = 0; i < 4; i++) {
		printf(" ID %d   ANGLE %lf  VTYPE %d \n", v[i], angles[i], qm -> vType[v[i] - 1]);
		if ( angles[i] >= ANGLE_ERR && qm -> vType[v[i] -1] == -1 ) {
			uvxyz[0] = qm -> uvs [2*(v[i] - 1)    ] ;
			uvxyz[1] = qm -> uvs [2*(v[i] - 1) + 1] ;
			uvxyz[2] = qm -> xyzs[3*(v[i] - 1)    ] ;
			uvxyz[3] = qm -> xyzs[3*(v[i] - 1) + 1] ;
			uvxyz[4] = qm -> xyzs[3*(v[i] - 1) + 2] ;
			//	printf(" ORI %lf  %lf  -> %lf  %lf  %lf \n", qm -> uvs [2*(v[i] - 1)    ], qm -> uvs [2*(v[i] - 1) +1   ],
			//		qm -> xyzs[3*(v[i] - 1)    ] , qm -> xyzs[3*(v[i] - 1)  +1  ] ,qm -> xyzs[3*(v[i] - 1)  +2  ] 		);
			stat = EG_averageCoords(qm, v[i]);
			//printf(" MOD %lf  %lf  -> %lf  %lf  %lf \n", qm -> uvs [2*(v[i] - 1)    ], qm -> uvs [2*(v[i] - 1) +1   ],
			//	qm -> xyzs[3*(v[i] - 1)    ] , qm -> xyzs[3*(v[i] - 1)  +1  ] ,qm -> xyzs[3*(v[i] - 1)  +2  ] 		);
			snprintf(buffer,33, "AVCOORD_%d",AV);
			++AV;
			print_mesh(qm, buffer);
			if ( stat != EGADS_SUCCESS) return stat;
			printf(" CALL AROUND VERTEX %d \n", v[i]);
			stat = quad_algebraic_area(qm, qID, v[i], &area, angles, &quadArea);
			printf(" AREA %d  stat %d \n",area,stat);
			if ( area == 1 ) break;
			else {
				qm -> uvs [2*(v[i] - 1)    ] = uvxyz[0];
				qm -> uvs [2*(v[i] - 1) + 1] = uvxyz[1];
				qm -> xyzs[3*(v[i] - 1)    ] = uvxyz[2];
				qm -> xyzs[3*(v[i] - 1) + 1] = uvxyz[3];
				qm -> xyzs[3*(v[i] - 1) + 2] = uvxyz[4];
			}
		}
	}
	PRINTAREA = 1;
	stat = quad_algebraic_area(qm, qID, v[0], &area, angles, &quadArea);
	PRINTAREA = 0 ;
	snprintf(buffer,33, "AVCOORD_%d",AV);
	++AV;
	print_mesh(qm, buffer);
	printf(" CHECK INVALID AREA FINALL RESULT %d \n", area);
	if ( area != 1 ) return EGADS_GEOMERR;
	return EGADS_SUCCESS;
}

static int
EG_projectToTangentPlane(double normal[], double O[], double p[], double *proj) {
	double c, dotNN = 0.0, dotNP = 0.0, dist, lambda;
	dist    = (p[0]- O[0]) * (p[0]- O[0]) + (p[1]- O[1])*(p[1]- O[1]) + (p[2]- O[2])*(p[2]- O[2]);
	dist    = sqrt(dist);
	if (dist < DEPS ) {
		proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
		return EGADS_SUCCESS;
	}
	c       = normal[0] *      O[0] + normal[1] *      O[1] + normal[2] *      O[2]; // Equation plane: a*x + b*y + c*z = C
	dotNN   = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
	dotNP   = normal[0] *      p[0] + normal[1] *      p[1] + normal[2] *      p[2];
	if ( PRINTAREA) {
		printf(" ORIGIN %lf %lf  %lf  NORMAL %lf %lf  %lf  TARGET %lf %lf %lf\n",
				O[0], O[1], O[2], normal[0], normal[1], normal[2], p[0], p[1], p[2]);
		printf("TOL %lf <= %lf ]]]] DIST %lf  c  %lf  DOT NN %lf  DOT NP %lf\n", dotNP -c, DEPS, dist, c, dotNN, dotNP);
	}
	if ( fabs(dotNP - c) <= 1.e-07) {
		proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
		return EGADS_SUCCESS;
	}
	lambda  = (c - dotNP)/dotNN;
	if ( PRINTAREA) printf(" DOT NN %lf  c- np %lf  lambda %lf\n", dotNN, c - dotNP, lambda);
	proj[0] = p[0] + lambda*normal[0];
	proj[1] = p[1] + lambda*normal[1];
	proj[2] = p[2] + lambda*normal[2];
	// check that point belongs to plane.
	dist  = normal[0] * proj[0] + normal[1] * proj[1] +  normal[2] * proj[2];
	if( fabs(dist - c) < DEPS) {
		return EGADS_SUCCESS;
	}
	else{
		printf(" ORIGIN %lf %lf  %lf  NORMAL %lf %lf  %lf  TARGET %lf %lf %lf\n", O[0], O[1], O[2], normal[0], normal[1], normal[2], p[0], p[1], p[2]);
		printf(" POINT SHOULD BELONG TO PLANE!!!!! %lf ~= %lf\n",dist,c);
		printf(" DOT NN %lf PN %lf LAMBDA %lf  POINT %lf %lf %lf\n", dotNN, dotNP, lambda, proj[0], proj[1], proj[2]);
		return EGADS_GEOMERR;
	}
}

static void
cross_product(double *A, double *B, double *cross) {
	cross[0] = A[1] * B[2] - A[2] * B[1];
	cross[1] = A[2] * B[0] - A[0] * B[2];
	cross[2] = A[0] * B[1] - A[1] * B[0];
}


static double
triang_cross_area(double A[], double B[], double C[], double *cross){
	int i;
	double AB[3], AC[3], norm;
	for ( i = 0 ; i < 3; ++i){
		AB[i] = B[i] - A[i];
		AC[i] = C[i] - A[i];
	}
	cross_product(AB, AC, cross);
	norm = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
	norm = 0.5*sqrt(norm);
	return (norm);
}





static int
computePenalty(quadMap *qm,  int starID, int qID, double *penalty, double *angles, double *quadArea) {
	double dot, dotK[2], z, zerr, alpha, dev[3], nAux,   minmax[2], scale, curv[8], projV[3*2],  proj[3], norms[2], normal[3], thetaK, kSUM, hTheta[2], fact[2], areaParallelogram, areaOpti;
	int i, j, id0, stat, validArea, q, v, link,  sID, movingVert, invalid = 0, idx, aux;
	q    = qm -> star[starID] -> quads[qID] - 1;
	sID  = qm -> star[starID] -> verts[0  ] - 1;
	stat = quad_algebraic_area(qm, q + 1, sID + 1, &validArea, angles, quadArea);
	for ( i = 0 ; i < 4; ++i) {
		if (qm -> quadIdx[4*q + i] == sID + 1 ) id0 = i;
	}
	if ( stat != EGADS_SUCCESS) {
		printf(" ERROR COMPUTING PENALTY FUNCTION  %d",stat);
		return stat;
	}
	// USE CURVATURES TO MEASURE MOVEMENT + ELEMENT SIZE
	stat = EG_curvature ( qm -> face, &qm -> uvs[2*sID], curv);
	cross_product(&curv[1], &curv[5], normal);
	for ( j = 0 ; j < 2; ++j) {
		aux  = 2*(qID + j) + 1;
		if ( aux >= qm -> star[starID] -> nV) aux -= (qm -> star[starID] -> nV - 1);
		link = qm ->star[starID] -> verts[aux] - 1;
		stat = EG_projectToTangentPlane(normal, &qm -> xyzs[3*sID], &qm -> xyzs[3*link], proj );
		for ( i = 0 ; i < 3; ++i) {
			proj [i       ] -= qm ->xyzs[3*sID + i ];
			projV[3*j + i ]  = proj[i];
		}
		EG_unitVector (proj    , &norms[j]);
		fact[j] = 0.0;
		if (j == 0 &&  norms[j] > HMAX ) {
			fact[j] = norms[j] / HMAX;
			invalid = 1;
		}
		else if (j == 0 && norms[j] < HMIN ) {
			fact[j] = HMIN / norms[j];
			invalid = 1;
		}
	}
	cross_product(&projV[0], &projV[3], proj );
	alpha = angles[0];
	idx   = sID;
	for ( i = 1; i < 4; ++i) {
		v = qm ->quadIdx[4*q + (i + id0)%4] - 1;
		if ( (i != 2) && qm -> vType[v] == -1) {
			movingVert = 1;
			for ( j = 0 ; j < qm -> nS; ++j) {
				if(v + 1 == qm -> star[j]->verts[0]) movingVert = 0;
			}
			if( (movingVert == 1) && (angles[i] >= M_PI*0.85) && (angles[i] > alpha)) {
				alpha = angles[i];
				idx   = v;
				break;
			}
			else if( (movingVert == 1) && ( fabs( (2.0 * M_PI)/(double)qm ->valence[v][0] - angles[i] ) > fabs( 2.0*M_PI/(double)qm->valence[idx][0] - alpha ) ) ) {
				alpha = angles[i];
				idx   = v;
			}
		}
	}
	dot      = cos(alpha - 0.5*M_PI);
	if ( alpha > 2.0*M_PI - acos(ANGPASS)) dot = -1.0;
	scale    = alpha;
	if ( dot   > 0 && dot < 0.5*ANGPASS) scale = 2.0*M_PI-alpha;
	z        = ERRCTT*((dot-ANGPASS*0.5)/ANGPASS);
	zerr     = erfc(z);
	*penalty = zerr * exp(zerr *scale/M_PI);
	if ( invalid == 1) {
		scale = zerr + exp( fact[0] + fact[1]);
		if ( scale > *penalty) *penalty = scale;
	}
	//#ifdef DEBUG
	if (OPTCOUNT%NOUT == 0 ) {
		printf("\n :::::::::: PENALTY FUNCTION AT VERTEX %d in QUAD %d : AREA = %lf  ANGLE AT VERT %lf-> VALIDAREA %d TOO SMALL TOO BIG = %d: AREA LENGTH %lf %lf \n",
				qm ->quadIdx[4*(q) + (i + id0)%4], qID, *quadArea, angles[0], validArea, invalid, norms[0],  norms[1]);
		printf("======= PENALTY : alpha = %lf  dot = %lf  scale = %lf   z  %lf  zerr = %lf  penalty %lf\n ",
				alpha, dot, scale, z, zerr, *penalty);
	}
	//#endif
	angles[1] = alpha;
	return EGADS_SUCCESS;

}



static int
quad_algebraic_area(quadMap *qm, int qID, int vID,  int *validArea, double *quadAngles, double *quadArea) {
	int i, id, v, vPos, qV[4], k, tri, piv[4], stat, tooClose, degenerated[4], ftype;
	double vTri[18*3], evalVert[18], cross[4], dotCross[2*4], area[2*4], norm[2], theta[2*4], vNormal[4], uv_dir[2], uv_eps[2],  vPro[8], vEps[8];
	double naux,  z, zerr, alpha, dot;
	vPos = -1;
	//printf(" QUAD ALG AREA %d  ID %d\n", qID, vID);
	for ( i = 0 ; i < 4; ++i) {
		qV[i] = qm ->quadIdx[4*(qID - 1) + i] - 1;
		if (qV[i] == (vID -1) ) vPos = i;
	}
	if ( vPos == -1 ) {
		printf(" VERTEX %d is STAR CENTRE AND SHOULD BELONG TO QUAD %d\n",vID, qID);
		return EGADS_INDEXERR;
	}
	for ( k = 0 ; k < 4; ++k) { // Compute area starting from all vertices to capture invalid shapes
		for ( tri = 0 ; tri < 2; ++ tri) { // quad area as a sum of 2 triangles
			piv[0]  = qV[(k  + vPos          )%4];
			piv[1]  = qV[(k  + vPos + 1 + tri)%4];
			piv[2]  = qV[(k  + vPos + 2 + tri)%4];
			if ( PRINTAREA ) printf(" AREA AROUND %d %d %d\n", piv[0] + 1, piv[1] + 1, piv[2] + 1);
			vTri[0] = qm ->xyzs[3*piv[0]]; vTri[1] = qm ->xyzs[3*piv[0]  + 1]; vTri[2] = qm ->xyzs[3*piv[0] + 2];
			stat = EG_evaluate(qm->face, &qm -> uvs[2*piv[0]], evalVert);
			if ( stat != EGADS_SUCCESS) return stat;
			cross_product(&evalVert[3], &evalVert[6], vNormal);
			if ( qm -> face -> mtype == SREVERSE )
				cross_product(&evalVert[6], &evalVert[3], vNormal);
			else
				cross_product(&evalVert[3], &evalVert[6], vNormal);
			EG_unitVector(vNormal, &vNormal[3]);
			for ( id = 1 ; id <= 2; ++id) { // move away from vertex towards B and C respectively: we are computing the angle very close to the centre
				uv_dir[0]      = qm -> uvs[2*piv[id]    ] - qm -> uvs[2*piv[0]    ];
				uv_dir[1]      = qm -> uvs[2*piv[id] + 1] - qm -> uvs[2*piv[0] + 1];
				naux           = uv_dir[0] * uv_dir[0] + uv_dir[1] * uv_dir[1];
				naux           = sqrt(naux);
				degenerated[k] = 0;
				if ( naux < DEPS ) { // vector OB or OC is too small => O ~B  or O ~C INVALID
					degenerated[  k    ] =    1;
					theta      [2*k    ] =  0.0;
					theta      [2*k + 1] =  0.0;
					area       [2*k    ] =  0.0;
					area       [2*k + 1] =  0.0;
					dotCross   [2*k    ] = -1.0;
					dotCross   [2*k + 1] = -1.0;
					id                   = 3;
					tri                  = 2;
				}
				if(!degenerated[k]) {
					if (PRINTAREA)
						printf(" NORM %.16e  \t", naux);
					uv_dir[0] /= naux;                uv_dir[1] /= naux;
					if (PRINTAREA)
						printf(" DIR %lf    %lf  \n",uv_dir[0], uv_dir[1]);
					uv_eps[0]  = qm -> uvs[2*piv[0]]; uv_eps[1]  = qm -> uvs[2*piv[0] + 1];
					do {
						uv_eps[0] += EPSDIR * uv_dir[0];
						uv_eps[1] += EPSDIR * uv_dir[1];
						stat       = EG_evaluate(qm -> face, uv_eps, evalVert);
						if( stat  != EGADS_SUCCESS ) return stat;
						if ( PRINTAREA)
							printf("  POINT %lf  %lf  %lf\n ", evalVert[0], evalVert[1], evalVert[2]);
						stat       = EG_projectToTangentPlane(vNormal,vTri, evalVert, &vTri[18*id]);
						if ( PRINTAREA) {
							printf(" PROJECTED POINT %lf  %lf  %lf\n ", vTri[18*id],  vTri[18*id + 1], vTri[18*id + 2]);
							PRINTAREA = 0 ;
						}
						if( stat  != EGADS_SUCCESS ) return stat;
						naux       = 0.0;
						for ( i = 0 ; i < 3; ++i)
							naux += (vTri[18*id + i] - vTri[i]) * (vTri[18*id + i] - vTri[i]);
						naux       = sqrt(naux);

					} while(naux < EPS);
				}
				PRINTAREA = 1;
				if ( PRINTAREA)
					printf(" PROJECTED POINT %lf  %lf  %lf\n ", vTri[18*id],  vTri[18*id + 1], vTri[18*id + 2]);
			}
			if(!degenerated[k]) {
				// get vectors vA vB from vertex v
				for ( i = 0 ; i < 3; ++i) {
					vEps[i]     = vTri[18*1 + i] - vTri[i];
					vEps[4 + i] = vTri[18*2 + i] - vTri[i];
				}
				EG_unitVector (&vEps[0], &vEps[3] );
				EG_unitVector (&vEps[4], &vEps[7] );
				if (PRINTAREA)
					printf(" UNIT VECTOR %.16e    %.16e \n", vEps[3], vEps[7]);
				if ( vEps[3] < 1.E-7 || vEps[7] < 1.E-7) theta[2*k + tri] = 0.0;
				else {
					dot = 0.0;
					for     ( i = 0 ; i < 3; ++i)      dot             += vEps[i] * vEps[i + 4];
					if ( PRINTAREA)
						printf(" DOT %lf   %lf  %lf  %lf   %lf  %lf   %lf\n",dot, vEps[0], vEps[1], vEps[2],
								vEps[4], vEps[5], vEps[6] );

					if      ( fabs(dot - 1.0) < DTOL ) theta[2*k + tri] = 0.0;
					else if ( fabs(dot + 1.0) < DTOL ) theta[2*k + tri] = M_PI;
					else                               theta[2*k + tri] = acos(dot);
				}
				// Compute the area using chords
				for ( i = 0 ; i < 3; ++i){
					vTri[       i] = qm -> xyzs[3*piv[0] + i ];
					vTri[  18 + i] = qm -> xyzs[3*piv[1] + i ];
					vTri[2*18 + i] = qm -> xyzs[3*piv[2] + i ];
				}
				area[2*k + tri]     = triang_cross_area(&vTri[0], &vTri[18], &vTri[36], cross);
				EG_unitVector(cross, &cross[3]);
				dotCross[2*k + tri] = cross[0] * vNormal[0] + cross[1] * vNormal[1] + cross[2] * vNormal[2];
#ifdef DEBUGG
				if (PRINTAREA)
					printf(" *********  DOT CROSS PRODUCT %lf  -> THETA = %lf \t ", dotCross[2*k + tri], theta[2*k + tri]);
#endif
				if ( dotCross[2*k + tri] < - 0.8 ) theta[2*k + tri] = 2.0 * M_PI - theta[2*k + tri];  // we have anti-clockwise orientation. Postive area implies reversed.
#ifdef DEBUGG
				if (PRINTAREA)
					printf("  -> THETA = %lf \t \n", theta[2*k + tri]);
#endif
			}
		}
#ifdef DEBUGG
		if (PRINTAREA)
			printf(" *********  AREA OF QUAD   %lf \n", area[2*k] + area[2*k + 1]);
#endif
	}
	*validArea = 1;
	if ( fabs( (area[0] + area[1]) - (area[2] + area[3]) ) > EPSAREA ) {
#ifdef DEBUGG
		if (PRINTAREA)
			printf(" AREA %lf != %lf\n", area[0] + area[1], area[2] + area[3] );
#endif
		*validArea = -1;
	}
	*quadArea = area[0] + area[1];
	for ( k = 0 ; k < 4; ++k) {
		quadAngles[k] = theta[2*k] + theta[2*k+1];
		//if ( quadAngles[k] >  2.0*M_PI ) quadAngles[k] -= 2.0*M_PI;
		if ( (qm -> vType[qV[(k  + vPos )%4]] == -1) && (quadAngles[k] >= ANGLE_ERR))  *validArea = -1;
		if ( (qm -> vType[qV[(k  + vPos )%4]] == -1) && (degenerated[k])            ) *validArea = -1;
#ifdef DEBUGG
		if (PRINTAREA)
			printf(" ANGLES AROUND %d =  %lf   AND AREAS %lf  ANGLE WITH SURFACE NORMAL  %lf  %lf\n",
					qm ->quadIdx[4*(qID - 1) + (vPos + k)%4] , quadAngles[k], area[2*k] + area[2*k + 1], dotCross[2*k], dotCross[2*k +1]);
#endif
	}
	if (PRINTAREA)
		printf(" LEAVING WITH VALID AREA = %d  SIZE  %lf \n",*validArea,*quadArea);
	return EGADS_SUCCESS;
}

static void
EG_vertexRegularity(quadMap *qm, int list[], int *vIrr, int n)
{
	int i;
	for ( i = 0; i < n; ++i) {
		if (qm -> vType[ list[2*i] -1] == 0) vIrr[i] = 100; //edge vertices have "SPECIAL" regularity//{
		//if (list[2*i + 1] == 2) vIrr[i] = 0 ;
		//else
		//}
		else if (qm -> vType[ list[2*i] -1] > 0) {
			if      ( list[2*i+1] == 2 ) vIrr[i] = -1;
			else if ( list[2*i+1] == 3 ) vIrr[i] =  0;
			else                         vIrr[i] = list[2*i + 1] + 1;
		}
		else {
			if      ( list[2*i+1]  < 4 ) vIrr[i] = -1;
			else if ( list[2*i+1] == 4 ) vIrr[i] =  0;
			else                         vIrr[i] = list[2*i + 1];
		}
	}
}


static int
EG_buildStar(quadMap *qm, vStar **star, int vID ){
	int stat, i, j, id0 = -1, k = 0, itQ, vROSE[100], vQUADS[100];
	int auxm, auxV, auxQ,  r, quad0, quadID, prevQuad, nPeaks, adjQ;
	static int qLoop[5] = {0, 1, 2, 3, 0};
	if ( vID <= 0 || vID > qm -> totVerts) {
		printf(" TRYING TO BUILD STAR ON  VERTEX out of bounds %d  MAX %d!!!\n", vID, qm -> totVerts);
		return EGADS_INDEXERR;
	}
	if (qm -> vType[vID-1] == -2) {
		printf(" We are trying to build on a deleted vertex!!!!!\n");
		return EGADS_INDEXERR;
	}
	quad0         = qm -> valence[vID - 1][1] - 1;
	printf(" STAR CENTRE %d START AT ID %d \n", vID, quad0 + 1);
	print_quad_specs(qm, quad0 + 1);
	quadID        = quad0;
	prevQuad      = -5;
	r             =  0;
	itQ           =  0;
	vROSE[r++]    = vID;
	vQUADS[itQ++] = quad0 + 1; // Accumulates the quads that we are going to touch
	for ( j = 0 ; j < 4; ++j) {
		if ((qm->quadIdx[4*quad0 + j ] == vID )) {
			id0 = j;
			break;
		}
	}
	if ( id0 == -1) {
		printf("I couldn't find the star centre (Vertex %d) in Quad %d\n. quadMap build UNSUCCESSFUL. Check updateUVs function!\n", vROSE[0], quad0 + 1);
		print_quad_specs(qm, quad0 + 1);
		return (EGADS_INDEXERR);
	}
	vROSE[r] = -1;
	do {
		for ( j = 0 ; j < 3; ++j) {
			auxm = (id0 + j)%4;
			if ( qm -> quadAdj[4*quadID + auxm] != (prevQuad +1) && qm -> quadAdj[4*quadID + auxm] != (quad0+1)) {
				if ( vROSE[1] != qm -> quadIdx[4*quadID + qLoop[auxm + 1] ]) { // At last quad, this link will be repeated; vROSE[2] is the first link
					vROSE[r++] = qm -> quadIdx[4*quadID + qLoop[auxm + 1] ];
				}
			}
		}
		id0      = (int)(id0 + 3)%4;
		prevQuad = quadID;
		quadID   = qm -> quadAdj[4*prevQuad + id0] - 1;
		print_quad_specs(qm, quadID + 1);
		if ( quadID < 0) { // find the other boundary vertex. We will loop around the valences
			for ( k = 0; k < qm ->valence[vID - 1][0]; ++k) {
				auxV = qm -> valence[vID -1][k + 2] - 1;
				auxQ = qm -> valence[auxV]  [1    ] - 1;
				for ( i = 0 ; i < r; ++i) {
					if ( (auxV + 1) == vROSE[i] ) {
						auxQ = -2; i = r;
					}
				}
				if ( auxQ >= 0) {
					// Find the quad common to vID andauxV
					auxm   = -1;
					for ( i = 0 ; i < 4; ++i) {
						if ( qm -> quadIdx[4*auxQ + i ] == vID) auxm = auxQ;
					}
					if ( auxm == -1) { // Wrong quad. Look around its adjacents to find the one that contains vID
						for ( i = 0 ; i < 4; ++i) {
							adjQ = qm -> quadAdj[4*auxQ + i ] - 1;
							if ( adjQ >= 0) {
								for ( j = 0 ; j < 4; ++j) {
									if ( qm -> quadIdx[4*adjQ + j ] == vID) {
										auxm = adjQ;
									}
								}
							}
						}
					}
					if ( auxm != -1 ) {
						for ( i = 0 ; i < 4; ++i) {
							if ( qm ->quadAdj[4*auxm + i] == -1 ) goto getId0;
						}
					}
				}
			}
			auxm = quad0;
		}
		getId0:
		if ( quadID < 0) {
			vROSE[r++]    = -1;
			vQUADS[itQ++] = -1;
			quadID        = auxm;
		}
		id0 = -1;
		for ( j = 0 ; j < 4; ++j) {
			if ( ( qm->quadIdx[4*quadID + j ] == vID  ) ) {//&&  (qm ->quadAdj[4*quadID + j] == (prevQuad + 1) ) ) {
				id0  = j;
				break;
			}
		}
		if ( id0 == -1) {
			printf(" LOOKING FOR VID %d IN QUAD %d \n", vID, quadID + 1);
			print_quad_specs(qm, quadID+ 1);

			printf(" WE MESSED UP. BUG ON FINDING ADJ QUAD!!\n");
			for ( j = 0 ; j < qm ->nQ; ++j) print_quad_specs(qm, j+1);
			exit(1);

		}
		if (quadID != quad0 ) vQUADS[itQ++] = quadID+1;

	}
	while(quadID != quad0);
	nPeaks           = r;
	*star            = (vStar*)EG_alloc (    sizeof(vStar));
	if ((*star) == NULL ) return EGADS_MALLOC;
	(*star) -> verts = (int*)  EG_alloc (r  *sizeof(int  ));
	(*star) -> nV    = nPeaks;
	(*star) -> quads = (int*)  EG_alloc (itQ*sizeof(int  ));
	(*star) -> nQ    = itQ;
	if ( (*star) -> verts == NULL || (*star) -> quads == NULL) {
		EG_free((*star));
		return EGADS_MALLOC;
	}
	for ( r = 0 ; r < itQ; ++r)      (*star) -> quads[r] = vQUADS[r];
	for ( r = 0 ; r < nPeaks; ++r) {
		(*star) -> verts[r] = vROSE[r];
	}
	//print_star((*star));
	return EGADS_SUCCESS;
}



static int
EG_vertexSplit(quadMap *qm, int sID) {
	int poly[2*3], q, aux, id0, i, j, val, nV, nQ;
	int id, idAux, stat, dist, validSplit, splitOK = 0, activity;
	int  *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, *locVtype = NULL, **locVal = NULL, *locQrem = NULL, *locVrem = NULL;
	double xyz[18], *uv0 = NULL;
	vStar *star;
	if ( qm -> vType  [sID - 1] == -2  ) return EGADS_SUCCESS;
	if ( qm -> valence[ sID -1][0] < 5 ) return EGADS_SUCCESS;
	stat = EG_buildStar(qm, &star, sID);
	if ( stat != EGADS_SUCCESS) return stat;
	stat = backupQm ( qm, &nV, &nQ, &locIdx, &locAdj, &locVtype, &locVal, &locQrem, &locVrem, &uv0);
	if ( stat != EGADS_SUCCESS) {
		EG_free(star -> verts);
		EG_free(star -> quads);
		EG_free(star);
		return stat;
	}
	for ( i = 0 ; i < 6; ++i) poly[i] = 100;
	poly[2*0    ] = star -> verts  [0];
	poly[2*0 + 1] = qm   -> valence[star -> verts[0] -1][0];
	for ( q = 0 ; q < star -> nQ; ++q) {
		id0           = 2*q + 1;
		poly[2*1    ] = star -> verts  [id0];
		if ( qm -> vType[poly[2*1] -1] != -1) continue;
		poly[2*1 + 1] = qm   -> valence[poly[2*1] - 1][0];
		for ( i = 0 ; i <= poly[1] - 4; ++i) {
			aux  = id0 + 4  + 2 * i;
			dist = 4 + 2*i;
			if ( aux >= star -> nV) aux -= (star -> nV - 1) ;
			poly[2*2    ] = star -> verts  [aux];
			if ( qm -> vType[poly[2*2] -1] != -1) continue;
			poly[2*2 + 1] = qm   -> valence[poly[2*2] - 1][0];
		}
		splitOK = 0 ;
		if ( (qm -> vType[poly[2*2] - 1] == -1) && (qm -> vType[ poly[2*1] - 1] == -1) ) {
			validSplit = 0 ;
			if ( (poly[2*2 + 1] == 3) && (poly[2*1 + 1] == 3) ) validSplit = 1;
			else if ( (poly[1] >= 6) && (poly[2*2 + 1] <= 4) && (poly[2*1 + 1] <= 4) )
				if ( poly[2*2 + 1] == 3 ||  poly[2*1 + 1] == 3) validSplit = 1;
			if ( validSplit )
			{
#ifdef DEBUG
				printf("WE FOUND A CANDIDATE FOR SPLITTING: SPLIT VERTEX V %d = %d FROM LINKS %d =  %d  and %d = %d\n", poly[0], poly[1], poly[2], poly[3], poly[4], poly[5]);
#endif
				snprintf(buffer,33,"VERTEXSPLIT_%d",VS);
				print_mesh(qm, buffer);
				stat = EG_splittingOperation(qm, star, poly, id0, dist, &activity);
				print_mesh(qm, buffer);
				if ( stat == EGADS_SUCCESS && activity) splitOK = 1;
				else {
#ifdef DEBUG
					printf(" SPLITTING WENT %d \n",stat);
#endif
					EG_free(star -> verts);
					EG_free(star -> quads);
					EG_free(star);
					stat = resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal, locQrem, locVrem, uv0);
					for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
					EG_free(locVal);
					EG_free(locIdx);
					EG_free(locAdj);
					EG_free(locVtype);
					EG_free(uv0);
					return stat;
				}
			}
		}
		if ( splitOK ) break;
	}
	EG_free(star -> verts);
	EG_free(star -> quads);
	EG_free(star);
	for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
	EG_free(locVal);
	EG_free(locIdx);
	EG_free(locAdj);
	EG_free(locVtype);
	EG_free(locQrem);
	EG_free(locVrem);
	EG_free(uv0);
	return EGADS_SUCCESS;
}


static int
EG_splittingOperation(quadMap *qm, vStar *star, int poly[], int id0, int distanceToId0, int *activity) {
	int qID[2], qAux[2], edges[4],  movingVerts[2000],  valAux[100], nMove;
	int q, newV, newQ, val, id, i, j, k, aux, stat, idAux, vj[2],vPull[2], v0, area ;
	double xyz[18], dir[2*4], diropp[2], angles[4], quadArea, uvCentre[2];
	vStar *checkStar;
	*activity = 0;
	if ( qm -> remQuads[0] > 0 ) {
		newV = qm -> remVerts[qm -> remVerts[0]];
		newQ = qm -> remQuads[qm -> remQuads[0]];
		--qm -> remVerts[0]; --qm -> remQuads[0];
		qm -> vType  [newV - 1] = -1;
	}
	else {
		newV                     = qm -> totVerts + 1;
		newQ                     = qm -> nQ       + 1;
		qm   -> vType [newV - 1] = -1;
		++qm -> totVerts;
		++qm -> nQ;
		if ( qm -> totVerts > 2* qm -> v0) {
			printf(" We have duplicated the number of initial vertices. This is too much. \n");
			return EGADS_INDEXERR;
		}
	}
	printf(" I'm going to split vertex %d between %d and %d and insert VERTEX %d AND QUAD %d\n",poly[0], poly[2], poly[4], newV, newQ);
	*activity = 1;
	// Split  Valences from original vertex
	val = (int)(distanceToId0/2) + 1;   // we compute star vertices by v = 2*q + 1, q =0,..,nQ
	id  =  poly[0] - 1;
	qm -> valence[id][0] = val;
	qm -> valence[id][1] = newQ;
	for ( i = 0; i < val; ++i) {
		aux = id0 + 2 * i; if ( aux >= star -> nV) aux -= (star -> nV - 1) ;
		qm -> valence[id][2 + i] = star -> verts[aux];
	}
	print_star(star);
	// Create Quad, point and reallocate adjacents
	qID[0] = star -> quads[(int)(id0                 - 1)/2    ];
	val    =               (int)(id0 + distanceToId0 - 1)/2 - 1;
	qID[1] = star -> quads[val % star ->nQ];
	printf(" ID0 %d  -> VERTEX %d at  QUAD %d   dist %d  Quad 1  %d\n", id0, star ->verts[id0], qID[0], distanceToId0, qID[1]);
	qm -> quadIdx[4*(newQ - 1)    ] = poly[2];
	qm -> quadIdx[4*(newQ - 1) + 1] = poly[0];
	qm -> quadIdx[4*(newQ - 1) + 2] = poly[4];
	qm -> quadIdx[4*(newQ - 1) + 3] = newV;
	//
	qm -> quadAdj[4*(newQ - 1)    ] = qID[0];
	qm -> quadAdj[4*(newQ - 1) + 1] = qID[1];
	for ( q = 0; q < 2; ++q ){
		stat  = EG_getAdjacentToVertexPair(qm, qID[q], poly[0], poly[2* (q + 1)], qAux);
		stat += EG_getCommonEdge (qm, qID[q], qAux[1], edges);
		if ( qAux[0] != edges[0] || stat != EGADS_SUCCESS) {
			printf(" I Messed up finding common faces. We have aux pointing at %d and edge pointing at %d\n",   qm -> quadAdj[4*(qID[q] - 1) + qAux[0]], edges[0]);

			return EGADS_INDEXERR;
		}
		qm -> quadAdj[4*(qID [q] - 1) + edges[0] ] = newQ;
		qm -> quadAdj[4*(qAux[1] - 1) + edges[1] ] = newQ;
		qm -> quadAdj[4*(newQ    - 1) + 3 - q    ] = qAux[1];
		for ( i = 0 ; i < 4; ++i)
			if (qm -> quadIdx[4*(qAux[1] - 1) + i] == poly[0] ) qm -> quadIdx[4*(qAux[1] - 1) + i] = newV;
	}
	q = ((int)(id0 + distanceToId0 - 1)/2)% star ->nQ;
	while ( q !=  (int)(id0 - 1)/2 ) {
		for ( i = 0 ; i < 4; ++i)
			if (qm -> quadIdx[4*(star->quads[q] - 1) + i] == poly[0] ) qm -> quadIdx[4*(star->quads[q] - 1) + i] = newV;
		q++; if ( q == star ->nQ) q = 0 ;
	}

	// Add valences to splitting vertices
	for ( j = 0 ; j < 2; ++j) {
		id  = poly[2*(j + 1)    ] - 1;
		val = poly[2*(j + 1) + 1];
		for ( i = 0 ; i < val; ++i) valAux[i] = qm -> valence[id][2 + i];
		++val;
		qm -> valence[id][0] = val;
		qm -> valence[id][1] = newQ;
		for ( i = 0; i < val - 1; ++i)
			qm -> valence[id][2 + i] = valAux[i];
		qm -> valence[id][2 + val - 1] = newV;
	}
	// Create valence pointer for new vertex and update links
	val = poly[1] - qm -> valence[poly[0] - 1][0] + 2;
	id  = newV    - 1;
	qm -> valence[id][0] = val;
	qm -> valence[id][1] = newQ;
	for ( i = 0; i < val; ++i) {
		aux = id0 + distanceToId0 + 2*i;
		if ( aux >= star -> nV) aux -= (star -> nV - 1);
		idAux                    = star -> verts[aux] - 1;
		qm -> valence[id][2 + i] = idAux + 1;  // add valence to new vertex
		if ( (idAux + 1 == poly[2]) || (idAux + 1 == poly[4]) )
			qm -> valence[id][2 + i] = idAux + 1;  // add valence to new vertex
		else {
			for ( j = 0 ; j < qm ->valence[idAux][0]; ++j) {
				if ( qm -> valence[idAux][2 + j] == poly[0])  {  // update link to new vertex
					qm -> valence[idAux][2 + j] = newV;
					j = qm ->valence[idAux][0];
				}
			}
		}
	}
#ifdef DEBUG
	printf(" VALENCES FOR OLD VERTEX = %d \n", poly[0]);
	for ( i = 0; i < qm -> valence[poly[0] - 1][0]; ++i) {
		idAux = qm -> valence[poly[0] - 1][2 + i];
		printf(" V %d = %d\n",i, idAux);
		for ( j = 0 ; j < qm -> valence[idAux-1][0]; ++j)
			printf(" points back %d = %d\n",j, qm ->valence[idAux -1][2+j]);
	}
	printf(" VALENCES FOR NEW VERTEX  = %d\n", id + 1);
	for ( i = 0; i < qm -> valence[id][0]; ++i) {
		idAux = qm -> valence[id][2 + i];
		printf(" V %d = %d\n",i, qm ->valence[id][2+i]);
		for ( j = 0 ; j < qm -> valence[idAux -1][0]; ++j)
			printf(" points back %d = %d\n",j, qm ->valence[idAux -1][2+j]);
	}
#endif
	stat  = EG_averageCoords (qm, poly[0]);
	stat += EG_averageCoords (qm, newV);
	if ( stat != EGADS_SUCCESS) return stat;
	snprintf(buffer,33,"INSERT_%d_%d",newV,ADDP);
	printf(" CALLING MESH %s \n", buffer);
	print_mesh(qm, buffer);
	++ADDP;
	vj[0]                = poly[0];
	vj[1]                = newV   ;
	nMove                = 0;
	movingVerts[nMove++] = vj[0];
	movingVerts[nMove++] = vj[1];
	for ( j = 0 ; j < 2; ++j) {
		stat = EG_buildStar(qm, &checkStar, vj[j]);
		if ( stat != EGADS_SUCCESS) return stat;
		for ( i = 0 ; i < checkStar -> nQ; ++i){
			stat = checkInvalidElement(qm, checkStar -> quads[i] );
			if ( stat != EGADS_SUCCESS) {
				EG_free(checkStar -> verts);
				EG_free(checkStar -> quads);
				EG_free(checkStar);
				return stat;
			}
		}
		for ( i = 0 ; i < qm -> valence[vj[j] - 1][0]; ++i) {
			aux = 1;
			if( qm -> vType[qm -> valence[vj[j] - 1][2 + i] -1] != -1 ) continue;
			for ( k = 0 ; k < nMove; ++k) {
				if ( qm -> valence[vj[j] - 1][2 + i] == movingVerts[k] ) {
					aux = 0;
					break;
				}
			}
			if ( aux )
				movingVerts[nMove++] = qm -> valence[vj[j] - 1][2 + i];
		}
	}
	EG_free(checkStar -> verts);
	EG_free(checkStar -> quads);
	EG_free(checkStar);
	for ( i = 0 ; i < nMove; ++i) printf(" MOVING VERT %d = %d TYPE %d\n",i, movingVerts[i], qm -> vType[movingVerts[i] -1]);
	printf(" IN SPLITTING OPERATION . WE HAVE INTRODUCED VERTEX %d. CALL OPTI\n", newV);
	snprintf(buffer,33,"SBEFOREOPTI_%d_%d",newV, ADDP);
	printf(" CALLING MESH %s\n",buffer);
	print_mesh(qm, buffer);
	stat = optimize_angles(qm, nMove, movingVerts);
	snprintf(buffer,33,"SAFTEROPTI_%d_%d",newV, ADDP);
	printf(" CALLING MESH %s\n",buffer);
	print_mesh(qm, buffer);
	++ADDP;
	return stat;
}

static int
EG_mergeVertices(quadMap *qm, int qC, int poly[], int collapseToCentre, int *activity) {
	int vC, vC2,vID, stat, i, e, q, v, vCommon, j, adjQ, k, updatedQuad, vVAL, colVal, aux, auxVal[100], addValence;
	int edges[4], auxAdj, auxAdj2, nV, nQ;
	double uv[2], vPos[18];
	vStar *star;
	--qC;
	*activity = 0;
	if ( qC < 0 ) {
		printf(" WE ARE CALLING TO MERGE VERTICES WITH A VOID QUAD!! BUGGGGG\n");
		return EGADS_INDEXERR;
	}
	for ( i = 0 ; i < 4; ++i)
		printf(" POLY %d = %d\n", poly[2*i], poly[2*i + 1]);
#ifdef DEBUG
	printf(" EG MERGE VERTICES:: COLLAPSE v(%d) = %d   TO v(%d) = %d AFFECTED VS v(%d) = %d v(%d) = %d (each will be -1) QUAD REMOVE %d\n ",
			poly[0], poly[1], poly[2*2], poly[2*2 + 1],
			poly[2*1], poly[2*1 + 1],poly[2*3], poly[2*3 + 1], qC + 1);
#endif
	// Get Quad Centre to place collapsed vertex
	vC  = poly[2*2] - 1;
	vC2 = poly[0]   - 1;
	stat = EG_buildStar(qm, &star, vC2 + 1);
	if ( stat != EGADS_SUCCESS) {
		return stat;
	}
	*activity = 1;
	if ( collapseToCentre) {
		uv[0] = 0.0; uv[1] = 0.0;
		for ( i = 0 ; i < 4; ++i) {
			uv[0] += qm->uvs[2*(poly[2*i] - 1)    ];
			uv[1] += qm->uvs[2*(poly[2*i] - 1) + 1];
		}
		uv[0] *= 0.25;
		uv[1] *= 0.25;
		stat   = EG_evaluate(qm->face, uv, vPos);
		if ( stat != EGADS_SUCCESS ) goto cleanup;
		// Collapsing cID with cID +2
		qm->uvs  [2*vC ] = uv  [0]; qm->uvs  [2*vC  + 1] = uv  [1];
		qm->uvs  [2*vC2] = uv  [0]; qm->uvs  [2*vC2 + 1] = uv  [1];
		qm->xyzs [3*vC ] = vPos[0]; qm->xyzs [3*vC  + 1] = vPos[1]; qm->xyzs[3*vC  + 2] = vPos[2];
		qm->xyzs [3*vC2] = vPos[0]; qm->xyzs [3*vC2 + 1] = vPos[1]; qm->xyzs[3*vC2 + 2] = vPos[2];
	}
	// Look at adjacent quads and point them to opposite adjacents
	int count = 0;
	for ( q = 0; q < 4; ++ q) {
		adjQ        = qm -> quadAdj[4*qC + q] - 1;
		updatedQuad = 0;
		for ( k = 0 ; k < 4; ++k)
			if ( qm -> quadIdx[4*adjQ + k] == vC2 + 1)
				updatedQuad = 1;
		if ( !updatedQuad ) continue;
		++count;
		stat = EG_getCommonEdge(qm, qC + 1,  adjQ + 1, edges);
		printf(" COMMON EDGES TO %d %d = %d %d  %d %d \n",qC + 1, adjQ + 1, edges[0], edges[1],edges[2], edges[3]);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		for ( e = 0 ; e < 2; ++ e) {
			if ( edges[2*e] != -1 && edges[2*e + 1] != -1) {
				vCommon = qm -> quadIdx[4*qC + edges[2*e]];
				printf(" edge %d  COMMON VERTEX %d \t",edges[2*e], vCommon);
				if (vCommon == vC2 + 1)
					vCommon = qm -> quadIdx[4*qC + (edges[2*e] + 1)%4];
				printf(" COMMON VERTEX %d \n",vCommon);
				for ( i = 0; i < 4; ++i) { // loop around quads and take the other quad that shares vertex vCommon
					auxAdj = qm -> quadAdj[4*qC + i ] - 1;
					//	printf(" AUX ADJ %d \n",auxAdj);
					if ( auxAdj >= 0 ) {
						for ( k = 0 ; k < 4; ++ k) {
							if ( auxAdj != adjQ  && qm -> quadIdx[4*auxAdj + k] == vCommon) {
								//k = 4;
								printf(" FOUND AN ADJACENCT to %d THAT HAS VCOMMON %d AND IT IS NOT %d\n", qC + 1, vCommon, adjQ + 1);
								print_quad_specs(qm, auxAdj + 1);
								print_quad_specs(qm, adjQ   + 1);
								for ( j = 0 ; j < 4; ++j) {
									printf(" j % d  ---------  %d \n", j, qm -> quadAdj[4*auxAdj + j] );
									if ( qm -> quadAdj[4*auxAdj + j] == qC + 1) { // found the other quad pointing at vCommon
										printf(" POINT %d TO %d", qm -> quadAdj[(4*adjQ + edges[2*e + 1] )], auxAdj + 1);
										qm -> quadAdj[(4*adjQ + edges[2*e + 1] )] = auxAdj + 1;
										printf(" POINT %d TO %d", qm -> quadAdj[4*auxAdj + j], adjQ + 1);
										qm -> quadAdj[4*auxAdj + j]               = adjQ   + 1;
									}
								}
							}
						}
					}
				}
			}
			for ( i = 0 ; i < 4; ++i)  // QC is a boundary quad so we need to point to the outside
				if ( qm -> quadAdj[4*adjQ + i] == qC + 1) qm -> quadAdj[4*adjQ + i] = -1;
		}
		// if ( count == 2) break;
	}
	for ( i = 0 ; i < 4; ++i) print_quad_specs(qm, qm ->quadAdj[4*qC + i]);
	// Transfer valencies from Collapsed vertex
	// Increase valence of collapsed vertex: 1-> copy old ones
	vVAL   = qm -> valence[vC][0];
	colVal = qm -> valence[vC][0] + qm -> valence[vC2][0] - 2;
	for ( i = 0 ; i < vVAL; ++i) auxVal[i] = qm -> valence[vC][2 + i];
	qm -> valence[vC][0] = colVal;
	// Assing a existing quad to vertex vC
	qm -> valence[vC][1] = - 1;
	for ( j = 0 ; j < 4; ++j) {
		for ( i = 0 ; i < 4; ++i) {
			q = qm->quadAdj[4*qC + j] - 1;
			if ( q >= 0) {
				if( qm -> quadIdx[ 4 * q + i ] == vC + 1 ) {
					qm -> valence[vC][1] = q + 1;
					break;
				}
			}
			if ( qm -> valence[vC][1] != -1 ) break;
		}
		if ( qm -> valence[vC][1] != -1 ) break;
	}
	for ( i = 0 ; i < vVAL; ++i) {
		qm -> valence[vC][2 + i] = auxVal[i];
		vID = qm -> valence[vC][2 + i] - 1;
		printf(" ALTER ID %d\n",vID + 1);
		if ( qm -> valence[vID][1] == qC + 1) {
			EG_getAdjacentToVertexPair(qm, qC + 1, vC + 1, vID + 1, edges);
			printf(" QUAD POINTED %d  \t", qm -> valence[vID][1]);
			qm -> valence[vID][1] = edges[1];
			printf(" NOW P OINTS  %d  \n", qm -> valence[vID][1]);
		}
	}
	// Get remaining valences from vertex 2 and elliminate vertex 2 as link from its neighbors
	k = vVAL + 2;
	for ( j = 0 ; j < qm -> valence[vC2][0]; ++j) {
		vID = qm -> valence[vC2][2 + j] - 1;
		printf(" ALTER ID %d\n",vID + 1);
		if ( qm -> valence[vID][1] == qC + 1) {
			printf(" QUAD POINTED %d  \t", qm -> valence[vID][1]);
			EG_getAdjacentToVertexPair(qm, qC + 1, vC2 + 1, vID + 1, edges);
			printf(" QUAD POINTED %d  \t", qm -> valence[vID][1]);
			qm -> valence[vID][1] = edges[1];
			printf(" NOW P OINTS  %d  \n", qm -> valence[vID][1]);
		}
		addValence = 1;
		for ( i = 0 ; i < vVAL; ++i) {
			if ( qm -> valence[vC][2 + i] == qm -> valence[vC2][2 + j]) {
				addValence = 0;
				i = vVAL;
			}
		}
		vID = qm -> valence[vC2][2 + j] - 1;
		aux = qm -> valence[vID][0];
		if ( addValence ) {
			qm -> valence[vC][k++] = vID + 1;  // increase k
			for ( v = 0 ; v < aux; ++v)
				if ( qm -> valence[vID][2 + v] == vC2 + 1) qm -> valence[vID][2 + v] = vC + 1;
		}
		else {
			for ( v = 0 ; v < aux; ++v) auxVal[v] = qm -> valence[vID][2 + v];
			--aux;
			qm -> valence[vID][0] = aux;
			for ( i = v = 0 ; v < aux + 1; ++v)
				if ( auxVal[v] != vC2 + 1)   qm -> valence[vID][2 + (i++)] = auxVal[v];
			// ensure that valence is pointing at valid quad
		}
		if ( qm -> valence[vID][1] == qC + 1) {
			printf(" ALTER ID %d\n",vID + 1);
			stat = EG_getAdjacentToVertexPair(qm , qC + 1, vID + 1, vC2 + 1, edges);
			printf(" QUAD POINTED %d  \t", qm -> valence[vID][1]);
			qm -> valence[vID][1] = edges[1];
			printf(" NOW POINTS  %d  \n", qm -> valence[vID][1]);
			if ( stat != EGADS_SUCCESS){
				printf(" I Can't find the other quad linking %d and %d near quad %d\n",vID+1, vC2 + 1, qC +1);
				for ( i = 0 ; i < 4; ++i) print_quad_specs(qm, qm->quadAdj[4*qC + i]);
				goto cleanup;
			}
		}
	}
	if ( k != qm -> valence[vC][0] + 2 ) {
		printf(" I MESSED UP TRANSFERING VALENCIES!!!\n");
		printf(" k %d  LENGTH VAL %d\n",k,qm -> valence[vC][0] + 2);
		stat = EGADS_INDEXERR;
		goto cleanup;
	}
	// Eliminate vertex vC2 from all the quads to which it belongs to
	for ( q = 0 ; q < star -> nQ; ++q) {
		for ( i = 0 ; i < 4; ++ i)
			if ( qm -> quadIdx[4*(star -> quads[q] -1) + i] == vC2 + 1) qm -> quadIdx[4*(star -> quads[q] -1) + i] = vC +  1;
	}
#ifdef DEBUG
	int pointsBack = 0;
	printf(" VC %d \n",vC + 1);
	for ( i = 0 ; i < qm -> valence[vC][0]; ++i) {
		printf(" VALENCE %d = %d  LETS CHECK THAT IT POINTS BACK \t",i, qm -> valence[vC][2 + i ]);
		vID = qm -> valence[vC][2 + i] - 1;
		for ( j = 0 ; j < qm -> valence[vID][0]; ++j) {
			if ( qm -> valence[vID][2 + j] == vC2 + 1) {
				printf(" VERTEX %d is still pointing at %d!!! We have just removed it\n",vID + 1, vC2 + 1);
			}
			else if ( qm -> valence[vID][2 + j] == vC + 1) {
				printf(" correct, valence (%d) = %d\n", j, qm -> valence[vID][2 + j]);
				pointsBack = 1;
			}
		}
		printf(" POINTS BACK %d \n",pointsBack );
		if (pointsBack != 1 ) {
			printf(" I MESSED UP THE VALENCES !@!@!@!!@!\n");
			exit(1);
		}

	}
#endif
	// Eliminate vertex quadID:
	printf(" I AM DESTROYING QUAD %d \n", qC + 1);
	for ( i = 0 ; i < 4; ++ i) {
		qm -> quadIdx[4*qC + i ] = -2;
		qm -> quadAdj[4*qC + i ] = -2;
	}
	print_quad_specs(qm, 2);
	print_quad_specs(qm, 12);
	// delete vertex vC2
	qm -> valence[  vC2    ][0] = -1;
	qm -> vType  [  vC2    ]    = -2; // -2 = removed
	qm -> uvs    [2*vC2    ]    = qm -> uvs[2*vC    ];
	qm -> uvs    [2*vC2 + 1]    = qm -> uvs[2*vC + 1];
    qm -> xyzs   [3*vC2    ]    = qm -> xyzs[3*vC    ];
	qm -> xyzs   [3*vC2 + 1]    = qm -> xyzs[3*vC + 1];
	qm -> xyzs   [3*vC2 + 2]    = qm -> xyzs[3*vC + 2];
	qm -> remQuads[0]++;
	qm -> remVerts[0]++;
	qm -> remQuads[qm -> remQuads[0]] = qC  + 1;
	qm -> remVerts[qm -> remVerts[0]] = vC2 + 1;
	print_quad_specs(qm, 2);
	print_quad_specs(qm, 12);
	cleanup:
	EG_free(star -> verts);
	EG_free(star -> quads);
	EG_free(star);
	//stat = resizeQm(qm);
	return stat;
}


static int
EG_removeDoublet(quadMap *qm, int vID) {
	int poly[2*4], polyOrd[2*4], i, j, qC, collId, stat, activity, link[2];
	vStar *star = NULL;
	if (qm -> vType[vID - 1] == -2) return EGADS_SUCCESS;
	if (qm -> vType[vID - 1] == 0 && qm -> valence[vID - 1][0] == 2) return EGADS_SUCCESS;
	if (qm -> vType[vID - 1]  > 0 && qm -> valence[vID - 1][0] == 2) return EGADS_GEOMERR;
	if (                             qm -> valence[vID - 1][0] != 2) return EGADS_SUCCESS;
	qC = qm -> valence[vID - 1][1] - 1;
	link[0] = qm -> valence[vID - 1][2];
	link[1] = qm -> valence[vID - 1][3];
	//print_quad_specs(qm, qC + 1);
	for ( i = 0 ; i < 4; ++i) {
		poly[2*i    ] = qm -> quadIdx[4*qC + i];
		poly[2*i + 1] = qm -> valence[poly[2*i] - 1][0];
		if (poly[2*i] == vID ) collId = i;
	}
	for ( i = 0 ; i < 4; ++i) {
		polyOrd[2*i    ] = poly[2*( (collId + i)%4)    ];
		polyOrd[2*i + 1] = poly[2*( (collId + i)%4) + 1];
	}
	printf(" REMOVE DOUBLET: DELETE V %d AND Q %d\n",polyOrd[0], qC + 1);
	snprintf(buffer,33,"DOUBLET_%d",DB);
	print_mesh(qm, buffer);
	stat = EG_mergeVertices(qm, qC + 1, polyOrd, 0, &activity);
	snprintf(buffer,33,"POSTDOUBLET_%d",DB);
	print_mesh(qm, buffer);
	DB++;
	printf(" EG_mergeVerticesSTAT %d  ACITVITY %d \n", stat, activity);
	if ( stat != EGADS_SUCCESS || activity == 0) {
		printf(" --->> IT CAN'T FAIL AT REMOVING DOUBLET! ");
		exit(1);
	}
#ifdef DEBUG
	printf(" IN DOUBLET WE HAVE COLLAPSED VERTEX %d to %d AND REMOVED QUAD %d\n", polyOrd[0], polyOrd[2*2], qC + 1);
#endif
	// Check that it has left valid elements at each old link
	for ( j = 0; j < 2; ++j )
		stat = EG_removeDoublet (qm, link[j]);
	for ( j = 0 ; j < 2; ++j) {
		if ( qm -> vType [link[j] - 1] != -2 ) {
			stat = EG_buildStar(qm, &star, link[j]);
			if ( stat != EGADS_SUCCESS)  return stat;
			for ( i = 0 ; i  < star->nQ;++i) {
				if ( star -> quads[i] < 0 ) continue;
				printf(" CHECK INVALID FOR COLLAPSED QUAD %d\n",star ->quads[i]);
				//stat = makeQuadValid(qm, star ->quads[i]);
				stat = checkInvalidElement(qm, star ->quads[i]);
				if ( stat != EGADS_SUCCESS) {
					printf(" IT CAN'T FAIL AT REMOVING DOUBLET! \n");
					exit(1);
					EG_free(star ->verts);
					EG_free(star ->quads);
					EG_free(star);
					return stat;
				}
			}
		}
	}
	EG_free (star -> verts);
	EG_free (star -> quads);
	EG_free (star         );
	printf(" IN REMOVE DOUBLET DELETE %d AND QUAD %d\n",polyOrd[0], qC + 1);
	++TOTCOLLAPSES;
	return EGADS_SUCCESS;
}






static int EG_collapseAndSplit(quadMap *qm, int qID) {
	int *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, *locVtype, **locVal = NULL, *locQrem = NULL, *locVrem = NULL ;
	int i, j, poly[2*4], polyOrd[2*4], vID, idCollapse, idSplit, stat, nV, nQ, triple = 1, activity;
	double *uv0;
	vStar *star;
	if ( qm -> quadIdx[4*(qID - 1)] < 0) return EGADS_SUCCESS;
	for ( i = 0 ;  i < 4; ++i) {
		printf(" Q %d  I %d  ID %d \n", qID, i,qm -> quadIdx[4*(qID - 1) + i ]);
		vID = qm -> quadIdx[4*(qID - 1) + i ] - 1;
		poly[2*i]     = vID + 1;
		poly[2*i + 1] = qm -> valence[vID][0];
		if (qm -> vType[vID] != -1 && poly[2*i + 1] <= 3 ) return EGADS_SUCCESS; //Operation valid only at interior quads
		triple *= poly[2*i + 1];
	}
	if ( triple != 3*5*5*4 ) return EGADS_SUCCESS;
	stat = backupQm ( qm, &nV, &nQ, &locIdx, &locAdj, &locVtype, &locVal, &locQrem, &locVrem, &uv0);
	for ( i = 0 ; i < 4; ++i)
		if ( poly[2*i + 1] == 3) break;
	for ( j = 0 ; j < 4; ++j) {
		polyOrd[2*j    ] = poly[2*( (j + i)%4)    ];
		polyOrd[2*j + 1] = poly[2*( (j + i)%4) + 1];
	}
	// Look for pair (3,5)-(4,4)  (or (3,4) - (5,5) )
	if ( qm -> vType[ polyOrd[2*2] -1 ] != -1 )
		stat = EG_mergeVertices(qm, qID, polyOrd, 0, &activity); // 0= pull vertex 0 to position 2 (v2 is fixed)
	else
		stat = EG_mergeVertices(qm, qID, polyOrd, 1, &activity); // 1 = reallocate collapsed at centre ( both 0,2) move
	if ( stat != EGADS_SUCCESS || !activity ) {
#ifdef DEBUG
		printf(" Inside EG_collapseAndSplit, EG_mergeVertices = %d\n", stat);
#endif
		goto cleanup;
	}
	if ( qm -> valence[ polyOrd[2*2] - 1][0] == 5 ) {
		stat = EGADS_SUCCESS;
		goto cleanup;
	}
	else if ( qm -> valence[ polyOrd[2*2] - 1][0] != 6 ) {
		printf(" I DONT GET IT. THIS SHOOULD BE VAL 6 BUT IT IS %d\n", qm -> valence[ polyOrd[2*2] - 1][0]);
		stat = resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal, locQrem, locVrem,  uv0);
		goto cleanup;
	}
	stat = EG_vertexSplit(qm, polyOrd[2*2]);
	if ( stat != EGADS_SUCCESS || !activity ) {
#ifdef DEBUG
		printf(" Inside EG_collapseAndSplit, EG_vertexSplit = %d\n", stat);
#endif
		stat += resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal,  locQrem, locVrem, uv0);
	}
	cleanup:
	for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
	EG_free(locVal);
	EG_free(locIdx);
	EG_free(locAdj);
	EG_free(locVtype);
	EG_free(uv0);
	EG_free(swapPtr);
	return stat;
}



static int
EG_vertexCollapse(quadMap *qm, int sID)  {
	int        i, j, k, q, v, qC, cID, vC, vM, vP, rV, id0, cp, vID, vVAL, addVal, nV, nQ;
	int        aux, auxID, quadID, stat, val, vCval,  qCount = 0, vCommon, updatedQuad, collapse, activity;
	int        adjQ, nMove, movingVerts[4], auxVal[100];  // local information about vertices and valence
	static int coll[2][4] = {{0, 2, 1, 3},{1, 3, 2, 0}}; // for collapsing edges
	double     uv[2], vPos[18], *uv0 = NULL;
	int        *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, *locVtype = NULL, **locVal = NULL, *locQrem = NULL, *locVrem = NULL, poly[2*4], polyOrd[2*4], vIrr[4];
	vStar      *star = NULL;
	if ( qm -> vType[sID -1] != -1 ) return EGADS_SUCCESS;
	// BACKUP: IN CASE OF INVALID COLLAPSE, RETRIEVE ORIGINAL MESH
	stat = EG_buildStar(qm, &star, sID);
	if ( stat != EGADS_SUCCESS) return stat;
	stat = backupQm ( qm, &nV, &nQ, &locIdx, &locAdj, &locVtype, &locVal, &locQrem, &locVrem, &uv0);
	if ( stat != EGADS_SUCCESS) {
		EG_free(star -> verts);
		EG_free(star -> quads);
		EG_free(star);
		return stat;
	}
	poly[0] = star -> verts  [0];
	poly[1] = qm   -> valence[poly[0] - 1][0];
	//print_star(star);
	for ( q = 0 ; q < star -> nQ; ++q) {
		collapse = 0;
		for ( i = 1; i < 4; ++i)
		{
			if ( (2 * q + i) == star -> nV ) aux = 1;
			else                             aux = 2 * q + i ;
			poly[2*i   ] = star -> verts  [aux   ];
			poly[2*i +1] = qm   -> valence[poly[2*i] - 1][0];
		}
		EG_vertexRegularity(qm, poly, vIrr, 4);
		for ( cp = 0 ; cp < 2; ++cp) {
			if ( (vIrr[coll[cp][0]] + vIrr[coll[cp][1]] == -2 ) && (vIrr[coll[cp][2]] + vIrr[coll[cp][3]] > 0 ) )
			{
				qC       = star -> quads[q] - 1;
				collapse = 1;
			}
			else if ( (vIrr[coll[cp][0]] + vIrr[coll[cp][1]] < 0 ) ) {
				if (     vIrr[coll[cp][2]] * vIrr[coll[cp][3]] > 0   )
				{
					qC       = star -> quads[q] - 1;
					collapse = 1;
				}
				else if ( (vIrr[coll[cp][2]] == -1 && vIrr[coll[cp][3]] > 0)  ||
						(      vIrr[coll[cp][3]] == -1 && vIrr[coll[cp][2]] > 0) )
				{
					qC       = star -> quads[q] - 1;
					collapse = 1;
				}
			}
			if (collapse) break;
		}
		if (collapse ) break;
	}
	if ( collapse != 1) goto cleanup;
	for ( i = 0 ; i < 4; ++i) {
		polyOrd[2*i    ] = poly[2*( (coll[cp][0] + i)%4)    ];
		polyOrd[2*i + 1] = poly[2*( (coll[cp][0] + i)%4) + 1];
	}
	if ( qm -> vType[ polyOrd[2*2] - 1 ] == -1 )
		stat   = EG_mergeVertices(qm, qC + 1, polyOrd,  1, &activity );
	else
		stat = EG_mergeVertices(qm, qC + 1, polyOrd,  0 , &activity);
	printf(" EG MERGE VERTICES %d   \n",stat);
	if (stat != EGADS_SUCCESS || !activity ) goto cleanup;
	stat = EG_buildStar(qm, &star, polyOrd[2*2]);
	if ( stat!= EGADS_SUCCESS) goto cleanup;
#ifdef DEBUG
	printf(" IN VERTEX COLLAPSE WE HAVE COLLAPSED VERTEX %d to %d AND REMOVED QUAD %d\n",polyOrd[0], polyOrd[2*2], qC + 1);
#endif
	for ( j = 0 ; j < star -> nV; ++j) {
		auxID =  - 1;
		if ( star -> verts[j] > 0) {
			auxID = star -> verts[j] - 1;
			if (qm -> valence[auxID][0] == 2 && qm -> vType[auxID] == -1) {
#ifdef DEBUG
				printf(" WE HAVE PRODUCED A DOUBLET AT %d\n",auxID + 1);
#endif
				stat = EG_removeDoublet(qm, auxID +1);
				if ( stat != EGADS_SUCCESS) goto cleanup;
			}
		}
	}
	// Optimize around quads
	nMove = 0 ;
	for ( i = 1 ; i < 4; ++i)
		if (  qm -> vType[polyOrd[2*i] -1] == -1 ) movingVerts[nMove++] = polyOrd[2*i];
	stat = optimize_angles(qm, nMove, movingVerts);
	if ( stat == EGADS_SUCCESS)   ++TOTCOLLAPSES;
	cleanup:
	if ( stat != EGADS_SUCCESS)
		stat = resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal,  locQrem, locVrem,  uv0);
	EG_free(star -> quads);
	EG_free(star -> verts);
	EG_free(star);
	for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
	EG_free(locVal);
	EG_free(locIdx);
	EG_free(locAdj);
	EG_free(locVtype);
	EG_free(uv0);
	EG_free(swapPtr);
	return stat;
}

static int
swappingOperation(quadMap *qm, vStar *star, int poly[], int swapIdx, int qID[], int pullVertex, int *activity) {
	int OK, stat, i, j, k, aux, auxVAL, auxID, auxQ, pos, it, adjQ, areaType, area, centreValency;
	static int perm[6] = {0, 1, 2, 3, 4, 5}, remV[4] = {0, 3, 3, 0};
	int Q1[4], Q2[4], edges[6], q0ID[6], q1ID[6], valence[100], v[4], addV[4];
	*activity = 0;
	addV[0] = swapIdx; addV[1] = swapIdx + 3;
	addV[2] = addV[1]; addV[3] = addV[0];
	printf(" TRY SWAP ========   %d, %d  <->  %d %d\n", remV[0], remV[1], addV[0], addV[1]);
	printf(" VERTICES ========   %d, %d  <->  %d %d\n", poly[2*remV[0]], poly[2*remV[1]], poly[2*addV[0]], poly[2*addV[1]]);
	snprintf(buffer,33,"swappOp%d", SO);
	print_mesh(qm, buffer);
	print_quad_specs(qm, qID[0]);
	print_quad_specs(qm, qID[1]);
	// +1 bias
	qID[0]--; qID[1]--;
	// We store the adjacent quad to each edge (ordered). LINKED EDGE 0-3
	for (j = 0 ; j < 4; ++j) if ( qm -> quadIdx[4*qID[0] + j] == poly[2*0] ) break;
	for (k = 0 ; k < 4; ++k) if ( qm -> quadIdx[4*qID[1] + k] == poly[2*3] ) break;
	if ( j >= 4 || k >= 4) {
		print_quad_specs(qm, qID[0] + 1);
		print_quad_specs(qm, qID[1] + 1);
		printf(" polys %d  =%d  %d =%d \n", poly[2*0],poly[1], poly[2*3],poly[2*2 + 1]);
		exit (1);
	}
	for ( i = 0 ; i <= 2 ; ++i) {
		edges[i    ] = qm -> quadAdj[4*qID[0] + (j + i)%4];
		edges[i + 3] = qm -> quadAdj[4*qID[1] + (k + i)%4];
	}
	for ( i = 0 ; i < 6 ; ++i) {
		printf(" adj %d = %d\n",i, edges[i]);
	}
	*activity = 1;
	// MODIFY CORRESPONDING QUADS IN MAP
	for ( i = 0 ; i < 4; ++ i) {
		Q1[i] = perm[(i + addV[0])%6];
		Q2[i] = perm[(i + addV[1])%6];
		printf(" Q1 %d  Q2 %d   ----> V1 %d V2  %d  \n",Q1[i], Q2[i], poly[2*Q1[i]], poly[2*Q2[i]]);
	}
	for ( i = 0; i <= 2; ++i) {
		q0ID[i    ] = qID[0] + 1;
		q0ID[i + 3] = qID[1] + 1;
		q1ID[Q1[i]] = qID[0] + 1;
		q1ID[Q2[i]] = qID[1] + 1;
	}
	for ( j = 0 ; j < 4; ++j) {
		qm -> quadIdx [4*qID[0] + j    ]    = poly[2*Q1[j]];
		qm -> quadIdx [4*qID[1] + j    ]    = poly[2*Q2[j]];
		qm -> valence [poly[2*Q1[j]] -1][1] = qID[0] + 1;  // Point valence at the appropriate quad
		qm -> valence [poly[2*Q2[j]] -1][1] = qID[1] + 1; // Point valence at the appropriate quad
		qm -> quadAdj [4*qID[0] + j    ]    = edges[Q1[j]];
		qm -> quadAdj [4*qID[1] + j    ]    = edges[Q2[j]];
	}
	qm -> quadAdj[4*qID[0] + 3] = qID[1] + 1;
	qm -> quadAdj[4*qID[1] + 3] = qID[0] + 1;
	print_quad_specs(qm, qID[0] +1);
	print_quad_specs(qm, qID[1] +1);
	for ( j = 0; j < 6; ++j) {
		for ( i = 0; i < 4; ++i) {
			if ( edges[j] != -1 ) {
				printf(" QUAD %d pointed AT %d \t ", edges[j], qm -> quadAdj[4*(edges[j] - 1) + i ]);
				if ( qm -> quadAdj[4*(edges[j] - 1) + i ] == q0ID[j]) qm -> quadAdj[4*(edges[j] - 1) + i ] = q1ID[j];  // Modify adjacent quads
				printf(" NOW   AT %d \n ", qm -> quadAdj[4*(edges[j] - 1) + i ]);
			}
		}
	}
#ifdef DEBUG
	for ( i = 0; i < 6; ++ i) {
		OK = 0;
		if (edges[i] >= 0) {
			for ( k = 0; k < 4; ++k) {
				printf(" EDGE %d \t COMP TO %d  \n", edges[i], q1ID[i]);
				if (qm ->quadAdj[4*(edges[i] - 1)  + k ] == q1ID[i] ) {
					OK = 1; k  = 5;
				}
			}
			printf(" EDGE %d \t COMP TO %d  \n", edges[i], q1ID[i]);
			if (OK == 0 ) {
				printf("Quad at edge %d = %d \n", i, edges[i]);
				print_quad_specs(qm, edges[i]);
				printf(" ONE ADJACENT SHOULD POINT AT %d !!!\n",q1ID[i]);
				exit(1);
			}
		}
	}
#endif
	// ASSOCIATE THE CORRECT QUAD TO EACH VERTEX
	for ( k = 0 ; k < 4; ++k){
		if (qm -> quadIdx[4*qID[0] + k] == poly[0]) break;
		else if (qm -> quadIdx[4*qID[0] + k] == poly[2*3]) {
			aux    = qID[0];
			qID[0] = qID[1];
			qID[1] = aux;
			break;
		}
	}
	// Now update the valencies : REMOVE PAIR 0-3 AND ADD  PAIR 1-4 ( 2-5)
	for ( j = 0 ; j < 2; ++j) {
		auxID  = poly[2*remV[2*j]   ] - 1; //central vertex   remV is a dumb array {0,3,3,0} to make my life easier
		auxVAL = poly[2*remV[2*j] +1]    ;
		if ( auxID < 0) {
			printf(" Found a problem whilst swapping. We are not going to swap that edge\n");
			return EGADS_INDEXERR;
		}
		for ( i = 0 ; i < auxVAL; ++i) valence[i] = qm -> valence[auxID][2 + i];
		--auxVAL;
		qm -> valence[auxID][0] = auxVAL;  // TOTAL VALENCY
		qm -> valence[auxID][1] = qID[j] + 1;
#ifdef DEBUG
		OK = 0;
		for ( i = 0 ; i < 4; ++i)
			if ( qm -> quadIdx[4*(qm -> valence[auxID][1] -1 ) + i ] == auxID + 1) OK = 1;
		if (!OK) {
			printf(" I have messed up assigining valence to quad. Vertex %d is supposed to belong to quad %d \t",auxID, qm -> valence[auxID][1]);
			print_quad_specs(qm, qm -> valence[auxID][1]);
		}
#endif
		for ( k = i = 0 ; i <= auxVAL; ++i) {
			if ( valence[i] != poly[2*remV[2*j+1]]){
				qm -> valence[auxID][2 + k] = valence[i];
				++k;
			}
		}
		auxID  = poly[2*addV[2*j]    ] - 1; //central vertex   addV is a dumb array for the same reason
		auxVAL = poly[2*addV[2*j] + 1];
		auxQ   = qm -> valence[auxID][1];
		if ( auxVAL != qm -> valence[auxID][0]) {
			printf(" ALERT! POLY DOESNT COINCIDE WITH VALENCE POLY%d VAL %d \n", auxVAL,qm -> valence[auxID][0]);
		}
		for ( i = 0 ; i < auxVAL; ++i) valence[i] = qm -> valence[auxID][2 + i];
		++auxVAL;
		qm -> valence[auxID][0] = auxVAL;
		qm -> valence[auxID][1] = auxQ;
		for ( i = 0 ; i < auxVAL ; ++i) qm -> valence[auxID][i + 2] = valence[i];
		qm ->valence[auxID ][2 + auxVAL - 1] = poly[2*addV[2*j + 1]];
	}
	// Finally: Check if we have created a doublet
	for ( j = 0 ; j < star -> nV; ++j) {
		auxID =  - 1;
		if ( star -> verts[j] > 0) {
			auxID = star -> verts[j] - 1;
			if (qm -> valence[auxID][0] == 2 && qm -> vType[auxID] == -1) {
				//#ifdef DEBUG
				printf(" IN EDGE SWAP. WE HAVE PRODUCED A DOUBLET AT %d\n",auxID);
				//#endif
				stat = EG_removeDoublet(qm, auxID +1);
				if ( stat != EGADS_SUCCESS) return stat;
			}
		}
	}
	snprintf(buffer,33,"swappOpAFF%d", SO);
	print_mesh(qm, buffer);
	++SO;
	printf(" WE HAVE SWAPPED ========   %d, %d  <->  %d %d\n", remV[0], remV[1], addV[0], addV[1]);
	return EGADS_SUCCESS;
}


static int
EG_edgeSwap(quadMap *qm, int sID) {
	int stat, i, j, k, aux, vL, movingVerts[100], nMove, v[4*2], vIrr[6], qID[2], link, swapped, nV, nQ, activity, invalid;
	int    *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, *locVtype = NULL, **locVal = NULL, *locQrem = NULL, *locVrem = NULL, poly[2*6];
	double fOpti, *uv0 = NULL;
	vStar *star = NULL, *optiStar = NULL;
	if ( sID > qm -> totVerts || sID <= 0)
		printf(" EIN??  WE ARE USING ID %d  \n",sID);
	if (  qm -> vType[sID - 1] == -2) return EGADS_SUCCESS;
	if ( (qm -> vType[sID - 1] ==  0) && qm -> valence[sID - 1][0] == 2) return EGADS_SUCCESS;
	if ( (qm -> vType[sID - 1]  >  0) && qm -> valence[sID - 1][0] == 3) return EGADS_SUCCESS;
	stat = EG_buildStar(qm, &star, sID);
#ifdef DEBUG
	printf(" EG_edgeSwap around vertex %d\n",sID);
#endif
	if ( stat != EGADS_SUCCESS) return stat;
	stat = backupQm ( qm, &nV, &nQ, &locIdx, &locAdj, &locVtype, &locVal, &locQrem, &locVrem, &uv0);
	if ( stat != EGADS_SUCCESS) goto cleanup;
	poly[0] = star-> verts[0];
	poly[1] = star-> nQ;
	swapPtr = (int *) EG_alloc(poly[1] *sizeof(int));
	if (swapPtr == NULL) {
		stat = EGADS_MALLOC;
		goto cleanup;
	}
	swapped = 0;
	for ( link = 0; link < poly[1]; ++link) {
		qID[0]        = star -> quads[link];
		swapPtr[link] = -1;
		vL            = 2*link + 1;
		if ( link + 1 == star -> nQ ) qID[1] = star -> quads[0];
		else    qID[1] = star -> quads[link + 1];
		if ( qID[0] == -1 || qID[1] == -1 ) continue;
		for ( i = 0 ; i  < 5; ++i) {
			aux = vL + i;
			if ( aux >= star ->nV) aux -= (star->nV -1);
			poly[2*(i+1)   ] = star ->verts[aux];
			if ( poly[2*(i+1)] < 0) { //This is a virtual quad used when the star centre is at the boundary. Skip
				poly[2*(i+1) +1] = -1;
				vL = -1;
				break;
			}
			else poly[2*(i+1) +1] = qm -> valence[poly[2*(i+1)] -1][0];
		}
		if (vL == -1) continue;
		EG_vertexRegularity(qm, poly, vIrr, 6);
		//#ifdef DEBUGG
		for ( i = 0 ; i < 6; ++i) printf(" VERT %d = %d \n",poly[2*i], poly[2*i+1]);
		//#endif
		// CHECK IF WE SHOULD SWAP
		if (qm -> vType[ poly[0] -1] == 0) { // this is a domain corner. SWAP FORCED
			if ( (vIrr[1] + vIrr[4]) < (vIrr[2] + vIrr[5]) )  swapPtr[link] = 1;
			else swapPtr[link] = 2;
		}
		else if (qm -> vType[ poly[0] -1] > 0 )
		{
			if ( (vIrr[1] != 100 ) && (vIrr[4] != 100) )
			{
				// not vertex
				if ( (qm -> vType[ poly[2*1] -1] > 0  && poly[2*1 + 1] < 5) ||
						( qm -> vType[ poly[2*4] -1] > 0  && poly[2*4 + 1] < 5) )
				{
					if ( (vIrr[1] + vIrr[4] <  vIrr[2] + vIrr[5] )  //higher valence that the other swapping option
							&& (vIrr[1] + vIrr[4] <= vIrr[0] + vIrr[3]) ) // improves or preserves current valence
						swapPtr[link] = 1;
				}
				else if ( (vIrr[1] + vIrr[4] < vIrr[2] + vIrr[5] )  //higher valence that the other swapping option
						&&    (vIrr[1] + vIrr[4] <= vIrr[0] + vIrr[3]) ) // improves or preserves current valence
					swapPtr[link] = 1;
			}
			if ( (vIrr[2] != 100 ) && (vIrr[5] != 100 ) )
			{
				// not vertex
				if ( (qm -> vType[ poly[2*2] -1] > 0  && poly[2*2 + 1] < 5) ||
						(qm -> vType[ poly[2*5] -1] > 0  && poly[2*5 + 1] < 5) ) {
					if ( (vIrr[2] + vIrr[5]) <= (vIrr[0] + vIrr[3]) ) // improves or preserves current valence
						swapPtr[link] = 2;
				}
				else if ( (vIrr[1] + vIrr[4] <= vIrr[0] + vIrr[3]) ) // improves or preserves current valence
					swapPtr[link] = 2;
			}
		}
		else {
			if   ( (vIrr[0] + vIrr[3]) >  0) {
				if      ( (vIrr[1] + vIrr[4]) == -2)  swapPtr[link] = 1;
				else if ( (vIrr[2] + vIrr[5]) == -2)  swapPtr[link] = 2;
				else if (( vIrr[0] + vIrr[3]) >  5)
				{ // two high valencies
					if      ( (qm -> vType[poly[2*1] - 1] == -1 && qm -> vType[poly[2*4] - 1] == -1 ) &&
							( (vIrr[1] == -1 && vIrr[4] == 0) || ( vIrr[4] == -1 && vIrr[1] == 0) ) ) swapPtr[link] = 1;
					else if ( (qm -> vType[poly[2*2] - 1] == -1 && qm -> vType[poly[2*5] - 1] == -1) &&
							( (vIrr[2] == -1 && vIrr[5] == 0) || ( vIrr[5] == -1 && vIrr[2] == 0) ) ) swapPtr[link] = 2;
				}
			}
		}
		if ( swapPtr[link] == -1) continue;
		stat = swappingOperation(qm, star, poly, swapPtr[link], qID, 1, &activity);
		if ( stat != EGADS_SUCCESS && activity ) {
			printf(" INVALID AREA. RESET STUFF \n");
			stat = resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal,  locQrem, locVrem,  uv0);
			goto cleanup;
		}
		if ( activity ) {
			stat = EGADS_SUCCESS;
			v[0] = poly[0];
			v[1] = poly[2*3];
			for ( j = 0 ; j < 2; ++j) {
				if (qm -> vType[v[j] - 1] != -2 ) {
					stat = EG_buildStar(qm, &optiStar, v[j]);
					printf(" BUILD STAR STAT %d \n",stat);
					if ( stat != EGADS_SUCCESS) {
						stat = resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal,  locQrem, locVrem,  uv0);
						goto cleanup;
					}
					for ( i = 0 ; i < optiStar -> nQ; ++i) {
						printf(" CALL INVALID %d\n", optiStar -> quads[i]);
						if ( optiStar -> quads[i] > 0)
							stat  = checkInvalidElement(qm, optiStar -> quads[i]);
						printf( " STAT %d \n",stat);
						if ( stat != EGADS_SUCCESS) {
							stat = resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal,  locQrem, locVrem,  uv0);
							goto cleanup;
						}
					}
				}
			}
			snprintf(buffer,33,"SWAP_NO_%d",NSWAP);
			print_mesh(qm, buffer);
			++NSWAP;
			swapped = 1;
		}
		printf(" OUT \n");
		snprintf(buffer,33,"SWAP_NO_%d",NSWAP);
		print_mesh(qm, buffer);
		++NSWAP;
		if ( swapped == 1 ) break;
		printf(" ? ddd\n");
	}
	snprintf(buffer,33,"SWAP_NO_%d",NSWAP);
	print_mesh(qm, buffer);
	++NSWAP;
	if ( swapped ) {
		snprintf(buffer,33,"SWAP_NO_%d",NSWAP);
		++NSWAP;
		print_mesh(qm, buffer);
		nMove = 0 ;
		for ( j = 0 ; j < 6; ++j) {
			printf(" _--------------   j %d  POLY %d  type %d ----------\n",
					j, poly[2*j], qm -> vType[poly[2*j] - 1]);
			if( qm -> vType[poly[2*j] - 1] == -1 ) {
				movingVerts[nMove++] = poly[2*j];
				stat = EG_buildStar(qm, &optiStar, poly[2*j]);
				print_star(star);
			}
		}
		for ( j = 0 ; j < 6; ++j) {
			for ( i = 0 ; i < qm -> valence[poly[2*j] - 1][0]; ++i) {
				aux = 1;
				if( qm -> vType[qm -> valence[poly[2*j] - 1][2 + i] -1] != -1 ) continue;
				for ( k = 0 ; k < nMove; ++k) {
					if ( qm -> valence[poly[2*j] - 1][2 + i] == movingVerts[k] ) {
						aux = 0;
						break;
					}
				}
				if ( aux )
					movingVerts[nMove++] = qm -> valence[poly[2*j] - 1][2 + i];
			}
		}
		stat = optimize_angles(qm, nMove, movingVerts);
		if ( stat != EGADS_SUCCESS) {
			stat =  resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal, locQrem, locVrem,  uv0);
		}
		else {
			snprintf(buffer,33,"OPTISWAP_NO_%d",NSWAP);
			print_mesh(qm, buffer);
			++NSWAP;
			++TOTSWAPS;
		}
	}
	cleanup:
	EG_free(star->verts);
	EG_free(star->quads);
	EG_free(star);
	for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
	EG_free(locVal);
	EG_free(locIdx);
	EG_free(locAdj);
	EG_free(locVtype);
	EG_free(locQrem);
	EG_free(locVrem);
	EG_free(uv0);
	EG_free(swapPtr);
	if ( optiStar ) {
		EG_free(optiStar -> quads);
		EG_free(optiStar -> verts);
		EG_free(optiStar);
	}
	//if ( stat != EGADS_SUCCESS )
	printf(" STAT  SWAPPING OPERATION %d\n",stat);
	return stat;
}



static int EG_doubleSplit(quadMap *qm, int sID ) {
	int    aux, v1, v2, s, i, id0, j, q, qq, sVAL, dist, d, stat, validSwap, qIDs[2],  nV, nQ, activity;
	int    swapID,  *locIdx = NULL, *locAdj = NULL, *locVtype = NULL, **locVal = NULL, *locQrem = NULL, *locVrem = NULL,  poly1[2*3], poly2[2*3];
	double *uv0 = NULL;
	char buff[33];
	vStar  *splitStar = NULL, *star = NULL;
	if ( qm -> vType[sID -1 ] != -1 || qm   -> valence[sID - 1][0] != 5 ) return EGADS_SUCCESS;
	stat = EG_buildStar(qm, &star, sID);
	if ( stat != EGADS_SUCCESS) return stat;
	sVAL = qm -> valence[sID - 1][0];
	// FIRST: LOOK AT POSSIBLE SWAP
	stat = backupQm ( qm, &nV, &nQ, &locIdx, &locAdj, &locVtype, &locVal, &locQrem, &locVrem, &uv0);
	for ( s = 0 ; s < sVAL; ++s) {
		id0 = 2*s + 1; // For splitting, need a few quads in between
		v1  = star -> verts[2*s + 1];
		if ( qm -> valence[v1 -1][0] != 3 || qm ->vType[v1 -1] != -1) continue;
		printf(" STAR WITH CENTRE %d  AND VAL %d HAS LINK %d  WITH VAL %d\n", sID, sVAL, v1, qm -> valence[v1 -1][0]);
		snprintf(buffer, 33, "DOUBLESPLIT_%d",DS);
		print_mesh(qm, buffer);
		// LOOK at quad +1:
		for ( d = 0 ; d <= sVAL - 4 ; ++d) {
			dist =   4 + 2*d;
			aux  = id0 + dist;
			if ( aux >= star -> nV) aux -= (star -> nV - 1) ;
			v2   = star -> verts[aux];
			if ( qm -> valence[v2 -1][0] != 5) continue;
			printf(" STAR WITH CENTRE %d  AND VAL %d HAS LINK %d WITH VAL %d\n", sID, sVAL, v2, qm -> valence[v2 -1][0]);
			stat = EG_buildStar(qm, &splitStar, v2);
			if (stat != EGADS_SUCCESS) {
				EG_free(star -> verts);
				EG_free(star -> quads);
				EG_free(star);
				return EGADS_SUCCESS;
			}
			// FIND STAR OUR VERTEX IN NEW STAR
			q = - 1;
			for ( i = 0 ; i < splitStar -> nQ; ++i)
				if ( splitStar -> verts[2*i + 1] == sID) {
					q = 2*i + 1; break;
				}
			if ( q == -1 ) continue;
			poly1[0] = v2;
			poly1[1] = qm -> valence[v2 - 1][0];
			poly1[2] = sID;
			poly1[3] = qm -> valence[sID - 1][0];
			dist = 6 - 2*d;
			aux  = q + dist;
			if ( aux >= splitStar -> nV) aux -= (splitStar -> nV - 1) ;
			poly1[4] = splitStar -> verts  [aux        ];
			poly1[5] = qm        -> valence[poly1[4] -1][0];
			if ( poly1[5] > 4 ) continue;  // this vertex will already be +1, allow maximum valence = 5
			printf(" FOUND ONE!!!!!!!!!!!! LETS SPLIT FIRST \n");
			snprintf(buffer, 33, "BEFOREDOUBLESPLIT_1DONE_%d",DS);
			print_mesh(qm, buffer);
			stat = EG_splittingOperation(qm, splitStar, poly1, q, dist, &activity);
			snprintf(buffer, 33, "DOUBLESPLIT_1DONE_%d",DS);
			print_mesh(qm, buffer);
			if ( stat != EGADS_SUCCESS && activity ) {
				stat =  resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal, locQrem, locVrem,  uv0);
				goto cleanup;
			}
			if ( stat == EGADS_SUCCESS && activity ) {
				printf(" GREAT!!! WE HAVE DONE THE FIRST SUCCESSFULLY\n");
				poly2[0] = sID;
				poly2[1] = qm -> valence[sID - 1][0];
				if ( poly2[1] != 6){
					printf(" I HAVE NO IDEA WHAT IVE DONT!!!!!!!!  VALENCE OF CENTRE  %d \n", poly2[1]);
					exit(1);
				}
				poly2[2] = v1;
				poly2[3] = qm -> valence[v1 - 1][0];
				if ( poly2[3] != 3){
					printf(" NOT VALENCE 3 I HAVE NO IDEA WHAT IVE DONT!!!!!!!!\n");
					exit(1);
				}
				stat = EG_buildStar(qm, &splitStar, sID);
				q = - 1;
				for ( i = 0 ; i < splitStar -> nQ; ++i)
					if ( splitStar -> verts[2*i + 1] == v1) {
						q = 2*i + 1; break;
					}
				aux  = q + 6;
				if ( aux >= splitStar -> nV) aux -= (splitStar -> nV - 1) ;
				poly2[4] = splitStar -> verts  [aux        ];
				poly2[5] = qm        -> valence[poly2[4] -1][0];
				if ( poly2[5] != 3 ) {
					printf(" I MESSEED UP WITH WHERE TO SWAP\n");
					exit(1);
				}
				stat = EG_splittingOperation(qm, splitStar, poly2, q, 6, &activity );
				snprintf(buffer, 33, "DOUBLESPLIT_2DONE_%d",DS);
				print_mesh(qm, buffer);
				++DS;
				if ( stat == EGADS_SUCCESS && activity ) goto cleanup;
				else{
#ifdef DEBUG
					printf(" SPLITTING WENT %d \n",stat);
#endif
					EG_free(star -> verts);
					EG_free(star -> quads);
					EG_free(star);
					stat =  resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal, locQrem, locVrem,  uv0);
					for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
					EG_free(locVal);
					EG_free(locIdx);
					EG_free(locAdj);
					EG_free(locVtype);
					EG_free(locQrem);
					EG_free(locVrem);
					EG_free(uv0);
					return stat;
				}
			}
		}
	}
	cleanup:
	if ( locVal)
		for ( i = 0 ; i < nV; ++i) EG_free(locVal[i]);
	EG_free(locVal);
	EG_free(locIdx);
	EG_free(locAdj);
	EG_free(locVtype);
	EG_free(uv0);
	EG_free(locQrem);
	EG_free(locVrem);
	EG_free(star->verts);
	EG_free(star->quads);
	EG_free(star);
	if ( splitStar ) {
		EG_free(splitStar->verts);
		EG_free(splitStar->quads);
	}
	EG_free(splitStar);
	return EGADS_SUCCESS;
}






static int EG_swapAndSplit(quadMap *qm, int sID ) {
	int    auxj, v1, v2, s, i, id0, j, q, qq, sVAL, dist, stat, validSwap, qIDs[2],  nV, nQ, activity;
	int    swapID,  *locIdx = NULL, *locAdj = NULL, *locVtype = NULL, **locVal = NULL, *locQrem = NULL, *locVrem = NULL, vSwap[2*6];
	double *uv0 = NULL;
	char buff[33];
	vStar  *swapStar = NULL, *star = NULL;
	if ( qm -> vType[sID -1 ] != -1 || qm   -> valence[sID - 1][0] < 5 ) return EGADS_SUCCESS;
	stat = EG_buildStar(qm, &star, sID);
	if ( stat != EGADS_SUCCESS) return stat;
	sVAL = qm -> valence[sID - 1][0];
	// FIRST: LOOK AT POSSIBLE SWAP
	for ( s = 0 ; s < sVAL; ++s){
		v1   = star -> verts[2*s + 1];
		if ( qm ->valence[v1 -1][0] != 3 || qm ->vType[v1 -1] != -1) continue;
		print_star(star);
		// LOOK at quad +1:
		for ( dist = 0 ; dist < 2; ++dist) {
			q    = (s + 2 + dist ) % star ->nQ;
			v2   = star -> verts[2*q + 1];
			if ( qm -> valence[v2 -1][0] != 4) continue;
			stat = EG_buildStar(qm, &swapStar, v2);
			print_star(swapStar);
			if (stat != EGADS_SUCCESS) {
				EG_free(star -> verts);
				EG_free(star -> quads);
				EG_free(star);
				return EGADS_SUCCESS;
			}
			q = - 1;
			for ( i = 0 ; i < swapStar -> nQ; ++i)
				if ( swapStar -> verts[2*i + 1] == sID) {
					q = i; break;
				}
			if ( q == -1 ) continue;
			if ( dist == 0 ) {
				qIDs[0] = swapStar -> quads[q];
				qq      = (q  + 1) %swapStar -> nQ;  //FWD
				qIDs[1] = swapStar -> quads[qq];
				swapID  = 1;
			} else  {
				q--;
				if ( q < 0 ) q = swapStar->nQ - 1;
				qIDs[1] = swapStar -> quads[q];
				q--;
				if ( q < 0 ) q = swapStar->nQ - 1;
				qIDs[0] = swapStar -> quads[q];
				swapID  = 2;
			}
			vSwap[2*0    ] = swapStar -> verts[0];
			vSwap[2*0 + 1] = qm       -> valence[ vSwap[0] -1][0];
			for ( j = 0; j < 5; ++j) {
				auxj                     = (2 * q + 1 + j);
				if ( auxj >= swapStar -> nV ) auxj -= (swapStar -> nV - 1);
				vSwap[2 * ( j + 1)     ] = swapStar -> verts  [auxj];
				vSwap[2 * ( j + 1)  + 1] =       qm -> valence[vSwap[2*(j + 1)] - 1 ][0];
			}
			validSwap = 0 ;
			if (     dist == 0 && ( vSwap[2*3 + 1] > 4 || vSwap[2*4 + 1] == 3) )
				validSwap = 1;
			else if (dist == 1 && ( vSwap[2*3 + 1] > 4 || vSwap[2*2 + 1] == 3) )
				validSwap = 1;
			if ( validSwap) break;
		}
		if ( validSwap) break;
	}
	if (validSwap) {
#ifdef DEBUG
		printf(" WE FOUND A VALID SWAP AND SPLIT! \n");
		printf(" FIRS SWAP QUADS %d %d and create edge %d %d\n", qIDs[0], qIDs[1], vSwap[2*1], vSwap[2*5] );
#endif
		// NOW LOOK AT POSSIBLE SPLIT IN STAR
		stat = backupQm ( qm, &nV, &nQ, &locIdx, &locAdj, &locVtype, &locVal, &locQrem, &locVrem, &uv0);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		stat = swappingOperation(qm, swapStar, vSwap, swapID, qIDs, 1, &activity);
		if (stat != EGADS_SUCCESS || !activity) {
			stat = resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal,  locQrem, locVrem, uv0);
			goto cleanup;
		}
		stat  = checkInvalidElement(qm, qm -> valence[vSwap[0  ] -1][1]);
		stat += checkInvalidElement(qm, qm -> valence[vSwap[2*3] -1][1]);
		if ( stat != EGADS_SUCCESS) {
			stat = resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal,  locQrem, locVrem,  uv0);
			goto cleanup;

		}
		stat = EG_vertexSplit(qm, sID);
		if ( stat != EGADS_SUCCESS || !activity) {
			stat = resetQm(qm, nV, nQ, locIdx, locAdj, locVtype, locVal,  locQrem, locVrem, uv0);
		}
		snprintf(buff,33,"ss_%i.txt",sscount);
		print_mesh(qm,buff);
		++sscount;
	}
	cleanup:
	if ( locVal)
		for ( i = 0 ; i < nV; ++i) EG_free(locVal[i]);
	EG_free(locVal);
	EG_free(locIdx);
	EG_free(locAdj);
	EG_free(locVtype);
	EG_free(locQrem);
	EG_free(locVrem);
	EG_free(uv0);
	EG_free(star->verts);
	EG_free(star->quads);
	EG_free(star);
	if ( swapStar ) {
		EG_free(swapStar->verts);
		EG_free(swapStar->quads);
	}
	EG_free(swapStar);
	return EGADS_SUCCESS;
}




static int
EG_vertexOperator(quadMap *qm, int vID, int opType, int *activity) {
	int stat, count, i;
	*activity = 0;
	if ( vID <=0 || vID > qm -> totVerts) {
		printf(" WHAT!   VERTEX %d \n", vID);
	}
	for ( i = 0 ; i < qm -> nQ; ++i) print_quad_specs(qm, i+1);
	switch(opType) {
	case TRYSWAPPING:
		count = TOTSWAPS;
		stat  = EG_edgeSwap(qm, vID);
		if ( stat != EGADS_SUCCESS) return stat;
		if ( TOTSWAPS > count) {
			printf(" \n\n ======  EG_vertexOperator SWAPPED VERTEX COUNT %d \n", TOTSWAPS);
			*activity = 1;
		}
		break;
	case TRYCOLLAPSING:
		count = TOTCOLLAPSES;
		//if ( qm -> vType[vID - 1] != -1 )return EGADS_SUCCESS;
		stat  = EG_vertexCollapse(qm,  vID);
		if ( stat != EGADS_SUCCESS) return stat;
		if ( TOTCOLLAPSES > count) {
			printf(" \n\n ======  EG_vertexOperator   COLLAPSED VERTEX COUNT %d: TRY SWAPPING WITH NEW MESH \n", TOTCOLLAPSES);
			*activity = 1;
		}
		break;
	case TRYSPLITTING:
		count = TOTSPLITS;
		//	if ( qm -> vType[vID - 1] != -1 )return EGADS_SUCCESS;
		stat = EG_vertexSplit(qm,  vID);
		if ( stat != EGADS_SUCCESS) return stat;
		if ( TOTSPLITS > count) {
			printf(" \n\n ======  EG_vertexOperator   SPLIT VERTEX COUNT %d: TRY SWAPPING WITH NEW MESH \n", TOTSPLITS);
			*activity = 1;
		}
		break;
	case TRYSWAPANDSPLIT:
		count = TOTSPLITS;
		stat  = EG_swapAndSplit(qm, vID);
		printf(" STAT IN SWAP AND PLIT %d \n",stat);
		if ( TOTSPLITS > count) {
			printf(" \n\n ======  EG_vertexOperator   SWAP + SPLIT VERTEX COUNT %d: TRY SWAPPING WITH NEW MESH \n", TOTSPLITS);
			*activity = 1;
		}
		break;
	case TRYDOUBLESPLIT:
		count = TOTSPLITS;
		stat  = EG_doubleSplit(qm, vID);
		printf(" STAT IN SWAP AND PLIT %d \n",stat);
		if ( TOTSPLITS > count) {
			printf(" \n\n ======  EG_vertexOperator   SWAP + SPLIT VERTEX COUNT %d: TRY SWAPPING WITH NEW MESH \n", TOTSPLITS);
			*activity = 1;
		}
		break;
	case TRYCOLLAPSESPLIT:
		count = TOTSPLITS;
		if ( vID <= qm -> nQ)
			stat  = EG_collapseAndSplit(qm, vID);
		printf(" STAT IN SWAP AND PLIT %d \n",stat);
		if ( TOTSPLITS > count) {
			printf(" \n\n ======  EG_vertexOperator   SWAP + SPLIT VERTEX COUNT %d: TRY SWAPPING WITH NEW MESH \n", TOTSPLITS);
			*activity = 1;
		}
		break;
	}
	//printf(" LEAVE VERTEX OPERATOR WITH ACVITIVY %d\n",*activity);
	return EGADS_SUCCESS;
}

static int
EG_meshOperator(quadMap *qm, int opType, int vertType, int *totActivity) {
	int i,  stat, activity = 0;
	*totActivity = 0;
	if ( opType == TRYCOLLAPSESPLIT ) {
		for (i = 0 ; i < qm -> nQ; ++i) {
			stat = EG_vertexOperator(qm, i + 1, opType, &activity) ;
			*totActivity += activity;
			if ( stat != EGADS_SUCCESS) return stat;
			if ( (opType != TRYSWAPPING) && activity ) {
				opType      = TRYSWAPPING;
				printf(" I HAVE COLLAPSED/SPLIT SOMETHING. RESET AND LOOK FOR SWAPS\n");
				break;
			}
		}
	}
	if ( opType == TRYCOLLAPSESPLIT) return stat;
	for (i = 0 ; i < qm -> totVerts; ++i) {
		if ( (qm -> vType[i] > 0 && vertType > 0) ||  qm -> vType[i] == vertType  ) {
			stat = EG_vertexOperator(qm, i + 1, opType, &activity) ;
			*totActivity += activity;
			if ( stat != EGADS_SUCCESS) break;
			if ( (opType != TRYSWAPPING) && activity ) {
				opType      = TRYSWAPPING;
				i           = -1;
				printf(" I HAVE COLLAPSED/SPLIT SOMETHING. RESET AND LOOK FOR SWAPS\n");
			}
		}
	}
	return stat;
}

static int EG_fullMeshRegularization(quadMap *qm, int ncall, int *meshImprovement, int globalOpti)
{
	int     stat, auxID, activity, aux, totActivity,i, j, k, nOpt, *optVec = NULL;
	double  uvbox[4], *lowerBounds = NULL, *upperBounds = NULL, fOpti, *uvOpt=NULL;
	static int vtype[3] = {0, 1, -1}; // corner, edge, interior
	FILE   *fop;
	// GET RANGE FOR EACH POINT
	printMeshStats(qm, ncall + 1);
	totActivity = 0;
	*meshImprovement = 0 ;
	printf(" PERFORM FULL MESH REGULARIZATION \n");
	printf(" START SWAPPING \n");
	for ( i = 0; i < 3; ++i) {
		printf(" SWAPPING      VERTEX TYPE %d \n",vtype[i]);
		do
		{
			stat              = EG_meshOperator(qm, TRYSWAPPING, vtype[i], &activity);
			if ( stat        != EGADS_SUCCESS) return stat;
			*meshImprovement += activity;
		}
		while (activity > 0);
	}
	printf(" DONE SWAPPING ROUND    TOT IMPROVEMENT %d \n", *meshImprovement);
	printf(" SPLIT + COLLAPSE VERTEX \n");
	do {
		stat              = EG_meshOperator (qm, TRYSPLITTING, -1, &activity);
		if ( stat        != EGADS_SUCCESS) return stat;
		*meshImprovement += activity;
	} while (activity > 0);
	do {
		stat              = EG_meshOperator (qm, TRYCOLLAPSESPLIT , -1, &activity);  //embedded swap + split
		if ( stat        != EGADS_SUCCESS) return stat;
		*meshImprovement += activity;
	} while (activity > 0);
	do {
		stat              = EG_meshOperator (qm, TRYDOUBLESPLIT,  -1, &activity);
		if ( stat        != EGADS_SUCCESS) return stat;
		*meshImprovement += activity;
	} while (activity > 0);
	do {
		stat              = EG_meshOperator (qm, TRYCOLLAPSING,  -1, &activity);
		if ( stat        != EGADS_SUCCESS) return stat;
		*meshImprovement += activity;
	} while (activity > 0);
	// SPLIT AND SWAP
	printf(" DONE =>  TOT IMPROVEMENT %d \n", *meshImprovement);
	snprintf (buffer, sizeof(char) * 32, "TOTAL_OPERATIONS_CALL%i.txt", ncall + 1);
	printf("BUFFER   %s\n",buffer);
	fop = fopen(buffer,"w");
	fprintf(fop," ---------------------------------------------\n");
	fprintf(fop," ROUND %d   TOTAL SWAPS IN   %d, COLLAPSES %d  SPLITS %d STAT = %d\n", k, TOTSWAPS, TOTCOLLAPSES,TOTSPLITS, stat);
	fprintf(fop," ---------------------------------------------\n");
	fclose(fop);
	printf(" ACTIVITY %d  TOTAL SWAPS IN THIS ROUND   %d, COLLAPSES %d  SPLITS %d STAT = %d\n", totActivity, TOTSWAPS, TOTCOLLAPSES,TOTSPLITS, stat);
	OVERALLSWAPS     += TOTSWAPS;
	TOTSWAPS          = 0 ;
	OVERALLCOLLAPSES += TOTCOLLAPSES;
	TOTCOLLAPSES      = 0 ;
	OVERALLSPLITS    += TOTSPLITS;
	TOTSPLITS         = 0 ;
	printf("-----------------------------------------------\n");
	printf(" ROUND %d OVERALL SWAPS  %d, COLLAPSES %d  SPLITS %d\n", k + 1, OVERALLSWAPS, OVERALLCOLLAPSES, OVERALLSPLITS);
	printf("-----------------------------------------------\n");
	if (globalOpti) {
		printf(" PERFORM GLOBAL OPTIMIZATION\n");
		for ( nOpt = i = 0 ; i < qm -> totVerts; ++i)
			if ( qm -> vType[i] == -1)  nOpt++;
		optVec = (int*) EG_alloc(nOpt*sizeof(int));
		if ( optVec == NULL )
			return EGADS_MALLOC;
		printf(" TOTAL VERTICES TO OPTIMIZE %d\n", nOpt);
		for ( nOpt = i = 0 ; i < qm -> totVerts; ++i)
			if ( qm -> vType[i] == -1) optVec[nOpt++] = i+1;
		printf(" TOTAL VERTICES TO OPTIMIZE %d\n", nOpt);
		print_mesh(qm, "BEFOREGLOBAL");
		stat = optimize_angles(qm, nOpt, optVec);
		print_mesh(qm, "AFTERGLOBAL");
		printf(" STAT FROM GLOBAL OPTIMIZER %d \n",stat);
	}
	printf("  MESH IMPROVEMENT IN THIS ROUND  %d\n",*meshImprovement);
	return stat;
}



static int
getCurvatureAngle(quadMap *qm, int sID, const double *uvs, double *err ) {
	int    stat, i, iMod, iMax, d, dMod, vID, vs[4];
	double dot[2*4], dMax, aux,  curv[8], proj[4*4], normal[3], theta, hOpt, k1, k2, kTheta, hTheta;
	vID  = qm -> star[sID] -> verts[0] - 1;
	if ( qm -> valence[vID][0] != 4) {
		*err = 100.0;
		return EGADS_INDEXERR;
	}
	stat = EG_curvature(qm -> face, uvs, curv);
	EG_unitVector(&curv[1], &aux);
	EG_unitVector(&curv[5], &aux);
	// FORCE CARTESIAN AXIS
	if ( PLANAR = 1) {
		curv[1] = 1.0; curv[2] = 0.0; curv[3] = 0.0;
		curv[5] = 0;   curv[6] = 1.0; curv[7] = 0.0;
	}
	//printf("\n PPAL DIRECTIONS %lf  %lf  %lf  %lf  %lf  %lf \n", curv[1], curv[2], curv[3], curv[5], curv[6], curv[7]);
	if ( stat != EGADS_SUCCESS) return stat;
	cross_product(&curv[1], &curv[5], normal);
	iMax = 0 ; dMax = 0.0;
	for ( i = 0 ; i < 4; i++) {
		vs[i]      = qm -> star[sID] -> verts[2*i + 1] - 1;
		//printf(" VERTEX %d  %d  TYPE %d \n", qm -> star[sID] -> verts[0], qm -> star[sID] -> verts[2*i + 1], qm -> vType[vs[i]]);
		stat       = EG_projectToTangentPlane(normal, &qm -> xyzs[3*vID], &qm -> xyzs[ 3*vs[i] ], &proj[4*i] );
		proj[4*i] -= qm -> xyzs[ 3*vID ]; proj[4*i + 1] -= qm -> xyzs[ 3*vID + 1]; proj[4*i + 2] -= qm -> xyzs[ 3*vID + 2];
		EG_unitVector(&proj[4*i], &proj[4*i + 3]);
		dot[2*i    ] = proj[4*i] * curv[1] + proj[4*i + 1] * curv[2] + proj[4*i + 2] * curv[3];
		dot[2*i + 1] = proj[4*i] * curv[5] + proj[4*i + 1] * curv[6] + proj[4*i + 2] * curv[7];
		//printf(" DOT %lf  %lf  NORM %lf \n",dot[2*i], dot[2*i + 1], proj[4*i + 3]);
		if ( fabs(dot[2*i]) > fabs(dMax) ) {
			iMax = i; d = 0; dMax = dot[2*i];
		}
		if ( fabs(dot[2*i + 1]) > fabs(dMax) ) {
			iMax = i; d = 1; dMax = dot[2*i +1];
		}
	}
	for ( i = 0 ; i < 4; i++) {
		if ( qm -> vType [vs[i]] >0  && qm -> valence[vs[i]][0] == 3 ) {
			iMax = i; d = 0; dMax = dot[2*i];
			break;
		}
	}
	*err = 0.0;
	if      ( fabs(curv[0] ) < KMIN ) curv[0] = KMIN;
	else if ( fabs(curv[0] ) > KMAX ) curv[0] = KMAX;
	if      ( fabs(curv[4] ) < KMIN ) curv[4] = KMIN;
	else if ( fabs(curv[4] ) > KMAX ) curv[4] = KMAX;
	for ( i = 0 ; i < 4; ++i)
	{
		iMod                   = (iMax + i ) % 4;
		// if ( qm -> vType [vs[iMod]] != -1 ) continue;
		if ( i % 2 == 0 ) dMod =  d;
		else              dMod = (d + 1)     % 2;
		aux                    = fabs(dot[2*iMod + dMod]);
		if   ( fabs(aux - 1.0) < DEPS ) theta = 0.0;
		else                            theta = fabs(acos(aux));
		// printf(" VERTEX %d = %d  dots % lf %lf  ANGLE %lf \n", iMod, vs[iMod], dot[2*iMod], dot[2*iMod + 1], theta);
		// NORMAL CURVATURE
		k1     = curv[ 4*   dMod           ];
		k2     = curv[ 4*( (dMod + 1 ) %2) ];
		hOpt   = GEPS / k1;
		kTheta = k1 * cos(theta) * cos(theta) + k2 * sin(theta) * sin(theta);
		hTheta = GEPS/kTheta;
		//printf(" HOPT %lf  HTHETA %lf  PROJ %lf error SIZE %lf ANGLE %lf \n",hOpt, hTheta, proj[4*iMod + 3], fabs(hOpt - proj[4*iMod +3]), theta);
		//*err  += (hOpt - hTheta ) * (hOpt - hTheta );
		*err  += (proj[4*iMod + 3] - hTheta ) * (hTheta - proj[4*iMod + 3]) + theta*theta;
	}
	//printf(" ACC ERROR %lf\n", *err);
	return EGADS_SUCCESS;
}

int CURVATURE = 0 ;
double NLOPT_objEqualAnglesTess(unsigned n, const double *uvs, double *grad, void *inQM)
{
	int    stat, i, j, k, vFix, vA, vB, aux, oID, peaks;
	double thetaOpt, thetas[4], dotErr, accPenalty, anglAcc;
	double vXYZ[18], oldXYZ[18], fQuad, fStar, fTot,  penalty, quadArea, curv[8];
	quadMap *qm = (quadMap*)inQM;
	if (OPTCOUNT == 0 ) SCOUNT = 0 ;
	++OPTCOUNT;
	for ( i = 0 ; i < qm->nS; ++i) {
		oID  = qm -> star[i] -> verts[0] -1 ;
		stat = EG_evaluate( qm->face, &uvs[2*i], vXYZ);
		if ( stat != EGADS_SUCCESS) return 1000000;
		qm-> xyzs[3*oID    ] = vXYZ[0];
		qm-> xyzs[3*oID + 1] = vXYZ[1];
		qm-> xyzs[3*oID + 2] = vXYZ[2];
		qm-> uvs [2*oID]     = uvs[2*i];
		qm-> uvs [2*oID + 1] = uvs[2*i + 1];
	}
	IOSTATUS = 0;
	if (OPTCOUNT%NOUT == 0 ) {
		IOSTATUS = 1;
		printf("#   CALL UPDATE OPTIMIZATION ITERATION %d\n",OPTCOUNT);
	}
	//  fprintf(filOpti,"#   CALL UPDATE OPTIMIZATION ITERATION %d\n",OPTCOUNT);
	fTot = 0.0;
	for ( i = 0 ; i < qm->nS; ++i) {
		//if ( CONVC == 30 ) print_starFile(qm ->star[i],qm);
		oID      = qm -> star[i] -> verts[0] -1;
		peaks    = qm -> valence[oID][0];
		thetaOpt = 2.0*PI/(double)(peaks);
		stat     = EG_evaluate( qm->face, &uvs[2*i], vXYZ);
		if ( stat != EGADS_SUCCESS) {
			printf(" NLOPT FAILED TO DO EG_Evaluate. EXIT PROGRAM \n");
			exit (1);
		}
		fStar      = 0.0;
		accPenalty = 0.0;
		anglAcc    = 0.0;
		for (j = 0 ; j < qm->star[i] -> nQ; ++j) {
			stat        = computePenalty(qm, i, j, &penalty, thetas, &quadArea);
			accPenalty += penalty;
			if ( stat != EGADS_SUCCESS)
				printf(" NLOPT FAILED TO DO computePenalty. EXIT PROGRAM \n");
			if(isnan(thetas[0])) {
				printf(" ROSE AT CENTRE %d  HAS A NAN VALUE FOR ANGLE AT QUAD %d theta %lf\n",qm ->star[i] -> verts[0], qm ->star[i]->quads[j], thetas[0]);
				exit(1);
			}
			//   fprintf(filOpti,"%d %d \t %lf %lf %lf\n", qm -> star[i]->verts[0] , qm -> star[i]->quads[j], penalty, quadArea, thetas[0]);
		}
		if ( peaks == 4) {
			stat   = getCurvatureAngle(qm, i, &uvs[2*i], &fStar);
			if(filOpti != NULL )
				fprintf(filOpti,"%d \t %lf %lf \n", qm -> star[i]->verts[0] , fStar,accPenalty);
		}
		if ( OPTCOUNT%NOUTMESH == 0 && OPTCOUNT != 0)
			printf("\n\n ------- ACCUMULATED PENALTY %lf ---- STAR ERROR %lf \n\n",accPenalty, fStar);
		fStar += accPenalty;
		if ( fStar < OPTITOL ) fStar  = 0.0;
		fTot += fStar;
		if(filOpti != NULL )
			fprintf(filOpti,"# TOTAL ITERATION FUNCTION %lf  %lf\n",fTot,fStar);
	}
	if ( OPTCOUNT%NOUTMESH == 0 && OPTCOUNT != 0){
		snprintf(buffer,33, "MESH_%i_ITERATION_%i.txt", CONVC, OPTCOUNT);
		print_mesh(qm, buffer);
	}
	return fTot;
}


static int optimize_angles(quadMap *qm, int nV, int *vID)
{
	int    regularStar, i, j, k, stat, iper, peaks, aux, area, stat2;
	double XYZ[18], uvbox[4], *lowerBounds = NULL, *upperBounds = NULL, fOpti, *uvOpt=NULL;
	double epsREL = 1.e-08, *epsABS = NULL, angles[4], vec, quadArea;
	double vals[18], normal[4];
	clock_t t_start, t_end;
	qm -> nS    = nV;
	//qm -> star  = EG_reall ( qm -> star, nV*sizeof(vStar*));
	uvOpt       = (double*) EG_alloc(2*nV*sizeof(double));
	lowerBounds = (double*) EG_alloc(2*nV*sizeof(double));
	upperBounds = (double*) EG_alloc(2*nV*sizeof(double));
	stat = EG_getRange(qm->face, uvbox, &iper);
	if (lowerBounds == NULL || upperBounds == NULL || uvOpt == NULL || qm ->star == NULL)  return EGADS_MALLOC;
	//#ifdef DEBUG
	for ( i = 0; i < nV; ++i)
		printf(" OPTIMIZE ANGLES %d  %d = %d \n",i,vID[i], qm -> vType[vID[i] -1]);
	//#endif
	IOSTATUS = 0;
	for ( i = 0; i < nV; ++i){
		if ( qm -> vType[vID[i] - 1] == -2 ){
			printf(" YOU ARE TRYING TO OPTIMIZE A DELETED VERTEX!!??\n");
			return EGADS_INDEXERR;
		}
		stat = EG_buildStar(qm, &qm ->star[i], vID[i]);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		uvOpt      [2*i    ] = qm -> uvs[2*(vID[i] -1)    ];
		uvOpt      [2*i + 1] = qm -> uvs[2*(vID[i] -1) + 1];
		lowerBounds[2*i    ] = uvbox[0];
		lowerBounds[2*i + 1] = uvbox[2];
		upperBounds[2*i    ] = uvbox[1];
		upperBounds[2*i + 1] = uvbox[3];
		for ( j = 0 ; j < qm -> star[i] -> nQ; ++j) {
			if ( qm -> star[i] -> quads[j] == -1 )
				printf(" I SHOULDNT BE OPTIMIZING A BOUNDARY VERTEX. STAR CENTRE %d \n", qm -> star[i] -> verts[0]);
			else
				printf(" CHECK INVALID FOR %d\n",qm -> star[i] -> quads[j]);
			stat = checkInvalidElement(qm, qm -> star[i]-> quads[j]);
			if ( stat != EGADS_SUCCESS) {
				//	printf("------------------------------------------ QUAD %d STAT %d\n", qm -> star[i] -> quads[j], stat);
				for ( j = 0 ; j < i; ++j ) {
					EG_free(qm -> star[j] -> verts);
					EG_free(qm -> star[j] -> quads);
					EG_free(qm -> star[j]);
				}
				EG_free(lowerBounds);
				EG_free(upperBounds);
				EG_free(uvOpt);
				return stat;
			}
		}
		//printf("\n");
	}
	OPTCOUNT = 0;
	// constraints for optimizer, u,v in face
	if (stat != EGADS_SUCCESS) {
		printf(" EG_getRange in optimize_angles = %d\n",stat);
		goto cleanup;
	}
	regularStar = 1;
	for (i = 0 ; i < nV; ++i ) {
		for ( j = 0 ; j < qm -> star[i] -> nQ; ++j) {
			aux = qm -> star[i] -> verts[2*j + 1] - 1;
			if (qm -> valence[aux][0] != 4 )regularStar = 0;
		}
	}
	if ( regularStar) {
		printf(" FOUND WHERE TO LOOK FOR OPTIMZER. EXIT PROGRAM AFTER opti\n");
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
	nlopt_set_stopval(opt,  OPTITOL);
	nlopt_set_maxeval(opt, NLOPTMAXEVAL);
	//nlopt_set_ftol_rel(opt, epsREL);
	nlopt_set_xtol_rel(opt, epsREL);
	//nlopt_set_ftol_abs(opt, epsABS);
	nlopt_set_xtol_abs(opt, epsABS);
	//nlopt_set_maxtime(opt, 120.0*2.0*(double)nV);
	//snprintf(buffer, sizeof(char) * 32, "OPTIMIZERMESHE_%d.txt",CONVC);
	//print_mesh(qm, buffer);
	//snprintf(optiName, sizeof(char) * 32, "OPTIMIZERCONVERGENCE_%d.txt",CONVC);
	OPTCOUNT = 0 ;
	/*printf(" OPEN FILE %s\n",optiName);
	 filOpti = fopen (optiName,"w");
	if(filOpti == NULL ) {
		printf(" I COULDNT OPEN FIL \n");
	}*/
	t_start = clock();
	CURVATURE = 1 ;
	//snprintf(buffer,33, "BEF_%d.txt",CONVC);
	//print_mesh(qm,buffer);
	// SET AXIS FOR DRIVING CURVATURE
	qm -> setIndex = 1;
	fOpti = NLOPT_objEqualAnglesTess(2*nV, uvOpt, NULL, qm);
	qm -> setIndex = 0;
	nlopt_result res = nlopt_optimize(opt, uvOpt, &fOpti);
	//snprintf(buffer,33, "AFT_%d.txt",CONVC);
	//print_mesh(qm,buffer);
	// if ( regularStar) exit(1);
	t_end = clock();
	//printf(" *-*-*-*-*-*-*- TIME REQUIRED TO OPTIMIZE %d POINTS = %lf  \n",nV, (double)(t_end - t_start) / CLOCKS_PER_SEC);
	NLOPTtermination(res);
	printf(" TERMINATION WITH OBJECTIVE %lf  TOTAL ITERATIONS %d\t",fOpti, OPTCOUNT);
	if ( res < 0 ) {
		stat = EGADS_GEOMERR;
		//	fclose(filOpti);
		goto cleanup;
	}
	else {
		stat = EGADS_SUCCESS;
		OPTCOUNT = -1 ;
		fOpti = NLOPT_objEqualAnglesTess(2*nV, uvOpt, NULL, qm);
	}
	// fclose(filOpti);
	CONVC++;
	//if (CONVC >= 10 ) exit(1);
	nlopt_destroy(opt);
	for (j = 0; j < nV; ++j) {
		qm ->uvs[2*(vID[j] -1 )    ] = uvOpt[2*j    ];
		qm ->uvs[2*(vID[j] -1 ) + 1] = uvOpt[2*j + 1];
		EG_evaluate(qm ->face, &qm ->uvs[2*(vID[j]-1)], XYZ);
		qm -> xyzs[3*(vID[j] -1 )   ] = XYZ[0];
		qm -> xyzs[3*(vID[j] -1 ) +1] = XYZ[1];
		qm -> xyzs[3*(vID[j] -1 ) +2] = XYZ[2];
		for ( i = 0 ; i < qm->star[j] -> nQ; ++i) {
			stat2 = quad_algebraic_area(qm, qm->star[j] ->quads[i], qm->star[j] ->verts[0], &area, angles, &quadArea);
			if (stat2 != EGADS_SUCCESS || area != 1) {
				stat = EGADS_GEOMERR;
				goto cleanup;
				continue;
				printf("OPTI  STAT2  %d  AREA %d \n",stat, area);
				FILE *file;
				snprintf(buffer, sizeof(char) * 32, "INVALID_ELEMENT%i.txt", INVALIDCOUNT);
				++INVALIDCOUNT;
				int qID =  qm->star[j] ->quads[i];
				file = fopen(buffer,"w");
				if (file == NULL ) {
					stat = EGADS_MALLOC;
					goto cleanup;
				}
				for ( k = 0 ; k < 4; ++k) fprintf(file,"%lf %lf %lf\n", qm -> xyzs[3*(qm ->quadIdx[4*(qID - 1) + k] - 1)], qm -> xyzs[3*(qm ->quadIdx[4*(qID - 1) + k] - 1) +1],qm -> xyzs[3*(qm ->quadIdx[4*(qID - 1) + k] - 1) +2]);
				fprintf(file,"%lf %lf %lf\n", qm -> xyzs[3*(qm ->quadIdx[4*(qID - 1)] - 1)], qm -> xyzs[3*(qm ->quadIdx[4*(qID - 1)] - 1) +1],qm -> xyzs[3*(qm ->quadIdx[4*(qID - 1)] - 1) +2]);
				fprintf(file,"\n\n");
				fclose(file);

			}
		}
	}
	cleanup:
	EG_free(lowerBounds);
	EG_free(upperBounds);
	if ( qm -> star) {
		for ( i = 0 ; i < nV; ++i){
			if ( qm -> star[i]) {
				EG_free(qm->star[i]->verts);
				qm->star[i]->verts = NULL;
				EG_free(qm->star[i]->quads);
				qm->star[i]->quads = NULL;
				EG_free(qm->star[i]);
				qm->star[i] = NULL;
			}
		}
		//EG_free(qm->star);
		//qm -> star = NULL;
	}
	EG_free(uvOpt);
	EG_free(epsABS);
	//snprintf(buffer, sizeof(char) * 32, "FIOPTIM_%d.txt",CONVC);
	//print_mesh(qm, buffer);
	//printf("\n");
	//exit(1);
	return stat;

}








#ifdef STANDALONE
int main(int argc, char *argv[])
{
	int          lastRound, stat, vt, v, otype, totActivity, activity, verts, f , e, i, j, k, q, m, n, iBody, oclass, mtype, len, nquad, ntri, sum;
	int          nseg, ngp, atype, alen, quad, *segs, *senses, *tris, *heads, nledges, nh, nloops,
	*lsenses, *esenses, meshActivity, globalOpti, *optVec = NULL, nt = 3, nOpt;
	const int    *ints;
	float        arg, color[3];
	double       box[6], size, params[3], *normals, eval[18];
	const double *xyzs, *uvs, *ts, *reals ;
	double *modxyzs;
	char         gpname[33], *startapp;
	const char   *OCCrev, *string;
	ego          context, tess, model, geom, *bodies, *dum, *loops, *ledges;
	wvData       items[6];
	float        eye[3]      = {0.0, 0.0, 7.0};
	float        center[3]   = {0.0, 0.0, 0.0};
	float        up[3]       = {0.0, 1.0, 0.0};
	static int   sides[3][2] = {{1,2}, {2,0}, {0,1}};
	static int   sideq[4][2] = {{1,2}, {2,5}, {5,0}, {0,1}};
	static int   neigq[4]    = {  0,     3,     4,     2};
	static int   vTri[6]     = {0, 1, 2, 0, 2, 3};
	static int region[3]     = {0, 1, -1};
	float         *lsegs;

	/* get our starting application line
	 *
	 * for example on a Mac:
	 * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
	 */
	startapp = getenv("WV_START");

	if ((argc != 2) && (argc != 5) && (argc != 6)) {
		printf("\n Usage: vAttr filename [angle maxlen sag] [deb]\n\n");
		return 1;
	}

	/* look at EGADS revision */
	EG_revision(&i, &j, &OCCrev);
	printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

	/* initialize */
	stat = EG_open(&context);
	if ( stat != EGADS_SUCCESS) return 1;

	stat = EG_loadModel(context, 0, argv[1], &model);
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
	// default granularity
	params[0] =  0.1 *size;
	params[1] =  0.01*size;
	params[2] =  20.0;
	// user defined triangulation
	if (argc >= 5) {
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
	for (iBody = 0; iBody < nbody; iBody++) {
#ifdef FORCETRIANGULATION   // to defeat a quad triangulation
		stat = EG_attributeAdd(bodies[iBody], ".qParams",
				ATTRSTRING, 4, NULL, NULL, "off");
		if (stat != EGADS_SUCCESS)
			printf(" Body %d: attributeAdd = %d\n", iBody, stat);
#endif
		mtype = 0;
		EG_getTopology(bodies[iBody], &geom, &oclass,
				&mtype, NULL, &j, &dum, &senses);
		bodydata[iBody].body  = bodies[iBody];
		bodydata[iBody].mtype = mtype;
		if (mtype == WIREBODY) {
			printf(" Body %d: Type = WireBody\n", iBody+1);
		} else if (mtype == FACEBODY) {
			printf(" Body %d: Type = FaceBody\n", iBody+1);
		} else if (mtype == SHEETBODY) {
			printf(" Body %d: Type = SheetBody\n", iBody+1);
		} else {
			printf(" Body %d: Type = SolidBody\n", iBody+1);
		}

		stat = EG_getBodyTopos(bodies[iBody], NULL, FACE,
				&bodydata[iBody].nfaces, &bodydata[iBody].faces);
		i    = EG_getBodyTopos(bodies[iBody], NULL, EDGE,
				&bodydata[iBody].nedges, &bodydata[iBody].edges);
		if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
			printf(" EG_getBodyTopos Face = %d\n", stat);
			printf(" EG_getBodyTopos Edge = %d\n", i);
			return 1;
		}
		stat = EG_makeTessBody(bodies[iBody], params, &bodydata[iBody].tess);
		if (stat != EGADS_SUCCESS) {
			printf(" EG_makeTessBody %d = %d\n", iBody, stat);
			continue;
		}
		stat = EG_makeTessBody(bodies[iBody], params, &bodydata[iBody].tess);
#ifdef FORCEQUADS
		tess = bodydata[iBody].tess;
		stat = EG_quadTess(tess, &bodydata[iBody].tess);
		if (stat != EGADS_SUCCESS) {
			printf(" EG_quadTess %d = %d  -- reverting...\n", iBody, stat);
			bodydata[iBody].tess = tess;
			continue;
		}
		EG_deleteObject(tess);
		QUADTESS = 1;
#endif
	}
	printf(" \n");

	cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
	if (cntxt == NULL) {
		printf(" failed to create wvContext!\n");
		for (iBody = 0; iBody < nbody; iBody++) {
			EG_deleteObject(bodydata[iBody].tess);
			EG_free(bodydata[iBody].edges);
			EG_free(bodydata[iBody].faces);
		}
		free(bodydata);

		printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
		printf(" EG_close          = %d\n", EG_close(context));
		return 1;
	}
	for (ngp = sum = stat = iBody = 0; iBody < nbody; iBody++) {
		quad = 0;
		stat = EG_attributeRet(bodydata[iBody].tess, ".tessType", &atype,
				&alen, &ints, &reals, &string);
		if (stat == EGADS_SUCCESS)
			if (atype == ATTRSTRING)
				if (strcmp(string, "Quad") == 0) quad = 1;

		/* get faces */
		bodydata[iBody].qm = (quadMap**) EG_alloc(bodydata[iBody].nfaces*sizeof(quadMap*));
		if (bodydata[iBody].qm == NULL ) {
			stat = EGADS_MALLOC;
			goto cleanup;
		}
		stat = EG_createQuadMap(&bodydata[iBody]);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		for ( f = 0 ; f < bodydata[iBody].nfaces; ++f) {
			sprintf(gpname, "Body %d Face %d", iBody+1, f+1);  // Preserve this name
			// set vertices
			len    = bodydata[iBody].qm[f]-> totVerts;
			nquad  = bodydata[iBody].qm[f]-> nQ;
			ntri   = bodydata[iBody].qm[f]-> nQ * 2;
			stat   = wv_setData(WV_REAL64, len, (void *) bodydata[iBody].qm[f]-> xyzs,  WV_VERTICES, &items[0]);
			if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", f, gpname);
			wv_adjustVerts(&items[0], focus);
			tris = (int *) malloc(3*ntri*sizeof(int));
			for ( i = j = 0 ; j < nquad; ++j) {
				for ( e = 0 ; e < 6; ++e )
					tris[i++] = bodydata[iBody].qm[f]-> quadIdx[4*j + vTri[e]];
			}
			if ( i != 3*ntri) {
				printf(" WHAT??? i = %d != %d \n",i, 3*ntri);
			}
			stat = wv_setData(WV_INT32, 3*ntri, (void *) tris, WV_INDICES, &items[1]);
			free(tris);
			if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", f, gpname);
			color[0]  = 1.0;
			color[1]  = iBody;
			color[1] /= nbody;
			color[2]  = 0.0;
			stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[2]);
			if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", f, gpname);
			// For now, draw a line every two vertices even though it produces duplicates from adjacent quads
			segs = (int *) malloc(2*4*nquad *sizeof(int));
			if (segs == NULL) {
				printf(" Can not allocate %d Sides!\n", nseg);
				stat = EGADS_MALLOC;
				goto cleanup;
			}
			for ( nseg = j = 0 ; j < nquad; ++j) {
				if (bodydata[iBody].qm[f] -> quadIdx[4*j] > 0) {
					for( k = 0 ; k < 4; ++k ) {
						segs[2*nseg     ] = bodydata[iBody].qm[f] -> quadIdx[4*j +  k     ];
						segs[2*nseg  + 1] = bodydata[iBody].qm[f] -> quadIdx[4*j + (k+1)%4];
						nseg++;
					}
				}
			}
			stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_LINDICES, &items[3]);
			if (stat < 0) printf(" wv_setData = %d for %s/item 3!\n", f, gpname);
			free(segs);
			//    color[0] = color[1] = color[2] = 0.8;
			color[0] = color[1] = color[2] = 0.0;
			stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_LCOLOR, &items[4]);
			if (stat < 0) printf(" wv_setData = %d for %s/item 4!\n", f, gpname);
			normals = (double *) malloc(3*len*sizeof(double));
			for ( i = 0 ; i < len; ++i) {
				stat = EG_evaluate(bodydata[iBody].qm[f]->face, &bodydata[iBody].qm[f]->uvs[2*i], eval);
				cross_product(&eval[3], &eval[6], &normals[3*i]);
			}
			stat = wv_setData  (WV_REAL64, len, (void *) normals, WV_NORMALS, &items[5]);
			stat = wv_addGPrim (cntxt, gpname, WV_TRIANGLE,      WV_ON|WV_ORIENTATION, 6, items);
			printf(" wv_addGPrim = %d for %s \n", stat, gpname);
			if (stat < 0)
				printf(" !!!!!!!!!!!!\n");
			if (stat > 0) ngp = stat+1;
			sum += ntri;
			free(normals);
			// get SURFACE edges //
			/* get loops */
			color[0] = color[1] = 0.0;
			color[2] = 1.0;
			for (i = 0; i < bodydata[iBody].nfaces; i++) {
				stat = EG_getTopology(bodydata[iBody].faces[i], &geom, &oclass,
						&mtype, NULL, &nloops, &loops, &lsenses);
				if (stat != EGADS_SUCCESS) continue;
				for (nh = j = 0; j < nloops; j++) {
					stat = EG_getTopology(loops[j], &geom, &oclass, &mtype, NULL,
							&nledges, &ledges, &esenses);
					if (stat != EGADS_SUCCESS) continue;

					/* count */
					for (ntri = nseg = k = 0; k < nledges; k++) {
						m = 0;
						while (ledges[k] != bodydata[iBody].edges[m]) {
							m++;
							if (m == bodydata[iBody].nedges) break;
						}
						/* assume that the edge is degenerate and removed */
						if (m == bodydata[iBody].nedges) continue;
						stat = EG_getTessEdge(bodydata[iBody].tess, m+1, &len,
								&xyzs, &ts);
						if (stat != EGADS_SUCCESS) {
							printf(" EG_getTessEdge %d = %d!\n", m+1, stat);
							nseg = 0;
							break;
						}
						if (len == 2)
							if ((xyzs[0] == xyzs[3]) && (xyzs[1] == xyzs[4]) &&
									(xyzs[2] == xyzs[5])) continue;
						nh++;
						nseg += len;
						ntri += 2*(len-1);
					}
					if (nseg == 0) continue;
					lsegs = (float *) malloc(3*nseg*sizeof(float));
					if (lsegs == NULL) {
						printf(" Can not allocate %d Segments!\n", nseg);
						continue;
					}
					segs = (int *) malloc(ntri*sizeof(int));
					if (segs == NULL) {
						printf(" Can not allocate %d Line Segments!\n", ntri);
						free(lsegs);
						continue;
					}
					heads = (int *) malloc(nh*sizeof(int));
					if (heads == NULL) {
						printf(" Can not allocate %d Heads!\n", nh);
						free(segs);
						free(lsegs);
						continue;
					}

					/* fill */
					for (nh = ntri = nseg = k = 0; k < nledges; k++) {
						m = 0;
						while (ledges[k] != bodydata[iBody].edges[m]) {
							m++;
							if (m == bodydata[iBody].nedges) break;
						}
						/* assume that the edge is degenerate and removed */
						if (m == bodydata[iBody].nedges) continue;
						EG_getTessEdge(bodydata[iBody].tess, m+1, &len, &xyzs, &ts);
						if (len == 2)
							if ((xyzs[0] == xyzs[3]) && (xyzs[1] == xyzs[4]) &&
									(xyzs[2] == xyzs[5])) continue;
						if (esenses[k] == -1) heads[nh] = -ntri/2 - 1;
						for (n = 0; n < len-1; n++) {
							segs[ntri] = n+nseg+1;
							ntri++;
							segs[ntri] = n+nseg+2;
							ntri++;
						}
						if (esenses[k] ==  1) heads[nh] = ntri/2;
						for (n = 0; n < len; n++) {
							lsegs[3*nseg  ] = xyzs[3*n  ];
							lsegs[3*nseg+1] = xyzs[3*n+1];
							lsegs[3*nseg+2] = xyzs[3*n+2];
							nseg++;
						}
						nh++;
					}
					sprintf(gpname, "Body %d Loop %d/%d", iBody+1, i+1, j+1);
					stat = wv_setData(WV_REAL32, nseg, (void *) lsegs,  WV_VERTICES, &items[0]);
					if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
					wv_adjustVerts(&items[0], focus);
					free(lsegs);
					stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[1]);
					if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
					stat = wv_setData(WV_INT32, ntri, (void *) segs, WV_INDICES, &items[2]);
					if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
					free(segs);
					stat = wv_addGPrim(cntxt, gpname, WV_LINE, WV_ON, 3, items);
					if (stat < 0) {
						printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
					} else {
						if (cntxt != NULL)
							if (cntxt->gPrims != NULL) {
								cntxt->gPrims[stat].lWidth = 1.0;
								n = wv_addArrowHeads(cntxt, stat, 0.05, nh, heads);
								if (n != 0) printf(" wv_addArrowHeads = %d\n", n);
							}
					}
					free(heads);
					/*      wv_printGPrim(cntxt, stat);  */
				}
			}



		}
		printf(" ** %d gPrims with %d triangles **\n", ngp, sum);
		// start the server code //
		if (argv[5]) {  //debugging-- no server call
			for (iBody = 0; iBody < nbody; iBody++) {
				for (f = 0; f < bodydata[iBody].nfaces; ++f) {
					if (bodydata[iBody].qm[f]-> uvs == NULL) continue;
					print_mesh(bodydata[iBody].qm[f],"INITMESH.txt");
					int it = 0;
					do  {
						stat = EG_fullMeshRegularization(bodydata[iBody].qm[f], 0, &meshActivity, 0);
						printf("\n\n\n-------- MESH IMPROVEMENT %d\n",meshActivity);
						snprintf(buffer,33,"FULLMESH_%d.txt",it);
						print_mesh(bodydata[iBody].qm[f],buffer);
						++it;
					}
					while(meshActivity > 0);
					stat = EG_fullMeshRegularization(bodydata[iBody].qm[f], 0, &meshActivity, 1);
					print_mesh(bodydata[iBody].qm[f],"finalMESH.txt");
				}
			}
		} else {
			meshActivity = 1;
			globalOpti   = 0;
			nOpt = 0;
			otype = 0; vt = 0 ; v = 0;  //initialize
			stat = 0 ;
			int count = 0;
			lastRound = 0 ;
			int iter = 0;
			print_mesh(bodydata[0].qm[0], "INITMESH.txt");
			if (wv_startServer(7681, NULL, NULL, NULL, 0, cntxt) == 0) {
				while (wv_statusServer(0)) {
					usleep(100000);
					if (stat == 0) {
						if (startapp != NULL) system(startapp);
						stat++;
					}
					if ( meshActivity ) {
						int vwcount = 0;
						//printf(" ENTER: TRY SWAPPING 0 COLLAPSING 1 SPLITTING 2\n");
						//scanf("%d",&otype);
						for (iBody = 0; iBody < nbody; iBody++) {
							for (f = 0; f < bodydata[iBody].nfaces; ++f) {
								totActivity = 0;
								if ( otype != TRYSWAPPING ) vt = 2;
								if ( globalOpti ) {
									activity = 1 ;
									for ( nOpt = i = 0 ; i < bodydata[iBody].qm[f] -> totVerts; ++i) {
										if ( bodydata[iBody].qm[f] -> vType[i] == -1)  ++nOpt;
									}
									optVec = (int*)EG_alloc (nOpt*sizeof(int));
									if ( optVec == NULL ) {
										printf(" MALLOC PROBLEM IN OPTVEC\n");
										return 1;
									}
									printf(" NOPT %d \n",nOpt);
									for ( nOpt = i = 0 ; i < bodydata[iBody].qm[f] -> totVerts; ++i)
									{
										printf(" I %d = %d\n",i,bodydata[iBody].qm[f]->vType[i]);
										if ( bodydata[iBody].qm[f] -> vType[i] == -1)  {
											printf(" AT POS %d ADD %d \n",nOpt,i+1);
											optVec[nOpt++] = i + 1;
										}
									}
									printf(" GLOBAL OPTIMIZATION %d\n", nOpt);
									stat = optimize_angles(bodydata[iBody].qm[f], nOpt, optVec);
									printf(" STAT FROM GLOBAL OPTIMIZER %d \n",stat);
									printf("FINAL ROUND: GLOBAL OPTIMIZATION");
									meshActivity = 1;
									snprintf(buffer, 33,"OPTIFINALMESH");
									print_mesh(bodydata[0].qm[0], buffer);
									EG_free(optVec);
									otype = TRYSWAPPING;
								}
								printf("\n\n ++++++++++++++++++++++++++ OTYPE %d   VERTEX = %d TYPE %d ++++++++++++++++++++++++++\n", otype, v+1,region[vt]);
								stat  = EG_vertexOperator(bodydata[iBody].qm[f], v + 1, otype, &activity) ;
								printf(" ACTIVITY  HERE %d\n",activity);
								if ( activity ) {
									// Visualize effect
									printf(" WV VIEWER SWAP STATUS   %d  ACTIVITY %d\n", stat, activity);
									printf("  TOTAL SWAPS IN   %d, COLLAPSES %d  SPLITS %d STAT = %d\n",  TOTSWAPS, TOTCOLLAPSES,TOTSPLITS, stat);
									sprintf(gpname, "Body %d Face %d", iBody+1, f + 1);
									//print_mesh(bodydata[iBody].qm[f], &vwcount, 1);
									printf(" TOT VERTS %d  REM VERTS %d  TOT QUADS %d reM QUADS %d\n",bodydata[iBody].qm[f]-> totVerts, bodydata[iBody].qm[f]-> remVerts[0],
											bodydata[iBody].qm[f]-> nQ, bodydata[iBody].qm[f]-> remQuads[0]										 );
									verts   = bodydata[iBody].qm[f]-> totVerts;
									modxyzs = (double*) EG_alloc (3*verts*sizeof(double));
									normals = (double*) EG_alloc (3*verts*sizeof(double));
									if (normals == NULL ||  xyzs == NULL) return EGADS_MALLOC;
									for (  i = 0 ; i < verts; ++i) {
										modxyzs[3*i    ] = bodydata[iBody].qm[f]-> xyzs[3*i    ];
										modxyzs[3*i + 1] = bodydata[iBody].qm[f]-> xyzs[3*i + 1];
										modxyzs[3*i + 2] = bodydata[iBody].qm[f]-> xyzs[3*i + 2];
										stat = EG_evaluate(bodydata[iBody].qm[f]-> face, &bodydata[iBody].qm[f]->uvs[2*i], eval);
										if (stat!= EGADS_SUCCESS){
											printf(" EG_evaluate after modifying mesh is %d ", stat);
											break;
										}
										if ( bodydata[iBody].qm[f] -> face -> mtype == SREVERSE )
											cross_product(&eval[6], &eval[3], &normals[3*i]);
										else
											cross_product(&eval[3], &eval[6], &normals[3*i]);
									}
									stat = wv_setData(WV_REAL64, verts, (void *) modxyzs, WV_VERTICES, &items[0]);
									free(modxyzs);
									if (stat < 0) printf(" wv_setData = %d for %s!\n", iBody, gpname);
									wv_adjustVerts(&items[0], focus);
									printf(" TOT SURVIVING QUADS %d  REMOVED %d \n", bodydata[iBody].qm[f]->nQ, bodydata[iBody].qm[f]->remQuads[0]);
									nquad = bodydata[iBody].qm[f]->nQ;
									if (bodydata[iBody].qm[f]->remQuads[0] > 0 )
										nquad -= bodydata[iBody].qm[f]->remQuads[0];
									ntri = 2*nquad;
									tris = (int *) malloc(3*ntri*sizeof(int));
									for ( i = j = 0 ; j < bodydata[iBody].qm[f] -> nQ ; ++j) {
										if ( bodydata[iBody].qm[f]-> quadIdx[4*j] < 0 ) continue;
										for ( e = 0 ; e < 6; ++e )
											tris[i++] = bodydata[iBody].qm[f]-> quadIdx[4*j + vTri[e]];
									}
									if ( i != 3* ntri) {
										printf(" ???? i %d   ntri  %d  \n", i, 3* ntri);
										exit (1 );
									}

									stat = wv_setData(WV_INT32, 3*ntri, (void *) tris, WV_INDICES, &items[1]);
									free(tris);
									if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", f, gpname);
									color[0]  = 1.0;
									color[1]  = iBody+1 ;
									color[1] /= nbody;
									color[2]  = 0.0;
									stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[2]);
									if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", f, gpname);
									// Get normals (must be the same that vertices and we have removed some
									segs  = (int *) malloc(2*4*nquad *sizeof(int));
									if (segs == NULL) {
										printf(" Can not allocate %d Sides!\n", nseg);
										return EGADS_MALLOC;
									}
									for ( nseg = j = 0 ; j < bodydata[iBody].qm[f] -> nQ; ++j) {
										if (bodydata[iBody].qm[f] -> quadIdx[4*j] > 0) {
											for( k = 0 ; k < 4; ++k ) {
												segs[2*nseg     ] = bodydata[iBody].qm[f] -> quadIdx[4*j + k          ];
												segs[2*nseg  + 1] = bodydata[iBody].qm[f] -> quadIdx[4*j + (k + 1 ) %4];
												nseg++;
											}
										}
									}
									stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_LINDICES, &items[3]);
									if (stat < 0) printf(" wv_setData = %d for %s/item 3!\n", f, gpname);
									EG_free(segs);
									color[0] = color[1] = color[2] = 0.0;
									stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_LCOLOR, &items[4]);
									if (stat < 0) printf(" wv_setData = %d for %s/item 4!\n", f, gpname);
									stat = wv_setData(WV_REAL64, verts, (void *) normals,  WV_NORMALS, &items[5]);
									if (stat < 0) printf(" wv_setData = %d for %s/item 5!\n", f, gpname);
									free(normals);
									//
									int index = wv_indexGPrim(cntxt, gpname);
									stat = wv_modGPrim(cntxt, index, 6, items);
									if (stat < 0)
										printf(" wv_modGPrim = %d for %s (%d)!\n", stat, gpname, index);
									if (globalOpti) meshActivity = 0;
									//snprintf(buffer,33,"MESH%d",count);
									//printf(" BUFFER %s \n",buffer);
									//print_mesh(bodydata[iBody].qm[f], buffer);
									//++count;
								}
								++v;
								if ( v ==  bodydata[iBody].qm[f] -> totVerts && lastRound == 0) {
									v = 0;
									vt++;  //change domain region: corner, edge or interior
									meshActivity = 1;
									printf(" RESET V VT %d \n",vt);
									if ( vt == 3 ) {
										vt = 0 ;
										printf(" TRIED EVERYTHING-> CHANGE OPERATION  %d =>",otype);
										if ( otype != TRYSWAPPING && activity )
											otype = TRYSWAPPING;
										else ++otype;
										printf(" TRIED EVERYTHING-> CHANGE OPERATION  %d =>",otype);
										if ( otype >= OPTYPES) {
											otype      = TRYSWAPPING;
											printf("FINAL ROUND: GLOBAL OPTIMIZATION");
											meshActivity = 1;
											snprintf(buffer,33,"FINALMESH_BEFOREOPTI");
											print_mesh(bodydata[0].qm[0], buffer);
											++iter;
											if ( iter > 500 ) lastRound = 1;
										}
									}
								}
								else if ( v ==  bodydata[iBody].qm[f] -> totVerts && lastRound == 1) globalOpti = 1;
							}
						}
					}
				}
			}
		}
	}
	wv_cleanupServers();
	/* finish up */
	cleanup:
	for (iBody = 0; iBody < nbody; iBody++) {
		EG_deleteObject(bodydata[iBody].tess);
		EG_free(bodydata[iBody].edges);
		EG_free(bodydata[iBody].faces);
		EG_destroyQuadMap(bodydata[iBody].qm, bodydata[iBody].nfaces);
	}
	free(bodydata);
	printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
	printf(" EG_close          = %d\n", EG_close(context));
	return 0;
}
#endif








/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *wsi, char *text, /*@unused@*/ int lena)
{
	if ((strcmp(text,"finer") == 0)) {
		printf(" DO STUFF ON MESH\n");

	}
}

/********************   IO FUNCTIONS **********************/


void sampleNormalPlane(double normal[], double point[], int vID, quadMap *qm) {
	double min[3], max[3], c, p[3], dt[3], r = 1.0;
	int  i, j, k, nP;
	FILE *f;
	char name[32];
	snprintf(name, sizeof(char) * 32, "PLANE%i.txt", vID);
	f = fopen(name,"w");
	if ( f == NULL) return ;
	nP = 20;
	c = 0.0;
	vID--;
	for (i = 0 ; i < 3; ++i) {
		min  [i] = qm -> xyzs[3*vID + i] -1; //qm -> xyzs[3*vID + i] - 100 *qm -> xyzs[3*vID + i] ;
		max  [i] = qm -> xyzs[3*vID + i]  +1;//0.5;// qm -> xyzs[3*vID + i] + 100 *qm -> xyzs[3*vID + i] ;
		if ( min[i] > max[i]) {
			dt[0] = min[i];
			min[i] = max[i];
			max[i] = dt[0];
		}
		c       += normal[i] * point[i];
	}
	for (i = 0 ; i < 3; ++i)
		dt[i] = (max[i] - min[i] )/(double)(nP -1);
	fprintf(f,"%lf %lf %lf\n",point[0], point[1], point[2]);
	for (i = 0 ; i < nP; ++i)
	{
		for ( j = 0 ; j < nP; ++j)
		{
			for ( k = 0 ; k < nP; ++k)
			{
				p[0] = min[0] + (double)i*dt[0];
				p[1] = min[1] + (double)j*dt[1];
				double aux = c - p[0] * normal[0] - p[1] * normal[1];
				if ( fabs(normal[2]) > EPS) p[2] = aux / normal[2];
				else {
					p[2] = min[2] + (double)k*dt[2];
					if (fabs(normal[1]) < EPS) {
						p[1] = min[1] + (double)j*dt[1];
						p[0] = point[0];
					}
					else {
						if (fabs(normal[0]) < EPS) p[1] = point[1];
						else {
							aux  = c - p[0] * normal[0];
							p[1] = aux / normal[1];
						}
					}
				}
				fprintf(f,"%lf %lf %lf\n",p[0], p[1], p[2]);
			}
		}
	}
	fclose(f);
}

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

static void
print_quad_specs(quadMap *qm, int id) {
	--id;
	if (id < 0 ) return ;
	int i, v;
	printf(" QUAD %d HAS VERTICES ",id +1);
	for ( i = 0 ; i < 4; ++i ) {
		v = qm -> quadIdx[4*id + i];
		if ( v > qm -> totVerts) {
			printf(" OUT OF BOUNDS IN QUAD SPECS\n");
			exit(1);
		}
		printf(" %d ",v);
	}
	printf("\t AND ADJACENT QUADS ");
	for ( i = 0 ; i < 4; ++i ) printf(" %d ",qm -> quadAdj[4*id + i]);
	printf("\n");
}
static void
print_star(vStar *star) {
	int i, j;
	printf(" ===== STAR CENTRED AT VERTEX %d \n", star->verts[0]);
	printf(" ROSA  DE %d PICOS \n",star->nV);
	for ( i = 1 ; i < star->nV; ++ i) {
		if ( i%2 == 0)
			printf(" OPP  PICO %d\n",star->verts[i]);
		else
			printf(" LINK  PICO %d\n",star->verts[i]);
	}
	for ( i = 0 ; i < star ->nQ; ++i) printf(" QUAD %d = %d\n",i,star->quads[i]);
}
static void
print_starFile(vStar *star, quadMap *qm) {
	int i, vID, k;
	FILE *fout;
	char buffer[33];
	//snprintf(buffer, sizeof(char) * 32, "STAR%i_CENTRE_%i.txt", CONVC, star->verts[0]);
	snprintf(buffer, sizeof(char) * 32, "STAR_%i.txt", star->verts[0]);
	printf(" ------WRITING ON FILE %s\n",buffer);
	fout = fopen(buffer,"w");
	if (fout == NULL ) return;
	for ( i = 0 ; i < star->nQ; ++i) {
		fprintf(fout, "%lf  %lf  %lf %d\n",qm->xyzs[3*(star->verts[0] - 1)],qm->xyzs[3*(star->verts[0] - 1) + 1],qm->xyzs[3*(star->verts[0] - 1) + 2],star->verts[0]);
		for ( k = 0 ; k < 3; ++k){
			vID =  2*i + k + 1;
			if ( vID >= star->nV) vID = 1;
			vID = star->verts[vID] - 1;
			fprintf(fout, "%lf  %lf  %lf %d\n", qm->xyzs[3*vID],qm->xyzs[3*vID + 1],qm->xyzs[3*vID + 2],vID + 1);
		}
		fprintf(fout, "%lf  %lf  %lf %d \n",qm->xyzs[3*(star->verts[0] - 1)],qm->xyzs[3*(star->verts[0] - 1) + 1],qm->xyzs[3*(star->verts[0] - 1) + 2],star->verts[0]);
		fprintf(fout,"\n\n");
	}
	for ( i = 0 ; i < star->nV; ++i) {
		vID = star->verts[i];
		printf(" STARFILE VID %d i %d TOT NV %d\n",vID, i, star->nV);
		if (vID <= 0 || vID > qm ->totVerts) exit(1);
		fprintf(fout, "# VERT %d = %lf  %lf  %lf\n",star->verts[i], qm ->xyzs[3*(vID - 1)], qm ->xyzs[3*(vID - 1) +1], qm ->xyzs[3*(vID - 1) + 2]);
	}
	fclose(fout);
	snprintf(buffer, sizeof(char) * 32, "CENTRE_%i.txt",  star->verts[0]);
	//snprintf(buffer, sizeof(char) * 32, "CENTRE%i_SCOUNT_%i.txt", CONVC, star->verts[0]);
	fout = fopen(buffer,"w");
	if (fout == NULL ) return;
	vID = star->verts[0] - 1;
	fprintf(fout, "%lf  %lf  %lf %d\n", qm->xyzs[3*vID],qm->xyzs[3*vID + 1],qm->xyzs[3*vID + 2],vID + 1);
	fclose(fout);
}








static void
printMeshStats(quadMap *qm, int sweep) {
	int fix, move, i,len, val ;
	int intVal[100], boundVal[100];
	FILE *fout;
	snprintf(buffer, sizeof(char) * 32, "MESH_STATS_%d.txt", sweep);
	printf(" WRITING ON FILE %s\n",buffer);
	fout = fopen(buffer, "w");
	len = qm ->totVerts;
	for ( i = 0; i < 100; ++i) {
		intVal  [i] = 0;
		boundVal[i] = 0;
	}
	for ( i = 0; i < len; ++i) {
		if ( qm -> vType[i] !=  -2){
			if ( qm -> vType[i] == -1){
				val = qm -> valence[i][0];
				++intVal[val];
			} else {
				++boundVal[qm ->valence[i][0]];
			}
		}
	}
	fprintf(fout,"--------------------- MESH STATS AFTER %d SWEEP --------------\n",sweep);
	fprintf(fout,"---- TOTAL VERTICES %d TOTAL QUADS %d --------------\n",qm -> totVerts - qm -> remVerts[0],qm -> nQ - qm -> remQuads[0] );
	fprintf(fout," INTERIOR VERTICES\n");
	for ( i = 0 ; i < 100; ++i) {
		if ( intVal[i]  > 0 ) fprintf(fout," VALENCE %d = %d VERTICES\n", i, intVal[i]);
	}
	fprintf(fout," BOUNDARY VERTICES\n");
	for ( i = 0 ; i < 100; ++i) {
		if ( boundVal[i]  > 0 ) fprintf(fout," VALENCE %d = %d VERTICES\n", i, boundVal[i]);
	}
	fclose(fout);
	return ;
}


static  void
print_mesh(quadMap *qm, char name[]) {
	int i,k, j, v, id, count, d;
	double eval[18];
	FILE *fout;
	printf("WRITING ON %s\n",buffer);
	fout = fopen(name, "w");
	if (fout == NULL ) return;
	for ( i = 0 ; i < qm ->nQ; ++i) {
		count = 0;
		for ( k = 0; k < 4; ++k) {
			if (qm -> quadIdx[4*i + k] == -2 ) count += -2;
			else if ( qm -> quadIdx[4*i + k] >0 &&  qm -> quadIdx[4*i + k] <= qm -> totVerts ) {
                 v =   qm -> quadIdx[4*i + k] - 1;
				fprintf(fout, "%lf  %lf  %lf %d  %lf  %lf\n",qm->xyzs[3*(qm->quadIdx[4*i + k] - 1)], qm->xyzs[3*(qm->quadIdx[4*i + k] - 1) +1], qm->xyzs[3*(qm->quadIdx[4*i + k] - 1) + 2], qm -> quadIdx[4*i + k],
						qm -> uvs[2*(qm->quadIdx[4*i + k] - 1)] , qm -> uvs[2*(qm->quadIdx[4*i + k] - 1) + 1]);
				double dist = 0.0;
				EG_evaluate(qm -> face, &qm -> uvs[2*v ], eval);
				for ( d = 0 ; d < 3; ++d) dist += ( eval[d] - qm->xyzs[3*v + d]) * ( eval[d] - qm->xyzs[3*v + d]);
				dist = sqrt  (dist);
				if ( dist  > 1.e-07) {
					printf(" IN QUAD %d  VERTEX %d. UVs and xyzs are mismatched.  UV %lf  %lf \n", i+1, v+1, qm -> uvs[2*v], qm -> uvs[2*v + 1]);
					for ( d = 0 ; d < 3; ++d) printf( "%lf  != %lf \t", eval[d], qm->xyzs[3*v - 1 + d]);
					exit(1);
				}
			}
			else
				printf(" WHAT !!!!  QUAD %d CLAIMS TO HAVE VERTEX %d\n", i + 1, qm -> quadIdx[4*i + k]);
		}
		if ( count != -8 && count != 0) {
			printf(" WHAT KIND OF QUAD DO WE HAVE ?!\n");
			print_quad_specs (qm, i + 1);
			exit(1);
		}
		if (qm -> quadIdx[4*i] != -2 ) {
			fprintf(fout, "%lf  %lf  %lf %d\n",qm->xyzs[3*(qm->quadIdx[4*i] - 1)], qm->xyzs[3*(qm->quadIdx[4*i] - 1) +1], qm->xyzs[3*(qm->quadIdx[4*i] - 1) + 2], qm -> quadIdx[4*i]);
			fprintf(fout,"\n\n");
		}
	}
	fclose(fout);
}



