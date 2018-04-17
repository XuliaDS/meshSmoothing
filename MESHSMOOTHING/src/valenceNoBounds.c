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


#ifndef SRC_VNB_H_
#define SRC_VNB_H_

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
#include "IO.h"
#include "wsserver.h"
#define PI              3.1415926535897931159979635
#define ERRCTT  3.554147

#define TRYSWAPPING         0
#define TRYCOLLAPSING       1
#define TRYSPLITTING        2
#define TRYSWAPANDSPLIT     3
#define TRYCOLLAPSEANDSPLIT 4


int TOTSWAPS         = 0;
int OVERALLSWAPS     = 0;
int TOTCOLLAPSES     = 0;
int OVERALLCOLLAPSES = 0;
int TOTSPLITS        = 0;
int OVERALLSPLITS    = 0;


int INVALIDCOUNT = 0;
int PLANESURFACE = 0;

#define MINRATIO 0.2
#define MAXRATIO 2.0
int OPTCOUNT  = 0;
int FCOUNT = 1;
int CONVC = 0;

#define EPSAREA 1.E-02
#define DEPS 1.E-09
#define EPSDIR 1.e-06
#define OPTITOL 0.1
#define DTOL 1.e-9

int  IOSTATUS = 0;
#define DEBUG_AREA
int NOUT = 150;
int NOUTMESH = 1000;
#define NLOPTMAXEVAL 100000

#define MULTIOPTI 1

#define DEBUG
int GLO = 0 ;

int QUADTESS  = 0;
int TESSCOUNT = 0;
#define ANGPASS 0.3
#define  FRONTFACE 0
int      SCOUNT = 0;
#define FORCETRIANGULATION
#define FORCEQUADS

/* globals used in these functions */
static wvContext *cntxt;
static bodyData  *bodydata;
static int       nbody;
static float     focus[4];
FILE *filOpti;
char optiName[33];

static int EG_removeDoublet(quadMap *qm, int vID) ;
static int EG_swapAndSplit(quadMap *qm, int vID);
static int EG_edgeSwap(quadMap *qm, int vID);
static int EG_vertexCollapse(quadMap *qm, int vID);
static int EG_vertexSplit(quadMap *qm, int vID, int );
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
	print_quad_specs(qm, qID);
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
				print_quad_specs(qm, auxQ + 1);
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
resetQm(quadMap *qm, int *quadIdx, int *quadAdj, int *vtype, int **valences, double *uvs) {
	int i, j, k, stat;
	double eval[18];
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
			qm -> valence[j]    = EG_reall( qm ->valence[j], sizeof(int));
			qm -> valence[j][0] = -1;
		} else {
			qm -> valence[j]     = EG_reall(qm -> valence[j], (valences[j][0] + 2)*sizeof(int));
			for ( i = 0 ; i < valences[j][0] + 2; ++i)
				qm -> valence[j][i] = valences[j][i];
		}
	}
	return EGADS_SUCCESS;
}

static int
backupQm(quadMap *qm, int **quadIdx, int **quadAdj, int **vtype, int ***valences, double **uvs) {
	int i, j;
	*quadIdx   = (int*)     EG_alloc(4*qm -> nQ       *sizeof(int   ));
	*quadAdj   = (int*)     EG_alloc(4*qm -> nQ       *sizeof(int   ));
	*vtype     = (int*)     EG_alloc(  qm -> totVerts *sizeof(int   ));
	*valences  = (int **)   EG_alloc(  qm -> totVerts *sizeof(int*  ));
	*uvs       = (double *) EG_alloc(2*qm -> totVerts *sizeof(double));
	if ( (*quadIdx) == NULL || (*quadAdj) == NULL || (*valences) == NULL
			|| (*vtype) == NULL || (*uvs) == NULL) return EGADS_MALLOC;
	for ( j = 0; j< qm -> totVerts; ++j) {
		(*valences)[j] = (int *) EG_alloc((2 + 100)*sizeof(int));
		if ((*valences)[j] == NULL ) {
			EG_free(*quadIdx);
			EG_free(*quadAdj);
			EG_free(*vtype);
			EG_free(*valences);
			return EGADS_MALLOC;
		}
		for ( i = 0 ; i < qm -> valence[j][0] + 2; ++i) (*valences)[j][i] = qm -> valence[j][i];
	}
	for ( i = 0 ; i < qm ->nQ; ++i) {
		for ( j = 0 ; j < 4; ++j){
			(*quadIdx)[4*i + j] = qm -> quadIdx[4*i + j];
			(*quadAdj)[4*i + j] = qm -> quadAdj[4*i + j];
		}
	}
	for ( i = 0 ; i < qm ->totVerts; ++i) {
		(*uvs)[2*i] = qm -> uvs[2*i]; (*uvs)[2*i + 1] = qm -> uvs[2*i + 1];
		(*vtype)[i] = qm -> vType[i];
	}
	return EGADS_SUCCESS;
}


static int
EG_splittingOperation(quadMap *qm, vStar *star, int poly[], int id0, int distanceToId0);

static int
EG_meshOperator(quadMap *qm, int opType, int region, int *activity) {
	int i, j, stat, doStuff, auxActivity, count;
	*activity = 0;
	if ( opType == TRYCOLLAPSEANDSPLIT) {
		for ( j = 0 ; j < qm -> nQ; ++j)
			count = TOTCOLLAPSES;
		stat = EG_collapseAndSplit(qm, j + 1);
		if ( count < TOTCOLLAPSES) *activity = 1;
	}
	if (region >=0 && opType != TRYSWAPPING) return EGADS_SUCCESS; //we don't split or collapse at boundary1
	for (i = 0 ; i < qm -> totVerts; ++i) {
		doStuff = 0;
		if (region == 0 && qm -> vType[i] == 0  && qm -> valence[i][0] > 2) {
#ifdef DEBUG
			printf(" Checking that corner are regular and force swap otherwise\n");
#endif
			doStuff = 1;
		}
		else if (region == 1 && qm -> vType[i] > 0  && qm -> valence[i][0] > 3) {
#ifdef DEBUG
			printf(" Checking that edges are regular and try to force swap otherwise\n");
#endif
			doStuff = 1;
		}
		else if (region == -1 && qm -> vType[i] == -1 ) {
#ifdef DEBUG
			printf(" Interior vertices \n");
#endif
			doStuff = 1;
		}
		if (doStuff) {
			switch(opType) {
			case TRYSWAPPING:
				count = TOTSWAPS;
				stat  = EG_edgeSwap(qm, i + 1);
				printf(" STAT IN EDGE SWAP %d \n",stat);
				if ( stat != EGADS_SUCCESS) return stat;
				if ( TOTSWAPS > count) {
					printf(" \n\n ======  SWAPPED VERTEX COUNT %d \n", TOTSWAPS);
					print_mesh(qm, &FCOUNT);
					*activity = 1;
				}
				break;
			case TRYCOLLAPSING:
				count = TOTCOLLAPSES;
				stat  = EG_vertexCollapse(qm,  i + 1);
				printf(" STAT IN COLLAPSEP %d \n",stat);
				if ( TOTCOLLAPSES > count) {
					printf(" \n\n ======  COLLAPSED VERTEX COUNT %d: TRY SWAPPING WITH NEW MESH \n", TOTCOLLAPSES);
					print_mesh(qm, &FCOUNT);
					auxActivity = 1;
					stat        = EG_meshOperator(qm, TRYSWAPPING, -1, activity);
					*activity   = auxActivity;
				}
				break;
			case TRYSPLITTING:
				count = TOTSPLITS;
				stat = EG_vertexSplit(qm,  i + 1, 1);
				printf(" STAT IN EDGE SPLIT %d \n",stat);
				if ( TOTSPLITS > count) {
					printf(" \n\n ======  SPLIT VERTEX COUNT %d: TRY SWAPPING WITH NEW MESH \n", TOTSPLITS);
					print_mesh(qm, &FCOUNT);
					auxActivity = 1;
					stat = EG_meshOperator(qm, TRYSWAPPING, -1,  activity);
					*activity = auxActivity;
				}
				break;
				/*case TRYSWAPANDSPLIT:
				count = TOTSPLITS;
				stat  = EG_swapAndSplit(qm, star);
				if ( TOTSPLITS > count) {
					printf(" SWAP AND SPLIT VERTEX: TRY SWAPPING WITH NEW MESH \n");
					print_mesh(qm, &FCOUNT);
					auxActivity = 1;
					stat        = EG_meshOperator(qm, TRYSWAPPING, -1, activity);
				 *activity   = auxActivity;
				}
				break;*/
			}
		}
	}
	return EGADS_SUCCESS;
}

static int checkInvalidElement(quadMap *qm, int qID, int v0ID, int pullVertex) {
	int    maxIt,i, j, k, v[4], stat, aux, v1ID = -1, adjQ = -1, area, it = 0, invalid;
	double dir[3], dirOpp[3], uvOri[2], vXYZ[18], penalty, angles[4], quadArea;
	vStar *star = NULL;
	if ( qm -> quadIdx[4*(qID-1)] == -1 || qm -> vType[v0ID - 1] == -2) {
		printf(" IT DOESN't EXIST. DON'T CHECK\n");
		return EGADS_SUCCESS; // empty quad. probably has been removed elsewhere
	}
	print_quad_specs(qm, qID);
	stat = EG_buildStar(qm, &star, v0ID);
	if ( stat != EGADS_SUCCESS ) {
		goto cleanup;
	}
	if ( qm -> vType[v0ID -1 ] >= 0 ) {
		stat = quad_algebraic_area(qm, qID, v0ID, &area, angles, &quadArea);
		EG_free(star -> verts);
		EG_free(star -> quads);
		EG_free(star);
		if ( area == 1) return EGADS_SUCCESS;
		else            return EGADS_GEOMERR;
	}
	invalid = 0 ;
	for ( i = 0 ; i < star -> nQ; ++i) {
		stat = quad_algebraic_area(qm, star ->quads[i], v0ID, &area, angles, &quadArea);
		if (stat != EGADS_SUCCESS ) {
			goto cleanup;
		}
		if ( area != 1) {
			invalid = 1;
			break;
		}
	}
	if ( pullVertex == 0 ){
		EG_free(star -> verts);
		EG_free(star -> quads);
		EG_free(star);
		if (invalid) return EGADS_GEOMERR;
		else         return EGADS_SUCCESS;
	}
	// INTERIOR VERTEX: MOVE IT AROUND
	// GET ADJACENT QUAD AND APPROPRIATE EDGE TO PULL FROM
	for ( k = 0 ; k < 4; ++k)
		if (qm -> quadIdx[4*(qID -1) + k] == v0ID) break;
	if ( qm -> quadAdj[4*(qID -1) + k] == -1)
	{
		printf(" VERTEX %d IS BOUNDARY. YOU CAN'T MOVE IT\n",v0ID + 1);
		print_quad_specs(qm, qID);
		stat = EGADS_GEOMERR;
		goto cleanup;
	}
	adjQ = qm ->quadAdj[4*(qID -1) +  k];
	print_quad_specs(qm, adjQ);
	// fin v0ID in adjQ
	for ( k = 0 ; k < 4; ++k) {
		if ( qm->quadIdx[4*(adjQ - 1) + k  ] == v0ID ) {
			v1ID = qm->quadIdx[4*(adjQ - 1) + (int)(k +2)%4  ];  // pull towards opposite vertex
			break;
		}
	}
	//printf(" UV V0 %lf  %lf  UV V1  %lf  %lf \n",qm ->uvs[2*(v0ID - 1)    ], qm ->uvs[2*(v0ID - 1)  +1  ], qm ->uvs[2*(v1ID - 1)    ], qm ->uvs[2*(v1ID - 1)  +1  ]);
	dir  [0]   = 0.1*(qm ->uvs[2*(v1ID - 1)    ] - qm ->uvs[2*(v0ID - 1)    ]);
	dir  [1]   = 0.1*(qm ->uvs[2*(v1ID - 1) + 1] - qm ->uvs[2*(v0ID - 1) + 1]);
	//printf(" DIR %lf  %lf  OPO %lf  %lf\n",dir[0], dir[1], dirOpp[0], dirOpp[1]);
	uvOri[0]   = qm ->uvs[2*(v0ID - 1)    ];
	uvOri[1]   = qm ->uvs[2*(v0ID - 1) + 1];
	maxIt      = 10000;
	printf(" +++++++++++++++++++++ CHECK INVALID ELEMENT  PULL %d TO %d ++++++++++++++++++++++++++++++++++\n", v0ID, v1ID);
	stat = EGADS_SUCCESS;
	IOSTATUS = 1;
	while (invalid)
	{
		qm ->uvs[2*(v0ID - 1)    ] += dir[0];
		qm ->uvs[2*(v0ID - 1) + 1] += dir[1];
		stat                        = EG_evaluate(qm -> face, &qm ->uvs[2*(v0ID - 1)], vXYZ);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		qm -> xyzs[3*(v0ID - 1)    ] = vXYZ[0];
		qm -> xyzs[3*(v0ID - 1) + 1] = vXYZ[1];
		qm -> xyzs[3*(v0ID - 1) + 2] = vXYZ[2];
		for ( i = 0 ; i < 2; ++ i)
			dirOpp[i] = qm ->uvs[2*(v0ID - 1) + i] - qm ->uvs[2*(v1ID - 1) + i];
		if ( (dirOpp[0]*dir[0] + dirOpp[1]*dir[1]) > 0) {
			printf(" WE HAVE GONE BEYOND THE POINT\n");
			stat = EGADS_GEOMERR;
			goto cleanup;
		}
		invalid = 0 ;
		for ( i = 0 ; i < star -> nQ; ++i) {
			stat = quad_algebraic_area(qm, star ->quads[i], v0ID, &area, angles, &quadArea);
			if (stat != EGADS_SUCCESS ) {
				printf(" FAILED TO COMPUTE ALGEBRAIC AREA %d\n",stat);
				goto cleanup;
			}
			if ( area != 1) invalid = 1;
		}
		if ( it > maxIt) {
			printf(" WE HAVE REACHED THE MAX ITERATIONS \n");
			stat = EGADS_INDEXERR;
			goto cleanup;
		}
		++it;
	}
	cleanup:
	EG_free(star -> verts);
	EG_free(star -> quads);
	EG_free(star);
	printf(" LEAVING CHECKINVALID WITH STAT %d \n", stat);
	return stat;
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
	if ( fabs(dotNP - c) ) {
		proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
		return EGADS_SUCCESS;
	}
	lambda  = (c - dotNP)/dotNN;
	proj[0] = p[0] + lambda*normal[0];
	proj[1] = p[1] + lambda*normal[1];
	proj[2] = p[2] + lambda*normal[2];
	// check that point belongs to plane.
	dist  = normal[0] * proj[0] + normal[1] * proj[1] +  normal[2] * proj[2];
	if( fabs(dist - c) < DEPS) {
		return EGADS_SUCCESS;
	}
	else{
		printf(" POINT SHOULD BELONG TO PLANE!!!!! %lf ~= %lf\n",dist,c);
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




#define MINKCHANGE 0.1

static int
computePenalty(quadMap *qm,  int starID, int qID, double *penalty, double *angles, double *quadArea) {
	double dot, z, zerr, alpha, dev[3], basis[3*3], normOri = 0.0, norm[2], scale;
	double curvature[8], invalidSize = 0.0;
	int i, j, id0, stat, validArea, q, v, movingVert, tooSmall = 0, invalid = 0, vv, idx;
	q    = qm ->star[starID] -> quads[qID];
	v    = qm ->star[starID] -> verts[0  ];
	stat = quad_algebraic_area(qm, q, v, &validArea, angles, quadArea);
	for ( i = 0 ; i < 4; ++i) {
		vv = qm ->quadIdx[4*(q - 1) + i];
		if (vv == v ) id0 = i;
	}
	if ( stat != EGADS_SUCCESS) {
		printf(" ERROR COMPUTING PENALTY FUNCTION  %d",stat);
		return stat;
	}
	--v; // +1 bias
	// USE CURVATURES TO MEASURE MOVEMENT
	// Distance to linking vertices
	// first link
	vv = qm -> star[starID] -> verts[2*qID + 1] - 1;
	for ( i = 0 ; i < 3; ++i) dev[i] = qm ->xyzs[3*vv + i ] - qm ->xyzs[3*v + i ];
	norm[0]  = sqrt(dev[0] * dev[0] + dev[1] * dev[1] + dev[2] * dev[2]);
	if ( norm[0] > qm -> vLengths[2*v] || norm[0] < qm->vLengths[2*v + 1]) {
		invalidSize = norm[0];
		invalid = 1;
	}
	// second link
	if ( (2 * qID + 3 ) ==  qm ->star[starID] -> nV ) vv = qm ->star[starID] -> verts[1] - 1;
	else  vv = qm ->star[starID] -> verts[2*qID + 3] - 1;
	for ( i = 0 ; i < 3; ++i) dev[i] = qm ->xyzs[3*vv + i ] - qm ->xyzs[3*v + i ];
	norm[1] = sqrt(dev[0] * dev[0] + dev[1] * dev[1] + dev[2] * dev[2]);
	if ( norm[1] > qm -> vLengths[2*v] || norm[1] < qm -> vLengths[2*v + 1]) {
		invalidSize = norm[1];
		invalid = 1;
	}
	// ERRROR FUNCTION
	alpha = angles[0];
	idx = qm ->star[starID] -> verts[0] - 1;
	for ( i = 1; i < 4; ++i) {
#ifdef DEBUG
		if (OPTCOUNT%NOUT == 0 )
			printf(" COMPARING ANGLE %lf  TO  %lf\t", angles[i], alpha);
#endif
		vv = qm ->quadIdx[4*(q - 1) + (i + id0)%4];
		if ( (i != 2) && qm ->vType[vv -1] == -1) {
			movingVert = 1;
			for ( j = 0 ; j < qm -> nS; ++j) {
				if(vv == qm -> star[j]->verts[0]) movingVert = 0;
			}
			if( (movingVert == 1) && (angles[i] >= M_PI*0.99)  && (angles[i] > alpha)) {
				alpha = angles[i];
				idx   = vv - 1;
				i = 4;
			}
			else if( (movingVert == 1) && ( fabs( (2.0 * M_PI)/(double)qm ->valence[vv - 1][0] - angles[i] ) > fabs( 2.0*M_PI/(double)qm->valence[idx][0] - alpha ) ) ) {
				alpha = angles[i];
				idx = vv - 1;
			}
		}
#ifdef DEBUG
		if (OPTCOUNT%NOUT == 0 )
			printf(" ===> ALPHA %lf\n", alpha);
#endif
	}
	dot = cos(alpha - 0.5*M_PI);
	if ( alpha > 2.0*M_PI - acos(ANGPASS)) dot = -1.0;
	scale   = alpha;
	if ( dot >0 && dot < 0.5*ANGPASS) scale = 2*M_PI-alpha;
	z      = ERRCTT*((dot-ANGPASS*0.5)/ANGPASS);
	zerr   = erfc(z);
	*penalty = zerr * exp(zerr *scale/M_PI);
	if ( invalid == 1) {
		scale = zerr + exp(1 - (norm[0] + norm[1]));
		if ( scale > *penalty) *penalty = scale;
	}

	if (OPTCOUNT%NOUT == 0 ) {
		printf("\n :::::::::: PENALTY FUNCTION AT VERTEX %d in QUAD %d : AREA = %lf  ANGLE AT VERT %lf-> VALIDAREA %d TOO SMALL TOO BIG = %d\n",
				qm ->quadIdx[4*(q - 1) + (i + id0)%4], qID, *quadArea, angles[0], validArea, invalid);
		printf(" PENALTY : alpha = %lf  dot = %lf  scale = %lf   z  %lf  zerr = %lf  penalty %lf\n ",
				alpha, dot, scale, z, zerr, *penalty);
	}
	angles[1] = alpha;
	return EGADS_SUCCESS;

}

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


static int
quad_algebraic_area(quadMap *qm, int qID, int vID,  int *validArea, double *quadAngles, double *quadArea) {
	int i, id, v, vPos, qV[4], k, tri, piv[4], stat, tooClose, degenerated[4], ftype;
	double vTri[18*3], evalVert[18], cross[4], dotCross[2*4], area[2*4], norm[2], theta[2*4], vNormal[4], uv_dir[2], uv_eps[2],  vPro[8], vEps[8];
	double naux,  z, zerr, alpha, dot;
	for ( i = 0 ; i < 4; ++i) {
		qV[i] = qm ->quadIdx[4*(qID - 1) + i] - 1;
		if (qV[i] == (vID -1) ) vPos = i;
	}
	if ( vPos == -1 ) {
		printf(" VERTEX %d is STAR CENTRE AND SHOULD BELONG TO QUAD %d\n",vID, qID);
		return EGADS_INDEXERR;
	}
	for ( k = 0 ; k < 4; ++k) { // Check additional area to ensure that the quad is correct
		for ( tri = 0 ; tri < 2; ++ tri) {
			piv[0]  = qV[(k  + vPos          )%4];
			piv[1]  = qV[(k  + vPos + 1 + tri)%4];
			piv[2]  = qV[(k  + vPos + 2 + tri)%4];
			vTri[0] = qm ->xyzs[3*piv[0]]; vTri[1] = qm ->xyzs[3*piv[0]  + 1]; vTri[2] = qm ->xyzs[3*piv[0] + 2];
			stat = EG_evaluate(qm->face, &qm -> uvs[2*piv[0]], evalVert);
			if ( stat != EGADS_SUCCESS) return stat;
			cross_product(&evalVert[3], &evalVert[6], vNormal);
			if ( qm -> face -> mtype == SREVERSE )
				cross_product(&evalVert[6], &evalVert[3], vNormal);
			else
				cross_product(&evalVert[3], &evalVert[6], vNormal);
			EG_unitVector(vNormal, &vNormal[3]);
			for ( id = 1 ; id <= 2; ++id) { // move away from vertex towards B and C respectively
				uv_dir[0] = qm -> uvs[2*piv[id]    ] - qm -> uvs[2*piv[0]    ];
				uv_dir[1] = qm -> uvs[2*piv[id] + 1] - qm -> uvs[2*piv[0] + 1];
				naux = uv_dir[0] * uv_dir[0] + uv_dir[1] * uv_dir[1];
				naux = sqrt(naux);
				degenerated[k] = 0;
				if ( naux < DEPS ) {
					degenerated[  k    ] = 1;
					theta      [2*k    ] = 0.0;
					theta      [2*k + 1] = 0.0;
					area       [2*k    ] = 0.0;
					area       [2*k + 1] = 0.0;
					dotCross   [2*k    ] = -1.0;
					dotCross   [2*k + 1] = -1.0;
					id = 3;
					tri = 2;
				}
				if(!degenerated[k]) {
					uv_dir[0] /= naux;                uv_dir[1] /= naux;
					uv_eps[0]  = qm -> uvs[2*piv[0]]; uv_eps[1]  = qm -> uvs[2*piv[0] + 1];
					do {
						uv_eps[0] += EPSDIR * uv_dir[0];
						uv_eps[1] += EPSDIR * uv_dir[1];
						stat       = EG_evaluate(qm -> face, uv_eps, evalVert);
						if( stat  != EGADS_SUCCESS ) return stat;
						stat       = EG_projectToTangentPlane(vNormal,vTri, evalVert, &vTri[18*id]);
						if( stat  != EGADS_SUCCESS ) return stat;
						naux  = 0.0;
						for ( i = 0 ; i < 3; ++i)
							naux += (vTri[18*id + i] - vTri[i]) * (vTri[18*id + i] - vTri[i]);
						naux       = sqrt(naux);
					} while(naux < EPS);
				}
			}
			if(!degenerated[k]) {
				// get vectors vA vB from vertex v
				for ( i = 0 ; i < 3; ++i) {
					vEps[i]     = vTri[18*1 + i] - vTri[i];
					vEps[4 + i] = vTri[18*2 + i] - vTri[i];
				}
				EG_unitVector (&vEps[0], &vEps[3] );
				EG_unitVector (&vEps[4], &vEps[7] );
				if ( vEps[3] < 1.E-7 || vEps[7] < 1.E-7) theta[2*k + tri] = 0.0;
				else {
					dot = 0.0;
					for     ( i = 0 ; i < 3; ++i)
						dot += vEps[i] * vEps[i + 4];
					if      ( fabs(dot - 1.0) < DTOL ) theta[2*k + tri] = 0.0;
					else if ( fabs(dot + 1.0) < DTOL ) theta[2*k + tri] = M_PI;
					else                               theta[2*k + tri] = acos(dot);
				}

				for ( i = 0 ; i < 3; ++i){
					vTri[       i] = qm -> xyzs[3*piv[0] + i ];
					vTri[  18 + i] = qm -> xyzs[3*piv[1] + i ];
					vTri[2*18 + i] = qm -> xyzs[3*piv[2] + i ];
				}
				area[2*k + tri] = triang_cross_area(&vTri[0], &vTri[18], &vTri[36], cross);
				EG_unitVector(cross, &cross[3]);
				dotCross[2*k + tri] = cross[0] * vNormal[0] + cross[1] * vNormal[1] + cross[2] * vNormal[2];
#ifdef DEBUGG
				if (IOSTATUS)
					printf(" *********  DOT CROSS PRODUCT %lf  -> THETA = %lf \t ", dotCross[2*k + tri], theta[2*k + tri]);
#endif
				if ( dotCross[2*k + tri] < - 0.8 ) theta[2*k + tri] = 2.0 * M_PI - theta[2*k + tri];  // we have anti-clockwise orientation. Postive area implies reversed.
#ifdef DEBUGG
				if (IOSTATUS)
					printf("  -> THETA = %lf \t \n", theta[2*k + tri]);
#endif
			}
		}
		if (IOSTATUS)
			printf(" *********  AREA OF QUAD   %lf \n", area[2*k] + area[2*k + 1]);
	}
	*validArea = 1;
	if ( fabs( (area[0] + area[1]) - (area[2] + area[3]) ) > EPSAREA ) {
		if (IOSTATUS)
			printf(" AREA %lf != %lf\n", area[0] + area[1], area[2] + area[3] );
		*validArea = -1;
	}
	*quadArea = area[0] + area[1];
	for ( k = 0 ; k < 4; ++k) {
		quadAngles[k] = theta[2*k] + theta[2*k+1];
		//if ( quadAngles[k] >  2.0*M_PI ) quadAngles[k] -= 2.0*M_PI;
		if ( (qm -> vType[qV[(k  + vPos )%4]] == -1) && (quadAngles[k] >= M_PI*0.99)) *validArea = -1;
		if ( (qm -> vType[qV[(k  + vPos )%4]] == -1) && (degenerated[k]))             *validArea = -1;
		if (IOSTATUS)
			printf(" ANGLES AROUND %d =  %lf   AND AREAS %lf  ANGLE WITH SURFACE NORMAL  %lf  %lf\n",
					qm ->quadIdx[4*(qID - 1) + (vPos + k)%4] , quadAngles[k], area[2*k] + area[2*k + 1], dotCross[2*k], dotCross[2*k +1]);
	}
	if (IOSTATUS)
		printf(" LEAVING WITH VALID AREA = %d  SIZE  %lf \n",*validArea,*quadArea);
	return EGADS_SUCCESS;
}

static void
EG_vertexRegularity(quadMap *qm, int list[], int *vIrr, int n)
{
	int i;
	for ( i = 0; i < n; ++i) {
		if (qm -> vType[ list[2*i] -1] == 0) {
			if (list[2*i + 1] == 2) vIrr[i] = 0 ;
			else vIrr[i] = 100; //edge vertices have "SPECIAL" regularity
		}
		else if (qm -> vType[ list[2*i] -1] > 0) {
			if      ( list[2*i+1] == 2 ) vIrr[i] = -1;
			else if ( list[2*i+1] == 3 ) vIrr[i] =  0;
			else                         vIrr[i] = list[2*i + 1];
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
	quad0         = qm -> valence[vID - 1][1] - 1;
	//printf(" STAR CENTRE %d START AT ID %d \n", vID, quad0 + 1);
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
		printf("I couldn't find the star centre (Vertex %d) in Quad %d\n. quadMap build UNSUCCESSFUL. Check updateUVs function!\n", vROSE[0], quad0);
		print_quad_specs(qm, quad0 + 1);
		return (EGADS_INDEXERR);
	}
	vROSE[r] = -1;
	do {
		//printf(" QUAD %d PREV %d ORI%d\n",quadID + 1, prevQuad + 1, quad0+1);
		//printf(" ID %d -> adj %d\n",id0, qm -> quadAdj[4*quadID + id0]);
		for ( j = 0 ; j < 3; ++j) {
			auxm = (id0 + j)%4;
			if ( qm -> quadAdj[4*quadID + auxm] != (prevQuad +1) && qm -> quadAdj[4*quadID + auxm] != (quad0+1)) {
				if ( vROSE[1] != qm -> quadIdx[4*quadID + qLoop[auxm + 1] ]) { // At last quad, this link will be repeated; vROSE[2] is the first link
					vROSE[r++] = qm -> quadIdx[4*quadID + qLoop[auxm + 1] ];
				}
			}
		}
		id0      = (int)(id0 + 3)%4;
		//printf(" ID %d -> adj %d\n",id0, qm -> quadAdj[4*quadID + id0]);
		prevQuad = quadID;
		//print_quad_specs(qm, prevQuad +1);
		quadID   = qm -> quadAdj[4*prevQuad + id0] - 1;
		//printf(" QUAD ID %d \n",qm -> quadAdj[4*prevQuad + id0]);
		//print_quad_specs(qm, quadID +1);
		//for ( k = 0 ; k < r; ++k) printf(" STAR %d = %d\n",k,vROSE[k]);
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
		//printf(" GET ID %d IN QUAD %d\n", vID, quadID + 1);
		if ( quadID < 0) {
			vROSE[r++]    = -1;
			vQUADS[itQ++] = -1;
			quadID        = auxm;
		}
		id0 = -1;
		//printf(" GET ID %d IN QUAD %d\n", vID, quadID + 1);
		//print_quad_specs(qm, quadID+1);
		for ( j = 0 ; j < 4; ++j) {
			if ( ( qm->quadIdx[4*quadID + j ] == vID  ) ) {//&&  (qm ->quadAdj[4*quadID + j] == (prevQuad + 1) ) ) {
				id0  = j;
				break;
			}
		}
		if ( id0 == -1) {
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
EG_vertexSplit(quadMap *qm, int sID, int performSplitting) {
	int poly[2*3], q, aux, id0, i, j, val;
	int id, idAux, stat, dist, validSplit;
	int  *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, *locVtype, **locVal = NULL;
	double xyz[18], *uv0 = NULL;
	vStar *star;
	if ( qm -> valence[ sID -1][0] < 5 ) return EGADS_SUCCESS;
	if ( qm -> remQuads[0] == 0 ) {
		printf(" At the moment, we only insert if we have collapsed vertices, to maintain overall mesh numbers\n");
		return EGADS_SUCCESS;
	}
	stat = EG_buildStar(qm, &star, sID);
	if ( stat != EGADS_SUCCESS) return stat;
	stat = backupQm ( qm, &locIdx, &locAdj, &locVtype, &locVal, &uv0);
	if ( stat != EGADS_SUCCESS) {
		EG_free(star -> verts);
		EG_free(star -> quads);
		EG_free(star);
		return stat;
	}
	poly[2*0    ] = star -> verts  [0];
	poly[2*0 + 1] = qm   -> valence[star -> verts[0] -1][0];
	for ( q = 0 ; q < star -> nQ; ++q) {
		id0           = 2*q + 1;
		poly[2*1    ] = star -> verts  [id0];
		poly[2*1 + 1] = qm   -> valence[poly[2*1] - 1][0];
		for ( i = 0 ; i < poly[1] - 4; ++i) {
			aux  = id0 + 4  + 2 * i;
			dist = 4 + 2*i;
			if ( aux >= star -> nV) aux -= (star -> nV - 1) ;
			poly[2*2    ] = star -> verts[aux];
			poly[2*2 + 1] = qm   -> valence[poly[2*2] - 1][0];
		}
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
				if ( performSplitting) {
					stat = EG_splittingOperation(qm, star, poly, id0, dist);
#ifdef DEBUG
					printf(" SPLITTING WENT %d \n",stat);
#endif
					EG_free(star -> verts);
					EG_free(star -> quads);
					EG_free(star);
					if ( stat != EGADS_SUCCESS)
						return resetQm(qm, locIdx, locAdj, locVtype, locVal, uv0);
					else return EGADS_SUCCESS;
				}
				else  return EGADS_SUCCESS;
			}
		}
	}
	EG_free(star -> verts);
	EG_free(star -> quads);
	EG_free(star);
	for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
	EG_free(locVal);
	EG_free(locIdx);
	EG_free(locAdj);
	EG_free(locVtype);
	EG_free(uv0);

	if ( !performSplitting ) return EGADS_EMPTY;
	return EGADS_SUCCESS;
}


static int
EG_splittingOperation(quadMap *qm, vStar *star, int poly[], int id0, int distanceToId0) {
	int qID[2], qAux[2], edges[4],  movingVerts[2],  valAux[100];
	int q, newV, newQ, val, id, i, j, aux, stat, idAux, v0, dir[2];
	double xyz[18];
	newV = qm -> remVerts[qm -> remVerts[0]];
	newQ = qm -> remQuads[qm -> remQuads[0]];
	--qm ->remVerts[0]; --qm ->remQuads[0];
	qm -> vType  [newV - 1] = -1;
	if ( qm -> remQuads[0] != qm -> remVerts[0] ) {
		printf(" The number of removed vertices (%d) and quads (%d) miss matches!! this is wrong! I am leaving program\n", qm -> remQuads[0], qm -> remVerts[0]);
		return EGADS_INDEXERR;
	}
	// Split  Valences from original vertex
	val = (int)(distanceToId0/2) + 1;
	id  =  poly[0] - 1;
	qm -> valence[id]    = EG_reall(qm -> valence[id], ( val + 2 )*sizeof(int));
	qm -> valence[id][0] = val;
	qm -> valence[id][1] = newQ;
	for ( i = 0; i < val; ++i) {
		aux = id0 + 2 * i; if ( aux >= star -> nV) aux -= (star -> nV - 1) ;
		qm -> valence[id][2 + i] = star -> verts[aux];
	}
	// Create Quad, point and reallocate adjacents
	qID[0] = star -> quads[(int)(id0                 - 1)/2    ];
	qID[1] = star -> quads[(int)(id0 + distanceToId0 - 1)/2 - 1];
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
	q = (int)(id0 + distanceToId0 - 1)/2;
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
		qm -> valence[id]    = EG_reall(qm -> valence[id], ( val + 2 )*sizeof(int));
		qm -> valence[id][0] = val;
		qm -> valence[id][1] = newQ;
		for ( i = 0; i < val - 1; ++i)
			qm -> valence[id][2 + i] = valAux[i];
		qm -> valence[id][2 + val - 1] = newV;
	}
	// Create valence pointer for new vertex and update links
	val = poly[1] - qm -> valence[poly[0] - 1][0] + 2;
	id  = newV    - 1;
	qm -> valence[id]    = EG_reall(qm -> valence[id], ( val + 2 )*sizeof(int));
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
	printf(" VALENCES FOR OLD VERTEX \n");
	for ( i = 0; i < qm -> valence[poly[0] - 1][0]; ++i) {
		idAux = qm -> valence[poly[0] - 1][2 + i];
		printf(" V %d = %d\n",i, idAux);
		for ( j = 0 ; j < qm -> valence[idAux-1][0]; ++j)
			printf(" points back %d = %d\n",j, qm ->valence[idAux -1][2+j]);
	}
	printf(" VALENCES FOR NEW VERTEX \n");
	for ( i = 0; i < qm -> valence[id][0]; ++i) {
		idAux = qm -> valence[id][2 + i];
		printf(" V %d = %d\n",i, qm ->valence[id][2+i]);
		for ( j = 0 ; j < qm -> valence[idAux -1][0]; ++j)
			printf(" points back %d = %d\n",j, qm ->valence[idAux -1][2+j]);
	}
#endif
	// Shift "up" old vertex
	aux = id0 + 2;
	if ( aux >= star -> nV) aux -= (star -> nV + 1);
	id = star -> verts[aux] - 1;
	v0 = poly[0] - 1;
	dir[0]                     = qm -> uvs[2*id    ] - qm -> uvs[2*v0    ];
	dir[1]                     = qm -> uvs[2*id + 1] - qm -> uvs[2*v0 + 1];
	qm -> uvs[2*(newV -1)    ] = qm -> uvs[2*v0    ] - 0.1 * dir[0];
	qm -> uvs[2*(newV -1) + 1] = qm -> uvs[2*v0 + 1] - 0.1 * dir[1];
	qm -> uvs[2*v0    ]        = 0.5* ( qm -> uvs[2*id    ] + qm -> uvs[2*v0    ]);//0.25 * dir[0];
	qm -> uvs[2*v0 + 1]        = 0.5* ( qm -> uvs[2*id + 1] + qm -> uvs[2*v0 + 1]);//0.25 * dir[0];
	stat   = EG_evaluate(qm->face, &qm ->uvs[2*v0], xyz);
	qm -> xyzs [3*v0] = xyz[0]; qm->xyzs [3*v0 + 1] = xyz[1]; qm->xyzs [3*v0 + 2] = xyz[2];
	stat   = EG_evaluate(qm->face, &qm ->uvs[2*(newV - 1)], xyz);
	if ( stat != EGADS_SUCCESS) return stat;
	qm -> xyzs [3*(newV - 1)] = xyz[0]; qm->xyzs [3*(newV - 1) + 1] = xyz[1]; qm->xyzs [3*(newV - 1) + 2] = xyz[2];
	stat = checkInvalidElement(qm, newQ, newV, 1);
	if ( stat != EGADS_SUCCESS)
		return stat;
	movingVerts[0] = newV;
	movingVerts[1] = poly[0];
	stat = optimize_angles(qm, 2, movingVerts);
	if ( stat != EGADS_SUCCESS)
		return stat;
	++TOTSPLITS;
	return EGADS_SUCCESS;
}

static int
EG_mergeVertices(quadMap *qm, int qC, int poly[], int collapseToCentre) {
	int vC, vC2,vID, stat, i, e, q, v, vCommon, j, adjQ, k, updatedQuad, vVAL, colVal, aux, auxVal[100], addValence;
	int edges[4], auxAdj, auxAdj2;
	double uv[2], vPos[18];
	vStar *star;
	--qC;
	if ( qC < 0 ) {
		printf(" WE ARE CALLING TO MERGE VERTICES WITH A VOID QUAD!! BUGGGGG\n");
		return EGADS_INDEXERR;
	}
	for ( i = 0 ; i < 4; ++i)
		printf(" POLY %d = %d\n", poly[2*i], poly[2*i + 1]);
#ifdef DEBUG
	printf(" EG MERGE VERTICES:: COLLAPSE v(%d) = %d   TO v(%d) = %d AFFECTED VS v(%d) = %d v(%d) = %d (each will be -1\n  ",
			poly[0], poly[1], poly[2*2], poly[2*2 + 1],
			poly[2*1], poly[2*1 + 1],poly[2*3], poly[2*3 + 1]);
#endif
	// Get Quad Centre to place collapsed vertex
	vC  = poly[2*2] - 1;
	vC2 = poly[0]   - 1;
	stat = EG_buildStar(qm, &star, vC2 + 1);
	if ( stat != EGADS_SUCCESS) {
		return stat;
	}
	if ( collapseToCentre) {
		uv[0] = 0.0; uv[1] = 0.0;
		for ( i = 0 ; i < 4; ++i) {
			uv[0] += qm->uvs[2*(poly[2*i] - 1)    ];
			uv[1] += qm->uvs[2*(poly[2*i] - 1) + 1];
		}
		uv[0] *= 0.25;
		uv[1] *= 0.25;
		stat   = EG_evaluate(qm->face, uv, vPos);
		if ( stat != EGADS_SUCCESS ) {
			EG_free(star -> verts);
			EG_free(star -> quads);
			EG_free(star);
			return stat;
		}
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
		if ( updatedQuad ) {
			++count;
			stat = EG_getCommonEdge(qm, qC + 1,  adjQ + 1, edges);
			if ( stat != EGADS_SUCCESS) {
				printf(" I can't find a common edge to %d  and  %d\n",qC + 1, adjQ  +1);
				EG_free(star -> verts);
				EG_free(star -> quads);
				EG_free(star);
				return EGADS_INDEXERR;
			}
			print_quad_specs(qm, qC + 1);
			print_quad_specs(qm, adjQ + 1);
			printf(" COMMON EDGE TO QUAD %d %d -> %d %d = %d %d, %d %d = %d %d\n",qC + 1, adjQ + 1,
					edges[0], edges[1], qm -> quadAdj[4*qC + edges[0]], qm -> quadAdj[4*adjQ + edges[1]],
					edges[2], edges[3], qm -> quadAdj[4*qC + edges[2]], qm -> quadAdj[4*adjQ + edges[3]]);
			//if ( stat == EGADS_SUCCESS) {
			for ( e = 0 ; e < 2; ++ e) {
				if ( edges[2*e] != -1 && edges[2*e + 1] != -1) {
					printf(" ADJ %d = %d  AND ADJ %d = %d\n", edges[2*e], qm -> quadAdj[4*qC + edges[2*e]],
							edges[2*e + 1], qm -> quadAdj[4*adjQ + edges[2*e + 1]] );
					printf(" CORRESPONDING VERTEX %d %d , %d  %d\n", qm -> quadIdx[4*qC + edges[2*e]],qm -> quadIdx[4*qC + (edges[2*e] + 1)%4]  ,
							qm -> quadIdx[4*adjQ + edges[2*e + 1]],qm -> quadIdx[4*adjQ + (edges[2*e + 1] + 1)%4]  );
					vCommon = qm -> quadIdx[4*qC + edges[2*e]];
					if (qm -> quadIdx[4*qC + edges[2*e]] == vC2 + 1)
						vCommon = qm -> quadIdx[4*qC + (edges[2*e] + 1)%4];
					printf(" COMMON VERTEX TO %d AND %d  = %d \n",qC+1, adjQ + 1, vCommon);
					for ( i = 0; i < 4; ++i) { // loop around quads and take the other quad that shares vertex vCommon
						auxAdj = qm -> quadAdj[4*qC + i ] - 1;
						if ( auxAdj >= 0 ) {
							print_quad_specs(qm, auxAdj + 1);
							for ( k = 0 ; k < 4; ++ k) {
								if ( auxAdj != adjQ  && qm -> quadIdx[4*auxAdj + k] == vCommon) {
									k = 4;
									for ( j = 0 ; j < 4; ++j) {
										if ( qm -> quadAdj[4*auxAdj + j] == qC + 1) { // found the other quad pointing at vCommon
											printf("--------\n QUAD %d pointed at %d\t",adjQ + 1, qm -> quadAdj[(4*adjQ + edges[2*e + 1] )] );
											qm -> quadAdj[(4*adjQ + edges[2*e + 1] )] = auxAdj + 1;
											printf("  now points to %d \n",qm -> quadAdj[(4*adjQ + edges[2*e + 1] )]);
											printf(" QUAD %d pointed at %d\t", auxAdj + 1, qC + 1);
											qm -> quadAdj[4*auxAdj + j]  = adjQ + 1;
											printf("  now points to %d \n",qm -> quadAdj[(4*auxAdj + j )]);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		//}
		if ( count == 2) break;
	}
	for ( i = 0 ; i < 4; ++i) print_quad_specs(qm, qm ->quadAdj[4*qC + i]);
	// Transfer valencies from Collapsed vertex
	// Increase valence of collapsed vertex: 1-> copy old ones
	vVAL   = qm -> valence[vC][0];
	colVal = qm -> valence[vC][0] + qm -> valence[vC2][0] - 2;
	for ( i = 0 ; i < vVAL; ++i) auxVal[i] = qm -> valence[vC][2 + i];
	qm -> valence[vC]    = EG_reall(qm -> valence[vC], (colVal + 2)*sizeof(int));
	qm -> valence[vC][0] = colVal;
	// Assing a existing quad to vertex vC
	qm -> valence[vC][1] = - 1;
	for ( j = 0 ; j < 4; ++j) {
		for ( i = 0 ; i < 4; ++i) {
			q = qm->quadAdj[4*qC + j];
			if ( q > 0) {
				if( qm -> quadIdx[ 4 * (q - 1) + i ] == vC + 1 ) {
					qm -> valence[vC][1] = qm->quadAdj[4*qC + j] ;
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
		if ( qm -> valence[vID][1] == qC + 1) {
			EG_getAdjacentToVertexPair(qm, qC + 1, vC + 1, vID + 1, edges);
			qm -> valence[vID][1] = edges[1];
		}
	}
	// Get remaining valences from vertex 2 and elliminate vertex 2 as link from its neighbors
	k = vVAL + 2;
	for ( j = 0 ; j < qm -> valence[vC2][0]; ++j) {
		vID = qm -> valence[vC2][2 + j] - 1;
		if ( qm -> valence[vID][1] == qC + 1) {
			EG_getAdjacentToVertexPair(qm, qC + 1, vC2 + 1, vID + 1, edges);
			qm -> valence[vID][1] = edges[1];
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
			printf(" K %d  valence %d = %d  out of bounds", k,vC + 1, qm ->valence[vC][0] + 2);
			qm -> valence[vC][k++] = vID + 1;  // increase k
			for ( v = 0 ; v < aux; ++v)
				if ( qm -> valence[vID][2 + v] == vC2 + 1) qm -> valence[vID][2 + v] = vC + 1;
		}
		else {
			for ( v = 0 ; v < aux; ++v) auxVal[v] = qm -> valence[vID][2 + v];
			--aux;
			qm -> valence[vID]    = EG_reall(qm -> valence[vID],(aux + 2)*sizeof(int));
			qm -> valence[vID][0] = aux;
			for ( i = v = 0 ; v < aux + 1; ++v)
				if ( auxVal[v] != vC2 + 1)   qm -> valence[vID][2 + (i++)] = auxVal[v];
			// ensure that valence is pointing at valid quad
		}
		if ( qm -> valence[vID][1] == qC + 1) {
			stat = EG_getAdjacentToVertexPair(qm , qC + 1, vID + 1, vC2 + 1, edges);
			if ( stat != EGADS_SUCCESS){
				printf(" I Can't find the other quad linking %d and %d near quad %d\n",vID+1, vC2 + 1, qC +1);
				for ( i = 0 ; i < 4; ++i) print_quad_specs(qm, qm->quadAdj[4*qC + i]);
				EG_free(star -> verts);
				EG_free(star -> quads);
				EG_free(star);
				return EGADS_INDEXERR;
			}
		}
	}
	if ( k != qm -> valence[vC][0] + 2 ) {
		printf(" I MESSED UP TRANSFERING VALENCIES!!!\n");
		printf(" k %d  LENGTH VAL %d\n",k,qm -> valence[vC][0] + 2);
		EG_free(star -> verts);
		EG_free(star -> quads);
		EG_free(star);
		return EGADS_INDEXERR;
	}
	// Eliminate vertex vC2 from all the quads to which it belongs to
	for ( q = 0 ; q < star -> nQ; ++q) {
		for ( i = 0 ; i < 4; ++ i)
			if ( qm -> quadIdx[4*(star -> quads[q] -1) + i] == vC2 + 1) qm -> quadIdx[4*(star -> quads[q] -1) + i] = vC +  1;
	}
#ifdef DEBUG
	for ( i = 0 ; i < qm -> valence[vC][0]; ++i) {
		printf(" VALENCE %d = %d  LETS CHECK THAT IT POINTS BACK \t",i, qm -> valence[vC][2 + i ]);
		vID = qm -> valence[vC][2 + i] - 1;
		for ( j = 0 ; j < qm -> valence[vID][0]; ++j) {
			if ( qm -> valence[vID][2 + j] == vC2 + 1) {
				printf(" VERTEX %d is still pointing at %d!!! We have just removed it\n",vID + 1, vC2 + 1);
			}
			else if ( qm -> valence[vID][2 + j] == vC + 1) {
				printf(" correct, valence (%d) = %d\n", j, qm -> valence[vID][2 + j]);
			}
		}
	}
#endif
	// Eliminate vertex quadID:
	for ( i = 0 ; i < 4; ++ i) {
		qm -> quadIdx[4*qC + i ] = -1;
		qm -> quadAdj[4*qC + i ] = -2;
	}
	// delete vertex vC2
	qm -> valence[vC2]    = EG_reall( qm ->valence[vC2], sizeof(int));
	qm -> valence[vC2][0] = -1;
	qm -> vType  [vC2]    = -2; // -2 = removed
	EG_free(star -> verts);
	EG_free(star -> quads);
	EG_free(star);
	return EGADS_SUCCESS;
}


static int
EG_removeDoublet(quadMap *qm, int vID) {
	int poly[2*4],polyOrd[2*4], i, j, qC, collId, stat;
	vStar *star;
	if (qm -> vType[vID - 1] == 0 && qm -> valence[vID - 1][0] == 2) return EGADS_SUCCESS;
	if (qm -> vType[vID - 1]  > 0 && qm -> valence[vID - 1][0] == 2) return EGADS_GEOMERR;
	if (                             qm -> valence[vID - 1][0] != 2) return EGADS_SUCCESS;
	qC = qm -> valence[vID - 1][1] - 1;
	print_quad_specs(qm, qC + 1);
	for ( i = 0 ; i < 4; ++i) {
		poly[2*i    ] = qm ->quadIdx[4*qC + i];
		poly[2*i + 1] = qm ->valence[poly[2*i] - 1][0];
		if (poly[2*i] == vID ) collId = i;
	}
	for ( i = 0 ; i < 4; ++i) {
		polyOrd[2*i    ] = poly[2*( (collId + i)%4)    ];
		polyOrd[2*i + 1] = poly[2*( (collId + i)%4) + 1];
	}
	stat = EG_mergeVertices(qm, qC + 1, polyOrd, 0);
	if ( stat != EGADS_SUCCESS) return stat;
#ifdef DEBUG
	printf(" IN DOUBLET WE HAVE COLLAPSED VERTEX %d to %d AND REMOVED QUAD %d\n", polyOrd[0], polyOrd[2*2], qC + 1);
#endif
	// Check that it has left valid elements at each old link
	for ( j = 0 ; j < 2; ++j) {
		stat = EG_buildStar(qm, &star, polyOrd[2*(2*j + 1)]);
		if ( stat != EGADS_SUCCESS)  return stat;
		for ( i = 0 ; i  < star->nQ;++i) {
			stat = checkInvalidElement(qm, star ->quads[i], star ->verts[0], 1);
			if ( stat != EGADS_SUCCESS) {
				EG_free(star ->verts);
				EG_free(star ->quads);
				EG_free(star);
				return stat;
			}
		}
		stat = EG_removeDoublet(qm, polyOrd[2*(2*j + 1)]);
		if ( stat != EGADS_SUCCESS) {
			EG_free(star ->verts);
			EG_free(star ->quads);
			EG_free(star);
			return stat;
		}
	}
	EG_free(star ->verts);
	EG_free(star ->quads);
	EG_free(star);
	if (qm -> remQuads[0] != qm -> remVerts[0] ){
		printf(" I HAVE REMOVED QUADS AND %d REMOVED VERTICES %d !!! THIS IS A DISASTER!\n",qm -> remQuads[0], qm -> remVerts[0]);
		exit(1);
	}
	qm -> remQuads[0]++; qm -> remVerts[0]++;
	qm -> remQuads[qm -> remQuads[0]] = qC  + 1; ;
	qm -> remVerts[qm -> remVerts[0]] = polyOrd[0];
	return EGADS_SUCCESS;
}






static int EG_collapseAndSplit(quadMap *qm, int qID) {
	int *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, *locVtype, **locVal = NULL ;
	int i, poly[2*4], polyOrd[2*4], vID, idCollapse, idSplit, stat;
	double *uv0;
	for ( i = 0 ;  i < 4; ++i) {
		vID = qm -> quadIdx[4*(qID - 1) + i ] - 1;
		if (qm -> vType[vID] != -1) return EGADS_SUCCESS; //Operation valid only at interior quads
		poly[2*i]     = vID + 1;
		poly[2*i + 1] = qm -> valence[vID][0];
	}
	// Look for pair (3,5)-(4,4)
	if ( (poly[2*0 + 1] * poly[2*2 + 1] == 15) && (poly[2*1 + 1] * poly[2*3 + 1] == 16))
		idCollapse = 0;
	else if ( (poly[2*0 + 1] * poly[2*2 + 1] == 16) && (poly[2*1 + 1] * poly[2*3 + 1] == 15))
		idCollapse = 1;
	else return EGADS_SUCCESS;
	for ( i = 0 ; i < 4; ++i) {
		polyOrd[2*i    ] = poly[2*( (idCollapse + i)%4)    ];
		polyOrd[2*i + 1] = poly[2*( (idCollapse + i)%4) + 1];
	}
	stat = backupQm ( qm, &locIdx, &locAdj, &locVtype, &locVal, &uv0);
	if ( stat != EGADS_SUCCESS) return stat;
	stat = EG_mergeVertices(qm, qID, polyOrd, 1);
	if ( stat != EGADS_SUCCESS) {
#ifdef DEBUG
		printf(" Inside EG_collapseAndSplit, EG_mergeVertices = %d\n", stat);
#endif
		stat = resetQm(qm, locIdx, locAdj, locVtype, locVal, uv0);
		goto cleanup;
	}
	vID  = poly[2*(idCollapse + 2)];
	stat = EG_vertexSplit(qm, vID , 1);
	if ( stat != EGADS_SUCCESS) {
#ifdef DEBUG
		printf(" Inside EG_collapseAndSplit, EG_vertexSplit = %d\n", stat);
#endif
		stat = resetQm(qm, locIdx, locAdj, locVtype, locVal, uv0);
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
	int        i, j, k, q, v, qC, cID, vC, vM, vP, rV, id0, cp, vID, vVAL, addVal;
	int        aux, auxID, quadID, stat, val, vCval,  qCount = 0, vCommon, updatedQuad, collapse;
	int        adjQ, nMove, movingVerts[4], auxVal[100];  // local information about vertices and valence
	static int coll[2][4] = {{0, 2, 1, 3},{1, 3, 2, 0}}; // for collapsing edges
	double     uv[2], vPos[18], *uv0 = NULL;
	int        *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, *locVtype, **locVal = NULL, poly[2*4], polyOrd[2*4], vIrr[4];
	vStar  *star;
	if ( qm -> vType[sID -1] != -1 ) return EGADS_SUCCESS;
	// BACKUP: IN CASE OF INVALID COLLAPSE, RETRIEVE ORIGINAL MESH
	stat = EG_buildStar(qm, &star, sID);
	if ( stat != EGADS_SUCCESS) return stat;
	stat = backupQm ( qm, &locIdx, &locAdj, &locVtype, &locVal, &uv0);
	if ( stat != EGADS_SUCCESS) goto cleanup;
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
			if ( qm -> vType[poly[2*i] - 1] != -1 ) collapse = - 1;
		}
		if ( collapse != - 1) {
			EG_vertexRegularity(qm, poly, vIrr, 4);
			for ( cp = 0 ; cp < 2; ++cp) {
				if ( (vIrr[coll[cp][0]] + vIrr[coll[cp][1]] == -2 ) && (vIrr[coll[cp][2]] + vIrr[coll[cp][3]] > 0 ) )
				{
					qC = star -> quads[q] - 1;
					collapse = 1;
				}
				else if ( (vIrr[coll[cp][0]] + vIrr[coll[cp][1]] < 0 ) && (vIrr[coll[cp][2]] * vIrr[coll[cp][3]] > 0 ) )
				{
					qC = star -> quads[q] - 1;
					collapse = 1;
				}
				if (collapse) break;
			}
		}
		if (collapse) break;
	}
	if ( collapse != 1) goto cleanup;
	for ( i = 0 ; i < 4; ++i) {
		polyOrd[2*i    ] = poly[2*( (coll[cp][0] + i)%4)    ];
		polyOrd[2*i + 1] = poly[2*( (coll[cp][0] + i)%4) + 1];
	}
	stat      = EG_mergeVertices(qm, qC + 1, polyOrd,  1 );
	if (stat != EGADS_SUCCESS) goto cleanup;
	EG_free(star -> quads);
	EG_free(star -> verts);
	EG_free(star);
	star = NULL;
	EG_buildStar(qm, &star, polyOrd[2*2]);
#ifdef DEBUG
	printf(" IN DOUBLET WE HAVE COLLAPSED VERTEX %d to %d AND REMOVED QUAD %d\n",polyOrd[0], polyOrd[2*2], qC + 1);
#endif
	for ( j = 0 ; j < star -> nV; ++j) {
		auxID =  - 1;
		if ( star -> verts[j] > 0) {
			auxID = star -> verts[j] - 1;
			if (qm -> valence[auxID][0] == 2 && qm -> vType[auxID] == -1) {
#ifdef DEBUG
				printf(" IN EDGE SWAP. WE HAVE PRODUCED A DOUBLET AT %d\n",auxID + 1);
#endif
				stat = EG_removeDoublet(qm, auxID +1);
				if ( stat != EGADS_SUCCESS) return stat;
				print_mesh(qm, &FCOUNT);
			}
		}
	}
	// Optimize around quads
	nMove = 0 ;
	for ( i = 1 ; i < 4; ++i)
		if (  qm -> vType[poly[2*i] -1] == -1 ) movingVerts[nMove++] = polyOrd[2*i];
	stat = optimize_angles(qm, nMove, movingVerts);
	if ( stat == EGADS_SUCCESS) {
		++TOTCOLLAPSES;
		if (qm -> remQuads[0] != qm -> remVerts[0] ){
			printf(" I HAVE REMOVED QUADS AND %d REMOVED VERTICES %d !!! THIS IS A DISASTER!\n",qm -> remQuads[0], qm -> remVerts[0]);
		}
		qm -> remQuads[0]++; qm -> remVerts[0]++;
		qm -> remQuads[qm -> remQuads[0]] = qC  + 1; ;
		qm -> remVerts[qm -> remVerts[0]] = polyOrd[0];
		printf(" OVERALL %d  %d", qm -> remVerts[0], qm -> remQuads[0]);
	}
	else stat = resetQm(qm, locIdx, locAdj, locVtype, locVal,  uv0);
	cleanup:
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
swappingOperation(quadMap *qm, vStar *star, int poly[], int swapIdx, int qID[], int pullVertex) {
	int stat, i, j, k, aux, auxVAL, auxID, pos, it, adjQ, areaType, area, centreValency;
	int v[4], addV[4], v03[2], OK;
	int perm[6] = {0, 1, 2, 3, 4, 5}, remV[4] = {0, 3, 3, 0}, qL[6]   = { 3, 0, 1, 2, 3, 0};
	int Q1[4], Q2[4], idx1[4], idx2[4], edges[6], q0ID[6], q1ID[6], valence[100];
	addV[0] = swapIdx; addV[1] = swapIdx + 3;
	addV[2] = addV[1]; addV[3] = addV[0];
	qID[0]--; qID[1]--;
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
		if ( auxVAL != qm -> valence[auxID][0]) {
			printf(" ALERT! POLY DOESNT COINCIDE WITH VALENCE POLY%d VAL %d \n", auxVAL,qm -> valence[auxID][0]);
		}
		for ( i = 0 ; i < auxVAL; ++i) valence[i] = qm -> valence[auxID][2 + i];
		++auxVAL;
		qm -> valence[auxID]    = EG_reall(qm -> valence[auxID ], (auxVAL + 2)*sizeof(int));
		qm -> valence[auxID][0] = auxVAL;
		for ( i = 0 ; i < auxVAL ; ++i) qm -> valence[auxID][i + 2] = valence[i];
		qm ->valence[auxID ][2 + auxVAL - 1] = poly[2*addV[2*j + 1]];
	}
	// FIND SIDE FROM WICH TO PULL
	for ( k = 0 ; k < 4; ++k){
		if (qm -> quadIdx[4*qID[0] + k] == v03[0] ) break;
		else if (qm -> quadIdx[4*qID[0] + k] == v03[1] ) {
			aux    = qID[0];
			qID[0] = qID[1];
			qID[1] = aux;
			break;
		}
	}
	// Finally: Check if we have created a doublet
	for ( j = 0 ; j < star -> nV; ++j) {
		auxID =  - 1;
		if ( star -> verts[j] > 0) {
			auxID = star -> verts[j] - 1;
			if (qm -> valence[auxID][0] == 2 && qm -> vType[auxID] == -1) {
#ifdef DEBUG
				printf(" IN EDGE SWAP. WE HAVE PRODUCED A DOUBLET AT %d\n",auxID);
#endif
				stat = EG_removeDoublet(qm, auxID +1);
				if ( stat != EGADS_SUCCESS) return stat;
				print_mesh(qm, &FCOUNT);
			}
		}
	}
	for ( j = 0 ; j < 2 ; ++ j) {
		if ( (checkInvalidElement(qm, qID[j] + 1, v03[j], pullVertex) ) != EGADS_SUCCESS) return EGADS_GEOMERR;
	}
	return EGADS_SUCCESS;
}


static int
EG_edgeSwap(quadMap *qm, int sID) {
	int stat, i, j, k, aux, vL, movingVerts[4], nMove;
	int v[4*2], vIrr[6], qID[2], link, swapped;
	int *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, *locVtype, **locVal = NULL, poly[2*6];
	double fOpti, *uv0 = NULL;
	vStar *star;
	stat = EG_buildStar(qm, &star, sID);
#ifdef DEBUG
	printf(" EG_edgeSwap around vertex %d\n",sID);
#endif
	if ( stat != EGADS_SUCCESS) return stat;
	if ( (qm -> vType[star -> verts[0] - 1] == 0 ) && qm -> valence[star -> verts[0] - 1][0] == 2) return EGADS_SUCCESS;
	if ( (qm -> vType[star -> verts[0] - 1]  > 0 ) && qm -> valence[star -> verts[0] - 1][0] == 3) return EGADS_SUCCESS;
	stat = backupQm ( qm, &locIdx, &locAdj, &locVtype, &locVal, &uv0);
	if ( stat != EGADS_SUCCESS) goto cleanup;
	poly[0]       = star->verts[0];
	poly[1]       = qm -> valence[poly[0] - 1][0];
	swapPtr       = (int *) EG_alloc(poly[1] *sizeof(int));
	if (swapPtr == NULL) {
		stat = EGADS_MALLOC;
		goto cleanup;
	}
	swapped = 0;
	for ( link = 0; link < poly[1]; ++link) {
		vL     = (2*link + 1);
		qID[0] = star -> quads[link];
		if ( link + 1 == star -> nQ ) qID[1] = star -> quads[0];
		else qID[1] = star->quads[link + 1];
		for ( i = 0 ; i  < 5; ++i) {
			aux = vL + i;
			if ( aux >= star ->nV) aux -= (star->nV -1);
			poly[2*(i+1)   ] = star ->verts[aux];
			if ( poly[2*(i+1)] < 0) { //This is a virtual quad used when the star centre is at the boundary. Skip
				i  = 5;
				vL = -1;
				poly[2*(i+1) +1] = -1;
			}
			else poly[2*(i+1) +1] = qm   ->valence[poly[2*(i+1)] -1][0];
		}
		swapPtr[link] = -1;
		if ( vL != -1 ) {
			EG_vertexRegularity(qm, poly, vIrr, 6);
#ifdef DEBUGG
			for ( i = 0 ; i < 6; ++i) printf(" VERT %d = %d \n",poly[2*i], poly[2*i+1]);
#endif
			// CHECK IF WE SHOULD SWAP
			if (qm -> vType[ poly[0] -1] == 0) { // this is a domain corner. SWAP FORCED
				if ( (vIrr[1] + vIrr[4]) < (vIrr[2] + vIrr[5]) )  swapPtr[link] = 1;
				else swapPtr[link] = 2;
			}
			else if (qm -> vType[ poly[0] -1] != -1) {
				if ( (vIrr[1] != 100) && (vIrr[4] != 100)  //not domain vertices
						&& (vIrr[1] + vIrr[4] <  vIrr[2] + vIrr[5])  //higher valence that the other swapping option
						&& (vIrr[1] + vIrr[4] <= vIrr[0] + vIrr[3]) ) // improves or preserves current valence
					swapPtr[link] = 1;
				else if ( (vIrr[2] != 100) && (vIrr[5] != 100 )
						&& (vIrr[2] + vIrr[5] <= vIrr[0] + vIrr[3]) )
					swapPtr[link] = 2;
			} else {
				if   ( (vIrr[0] + vIrr[3]) >  0) {
					if      ( (vIrr[1] + vIrr[4]) == -2)  swapPtr[link] = 1;
					else if ( (vIrr[2] + vIrr[5]) == -2)  swapPtr[link] = 2;
					else if (( vIrr[0] + vIrr[3]) >  5)
					{ // two high valencies
						if      ( (vIrr[1] == -1 && vIrr[4] == 0) ||
								(  vIrr[4] == -1 && vIrr[1] == 0) ) swapPtr[link] = 1;
						else if ( (vIrr[2] == -1 && vIrr[5] == 0) ||
								(  vIrr[5] == -1 && vIrr[2] == 0) ) swapPtr[link] = 2;
					}
				}
			}
			if ( swapPtr[link] != -1 ) {
#ifdef DEBUGG
				printf(" WE ARE GOING TO SWAP QUADS %d %d  USING LINK %d %d POLY SPECS: \n",qID[0], qID[1],swapPtr[link], swapPtr[link] + 3);
				printf(" EXCHANGE VALENCE V(%d) = %d V(%d) = %d  AND V(%d) = %d V(%d) = %d \n", poly[0], poly[1], poly[2*3], poly[2*3 + 1],
						poly[2*swapPtr[link]], poly[2*swapPtr[link] + 1], poly[2*(swapPtr[link] + 3)], poly[2*(swapPtr[link] + 3) + 1]);
				for ( i = 0 ; i < 6; ++i) printf(" VERT %d = %d \n",poly[2*i], poly[2*i+1]);
#endif
				stat = swappingOperation(qm, star, poly, swapPtr[link], qID, 0);
#ifdef DEBUGG
				printf(" SWAPPING OPERATION RESULTED IN %d\n",stat);
				for ( i = 0 ; i < 6; ++i) printf(" VERT %d = %d \n",poly[2*i], poly[2*i+1]);
#endif
				if ( stat == EGADS_SUCCESS ) {
					swapped = 1;
					link    = 100;
					++TOTSWAPS;
					break;
				}
				else  {
					stat = resetQm(qm, locIdx, locAdj, locVtype, locVal,  uv0);
					if ( stat != EGADS_SUCCESS) goto cleanup;
				}
			}
		}
		if (swapped ) break;
	}
	if(!swapped) {
		for ( link = 0 ; link < poly[1]; ++link) {
			if (swapPtr[link] != -1) {
				vL     = (2*link + 1);
				for (i = 0 ; i < star -> nQ; ++i ) printf(" STAR Q %d = %d\n",i, star->quads[i]);
				qID[0] = star -> quads[link];
				if ( link + 1 == star -> nQ ) qID[1] = star -> quads[0];
				else qID[1] = star->quads[link + 1];
				for ( i = 0 ; i  < 5; ++i) {
					aux = vL + i;
					if ( aux >= star ->nV) aux -= (star->nV -1);
					poly[2*(i+1)   ] = star->verts[aux   ];
					if ( poly[2*(i+1)] < 0) { //This is a virtual quad used when the star centre is at the boundary. Skip
						printf(" THIS SHOULDN'T HAVE BEEN CONSIDERED A VALID SWAP. BUG\n");
						stat = EGADS_INDEXERR;
						goto cleanup;
					}
					poly[2*(i+1) +1] = qm -> valence[poly[2*(i+1)] - 1][0];
				}
#ifdef DEBUGG
				printf(" WE ARE GOING TO SWAP QUADS %d %d  USING LINK %d %d POLY SPECS: \n",qID[0], qID[1],swapPtr[link], swapPtr[link] + 3);
				printf(" EXCHANGE VALENCE V(%d) = %d V(%d) = %d  AND V(%d) = %d V(%d) = %d \n", poly[0], poly[1], poly[2*3], poly[2*3 + 1],
						poly[2*swapPtr[link]], poly[2*swapPtr[link] + 1], poly[2*(swapPtr[link] + 3)], poly[2*(swapPtr[link] + 3) + 1]);
				for ( i = 0 ; i < 6; ++i) printf(" VERT %d = %d \n",poly[2*i], poly[2*i+1]);
#endif

				stat = swappingOperation( qm, star, poly, swapPtr[link], qID, 1);
				if ( stat == EGADS_SUCCESS ) {
					swapped = 1;
					++TOTSWAPS;
					link = 100;
					break;
				}
				else {
					stat = resetQm(qm, locIdx, locAdj, locVtype, locVal, uv0);
					goto cleanup;
				}
			}
		}
	}
	if ( swapped ) {
#if MULTIOPTI
		nMove = 0 ;
		for ( i = 0 ; i < 6; ++i)
			if ( qm -> vType[poly[2*i] -1] == -1) movingVerts[nMove++] = poly[2*i];
		stat = optimize_angles(qm, nMove, movingVerts);
		if ( stat != EGADS_SUCCESS)
			stat = resetQm(qm, locIdx, locAdj, locVtype, locVal, uv0);
#else
		nMove = 0 ;
		for ( k = 0 ; k < star ->nV; ++k) {
			optiVerts[nMove] = star -> verts[k];
			if (qm -> vType[optiVerts[nMove] -1] == -1 ) {
				++nMove;
			}
		}
		stat = optimize_angles(qm, nMove, optiVerts);
#endif
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
	EG_free(uv0);
	EG_free(swapPtr);
	return stat;
}




static int EG_swapAndSplit(quadMap *qm, int sID ) {
	int    v0, v1, v2, s, i, j, q, qq, sVAL, link, stat, validSwap, validSplit, qIDs[2], adj[2], id;
	int    *swapPtr = NULL, *locIdx = NULL, *locAdj = NULL, *locVtype, **locVal = NULL, vSwap[2*6];
	double *uv0 = NULL;
	vStar  *swapStar, *star;
	stat = EG_buildStar(qm, &star, sID);
	if ( stat != EGADS_SUCCESS) return stat;
	v0   = star -> verts [0] - 1;
	sVAL = qm   -> valence[v0][0];
	if ( sVAL < 4 ) return EGADS_SUCCESS;
	// FIRST: LOOK AT POSSIBLE SWAP
	for ( s = 0 ; s < sVAL; ++s){
		v1      = star -> verts[2*s + 1] - 1;
		qIDs[0] = star -> quads[s]       - 1;
		if ( qm -> valence[v1][0] >= 5 && qm -> vType[v1] == -1) {
			// build a star around the high valence link in order to swap
			stat  = EG_buildStar(qm, &swapStar, v1 + 1);
			// Get adjacent to quad qID in swapStar and NOT in star
			validSwap = 0;
			for ( q = 0 ; q < swapStar -> nQ; q++)
			{
				qIDs[1]   = swapStar -> quads[q] - 1;
				validSwap = 0;
				for ( i = 0 ; i < 4; ++i) {
					if (qm -> quadAdj[4*qIDs[1] + i] == qIDs[0]) {
						validSwap = 1; // qID points at qID[0]
						for ( qq = 0 ; qq < star -> nQ; qq++)
							if ( star -> quads[qq] == qIDs[1] ) // It belongs to star. BAD
								validSwap = 0;
					}
					if (validSwap ) break;
				}
				if ( validSwap ) { // qID is in swapStar and not in star. lets see the valences
					id = 0;
					// Collect vertices from the two swapping quads: qIDs{0,1}
					for ( i = 0 ; i < 4; i++) if ( qm -> quadIdx[4*qIDs[0] + i] == v1 + 1) id = i;
					for ( i = 0 ; i < 4; i++) {
						vSwap[2*i]     = qm -> quadIdx[4*qIDs[0] + (id + i)%2];
						vSwap[2*i + 1] = qm -> valence[vSwap[2*i] - 1][0];
					}
					// Add the remainnig two vertices from second quad
					for ( i = 0 ; i < 4; i++) if ( qm -> quadIdx[4*qIDs[1] + i] == v1 + 1) id = i;
					for ( i = 0 ; i < 2; i++) {
						vSwap[2*i + 4*2    ] = qm -> quadIdx[4*qIDs[0] + (id + 2 + i)%2];
						vSwap[2*i + 4*2 + 1] = qm -> valence[vSwap[2*i + 4*2] - 1][0];
					}
					// try to remove edge 0,3 for edge 1,4
					validSwap = 0;
					if      ( vSwap[3*2 + 1]  > 4 && vSwap[4*2 + 1] <= 4 ) validSwap = 1;
					else if ( vSwap[3*2 + 1] == 4 && vSwap[4*2 + 1] == 3 ) validSwap = 1;
				}
			}
		}
		if (validSwap) break;
	}
	if (validSwap) {
#ifdef DEBUG
		printf(" WE FOUND A VALID SWAP AND SPLIT! \n");
		printf(" FIRS SWAP QUADS %d %d and create edge %d %d\n", qIDs[0], qIDs[1], vSwap[2*1], vSwap[2*5] );
#endif
		// NOW LOOK AT POSSIBLE SPLIT IN STAR
		if ( EG_vertexSplit(qm, star ->verts[0], 0)) {
#ifdef DEBUG
			printf(" WE ALSO FOUND A VALID SPLIT. GREAT!\n");
#endif
			stat = swappingOperation(qm, swapStar, vSwap, 1, qIDs, 1);
			if (stat != EGADS_SUCCESS) {
				EG_free(star->verts);
				EG_free(star->quads);
				EG_free(star);
				return resetQm(qm, locIdx, locAdj, locVtype, locVal, uv0);
			}
			stat = EG_vertexSplit(qm, star ->verts[0], 1);
			if ( stat != EGADS_SUCCESS) {
				EG_free(star->verts);
				EG_free(star->quads);
				EG_free(star);
				return resetQm(qm, locIdx, locAdj, locVtype, locVal, uv0);
			}
		}
	}
	for ( i = 0 ; i < qm -> totVerts; ++i) EG_free(locVal[i]);
	EG_free(locVal);
	EG_free(locIdx);
	EG_free(locAdj);
	EG_free(locVtype);
	EG_free(uv0);

	EG_free(star->verts);
	EG_free(star->quads);
	EG_free(star);
	return EGADS_SUCCESS;
}

static int updateUVs(ego tess, int iFace, ego face, double *uvs, double *xyzs)
{
	int          nSweeps = 100, activity, totActivity, itQ, nPeaks, i, j, k, count, kOK, id0, auxID, auxQ, nseg, newUV, stat, nVars, quadID, prevQuad;
	int          kk,  m, n, iper, len, nTri, optVec[100], nOpt;
	const int    *tris, *tric, *ptype, *pindex;
	double       uvbox[4], *lowerBounds = NULL, *upperBounds = NULL, fOpti, *uvOpt=NULL, zeroplane;
	const double *xyzTess, *uvTess;
	static int   qV[6]    = { 0, 1, 2, 5, 0, 1};
	static int   qLoop[5] = { 0, 1, 2, 3, 0};
	quadMap      *qm = NULL;
	FILE *fop;
	fop = fopen("OPERATIONS","w");
	if (fop == NULL) return EGADS_MALLOC;
	fclose(fop);
	stat = EG_getTessFace(tess, iFace, &len,
			&xyzTess, &uvTess, &ptype, &pindex, &nTri, &tris, &tric);
	if (stat != EGADS_SUCCESS) {
		return stat;
	}
	qm = (quadMap*) EG_alloc(sizeof(quadMap));
	if ( qm == NULL ) {
		stat =  EGADS_MALLOC;
		goto cleanup;
	}
	qm -> face     = face;
	qm -> totVerts = len;
	qm -> nQ       = nTri/2;
	qm -> xyzs     = (double*) EG_alloc(3*len    *sizeof(double));
	qm -> oriXYZ   = (double*) EG_alloc(3*len    *sizeof(double));
	qm -> uvs      = (double*) EG_alloc(2*len    *sizeof(double));
	qm -> vType    = (int*)    EG_alloc(  len    *sizeof(   int));
	qm -> remVerts = (int*)    EG_alloc(  len    *sizeof(   int));
	qm -> quadIdx  = (int*)    EG_alloc(4*qm ->nQ*sizeof(   int));
	qm -> quadAdj  = (int*)    EG_alloc(4*qm ->nQ*sizeof(   int));
	qm -> vLengths = (double*) EG_alloc(2*len    *sizeof(double));
	qm -> vCurvs   = (double*) EG_alloc(8*len    *sizeof(double));
	qm -> remQuads = (int*)    EG_alloc(qm -> nQ *sizeof(   int));
	qm -> star     = (vStar**) EG_alloc(     len *sizeof(vStar*));
	qm -> nS       = 0;
	if ( qm ->quadIdx == NULL || qm ->quadAdj == NULL || qm -> xyzs == NULL || qm -> oriXYZ == NULL ||
			qm ->uvs == NULL || qm ->vType == NULL || qm ->vLengths == NULL || qm -> vCurvs == NULL ||
			qm -> remQuads == NULL || qm -> remVerts == NULL || qm -> star == NULL) {
		stat = EGADS_MALLOC;
		goto cleanup;
	}
	qm -> remQuads[0] = 0 ; qm -> remVerts[0] = 0;
	zeroplane = 0.0;
	for (j = 0; j < len; j++) {
		qm -> xyzs[3*j] = xyzTess[3*j]; qm ->xyzs[3*j + 1] = xyzTess[3*j +1];qm ->xyzs[3*j +2] = xyzTess[3*j +2];
		qm -> uvs [2*j] = uvTess [2*j]; qm ->uvs  [2*j + 1] = uvTess [2*j +1];
		qm -> vType [j] = ptype[j];
		zeroplane += qm -> xyzs[3*j + 2];
	}
	if ( fabs(zeroplane) < DEPS){
		PLANESURFACE = 1;
		printf(" WE HAVE A 2D SURFACE ALONG PLANE z=0\n");
	}

	for (j = 0; j < qm -> nQ; j++)
		for ( k = 0; k < 4; ++k)  qm -> quadIdx[4*j + k ] = tris[6*j + qV[k+1]];
	qm -> valence = (int **) EG_alloc(len*sizeof(int*));
	if ( qm -> valence == NULL ) {
		stat = EGADS_MALLOC;
		goto cleanup;
	}
	for ( j = 0; j < len; ++j) {
		qm -> valence[j] = (int *) EG_alloc((2 + 100)*sizeof(int));
		if (qm -> valence[j] == NULL ) {
			stat = EGADS_MALLOC;
			goto cleanup;
		}
		qm -> valence[j][0] = 0;
	}
	for ( j = 0; j < qm -> nQ; j++) {
		kk  = 0;
		kOK = 0;
		for ( i = 0 ; i < qm -> nQ; ++i) {
			if ( i == j || i == -1) ++i;
			for ( k = 0 ; k < 4; ++k ) {
				if ( i == j || i == -1) ++i;
				if ( (qm -> quadIdx[4*j + qLoop[kk    ]] == qm -> quadIdx[4*i + qLoop[k]] || qm -> quadIdx[4*j + qLoop[kk    ]] == qm -> quadIdx[4*i + qLoop[k + 1]] )&&
						( qm -> quadIdx[4*j + qLoop[kk + 1]] == qm -> quadIdx[4*i + qLoop[k]] || qm -> quadIdx[4*j + qLoop[kk + 1]] == qm -> quadIdx[4*i + qLoop[k + 1]]) ) {
					qm -> quadAdj[4*j + kk] = i+1;
					++kk;
					i   = -1;
					kOK = 1;
					k   = 4;
					if (kk == 4)  i = qm -> nQ;
				}
			}
			if ( (kOK == 0) && (i == qm -> nQ -1) ){
				qm -> quadAdj[4*j + kk] = -1;
				i                       = -1;
				++kk;
				if (kk == 4) i = qm -> nQ;
			}
			else  kOK = 0 ;
		}
	}
	for ( j = 0 ; j < qm ->nQ; ++j) {
		for ( k = 0 ; k < 4; ++k) {
			if ( qm -> quadAdj[4*j +k] > j +1  || qm -> quadAdj[4*j +k] == -1) {
				auxID = qm -> quadIdx [4*j + qLoop[k]] -1;
				qm   -> valence[auxID][1] = j + 1;
				qm   -> valence[auxID][2 + qm -> valence[auxID][0] ] = qm -> quadIdx[4*j+ qLoop[k+1]];
				++qm -> valence[auxID][0];
				auxID = qm -> quadIdx [4*j + qLoop[k+1]] -1;
				qm   -> valence[auxID][1] = j + 1;
				qm   -> valence[auxID][2 + qm -> valence[auxID][0] ] = qm -> quadIdx[4*j+ qLoop[k]];
				++qm -> valence[auxID][0];
			}
		}
	}
	// GET RANGE FOR EACH POINT

	for ( j = 0 ; j < qm ->totVerts; ++j) {
		if ( qm -> vType[j] == -1) {
			// LINKS BY DIRECTION
			double norm, maxN = 0.0, minN = 100.0;;
			int maxID, minID;
			for ( i = 0 ; i < qm ->valence[j][0]; ++i) {
				auxID = qm -> valence[j][2 + i] - 1;
				norm  = 0.0;
				for ( k = 0 ; k < 3; ++ k) norm += (qm-> xyzs[3*j + k] - qm ->xyzs[3*auxID +k]) * (qm-> xyzs[3*j + k] - qm ->xyzs[3*auxID +k]);
				norm = sqrt(norm);
				if ( norm > maxN ){
					maxN = norm; maxID = auxID;
				}
				if ( norm < minN ){
					minN = norm; minID = auxID;
				}
			}
			qm -> vLengths[2*j]     = maxN*1.1;
			qm -> vLengths[2*j + 1] = minN*0.5;
		}
		stat = EG_curvature (qm -> face, &qm -> uvs[2*j], &qm -> vCurvs[8*j]);
		if ( stat != EGADS_SUCCESS) goto cleanup;

		/*(qm ->xyzs[3*maxID    ] - qm-> xyzs[3*j    ]);
      qm -> vLengths[6*j + 1] = (qm ->xyzs[3*maxID + 1] - qm-> xyzs[3*j + 1]);
      qm -> vLengths[6*j + 2] = (qm ->xyzs[3*maxID + 2] - qm-> xyzs[3*j + 2]);
      for ( i = 0 ; i < qm ->valence[j][0]; ++i) {
        auxID = qm -> valence[j][2 + i] - 1;
        if ( auxID != maxID) {
          double dot = 0.0;
        for ( k = 0 ; k < 3; ++ k) norm += (qm-> xyzs[3*j + k] - qm ->xyzs[3*auxID +k]) * (qm-> xyzs[3*j + k] - qm ->xyzs[3*auxID +k]);
        norm = sqrt(norm);
        if ( norm > maxN ){
          maxN = norm; maxID = auxID;
        }
        if ( norm < minN){
          minN = norm; minID = auxID;
        }
      }
      qm -> vLengths[6*j    ] = (qm ->xyzs[3*maxID    ] - qm-> xyzs[3*j    ]);
      qm -> vLengths[6*j + 2] = (qm ->xyzs[3*maxID + 2] - qm-> xyzs[3*j + 2]);
      qm -> vLengths[6*j + 3] = (qm ->xyzs[3*minID    ] - qm-> xyzs[3*j    ]);
      qm -> vLengths[6*j + 4] = (qm ->xyzs[3*minID + 1] - qm-> xyzs[3*j + 1]);
      qm -> vLengths[6*j + 5] = (qm ->xyzs[3*minID + 2] - qm-> xyzs[3*j + 2]);
    }*/
		qm -> oriXYZ[3*j    ] = qm -> xyzs[3*j    ];
		qm -> oriXYZ[3*j + 1] = qm -> xyzs[3*j + 1];
		qm -> oriXYZ[3*j + 2] = qm -> xyzs[3*j + 2];
	}
#ifdef DEBUG
	print_mesh(qm, &FCOUNT);
#endif
	count = 0;
	id0   = 0;
	printMeshStats(qm,0);
	totActivity = 1; k = 0 ;
	while ( k < nSweeps || totActivity > 0) {
		totActivity = 0;
		// DOMAIN CORNERS
		//
		printf("------------------------\n  DOMAIN CORNERS ROUND %d \n ----------------------\n",k);
		stat = EG_meshOperator(qm, TRYSWAPPING,      0, &activity);
		totActivity += activity;
		printf(" EG MESH OPERATOR FOR SWAPPING ACTIVITY %d  STAT %d\n", activity,stat);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		// DOMAIN EDGES
		//
		printf("------------------------\n  DOMAIN EDGES  ROUND %d \n ----------------------\n",k);
		stat = EG_meshOperator(qm, TRYSWAPPING,      1, &activity);
		totActivity += activity;
		printf(" EG MESH OPERATOR FOR SWAPPING ACTIVITY %d  STAT %d\n", activity,stat);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		// INTERIOR
		//
		printf("------------------------\n  INTERIOR VERTICES ROUND %d \n ----------------------\n",k);
		stat = EG_meshOperator(qm, TRYSWAPPING,     -1, &activity);
		totActivity += activity;
		printf(" EG MESH OPERATOR FOR SWAPPING ACTIVITY %d  STAT %d\n", activity,stat);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		stat = EG_meshOperator(qm, TRYCOLLAPSING,   -1, &activity);
		totActivity += activity;
		printf(" EG MESH OPERATOR FOR COLLAPSING ACTIVITY %d  STAT %d\n", activity,stat);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		stat = EG_meshOperator(qm, TRYSPLITTING, -1, &activity);
		totActivity += activity;
		printf(" EG MESH OPERATOR FOR SPLITTING ACTIVITY %d  STAT %d\n", activity,stat);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		fop = fopen("OPERATIONS","a");
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
		printf(" ROUND %d OVERALL SWAPS  %d, COLLAPSES %d  SPLITS %d\n", k + 1, OVERALLSWAPS, OVERALLCOLLAPSES,OVERALLSPLITS);
		printf("-----------------------------------------------\n");
		++k;
		printMeshStats(qm, k + 1);
	}
	printf(" PERFORM GLOBAL OPTIMIZATION\n");
	printMeshStats(qm, k + 1);
	for ( nOpt = i = 0 ; i < len; ++i)
		if ( qm -> vType[i] == -1)  optVec[nOpt++] = i+1;
	stat = optimize_angles(qm, nOpt, optVec);
	printf(" STAT FROM GLOBAL OPTIMIZER %d \n",stat);
	print_mesh(qm, &FCOUNT);
	cleanup:
	EG_free(qm -> quadIdx);
	EG_free(qm -> quadAdj);
	EG_free(qm -> uvs);
	EG_free(qm -> xyzs);
	for ( i = 0 ; i < len; ++i) {
		EG_free(qm ->valence[i]);
		EG_free(qm ->star[i]);
	}
	EG_free(qm -> valence);
	EG_free(qm -> vType);
	EG_free(qm -> star);
	EG_free(qm);
	EG_free(tess);
	return stat;
}



double NLOPT_objEqualAnglesTess(unsigned n, const double *uvs, double *grad, void *inQM)
{
	int    stat, i, j, k, vFix, vA, vB, aux, oID, peaks;
	double thetaOpt, thetas[4];
	double vXYZ[18], oldXYZ[18], fQuad, fStar, fTot,  penalty, quadArea;
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
		qm->uvs  [2*oID]     = uvs[2*i];
		qm->uvs  [2*oID + 1] = uvs[2*i + 1];
	}
	IOSTATUS = 0;
	if (OPTCOUNT%NOUT == 0 ) {
		printf(" ITERATION    %d  \n",OPTCOUNT);
		IOSTATUS = 1;
		printf(" IO STATUS %d\n",IOSTATUS);
	}
	fprintf(filOpti,"#   CALL UPDATE OPTIMIZATION ITERATION %d\n",OPTCOUNT);
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
		fStar = 0.0;
		for (j = 0 ; j < qm->star[i] ->nQ; ++j) {
			stat = computePenalty(qm, i, j, &penalty, thetas, &quadArea);
			if (OPTCOUNT%NOUT == 0 )
				printf(" #########   NLOPT FUNCTION FOR STAR %d = %d AND QUAD %d: PENALTY %lf AREA SIZE %lf  WITH MAX ANGLE QUAD %lf   ANGLE AT VERTEX THETA = %lf  OPTITHETA %lf   ACCUMULATED %lf \n",
						i+1, qm -> star[i]->verts[0], qm ->star[i]->quads[j],penalty, quadArea, thetas[1],thetas[0], thetaOpt,fTot);
			if ( stat != EGADS_SUCCESS)  {
				printf(" NLOPT FAILED TO DO computePenalty. EXIT PROGRAM \n");
				exit (1);
			}
			if(isnan(thetas[0])) {
				printf(" ROSE AT CENTRE %d  HAS A NAN VALUE FOR ANGLE AT QUAD %d theta %lf\n",qm ->star[i] -> verts[0], qm ->star[i]->quads[j], thetas[0]);
				return EGADS_EMPTY;
			}
			fQuad  = (thetaOpt - thetas[0])*(thetaOpt - thetas[0]);
			fStar += fQuad + penalty;
			fprintf(filOpti,"%d %d\t %lf %lf %lf\t %lf %lf %lf %lf %lf %lf\n",qm -> star[i] -> quads[j],  qm -> star[i]->verts[0] ,thetas[0], fabs(thetaOpt - thetas[0]),thetaOpt, fQuad,fStar,fTot, thetas[1], penalty,quadArea);
		}
		if ( fStar < OPTITOL ) fStar  = 0.0;
		fTot += fStar;
		fprintf(filOpti,"# TOTAL ITERATION FUNCTION %lf  %lf\n",fTot,fStar);
	}
	if ( OPTCOUNT%NOUTMESH == 0) print_mesh(qm, &FCOUNT);
	return fTot;
}


static int optimize_angles(quadMap *qm, int nV, int *vID)
{
	int    i, j, k, stat, iper, peaks, aux, area, stat2;
	double XYZ[18], uvbox[4], *lowerBounds = NULL, *upperBounds = NULL, fOpti, *uvOpt=NULL;
	double epsREL = 1.e-08, *epsABS = NULL, *uv0 = NULL, angles[4], vec, quadArea;
	double vals[18], normal[4];
	clock_t t_start, t_end;
	qm -> nS     = nV;
	//qm ->star   = (vStar**) EG_alloc(  nV*sizeof(vStar*));
	uvOpt       = (double*) EG_alloc(2*nV*sizeof(double));
	uv0         = (double*) EG_alloc(2*nV*sizeof(double));
	lowerBounds = (double*) EG_alloc(2*nV*sizeof(double));
	upperBounds = (double*) EG_alloc(2*nV*sizeof(double));
	stat = EG_getRange(qm->face, uvbox, &iper);
	if (lowerBounds == NULL || upperBounds == NULL || uvOpt == NULL || uv0 == NULL || qm ->star == NULL)  return EGADS_MALLOC;
#ifdef DEBUG
	for ( i = 0; i < nV; ++i)
		printf(" OPTIMIZE ANGLES %d  %d\n",i,vID[i]);
#endif
	for ( i = 0; i < nV; ++i){
		stat = EG_buildStar(qm, &qm ->star[i], vID[i]);
		if ( stat != EGADS_SUCCESS) goto cleanup;
		uvOpt      [2*i    ] = qm ->uvs[2*(vID[i] -1)    ];
		uvOpt      [2*i + 1] = qm ->uvs[2*(vID[i] -1) + 1];
		uv0        [2*i    ] = qm ->uvs[2*(vID[i] -1)    ];
		uv0        [2*i + 1] = qm ->uvs[2*(vID[i] -1) + 1];
		lowerBounds[2*i    ] = uvbox[0];
		lowerBounds[2*i + 1] = uvbox[2];
		upperBounds[2*i    ] = uvbox[1];
		upperBounds[2*i + 1] = uvbox[3];
		for ( j = 0 ; j < qm -> star[i] -> nQ; ++j) {
			if ( qm -> star[i] -> quads[j] == -1 ) printf(" I SHOULDNT BE OPTIMIZING A BOUNDARY VERTEX\n");
			else
				stat = checkInvalidElement(qm, qm -> star[i] -> quads[j], qm -> star[i] -> verts[0], 1);
			if ( stat != EGADS_SUCCESS) goto cleanup;
			printf("------------------------------------------ QUAD %d STAT %d\n", qm -> star[i] -> quads[j], stat);
		}
		printf("\n");
	}
	OPTCOUNT = 0;
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
	nlopt_set_stopval(opt,  OPTITOL);
	nlopt_set_maxeval(opt, NLOPTMAXEVAL);
	//nlopt_set_ftol_rel(opt, epsREL);
	nlopt_set_xtol_rel(opt, epsREL);
	//nlopt_set_ftol_abs(opt, epsABS);
	nlopt_set_xtol_abs(opt, epsABS);
	//nlopt_set_maxtime(opt, 120.0*2.0*(double)nV);
	snprintf(optiName, sizeof(char) * 32, "OPTIMIZERCONVERGENCE_%d.txt",CONVC);
	OPTCOUNT = 0 ;
	printf(" OPEN FILE %s\n",optiName);
	filOpti = fopen (optiName,"w");
	t_start = clock();
	nlopt_result res = nlopt_optimize(opt, uvOpt, &fOpti);
	for ( i = 0 ; i < qm->nS; ++i) {
		print_starFile(qm -> star[i], qm, CONVC, OPTCOUNT );
		for ( j = 0; j < qm -> star[i] ->nV; ++j ) {
			stat = EG_evaluate(qm->face, &qm -> uvs[2*(qm ->star[i] ->verts[j] - 1)], vals);
			if ( stat != EGADS_SUCCESS) return stat;
			cross_product(&vals[3], &vals[6], normal);
			EG_unitVector(normal, &normal[3]);
			sampleNormalPlane(normal, vals, qm ->star[i] ->verts[j],  qm, CONVC, OPTCOUNT);
		}
	}
	t_end = clock();
	//printf(" *-*-*-*-*-*-*- TIME REQUIRED TO OPTIMIZE %d POINTS = %lf  \n",nV, (double)(t_end - t_start) / CLOCKS_PER_SEC);
	NLOPTtermination(res);
	printf(" TERMINATION WITH OBJECTIVE %lf  TOTAL ITERATIONS %d\n",fOpti,OPTCOUNT);
	if ( res < 0 ) {
		stat = EGADS_GEOMERR;
		fclose(filOpti);
		goto cleanup;
	}
	else {
		stat = EGADS_SUCCESS;
		OPTCOUNT = -1 ;
		printf(" OUTPUT FINAL ITERATION\n");
		fOpti = NLOPT_objEqualAnglesTess(2*nV, uvOpt, NULL, qm);
	}
	fclose(filOpti);
	CONVC++;
	nlopt_destroy(opt);
	for (j = 0; j < nV; ++j) {
		qm ->uvs[2*(vID[j] -1 )    ] = uvOpt[2*j    ];
		qm ->uvs[2*(vID[j] -1 ) + 1] = uvOpt[2*j + 1];
		EG_evaluate(qm ->face, &qm ->uvs[2*(vID[j]-1)], XYZ);
		qm -> xyzs[3*(vID[j] -1 )   ] = XYZ[0];
		qm -> xyzs[3*(vID[j] -1 ) +1] = XYZ[1];
		qm -> xyzs[3*(vID[j] -1 ) +2] = XYZ[2];
		IOSTATUS = 1;
		for ( i = 0 ; i < qm->star[j] -> nQ; ++i) {
			stat2 = quad_algebraic_area(qm, qm->star[j] ->quads[i], qm->star[j] ->verts[0], &area, angles, &quadArea);
			if (stat2 != EGADS_SUCCESS || area == -1) {
				printf("OPTI  STAT2  %d  AREA %d \n",stat, area);
				FILE *file;
				char buffer[32];
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
				stat = EGADS_GEOMERR;
				goto cleanup;
			}
		}
	}
	cleanup:
	if ( stat != EGADS_SUCCESS) {
		printf(" NLOPT FAILED TO PRODUCE A VALID MESH. RESET VALUES AND RETURN ERROR %d\n", stat);
	}
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
	EG_free(uv0);
	EG_free(epsABS);
	return stat;
	//if ( stat != EGADS_SUCCESS) exit(1);
	//return  stat;
}









int main(int argc, char *argv[])
{
	int          i, j, k, ibody, stat, oclass, mtype, len, ntri, sum;
	int          nseg, ngp, atype, alen, quad, *segs, *senses;
	const int    *tris, *tric, *ptype, *pindex, *ints;
	float        arg, color[3];
	double       box[6], size, params[3];
	const double *xyzs, *uvs, *ts, *reals ;
	double *mod_xyzs;
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
			printf(" TESSELATION HAS %d POINTS\n",len);
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
	//double       *xyzs, result[18];
	stat = 0;
	stat = updateUVs(bodydata[0].tess, FRONTFACE +1, bodydata[0].faces[FRONTFACE],
			bodydata[0].uvs,mod_xyzs);
	/*  if (wv_startServer(7681, NULL, NULL, NULL, 0, cntxt) == 0) {

    // we have a single valid server -- stay alive a long as we have a client //
    while (wv_statusServer(0)) {
      usleep(500000);
      if (stat == 0) {
        if (startapp != NULL) system(startapp);
        stat++;
      }
    }
  }


  wv_cleanupServers();*/
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


#endif
