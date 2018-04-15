/*
 * IO.h
 *
 *  Created on: Apr 13, 2018
 *      Author: docampo
 */

#ifndef SRC_IO_H_
#define SRC_IO_H_

#include"egads.h"
#include <math.h>



#define MAX_LINKS 100
#define EPS  1.E-05


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




typedef struct{
  int    *verts, *quads;
  int  nV, nQ; // nV = origin(1) + peaks (n)
} vStar;

typedef struct {
  int     *quadIdx, *quadAdj, **valence, *vType, *remQuads, *remVerts;
  int      nQ, nS, totVerts;
  ego      face;
  double  *xyzs, *oriXYZ, *uvs,  *vLengths, *vCurvs;
  vStar   **star;
} quadMap;






void
print_quad_specs(quadMap *qm, int id) ;

void sampleNormalPlane(double normal[], double point[], int vID, quadMap *qm, int, int) ;
void NLOPTtermination(int n) ;

void print_starFile(vStar *star, quadMap *qm, int, int) ;
void print_star(vStar *star) ;
void
printMeshStats(quadMap *qm, int sweep);
void
print_mesh(quadMap *qm, int *fn);





#endif /* SRC_IO_H_ */
