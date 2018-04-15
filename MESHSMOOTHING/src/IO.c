/*
 * IO.c
 *
 *  Created on: Apr 13, 2018
 *      Author: docampo
 */

#include "IO.h"



void sampleNormalPlane(double normal[], double point[], int vID, quadMap *qm, int nopt,int ncount) {
  double min[3], max[3], c, p[3], dt[3], r = 1.0;
  int  i, j, k, nP;
  FILE *f;
  char name[32];
  snprintf(name, sizeof(char) * 32, "PLANE%i_OPTIN_%i_SCOUNT_%i.txt", nopt, ncount,vID);
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

 void
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
 void
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
 void
print_starFile(vStar *star, quadMap *qm, int call, int count) {
  int i, vID, k;
  FILE *fout;
  char buffer[33];
  //snprintf(buffer, sizeof(char) * 32, "STAR%i_CENTRE_%i.txt", CONVC, star->verts[0]);
  snprintf(buffer, sizeof(char) * 32, "STAR%i_OPTIN_%i_SCOUNT_%i.txt", call, count, star->verts[0]);
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
  snprintf(buffer, sizeof(char) * 32, "CENTRE%i_OPTIN_%i_SCOUNT_%i.txt", call, count, star->verts[0]);
  //snprintf(buffer, sizeof(char) * 32, "CENTRE%i_SCOUNT_%i.txt", CONVC, star->verts[0]);
  fout = fopen(buffer,"w");
  if (fout == NULL ) return;
  vID = star->verts[0] - 1;
  fprintf(fout, "%lf  %lf  %lf %d\n", qm->xyzs[3*vID],qm->xyzs[3*vID + 1],qm->xyzs[3*vID + 2],vID + 1);
  fclose(fout);
}








 void
printMeshStats(quadMap *qm, int sweep) {
  int fix, move, i,len, val ;
  int intVal[100], boundVal[100];
  FILE *fout;
  char buffer[33];
  snprintf(buffer, sizeof(char) * 32, "MESH_STATS_%d.txt", sweep);
  printf(" WRITING ON FILE %s\n",buffer);
  fout = fopen(buffer, "w");
  len = qm ->totVerts;
  for ( i = 0; i < 100; ++i) {
    intVal  [i] = 0;
    boundVal[i] = 0;
  }
  for ( i = 0; i < len; ++i) {
    if ( qm -> vType[i] == -1){
      printf(" VERTEX %d HAS VALENCE %d. ADD TO LIST\n", i + 1, qm -> valence[i][0]);
      val = qm -> valence[i][0];
      ++intVal[val];
      printf(" ADD VALENCE TO %d = %d \n",val, intVal[val]);
      if ( val == 99 ) printf(" ALERT !!!!!!!!!!!!!!!!! %d = %d\n",val, intVal[val]);
    } else {
      ++boundVal[qm ->valence[i][0]];
    }
  }
  fprintf(fout,"--------------------- MESH STATS AFTER %d SWEEP --------------\n",sweep);
  fprintf(fout,"---- TOTAL VERTICES %d TOTAL QUADS %d --------------\n",qm -> totVerts - qm -> remVerts[0],qm -> nQ - qm -> remQuads[0] );
  fprintf(fout," INTERIOR VERTICES\n");
  for ( i = 0 ; i < 100; ++i) {
    printf(" VALENCE %d - > %d\n", i, intVal[i]);
    if ( intVal[i]  > 0  && intVal[i] < 99) fprintf(fout," VALENCY %d = %d VERTICES\n", i, intVal[i]);
  }
  fprintf(fout," BOUNDARY VERTICES\n");
  for ( i = 0 ; i < 100; ++i) {
    if ( boundVal[i]  > 0 ) fprintf(fout," VALENCY %d = %d VERTICES\n", i, boundVal[i]);
  }
  fclose(fout);
  return ;
}


 void
print_mesh(quadMap *qm, int *fn) {
  int i,k, j, id;
  FILE *fout;
  char buffer[33];
  double eval[18];
  snprintf(buffer, sizeof(char) * 32, "MESH%i.txt", *fn);
  printf(" WRITING ON FILE %s\n",buffer);
  fout = fopen(buffer,"w");
  if (fout == NULL ) return;
  for ( i = 0 ; i < qm ->nQ; ++i) {
    if ( qm -> quadIdx[4*i] != -1 ) {
      for ( k = 0; k < 4; ++k)
        fprintf(fout, "%lf  %lf  %lf %d\n",qm->xyzs[3*(qm->quadIdx[4*i + k] - 1)], qm->xyzs[3*(qm->quadIdx[4*i + k] - 1) +1], qm->xyzs[3*(qm->quadIdx[4*i + k] - 1) + 2], qm -> quadIdx[4*i + k]);
      fprintf(fout, "%lf  %lf  %lf %d\n",qm->xyzs[3*(qm->quadIdx[4*i] - 1)], qm->xyzs[3*(qm->quadIdx[4*i] - 1) +1], qm->xyzs[3*(qm->quadIdx[4*i] - 1) + 2], qm -> quadIdx[4*i]);
      fprintf(fout,"\n\n");
    }
  }
  fclose(fout);
  (*fn)++;
}

