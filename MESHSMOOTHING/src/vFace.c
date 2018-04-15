/*
 * planarFace.c
 *
 *  Created on: Feb 16, 2018
 *      Author: docampo
 */






/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             An example of "bottom-up" construction
 *
 *      Copyright 2011-2017, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

//#define MKCYL

#include "egads.h"
#ifdef MKCYL


int main(/*@unused@*/int argc, /*@unused@*/ char *argv[])
{
  int    i, nface, nseg;
  int senses[4];
  double data[7], bbox[6];
  ego    context, body, bcut, model, *faces, *facEdg;
  senses[0] = 1;
  senses[1] = 1;
  senses[2] = 1;
  senses[3] = 1;
  if (argc != 1) {
    printf("\n Usage: scribe\n\n");
    return 1;
  }

  /* initialize */
  printf(" EG_open           = %d\n", EG_open(&context));
  data[0] = data[2] = data[3] = data[5] = 0.0;
  data[1] = -1.0;
  data[4] =  1.0;
  data[6] =  1;
  printf(" EG_makeSolidBody  = %d\n", EG_makeSolidBody(context, CYLINDER, data,
      &body));
  printf(" EG_getBodyTopos   = %d\n", EG_getBodyTopos(body, NULL, FACE,
      &nface, &faces));


  printf(" EG_makeTopology L  = %d\n", EG_makeTopology(context, NULL, BODY, FACEBODY,
      NULL, 1, &faces[0], senses, &body));

  printf(" EG_makeTopology M  = %d\n", EG_makeTopology(context, NULL, MODEL, 0,
      NULL, 1, &body,
      NULL, &model));
  printf(" EG_saveModel       = %d\n", EG_saveModel(model, "cylinderFace.egads"));
  printf("\n");
  EG_setOutLevel(context, 2);
  printf(" EG_close           = %d\n", EG_close(context));

  return 0;
}

#else








int main(/*@unused@*/int argc, /*@unused@*/ char *argv[])
{
  int    senses[4];
  double xyz[3], data[6], dum[3], range[2];
  ego    context, nodes[4], lines[4], edges[4], objs[2], loop, face, body;
  ego    model = NULL;

  printf(" EG_open            = %d\n", EG_open(&context));

  /* make the Nodes */
  xyz[0] = xyz[1] = xyz[2] = 0.0;
  printf(" EG_makeTopology N0 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
      xyz, 0, NULL, NULL,
      &nodes[0]));
  xyz[0] = 1.0;
  printf(" EG_makeTopology N1 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
      xyz, 0, NULL, NULL,
      &nodes[1]));
  xyz[1] = 1.0;
  printf(" EG_makeTopology N2 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
      xyz, 0, NULL, NULL,
      &nodes[2]));
  xyz[0] = 0.0;
  xyz[1] = 1.0;
  printf(" EG_makeTopology N2 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
      xyz, 0, NULL, NULL,
      &nodes[3]));
  /* make the Curves */
  data[0] = data[1] = data[2] = data[4] = data[5] = 0.0;
  data[3] = 1.0;
  printf(" EG_makeGeometry L0 = %d\n", EG_makeGeometry(context, CURVE, LINE,
      NULL, NULL, data,
      &lines[0]));
  data[0] = 1.0;
  data[3] = 0.0;
  data[4] = 1.0;
  printf(" EG_makeGeometry L1 = %d\n", EG_makeGeometry(context, CURVE, LINE,
      NULL, NULL, data,
      &lines[1]));
  data[1] =  1.0;
  data[3] = -1.0;
  data[4] = 0.0;
  printf(" EG_makeGeometry L2 = %d\n", EG_makeGeometry(context, CURVE, LINE,
      NULL, NULL, data,
      &lines[2]));
  data[0] =  0.0;
  data[3] = 0.0;
  data[4] = -1.0;
  printf(" EG_makeGeometry L2 = %d\n", EG_makeGeometry(context, CURVE, LINE,
      NULL, NULL, data,
      &lines[3]));

  /* construct the Edges */
  xyz[0] = xyz[1] = xyz[2] = 0.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[0], xyz, &range[0],
      dum));
  xyz[0] = 1.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[0], xyz, &range[1],
      dum));
  printf("                      range = %lf %lf\n", range[0], range[1]);
  objs[0] = nodes[0];
  objs[1] = nodes[1];
  printf(" EG_makeTopology E0 = %d\n", EG_makeTopology(context, lines[0], EDGE,
      TWONODE, range, 2,
      objs, NULL, &edges[0]));

  xyz[0] = 1.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[1], xyz, &range[0],
      dum));
  xyz[1] = 1.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[1], xyz, &range[1],
      dum));
  printf("                      range = %lf %lf\n", range[0], range[1]);
  objs[0] = nodes[1];
  objs[1] = nodes[2];
  printf(" EG_makeTopology E1 = %d\n", EG_makeTopology(context, lines[1], EDGE,
      TWONODE, range, 2,
      objs, NULL, &edges[1]));

  xyz[0] = 1.0;
  xyz[1] = 1.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[2], xyz, &range[0],
      dum));
  xyz[0] = 0.0;
  xyz[1] = 1.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[2], xyz, &range[1],
      dum));
  printf("                      range = %lf %lf\n", range[0], range[1]);
  objs[0] = nodes[2];
  objs[1] = nodes[3];
  printf(" EG_makeTopology E2 = %d\n", EG_makeTopology(context, lines[2], EDGE,
      TWONODE, range, 2,
      objs, NULL, &edges[2]));
  xyz[0] = 0.0;
  xyz[1] = 1.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[3], xyz, &range[0],
      dum));
  xyz[0] = 0.0;
  xyz[1] = 0.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[3], xyz, &range[1],
      dum));
  printf("                      range = %lf %lf\n", range[0], range[1]);
  objs[0] = nodes[3];
  objs[1] = nodes[0];
  printf(" EG_makeTopology E2 = %d\n", EG_makeTopology(context, lines[3], EDGE,
      TWONODE, range, 2,
      objs, NULL, &edges[3]));
  /* make the Loop and then the Face */
  senses[0] = 1;
  senses[1] = 1;
  senses[2] = 1;
  senses[3] = 1;
  printf(" EG_makeTopology L  = %d\n", EG_makeTopology(context, NULL, LOOP,
      CLOSED, NULL, 4, edges,
      senses, &loop));

  printf(" EG_makeFace        = %d\n", EG_makeFace(loop, SREVERSE, NULL,
      &face));
  printf(" EG_makeTopology L  = %d\n", EG_makeTopology(context, NULL, BODY, FACEBODY,
      NULL, 1, &face, senses, &body));

  printf(" EG_makeTopology M  = %d\n", EG_makeTopology(context, NULL, MODEL, 0,
      NULL, 1, &body,
      NULL, &model));
  printf(" EG_saveModel       = %d\n", EG_saveModel(model, "planarFace.egads"));
  printf("\n");



  EG_setOutLevel(context, 2);
  printf(" EG_close           = %d\n", EG_close(context));

  return 0;
}
#endif
