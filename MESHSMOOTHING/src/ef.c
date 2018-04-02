#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define SUPP 3.554147
#define EPS 0.1
int main (int argc, char *argv[]) {
	if (argc != 2) {
	printf(" INPUT DIVIDING FACTOR\n");
	exit(1);
	}
	double div, dt, z, zerr,zerrc, dot = 0.0, alpha = 0.0, penalty;
	div = atof(argv[1]);
	FILE *fout;
	char buff[33];
	snprintf(buff,sizeof(char)*32, "errof_%f",div);
	fout = fopen(buff,"w");
	if (fout == NULL) return -1;
	dt = 0.001;
	while(alpha <=2*M_PI)
	{
	  dot   = cos(alpha - M_PI*0.5);
	  if ( alpha > 2*M_PI - acos(div)) dot = -1.0;
/*          if ( alpha < M_PI - EPS)	z = SUPP*(1.0 - dot)/div;
	  else  		z = SUPP*(dot - 1)/div;
*/	  
	  z       = SUPP*((dot-div*0.5)/div);
	  zerr    = erf(z);
	  zerrc   = erfc(z);
          penalty = zerrc * exp(SUPP*(1.0 + dot));
          fprintf(fout,"%lf %lf %lf %lf %lf %lf\n",z,zerr,zerrc,alpha*180.0/M_PI, dot,penalty);
	  alpha += dt;
	}
	return 0;
}	


