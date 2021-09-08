/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"

extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double []);

double f1dim(double x)
/* Must accompany linmin.
See README-nr.txt for more information.
*/
{
	int j;
	double f,*xt;

	xt=vector_double(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_vector_double(xt,1,ncom);
	return f;
}
#undef NRANSI
