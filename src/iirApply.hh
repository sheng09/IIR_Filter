#ifndef __IIR_APPLY__
#define __IIR_APPLY__

int apply(float *data, int nsamps, bool zp, float *sn, float *sd, int nsects);

#endif // !__IIR_APPLY__