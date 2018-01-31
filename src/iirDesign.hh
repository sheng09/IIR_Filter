#ifndef __IIR_DESIGN_H__
#define __IIR_DESIGN_H__

//#include <string>
#include "complex.hh"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */ 


#define IIR_LP 0
#define IIR_HP 1
#define IIR_BP 2
#define IIR_BR 3

#define IIR_BU 0
#define IIR_BE 1
#define IIR_C1 2
#define IIR_C2 3

typedef int IIR_TYPE;
typedef int IIR_APROTO;

int design(int iord, IIR_TYPE type, int IIR_APROTO, double a, double trbndw,
            double fl, double fh, double ts, float *sn, float *sd, int *nsects);
void buroots(complexf *p, char *rtype, int rtype_s, float *dcvalue, int *nsects, int iord);
void beroots(complexf *p, char *rtype, int rtype_s, float *dcvalue, int *nsects, int iord);
void chebparm(double a, double trbndw, int iord, float *eps, float *ripple);
void c1roots(complexf *p, char *rtype, int rtype_s, float *dcvalue, int *nsects, int iord, double eps);
void c2roots(complexf *p, complexf *z, char *rtype, int rtype_s, float *dcvalue, int *nsects, int iord, double a, double omegar);
void lp(complexf *p, complexf *z, char *rtype, int rtype_s, double dcvalue, int nsects, float *sn, float *sd);
void lptbp(complexf *p, complexf *z, char *rtype, int rtype_s, double dcvalue, int *nsects, double fl, double fh, float *sn, float *sd);
void lptbr(complexf *p, complexf *z, char *rtype, int rtype_s, double dcvalue, int *nsects, double fl, double fh, float *sn, float *sd);
void lpthp(complexf *p, complexf *z, char *rtype, int rtype_s, double dcvalue, int nsects, float *sn, float *sd);

double warp(double f, double ts);
void cutoffs(float *sn, float *sd, int nsects, double f);

void bilin2(float *sn, float *sd, int nsects);

char *fstrncpy(char *to, int tolen, std::string s_from, int fromlen);


#ifdef __cplusplus
}
#endif /* __cplusplus */ 


#endif //!__IIR_DESIGN_H__
