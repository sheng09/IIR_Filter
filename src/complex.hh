#ifndef __COMPLEX_H__
#define __COMPLEX_H__

#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

typedef struct complexf_t
{
    float re; /** Real Part */
    float im; /** Imaginary Part */
} complexf;

double   aimag(complexf fc);
double   cmplxabs(complexf c);
double   cmplxang(complexf c);
double   cmplxtof(complexf c);

complexf cmplxadd(complexf c1, complexf c2);
complexf cmplxcj(complexf c);
complexf cmplxmul(complexf c1, complexf c2);
complexf flttocmplx(double d1, double d2);
complexf cmplxsub(complexf c1, complexf c2);

complexf cmplxsqrt(complexf c);
complexf cmplxdiv(complexf c1, complexf c2);
complexf cmplxlog(complexf c);
complexf cmplxexp(complexf c);
complexf cmplxpow(complexf c, double d);
complexf cmplxneg(complexf c);
double powi(double b, int x);


#endif // !__COMPLEX_H__
