#include <cstdio>
#include <cstring>
#include <cstdlib>
//#include <string>

#include <iostream>
#include <cmath>
#include "iirDesign.hh"
#include "complex.hh"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */ 


/** 
 * Design IIR Digital Filters from Analog Prototypes
 * 
 * @param iord 
 *    Filter Order, Maximum of 10
 * @param type 
 *    Filter Type, Character *2
 *      - 0 'LP'  Lowpass
 *      - 1 'HP'  Highpass
 *      - 2 'BP'  Bandpass
 *      - 4 'BR'  Bandreject
 * @param aproto 
 *    Analog Prototype, Character *2
 *      - 0 'BU'  Butterworth
 *      - 1 'BE'  Bessel
 *      - 2 'C1'  Cheyshev Type I
 *      - 4 'C2'  Cheyshev Type II
 * @param a 
 *    Chebyshev stopband Attenuation Factor
 * @param trbndw 
 *    Chebyshev transition bandwidth, fraction of lowpass prototype 
 *    passband width
 * @param fl 
 *    Low Frequency cutoff
 * @param fh 
 *    High Frequency cutoff
 * @param ts 
 *    Sampling Rate / Delta
 * @param sn 
 *    Array containing numerator coefficients of 2nd Order Sections
 *    Packed Head to Tail
 * @param sd 
 *    Array containing denominator coefficients of 2nd Order Sections
 *    Packed Head to Tail
 * @param nsects 
 *    Length of arrays \p sn and \p sd
 *
 *  @copyright 1990  Regents of the University of California                      
 *
 *  @author  Dave Harris
 *           Lawrence Livermore National Laboratory
 *           L-205
 *           P.O. Box 808
 *           Livermore, CA  94550
 *           USA
 * 
 *  @date Documented/Reviewed
 *
 */
int design(int iord, IIR_TYPE type, IIR_APROTO aproto, double a, double trbndw,
            double fl, double fh, double ts, float *sn, float *sd, int *nsects)
{
    char stype[10][4];
    float dcvalue, eps, fhw, flw, omegar, ripple;
    complexf p[10], z[10];

    //std::cout << "b desin\n";
    /*  Analog prototype selection                                                   
	 * */

    if (aproto == IIR_BU)  // BU
    {
        buroots(p, (char *)stype, 4, &dcvalue, nsects, iord);
    }
    else if (aproto == IIR_BE)
    {
        beroots(p, (char *)stype, 4, &dcvalue, nsects, iord);
    }
    else if (aproto == IIR_C1)
    {
        chebparm(a, trbndw, iord, &eps, &ripple);
        c1roots(p, (char *)stype, 4, &dcvalue, nsects, iord, eps);
    }
    else if (aproto == IIR_C2)
    {
        omegar = 1. + trbndw;
        c2roots(p, z, (char *)stype, 4, &dcvalue, nsects, iord, a, omegar);
    }
    else
	{
        fprintf(stderr, "filter: Unknown Analog filter prototype: '%d'\n", aproto);
        fprintf(stderr, "        Expected: BU, BESSEL, C1, C2\n");
        return -1;
    }

    /*  Analog mapping selection                                                     
	 * */
    if (type == IIR_BP) // BP
    {
        flw = warp(fl * ts / 2., 2.);
        fhw = warp(fh * ts / 2., 2.);
        lptbp(p, z, (char *)stype, 4, dcvalue, nsects, flw, fhw, sn, sd);
    }
    else if (type == IIR_BR) // BR
    {
        flw = warp(fl * ts / 2., 2.);
        fhw = warp(fh * ts / 2., 2.);
        lptbr(p, z, (char *)stype, 4, dcvalue, nsects, flw, fhw, sn, sd);
    }
    else if (type == IIR_LP) // LP
    {
        fhw = warp(fh * ts / 2., 2.);
        lp(p, z, (char *)stype, 4, dcvalue, *nsects, sn, sd);
        cutoffs(sn, sd, *nsects, fhw);
    }
    else if (type == IIR_HP) // HP
    {
        flw = warp(fl * ts / 2., 2.);
        lpthp(p, z, (char *)stype, 4, dcvalue, *nsects, sn, sd);
        cutoffs(sn, sd, *nsects, flw);
    }

    /*  Bilinear analog to digital transformation                                    
	 * */
    bilin2(sn, sd, *nsects);
    return 0;
}

/** 
 * Compute the Butterworth Poles for a Normalized Low Pass (LP) Filter
 * 
 * @param p 
 *    Complex Array containing Poles. Contains only one of from each
 *     - Complex Conjugate Pole-Zero Pair
 *     - Complex Conjugate Pole Pair
 *     - Single Real Pole
 * @param rtype 
 *    Character Array indicating 2nd Order Section Types
 *     - 'CPZ' Complex Conjugate Pole-Zero Pair
 *     - 'CP'  Complex Conjugate Pole Pair
 *     - 'SP'  Single Real Pole
 * @param rtype_s 
 *     Length of string \p rtype
 * @param dcvalue 
 *     Magnitude of the filter at Zero Frequency
 * @param nsects 
 *     Number of 2nd Order Sections
 * @param iord 
 *     Desired Filter Order, Must be between 1 and 8
 *
 * @return Nothing
 *
 * @copyright 1990  Regents of the University of California
 *
 *
 * @author  Dave Harris                                                         
 *          Lawrence Livermore National Laboratory                              
 *          L-205                                                               
 *          P.O. Box 808                                                        
 *          Livermore, CA  94550                                                
 *          USA                                                                 
 *          (415) 423-0617                                                      
 *
 * \note
 *   \f[
 *     \begin{array}{rcl}
 *        n &=& iord \\
 *     H(s) &=& k_0 / \displaystyle\prod_{k=1}^n  s - s_k \\
 *      s_k &=& exp ( i \pi [0.5 + 2(2 k - 1)/2n ] ); k = 1, 2, ..., n \\
 *      k_0 &=& 1.0
 *     \end{array}
 *   \f]
 *
 * \date 900907   LAST MODIFIED
 * 
 */
void buroots(complexf *p, char *rtype, int rtype_s, float *dcvalue, int *nsects, int iord)
{
    #define RTYPE(I_, J_) (rtype + (I_) * (rtype_s) + (J_))
    int half, k;
    float angle, pi;

    complexf *const P = &p[0] - 1;

    pi = 3.14159265;

    half = iord / 2;

    /* TEST FOR ODD ORDER, AND ADD POLE AT -1                                        
	 * */
    *nsects = 0;
    if (2 * half < iord)
    {
        P[1] = flttocmplx(-1., 0.);
        fstrncpy(RTYPE(0, 0), rtype_s - 1, "SP", 2);
        *nsects = 1;
    }
    for (k = 1; k <= half; k++)
    {
        angle = pi * (.5 + (float)(2 * k - 1) / (float)(2 * iord));
        *nsects = *nsects + 1;
        P[*nsects] = flttocmplx(cos(angle), sin(angle));
        fstrncpy(RTYPE(*nsects - 1, 0), rtype_s - 1, "CP", 2);
    }
    *dcvalue = 1.0;

    return ;
    #undef RTYPE
}

/** 
 * Compute Bessel Poles For a Normalized Low Pass (LP) Filter
 * 
 * @param p 
 *    Complex Array containing Poles. Contains only one of from each
 *     - Complex Conjugate Pole-Zero Pair
 *     - Complex Conjugate Pole Pair
 *     - Single Real Pole
 * @param rtype 
 *    Character Array indicating 2nd Order Section Types
 *     - 'CPZ' Complex Conjugate Pole-Zero Pair
 *     - 'CP'  Complex Conjugate Pole Pair
 *     - 'SP'  Single Real Pole
 * @param rtype_s 
 *     Length of string \p rtype
 * @param dcvalue 
 *     Magnitude of the filter at Zero Frequency
 * @param nsects 
 *     Number of 2nd Order Sections
 * @param iord 
 *     Desired Filter Order, Must be between 1 and 8
 *
 * @return Nothing
 *
 * @copyright 1990  Regents of the University of California                      
 *
 * \bug Poles are defined explicitly within the routine. 
 *        Pull them out and define and document them.
 *
 * @author  Dave Harris                                                         
 *           Lawrence Livermore National Laboratory                              
 *           L-205                                                               
 *           P.O. Box 808                                                        
 *           Livermore, CA  94550                                                
 *           USA                                                                 
 *           (415) 423-0617                                                      
 *
 *
 * \date 920415   Changed P and RTYPE to adjustable array by 
 *                     using an "*" rather than a "1".     
 *
 */
void beroots(complexf *p, char *rtype, int rtype_s, float *dcvalue, int *nsects, int iord)
{
    #define RTYPE(I_, J_) (rtype + (I_) * (rtype_s) + (J_))

    complexf *const P = &p[0] - 1;

    if (iord == 1)
    {

        P[1] = flttocmplx(-1.0, 0.0);
        fstrncpy(RTYPE(0, 0), rtype_s - 1, "SP", 2);
    }
    else if (iord == 2)
    {

        P[1] = flttocmplx(-1.1016013, 0.6360098);
        fstrncpy(RTYPE(0, 0), rtype_s - 1, "CP", 2);
    }
    else if (iord == 3)
    {

        P[1] = flttocmplx(-1.0474091, 0.9992645);
        fstrncpy(RTYPE(0, 0), rtype_s - 1, "CP", 2);
        P[2] = flttocmplx(-1.3226758, 0.0);
        fstrncpy(RTYPE(1, 0), rtype_s - 1, "SP", 2);
    }
    else if (iord == 4)
    {

        P[1] = flttocmplx(-0.9952088, 1.2571058);
        fstrncpy(RTYPE(0, 0), rtype_s - 1, "CP", 2);
        P[2] = flttocmplx(-1.3700679, 0.4102497);
        fstrncpy(RTYPE(1, 0), rtype_s - 1, "CP", 2);
    }
    else if (iord == 5)
    {

        P[1] = flttocmplx(-0.9576766, 1.4711244);
        fstrncpy(RTYPE(0, 0), rtype_s - 1, "CP", 2);
        P[2] = flttocmplx(-1.3808774, 0.7179096);
        fstrncpy(RTYPE(1, 0), rtype_s - 1, "CP", 2);
        P[3] = flttocmplx(-1.5023160, 0.0);
        fstrncpy(RTYPE(2, 0), rtype_s - 1, "SP", 2);
    }
    else if (iord == 6)
    {

        P[1] = flttocmplx(-0.9306565, 1.6618633);
        fstrncpy(RTYPE(0, 0), rtype_s - 1, "CP", 2);
        P[2] = flttocmplx(-1.3818581, 0.9714719);
        fstrncpy(RTYPE(1, 0), rtype_s - 1, "CP", 2);
        P[3] = flttocmplx(-1.5714904, 0.3208964);
        fstrncpy(RTYPE(2, 0), rtype_s - 1, "CP", 2);
    }
    else if (iord == 7)
    {

        P[1] = flttocmplx(-0.9098678, 1.8364514);
        fstrncpy(RTYPE(0, 0), rtype_s - 1, "CP", 2);
        P[2] = flttocmplx(-1.3789032, 1.1915667);
        fstrncpy(RTYPE(1, 0), rtype_s - 1, "CP", 2);
        P[3] = flttocmplx(-1.6120388, 0.5892445);
        fstrncpy(RTYPE(2, 0), rtype_s - 1, "CP", 2);
        P[4] = flttocmplx(-1.6843682, 0.0);
        fstrncpy(RTYPE(3, 0), rtype_s - 1, "SP", 2);
    }
    else if (iord == 8)
    {

        P[1] = flttocmplx(-0.8928710, 1.9983286);
        fstrncpy(RTYPE(0, 0), rtype_s - 1, "CP", 2);
        P[2] = flttocmplx(-1.3738431, 1.3883585);
        fstrncpy(RTYPE(1, 0), rtype_s - 1, "CP", 2);
        P[3] = flttocmplx(-1.6369417, 0.8227968);
        fstrncpy(RTYPE(2, 0), rtype_s - 1, "CP", 2);
        P[4] = flttocmplx(-1.7574108, 0.2728679);
        fstrncpy(RTYPE(3, 0), rtype_s - 1, "CP", 2);
    }

    *nsects = iord - iord / 2;

    *dcvalue = 1.0;

    /*  DONE                                                                         
	 * */
    return;
    #undef RTYPE
}

/** 
 * Calculate Chebyshev Type I and II Design Parameters
 * 
 * @param a 
 *    Desired Stopband Attenuation, i.e. max stopband 
 *    amplitude is 1/ATTEN
 * @param trbndw 
 *    Transition bandwidth between stop and passband as a 
 *    fraction of the passband width
 * @param iord 
 *    Filter Order (number of Poles)
 * @param eps 
 *    Output Chebyshev passband parameter
 * @param ripple 
 *    Passband ripple
 *
 * @return Nothing
 *
 * \note
 *   \f[
 *     \begin{array}{rcl}
 *  \omega &=& 1.0 + trbndw \\
 *  \alpha &=& (\omega + \sqrt{\omega^2 - 1.0} )^{iord}\\
 *       g &=& \alpha^2 + 1 / 2\alpha \\
 *     eps &=& \sqrt{a^2 - 1.0 } / g \\
 *  ripple &=& 1 / \sqrt{ 1.0 + eps^2 } \\
 *     \end{array}
 *   \f]
 *
 * \copyright 1990  Regents of the University of California
 *
 * \author   Dave Harris
 *           Lawrence Livermore National Laboratory
 *           L-205
 *           P.O. Box 808
 *           Livermore, CA  94550
 *           USA
 *           (415) 423-0617
 */
void chebparm(double a, double trbndw, int iord, float *eps, float *ripple)
{
    float alpha, g, omegar;

    omegar = 1. + trbndw;
    alpha = powi(omegar + sqrt(powi(omegar, 2) - 1.), iord);
    g = (powi(alpha, 2) + 1.) / (2. * alpha);
    *eps = sqrt(powi(a, 2) - 1.) / g;
    *ripple = 1. / sqrt(1. + powi(*eps, 2));

    return;
}

/** 
 * Compute Chebyshev Type I Poles for a Normalized Low Pass (LP) Filter
 * 
 * @param p 
 *    Complex Array containing Poles. Contains only one of from each
 *     - Complex Conjugate Pole-Zero Pair
 *     - Complex Conjugate Pole Pair
 *     - Single Real Pole
 * @param rtype 
 *    Character Array indicating 2nd Order Section Types
 *     - 'CPZ' Complex Conjugate Pole-Zero Pair
 *     - 'CP'  Complex Conjugate Pole Pair
 *     - 'SP'  Single Real Pole
 * @param rtype_s 
 *     Length of string \p rtype
 * @param dcvalue 
 *     Magnitude of the filter at Zero Frequency
 * @param nsects 
 *     Number of 2nd Order Sections
 * @param iord 
 *     Desired Filter Order, Must be between 1 and 8
 * @param eps
 *     Output Chebyshev Parameter Related to Passband Ripple
 *
 * @return Nothing
 *
 * @copyright 1990  Regents of the University of California
 *
 * @author  Dave Harris                                                         
 *          Lawrence Livermore National Laboratory                              
 *          L-205                                                               
 *          P.O. Box 808                                                        
 *          Livermore, CA  94550                                                
 *          USA                                                                 
 *          (415) 423-0617                                                      
 *
 * \todo Documentation for chebyshev Type I Filter.  
 *          Find roots and multiply them together.
 *
 *
 * \date 900907    LAST MODIFIED
 *
 */
void c1roots(complexf *p, char *rtype, int rtype_s, float *dcvalue, int *nsects, int iord, double eps)
{
    #define RTYPE(I_, J_) (rtype + (I_) * (rtype_s) + (J_))
    int half, i, i_;
    float angle, c, gamma, omega, pi, s, sigma;

    complexf *const P = &p[0] - 1;

    pi = 3.14159265;
    half = iord / 2;

    /*  INTERMEDIATE DESIGN PARAMETERS                                               
	 * */
    gamma = (1. + sqrt(1. + eps * eps)) / eps;
    gamma = log(gamma) / (float)(iord);
    gamma = exp(gamma);
    s = .5 * (gamma - 1. / gamma);
    c = .5 * (gamma + 1. / gamma);

    /*  CALCULATE POLES                                                              
	 * */
    *nsects = 0;
    for (i = 1; i <= half; i++)
    {
        i_ = i - 1;
        fstrncpy(RTYPE(i_, 0), rtype_s - 1, "CP", 2);
        angle = (float)(2 * i - 1) * pi / (float)(2 * iord);
        sigma = -s * sin(angle);
        omega = c * cos(angle);
        P[i] = flttocmplx(sigma, omega);
        *nsects = *nsects + 1;
    }
    if (2 * half < iord)
    {
        fstrncpy(RTYPE(half, 0), rtype_s - 1, "SP", 2);
        P[half + 1] = flttocmplx(-s, 0.0);
        *nsects = *nsects + 1;
        *dcvalue = 1.0;
    }
    else
    {
        *dcvalue = 1. / sqrt(1 + powi(eps, 2));
    }

    /*  DONE                                                                         
	 * */
    return;
    #undef RTYPE
}

/** 
 * Compute root for normailzed Low Pass Chebyshev Type II Filter
 * 
 * @param p 
 *    Complex Array containing Poles. Contains only one of from each
 *     - Complex Conjugate Pole-Zero Pair
 *     - Complex Conjugate Pole Pair
 *     - Single Real Pole
 * @param z 
 *    Complex Array containing Zeros Contains only one of from each
 *     - Complex Conjugate Pole-Zero Pair
 *     - Complex Conjugate Pole Pair
 *     - Single Real Pole
 * @param rtype 
 *    Character Array indicating 2nd Order Section Types
 *     - 'CPZ' Complex Conjugate Pole-Zero Pair
 *     - 'CP'  Complex Conjugate Pole Pair
 *     - 'SP'  Single Real Pole
 * @param rtype_s 
 *    Length of string \p rtype
 * @param dcvalue 
 *    Magnitude of filter at Zero Frequency
 * @param nsects 
 *    Number of 2nd order sections
 * @param iord 
 *   Input Desired Filter Order
 * @param a 
 *   Input Stopband attenuation factor
 * @param omegar 
 *   Input Cutoff frequency of stopband passband cutoff is 1.0 Hz
 *
 * @return Nothing
 *
 * \copyright Copyright 1990  Regents of the University of California
 *
 * \author   Dave Harris
 *           Lawrence Livermore National Laboratory
 *           L-205
 *           P.O. Box 808
 *           Livermore, CA  94550
 *           USA
 *           (415) 423-0617
 *
 * \date 900907:  LAST MODIFIED
 *
 */
void c2roots(complexf *p, complexf *z, char *rtype, int rtype_s, float *dcvalue, int *nsects, int iord, double a, double omegar)
{
    #define RTYPE(I_, J_) (rtype + (I_) * (rtype_s) + (J_))
    int half, i, i_;
    float alpha, angle, beta, c, denom, gamma, omega, pi, s, sigma;

    complexf *const P = &p[0] - 1;
    complexf *const Z = &z[0] - 1;

    pi = 3.14159265;
    half = iord / 2;

    /*  INTERMEDIATE DESIGN PARAMETERS                                               
	 * */
    gamma = a + sqrt(a * a - 1.);
    gamma = log(gamma) / (float)(iord);
    gamma = exp(gamma);
    s = .5 * (gamma - 1. / gamma);
    c = .5 * (gamma + 1. / gamma);

    *nsects = 0;
    for (i = 1; i <= half; i++)
    {
        i_ = i - 1;

        /*  CALCULATE POLES                                                              
		 * */
        fstrncpy(RTYPE(i_, 0), rtype_s - 1, "CPZ", 2);

        angle = (float)(2 * i - 1) * pi / (float)(2 * iord);
        alpha = -s * sin(angle);
        beta = c * cos(angle);
        denom = alpha * alpha + beta * beta;
        sigma = omegar * alpha / denom;
        omega = -omegar * beta / denom;
        P[i] = flttocmplx(sigma, omega);

        /*  CALCULATE ZEROS                                                              
		 * */
        omega = omegar / cos(angle);
        Z[i] = flttocmplx(0.0, omega);

        *nsects = *nsects + 1;
    }

    /*  ODD-ORDER FILTERS                                                            
	 * */
    if (2 * half < iord)
    {
        fstrncpy(RTYPE(half, 0), rtype_s - 1, "SP", 2);
        P[half + 1] = flttocmplx(-omegar / s, 0.0);
        *nsects = *nsects + 1;
    }

    /*  DC VALUE                                                                     
	 * */
    *dcvalue = 1.0;

    /*  DONE                                                                         
	 * */
    return;
    #undef RTYPE
}

/** 
 *  Subroutine to generate second order section parameterization
 *    from an pole-zero description for lowpass filters.
 * 
 * 
 * @param p 
 *    Array of Poles
 * @param z 
 *    Array of Zeros
 * @param rtype 
 *    Character array containing root type information
 *      - "SP"  Single real pole
 *      - "CP"  Complex conjugate pole pair
 *      - "CPZ" Complex conjugate pole and zero pairs
 * @param rtype_s
 *     Length of \p rtype  
 * @param dcvalue 
 *     Zero-frequency value of prototype filter
 * @param nsects 
 *     Number of second-order sections
 * @param sn 
 *     Output Numerator polynomials for second order sections.
 * @param sd 
 *     Output Denominator polynomials for second order sections.
 *
 * @copyright  Copyright 1990  Regents of the University of California
 *
 * @author:  Dave Harris
 *           Lawrence Livermore National Laboratory
 *           L-205
 *           P.O. Box 808
 *           Livermore, CA  94550
 *           USA
 *           (415) 423-0617
 *
 * \date   020416: Changed SN and SD adjustable arrays to use 
 *                 "*" rather than "1". - wct
 */
void lp(complexf *p, complexf *z, char *rtype, int rtype_s, double dcvalue, int nsects, float *sn, float *sd)
{

    #define RTYPE(I_, J_) (rtype + (I_) * (rtype_s) + (J_))

    int i, i_, iptr;
    float scale;

    complexf *const P = &p[0] - 1;
    float *const Sd = &sd[0] - 1;
    float *const Sn = &sn[0] - 1;
    complexf *const Z = &z[0] - 1;

    iptr = 1;
    for (i = 1; i <= nsects; i++)
    {
        i_ = i - 1;

        if (memcmp(RTYPE(i_, 0), "CPZ", 3) == 0)
        {

            scale = cmplxtof(cmplxmul(P[i], cmplxcj(P[i]))) / cmplxtof(cmplxmul(Z[i],
                                                                                cmplxcj(Z[i])));
            Sn[iptr] = cmplxtof(cmplxmul(Z[i], cmplxcj(Z[i]))) * scale;
            Sn[iptr + 1] = -2. * cmplxtof(Z[i]) * scale;
            Sn[iptr + 2] = 1. * scale;
            Sd[iptr] = cmplxtof(cmplxmul(P[i], cmplxcj(P[i])));
            Sd[iptr + 1] = -2. * cmplxtof(P[i]);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;
        }
        else if (memcmp(RTYPE(i_, 0), "CP", 2) == 0)
        {

            scale = cmplxtof(cmplxmul(P[i], cmplxcj(P[i])));
            Sn[iptr] = scale;
            Sn[iptr + 1] = 0.;
            Sn[iptr + 2] = 0.;
            Sd[iptr] = cmplxtof(cmplxmul(P[i], cmplxcj(P[i])));
            Sd[iptr + 1] = -2. * cmplxtof(P[i]);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;
        }
        else if (memcmp(RTYPE(i_, 0), "SP", 2) == 0)
        {

            scale = -cmplxtof(P[i]);
            Sn[iptr] = scale;
            Sn[iptr + 1] = 0.;
            Sn[iptr + 2] = 0.;
            Sd[iptr] = -cmplxtof(P[i]);
            Sd[iptr + 1] = 1.;
            Sd[iptr + 2] = 0.;
            iptr = iptr + 3;
        }
    }

    Sn[1] = dcvalue * Sn[1];
    Sn[2] = dcvalue * Sn[2];
    Sn[3] = dcvalue * Sn[3];

    return;
    #undef RTYPE
}

/** 
 *  Subroutine to convert an prototype lowpass filter to a bandpass filter via
 *    the analog polynomial transformation.  The lowpass filter is
 *    described in terms of its poles and zeros (as input to this routine).
 *    The output consists of the parameters for second order sections.
 * 
 * @param p 
 *    Array of Poles
 * @param z 
 *    Array of Zeros
 * @param rtype 
 *    Character array containing root type information
 *      - "SP"  Single real pole
 *      - "CP"  Complex conjugate pole pair
 *      - "CPZ" Complex conjugate pole and zero pairs
 * @param rtype_s 
 *     Length of \p rtype  
 * @param dcvalue 
 *     Zero-frequency value of prototype filter
 * @param nsects 
 *     Number of second-order sections.
 *     On output this subroutine doubles the number of            
 *     sections.                                        
 * @param fl 
 *     Low Frequency cutoff
 * @param fh 
 *     High Frequency cutoff
 * @param sn 
 *     Output Numerator polynomials for second order sections.
 * @param sd 
 *     Output Denominator polynomials for second order sections.
 *
 * \bug This routine defines PI and TWO PI, which are available in math.h
 *
 * @copyright  Copyright 1990  Regents of the University of California
 *
 * @author:  Dave Harris
 *           Lawrence Livermore National Laboratory
 *           L-205
 *           P.O. Box 808
 *           Livermore, CA  94550
 *           USA
 *           (415) 423-0617
 *
 * 
 */
void lptbp(complexf *p, complexf *z, char *rtype, int rtype_s, double dcvalue, int *nsects, double fl, double fh, float *sn, float *sd)
{

    #define RTYPE(I_, J_) (rtype + (I_) * (rtype_s) + (J_))

    int idx, idx_, iptr, n;
    float a, b, pi, scale, twopi;
    complexf ctemp, h, p1, p2, s, z1, z2;

    complexf *const P = &p[0] - 1;
    float *const Sd = &sd[0] - 1;
    float *const Sn = &sn[0] - 1;
    complexf *const Z = &z[0] - 1;

    pi = 3.14159265;
    twopi = 2. * pi;
    a = twopi * twopi * fl * fh;
    b = twopi * (fh - fl);

    n = *nsects;
    *nsects = 0;
    iptr = 1;
    for (idx = 1; idx <= n; idx++)
    {
        idx_ = idx - 1;

        if (memcmp(RTYPE(idx_, 0), "CPZ", 3) == 0)
        {

            ctemp = cmplxsub(cmplxpow(
                                 (cmplxmul(flttocmplx(b, 0.), Z[idx])),
                                 (double)2),
                             flttocmplx(4. * a, 0.));

            ctemp = cmplxsqrt(ctemp);
            z1 = cmplxmul(flttocmplx(0.5, 0.),
                          (cmplxadd(cmplxmul(flttocmplx(b, 0.), Z[idx]),
                                    ctemp)));
            z2 = cmplxmul(flttocmplx(0.5, 0.),
                          (cmplxsub(cmplxmul(flttocmplx(b, 0.), Z[idx]),
                                    ctemp)));
            ctemp = cmplxsub(cmplxpow((cmplxmul(flttocmplx(b, 0.),
                                                P[idx])),
                                      (double)2),
                             flttocmplx(4. * a, 0.));
            ctemp = cmplxsqrt(ctemp);
            p1 = cmplxmul(flttocmplx(0.5, 0.),
                          (cmplxadd(cmplxmul(flttocmplx(b, 0.), P[idx]),
                                    ctemp)));
            p2 = cmplxmul(flttocmplx(0.5, 0.),
                          (cmplxsub(cmplxmul(flttocmplx(b, 0.), P[idx]),
                                    ctemp)));
            Sn[iptr] = cmplxtof(cmplxmul(z1, cmplxcj(z1)));
            Sn[iptr + 1] = -2. * cmplxtof(z1);
            Sn[iptr + 2] = 1.;
            Sd[iptr] = cmplxtof(cmplxmul(p1, cmplxcj(p1)));
            Sd[iptr + 1] = -2. * cmplxtof(p1);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;
            Sn[iptr] = cmplxtof(cmplxmul(z2, cmplxcj(z2)));
            Sn[iptr + 1] = -2. * cmplxtof(z2);
            Sn[iptr + 2] = 1.;
            Sd[iptr] = cmplxtof(cmplxmul(p2, cmplxcj(p2)));
            Sd[iptr + 1] = -2. * cmplxtof(p2);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;

            *nsects = *nsects + 2;
        }
        else if (memcmp(RTYPE(idx_, 0), "CP", 2) == 0)
        {

            ctemp = cmplxsub(cmplxpow(
                                 (cmplxmul(flttocmplx(b, 0.), P[idx])),
                                 (double)2),
                             flttocmplx(4. * a, 0.));
            ctemp = cmplxsqrt(ctemp);
            p1 = cmplxmul(flttocmplx(0.5, 0.),
                          (cmplxadd(cmplxmul(flttocmplx(b, 0.), P[idx]),
                                    ctemp)));
            p2 = cmplxmul(flttocmplx(0.5, 0.),
                          (cmplxsub(cmplxmul(flttocmplx(b, 0.), P[idx]),
                                    ctemp)));
            Sn[iptr] = 0.;
            Sn[iptr + 1] = b;
            Sn[iptr + 2] = 0.;
            Sd[iptr] = cmplxtof(cmplxmul(p1, cmplxcj(p1)));
            Sd[iptr + 1] = -2. * cmplxtof(p1);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;
            Sn[iptr] = 0.;
            Sn[iptr + 1] = b;
            Sn[iptr + 2] = 0.;
            Sd[iptr] = cmplxtof(cmplxmul(p2, cmplxcj(p2)));
            Sd[iptr + 1] = -2. * cmplxtof(p2);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;

            *nsects = *nsects + 2;
        }
        else if (memcmp(RTYPE(idx_, 0), "SP", 2) == 0)
        {

            Sn[iptr] = 0.;
            Sn[iptr + 1] = b;
            Sn[iptr + 2] = 0.;
            Sd[iptr] = a;
            Sd[iptr + 1] = -b * cmplxtof(P[idx]);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;

            *nsects = *nsects + 1;
        }
    }

    /*  Scaling - use the fact that the bandpass filter amplitude 
     *  at sqrt( omega_l 
     *  equals the amplitude of the lowpass prototype at d.c.
     * */
    s = flttocmplx(0., sqrt(a));
    h = flttocmplx(1., 0.);

    iptr = 1;
    for (idx = 1; idx <= *nsects; idx++)
    {
        h = cmplxdiv(cmplxmul(h, (cmplxadd(cmplxmul(
                                               (cmplxadd(cmplxmul(flttocmplx(Sn[iptr + 2], 0.), s),
                                                         flttocmplx(Sn[iptr + 1], 0.))),
                                               s),
                                           flttocmplx(Sn[iptr], 0.)))),
                     (cmplxadd(cmplxmul(
                                   (cmplxadd(cmplxmul(flttocmplx(Sd[iptr + 2], 0.), s),
                                             flttocmplx(Sd[iptr + 1], 0.))),
                                   s),
                               flttocmplx(Sd[iptr], 0.))));
        iptr = iptr + 3;
    }
    scale = cmplxtof(cmplxdiv(flttocmplx(dcvalue, 0.),
                              cmplxsqrt(cmplxmul(flttocmplx(cmplxtof(h), 0.),
                                                 cmplxcj(h)))));

    Sn[1] = Sn[1] * scale;
    Sn[2] = Sn[2] * scale;
    Sn[3] = Sn[3] * scale;

    return;
    #undef RTYPE
}

/** 
 * 
 *  Subroutine to convert a lowpass filter to a band reject filter
 *    via an analog polynomial transformation.  The lowpass filter is
 *    described in terms of its poles and zeros (as input to this routine).
 *    The output consists of the parameters for second order sections.
 * 
 * @param p 
 *    Array of Poles
 * @param z 
 *    Array of Zeros
 * @param rtype 
 *    Character array containing root type information
 *      - "SP"  Single real pole
 *      - "CP"  Complex conjugate pole pair
 *      - "CPZ" Complex conjugate pole and zero pairs
 * @param rtype_s 
 *     Length of \p rtype  
 * @param dcvalue 
 *     Zero-frequency value of prototype filter
 * @param nsects 
 *     Number of second-order sections.
 *     On output this subroutine doubles the number of            
 *     sections.                                        
 * @param fl 
 *     Low Frequency cutoff
 * @param fh 
 *     High Frequency cutoff
 * @param sn 
 *     Output Numerator polynomials for second order sections.
 * @param sd 
 *     Output Denominator polynomials for second order sections.
 *
 * \bug This routine defines PI and TWO PI, which are available in math.h
 *
 * @copyright  Copyright 1990  Regents of the University of California
 *
 * @author:  Dave Harris
 *           Lawrence Livermore National Laboratory
 *           L-205
 *           P.O. Box 808
 *           Livermore, CA  94550
 *           USA
 *           (415) 423-0617
 *
 */
void lptbr(complexf *p, complexf *z, char *rtype, int rtype_s, double dcvalue, int *nsects, double fl, double fh, float *sn, float *sd)
{

    #define RTYPE(I_, J_) (rtype + (I_) * (rtype_s) + (J_))

    int i, i_, iptr, n;
    float a, b, h, pi, scale, twopi;
    complexf cinv, ctemp, p1, p2, z1, z2;

    complexf *const P = &p[0] - 1;
    float *const Sd = &sd[0] - 1;
    float *const Sn = &sn[0] - 1;
    complexf *const Z = &z[0] - 1;

    pi = 3.14159265;
    twopi = 2. * pi;
    a = twopi * twopi * fl * fh;
    b = twopi * (fh - fl);

    n = *nsects;
    *nsects = 0;
    iptr = 1;
    for (i = 1; i <= n; i++)
    {
        i_ = i - 1;

        if (memcmp(RTYPE(i_, 0), "CPZ", 3) == 0)
        {

            cinv = cmplxdiv(flttocmplx(1., 0.), Z[i]);
            ctemp = cmplxsub(cmplxpow((cmplxmul(flttocmplx(b, 0.), cinv)), (double)2),
                             flttocmplx(4. * a, 0.));
            ctemp = cmplxsqrt(ctemp);
            z1 = cmplxmul(flttocmplx(0.5, 0.), (cmplxadd(cmplxmul(flttocmplx(b, 0.), cinv),
                                                         ctemp)));
            z2 = cmplxmul(flttocmplx(0.5, 0.), (cmplxsub(cmplxmul(flttocmplx(b, 0.), cinv),
                                                         ctemp)));
            cinv = cmplxdiv(flttocmplx(1., 0.), P[i]);
            ctemp = cmplxsub(cmplxpow((cmplxmul(flttocmplx(b, 0.), cinv)), (double)2),
                             flttocmplx(4. * a, 0.));
            ctemp = cmplxsqrt(ctemp);
            p1 = cmplxmul(flttocmplx(0.5, 0.), (cmplxadd(cmplxmul(flttocmplx(b, 0.), cinv),
                                                         ctemp)));
            p2 = cmplxmul(flttocmplx(0.5, 0.), (cmplxsub(cmplxmul(flttocmplx(b, 0.), cinv),
                                                         ctemp)));
            Sn[iptr] = cmplxtof(cmplxmul(z1, cmplxcj(z1)));
            Sn[iptr + 1] = -2. * cmplxtof(z1);
            Sn[iptr + 2] = 1.;
            Sd[iptr] = cmplxtof(cmplxmul(p1, cmplxcj(p1)));
            Sd[iptr + 1] = -2. * cmplxtof(p1);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;
            Sn[iptr] = cmplxtof(cmplxmul(z2, cmplxcj(z2)));
            Sn[iptr + 1] = -2. * cmplxtof(z2);
            Sn[iptr + 2] = 1.;
            Sd[iptr] = cmplxtof(cmplxmul(p2, cmplxcj(p2)));
            Sd[iptr + 1] = -2. * cmplxtof(p2);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;

            *nsects = *nsects + 2;
        }
        else if (memcmp(RTYPE(i_, 0), "CP", 2) == 0)
        {

            cinv = cmplxdiv(flttocmplx(1., 0.), P[i]);
            ctemp = cmplxsub(cmplxpow((cmplxmul(flttocmplx(b, 0.), cinv)), (double)2),
                             flttocmplx(4. * a, 0.));
            ctemp = cmplxsqrt(ctemp);
            p1 = cmplxmul(flttocmplx(0.5, 0.), (cmplxadd(cmplxmul(flttocmplx(b, 0.), cinv),
                                                         ctemp)));
            p2 = cmplxmul(flttocmplx(0.5, 0.), (cmplxsub(cmplxmul(flttocmplx(b, 0.), cinv),
                                                         ctemp)));
            Sn[iptr] = a;
            Sn[iptr + 1] = 0.;
            Sn[iptr + 2] = 1.;
            Sd[iptr] = cmplxtof(cmplxmul(p1, cmplxcj(p1)));
            Sd[iptr + 1] = -2. * cmplxtof(p1);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;
            Sn[iptr] = a;
            Sn[iptr + 1] = 0.;
            Sn[iptr + 2] = 1.;
            Sd[iptr] = cmplxtof(cmplxmul(p2, cmplxcj(p2)));
            Sd[iptr + 1] = -2. * cmplxtof(p2);
            Sd[iptr + 2] = 1.;
            iptr = iptr + 3;

            *nsects = *nsects + 2;
        }
        else if (memcmp(RTYPE(i_, 0), "SP", 2) == 0)
        {

            Sn[iptr] = a;
            Sn[iptr + 1] = 0.;
            Sn[iptr + 2] = 1.;
            Sd[iptr] = -a * cmplxtof(P[i]);
            Sd[iptr + 1] = b;
            Sd[iptr + 2] = -cmplxtof(P[i]);
            iptr = iptr + 3;

            *nsects = *nsects + 1;
        }
    }

    /*  Scaling - use the fact that the bandreject filter amplitude at d.c.
	 *            equals the lowpass prototype amplitude at d.c.
	 * */
    h = 1.0;

    iptr = 1;
    for (i = 1; i <= *nsects; i++)
    {
        h = h * Sn[iptr] / Sd[iptr];
        iptr = iptr + 3;
    }
    scale = dcvalue / fabs(h);
    Sn[1] = Sn[1] * scale;
    Sn[2] = Sn[2] * scale;
    Sn[3] = Sn[3] * scale;

    return;
    #undef RTYPE
}

/*
 *  Subroutine to convert a lowpass filter to a highpass filter via
 *    an analog polynomial transformation.  The lowpass filter is
 *    described in terms of its poles and zeroes (as input to this routine).
 *    The output consists of the parameters for second order sections.
 *
 * @param p 
 *    Array of Poles
 * @param z 
 *    Array of Zeros
 * @param rtype 
 *    Character array containing root type information
 *      - "SP"  Single real pole
 *      - "CP"  Complex conjugate pole pair
 *      - "CPZ" Complex conjugate pole and zero pairs
 * @param rtype_s 
 *     Length of \p rtype  
 * @param dcvalue 
 *     Zero-frequency value of prototype filter
 * @param nsects 
 *     Number of second-order sections.
 * @param sn 
 *     Output Numerator polynomials for second order sections.
 * @param sd 
 *     Output Denominator polynomials for second order sections.
 *
 * \bug This routine defines PI and TWO PI, which are available in math.h
 *
 * @copyright  Copyright 1990  Regents of the University of California
 *
 * @author:  Dave Harris
 *           Lawrence Livermore National Laboratory
 *           L-205
 *           P.O. Box 808
 *           Livermore, CA  94550
 *           USA
 *           (415) 423-0617
 *
 * 
 */

void lpthp(complexf *p, complexf *z, char *rtype, int rtype_s, double dcvalue, int nsects, float *sn, float *sd)
{

    #define RTYPE(I_, J_) (rtype + (I_) * (rtype_s) + (J_))

    int i, i_, iptr;
    float scale;

    complexf *const P = &p[0] - 1;
    float *const Sd = &sd[0] - 1;
    float *const Sn = &sn[0] - 1;
    complexf *const Z = &z[0] - 1;

    iptr = 1;
    for (i = 1; i <= nsects; i++)
    {
        i_ = i - 1;

        if (memcmp(RTYPE(i_, 0), "CPZ", 3) == 0)
        {

            scale = cmplxtof(cmplxmul(P[i], cmplxcj(P[i]))) / cmplxtof(cmplxmul(Z[i],
                                                                                cmplxcj(Z[i])));
            Sn[iptr] = 1. * scale;
            Sn[iptr + 1] = -2. * cmplxtof(Z[i]) * scale;
            Sn[iptr + 2] = cmplxtof(cmplxmul(Z[i], cmplxcj(Z[i]))) * scale;
            Sd[iptr] = 1.;
            Sd[iptr + 1] = -2. * cmplxtof(P[i]);
            Sd[iptr + 2] = cmplxtof(cmplxmul(P[i], cmplxcj(P[i])));
            iptr = iptr + 3;
        }
        else if (memcmp(RTYPE(i_, 0), "CP", 2) == 0)
        {

            scale = cmplxtof(cmplxmul(P[i], cmplxcj(P[i])));
            Sn[iptr] = 0.;
            Sn[iptr + 1] = 0.;
            Sn[iptr + 2] = scale;
            Sd[iptr] = 1.;
            Sd[iptr + 1] = -2. * cmplxtof(P[i]);
            Sd[iptr + 2] = cmplxtof(cmplxmul(P[i], cmplxcj(P[i])));
            iptr = iptr + 3;
        }
        else if (memcmp(RTYPE(i_, 0), "SP", 2) == 0)
        {

            scale = -cmplxtof(P[i]);
            Sn[iptr] = 0.;
            Sn[iptr + 1] = scale;
            Sn[iptr + 2] = 0.;
            Sd[iptr] = 1.;
            Sd[iptr + 1] = -cmplxtof(P[i]);
            Sd[iptr + 2] = 0.;
            iptr = iptr + 3;
        }
    }

    Sn[1] = Sn[1] * dcvalue;
    Sn[2] = Sn[2] * dcvalue;
    Sn[3] = Sn[3] * dcvalue;

    return;
    #undef RTYPE
}

/** 
 * Transform an analog filter to a digital filter via the bilinear 
 *    transformation.  Assumes both filters are stored as 2nd Order
 *    sections and the transform is done in place.
 * 
 * @param sn 
 *    Array containing numerator polynomial coefficeients. Packed head to
 *      tail and using groups of 3.  Length is 3 * \p nsects
 * @param sd 
 *    Array containing demoninator polynomial coefficeients. Packed head to
 *      tail and using groups of 3.  Length is 3 * \p nsects
 * @param nsects 
 *    Number of 2nd order sections.
 * 
 * @return Nothing
 *
 * @copyright 1990  Regents of the University of California
 *
 * @author  Dave Harris                                                         
 *           Lawrence Livermore National Laboratory
 *           L-205
 *           P.O. Box 808
 *           Livermore, CA  94550
 *           USA
 *           (415) 423-0617
 *
 * \note
 *    - scale = a0 + a1 + a2
 *    - 
 *    - a_0 = 1.0
 *    - a_1 = 2.0 * (a_0 - a_2) / scale
 *    - a_2 = (a_2 - a_1 + a_0) / scale
 *    - 
 *    - b_0 = (b_2 + b_1 + b_0) / scale
 *    - b_1 = 2.0 * (b_0 - b_2) / scale
 *    - b_2 = (b_2 - b_1 + b_0) / scale
 *
 * \todo Further Documentation 
 */
void bilin2(float *sn, float *sd, int nsects)
{
    int i, iptr;
    double a0, a1, a2, scale;

    float *const Sd = &sd[0] - 1;
    float *const Sn = &sn[0] - 1;

    iptr = 1;
    for (i = 1; i <= nsects; i++)
    {

        a0 = Sd[iptr];
        a1 = Sd[iptr + 1];
        a2 = Sd[iptr + 2];

        scale = a2 + a1 + a0;
        Sd[iptr] = 1.;
        Sd[iptr + 1] = (2. * (a0 - a2)) / scale;
        Sd[iptr + 2] = (a2 - a1 + a0) / scale;

        a0 = Sn[iptr];
        a1 = Sn[iptr + 1];
        a2 = Sn[iptr + 2];

        Sn[iptr] = (a2 + a1 + a0) / scale;
        Sn[iptr + 1] = (2. * (a0 - a2)) / scale;
        Sn[iptr + 2] = (a2 - a1 + a0) / scale;

        iptr = iptr + 3;
    }

    return;
}



/**
* Applies tangent frequency warping to compensate
*   for bilinear analog -> digital transformation
*
* @param f
*   Original Design Frequency Specification (Hz)
* @param ts
*   Sampling Internal (seconds)
*
* @return
*
* @date 200990  Last Modified:  September 20, 1990
*
* @opyright   1990  Regents of the University of California
*
* @author   Dave Harris
*           Lawrence Livermore National Laboratory
*           L-205
*           P.O. Box 808
*           Livermore, CA  94550
*           USA
*           (415) 423-0617
*
*/
double warp(double f, double ts)
{
	float angle, twopi, warp_v;

	twopi = 2.0 * M_PI;
	angle = twopi*f*ts / 2.;
	warp_v = 2.*tan(angle) / ts;
	warp_v = warp_v / twopi;

	return(warp_v);
}


/**
* Alter the cutoff of a filter.  Assumed that the filter
* is structured as 2nd order sections.  Changes the cutoffs
* of a normalized lowpass or highpass filters through a
* simple polynomial transformation.
*
* @param sn
*   Numerator polynomials for 2nd order sections
* @param sd
*   Denominator polynomials for 2nd order sections
* @param nsects
*   Number of 2nd order sections
* @param f
*   New cutoff frequency
*
* @return Nothing
*
* \copyright 1990  Regents of the University of California
*
* \author   Dave Harris
*           Lawrence Livermore National Laboratory
*           L-205
*           P.O. Box 808
*           Livermore, CA  94550
*           USA
*           (415) 423-0617
*
*/
void cutoffs(float *sn, float *sd, int nsects, double f)
{
	int i, iptr;
	float scale;

	float *const Sd = &sd[0] - 1;
	float *const Sn = &sn[0] - 1;

	scale = 2.*3.14159265*f;

	iptr = 1;
	for (i = 1; i <= nsects; i++) {
		Sn[iptr + 1] = Sn[iptr + 1] / scale;
		Sn[iptr + 2] = Sn[iptr + 2] / (scale*scale);
		Sd[iptr + 1] = Sd[iptr + 1] / scale;
		Sd[iptr + 2] = Sd[iptr + 2] / (scale*scale);
		iptr = iptr + 3;
	}

	return;
}



/**
* Copy a string from \p from to \p to.  Spaces are padded at the end
*    of the to string
*
* @param to
*    Where to place the string
* @param tolen
*    Length of \p to
* @param from
*    Where to copy the string from
* @param fromlen
*    Length of \p from
*
* @return
*    The copied string
*
* @bug Should be replaced with strncpy() or equivalent
*
*/
char *fstrncpy(char *to, int tolen, std::string s_from,	int fromlen) {
	int cpylen;
    const char *from = s_from.c_str();

	if (to == NULL || from == NULL || tolen <= 0 || fromlen <= 0)
		return(NULL);

	cpylen = fromlen;
	if (fromlen > tolen) cpylen = tolen;

	memcpy(to, from, cpylen);
	if (cpylen < tolen)
		memset(to + cpylen, (int)' ', tolen - cpylen);
	to[tolen] = '\0';

	return(to);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */ 
