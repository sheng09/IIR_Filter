#include "iirFilter.hh"
#include "iirDesign.hh"
#include "complex.hh"
#include "iirApply.hh"

#include <cstring>





/** 
 *  IIR filter design and implementation
 * 
 * @param data 
 *    real array containing sequence to be filtered
 *    original data destroyed, replaced by filtered data
 * @param nsamps 
 *    number of samples in data
 * @param aproto 
 *    character*8 variable, contains type of analog prototype filter
 *      - '(BU)tter  ' -- butterworth filter
 *      - '(BE)ssel  ' -- bessel filter
 *      - 'C1      ' -- chebyshev type i
 *      - 'C2      ' -- chebyshev type ii
 * @param trbndw 
 *    transition bandwidth as fraction of lowpass
 *    prototype filter cutoff frequency.  used
 *    only by chebyshev filters.
 * @param a 
 *    attenuation factor.  equals amplitude
 *    reached at stopband edge.  used only by
 *    chebyshev filters.
 * @param iord 
 *    order (#poles) of analog prototype
 *    not to exceed 10 in this configuration.  4 - 5
 *    should be ample.
 * @param type 
 *    character*8 variable containing filter type
 *     - 'LP' -- low pass
 *     - 'HP' -- high pass
 *     - 'BP' -- band pass
 *     - 'BR' -- band reject
 * @param flo 
 *    low frequency cutoff of filter (hertz)
 *    ignored if type = 'lp'
 * @param fhi 
 *    high frequency cutoff of filter (hertz)
 *    ignored if type = 'hp'
 * @param ts 
 *    sampling interval (seconds)
 * @param passes 
 *    integer variable containing the number of passes
 *     - 1 -- forward filtering only
 *     - 2 -- forward and reverse (i.e. zero phase) filtering
 *
 * @author:  Dave B. Harris
 *
 * @date 120990 Last Modified:  September 12, 1990
 */
void xapiir(float *data, int nsamps, IIR_APROTO aproto,
			double trbndw, double a, int iord, IIR_TYPE type,
            double flo, double fhi, double ts, int passes) {
    bool zp;
    int nsects;
    float sd[30], sn[30];

    design(iord, type, aproto, a, trbndw, flo, fhi, ts, sn, sd, &nsects);

    /*  Filter data  */
    apply(data, nsamps, ((passes == 1) ? false : true), sn, sd, nsects);
    return;
}
