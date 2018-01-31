#ifndef __IIR_FILTER__
#define __IIR_FILTER__

#include <iostream>
#include "iirDesign.hh"
#include "iirApply.hh"
#include <string>

class iirfilter
{
public:
	iirfilter(IIR_TYPE type, IIR_APROTO aproto, int ord, double a, double trbndw, double f1, double f2, double dt)
		: d_type(type), d_aproto(aproto), d_ord(ord), d_a(a), d_trbndw(trbndw), d_f1(f1), d_f2(f2), d_dt(dt)
	{
		init(d_type, d_aproto, d_ord, d_a, d_trbndw, d_f1, d_f2, d_dt);
	}
	~iirfilter() {};

	int init(IIR_TYPE type, IIR_APROTO aproto, int ord, double a, double trbndw, double f1, double f2, double dt) {
		return design(ord, type, aproto, a, trbndw, f1, f2, dt, d_sn, d_sd, &d_nsects);
	}
	int filter(float *data, int len, int passes) {
		return apply(data, len, ((passes == 1) ? false : true), d_sn, d_sd, d_nsects);
	}
	int info() {
		static std::string ty_lst[4] = {"LP", "HP", "BP", "BR"};
		static std::string ap_lst[4] = {"Butter", "Bessel", "Chebyshev1", "Chebyshev2"};
		std::cout << "Type: " << ty_lst[d_type] << "\n";
		std::cout << "Apro: " << ap_lst[d_aproto] << "\n";
		std::cout << "Order:" << d_ord << "\n";
		std::cout << "Band: [ " << d_f1 << " , " << d_f2 << " ] (Hz)\n";
		return 0;
	}

private:
	IIR_APROTO d_aproto;  // BU, BE, C1, C2
	double d_trbndw;      // transition bandwidth as fraction of lowpass prototype filter cutoff frequency. used only by chebyshev filters.
	double d_a;           // attenuation factor.  equals amplitude reached at stopband edge. used only by chebyshev filters.
	int    d_ord;         // order, not to exceed 10
	IIR_TYPE d_type;      // LP, HP, BP, BR
	double d_f1;          // low frequency cutoff of filter (hertz), ignored if type = 'lp'
	double d_f2;          // high frequency cutoff of filter (hertz), ignored if type = 'hp'
	double d_dt;          // sampling interval (seconds)
	//int    d_passes;      // integer variable containing the number of passes
					      //   1 --forward filtering only
					      //   2 --forward and reverse(i.e.zero phase) filtering
private:
	float d_sn[30];
	float d_sd[30];
	int   d_nsects;
};
// BUTTER ////////////////////////////////////////////////////////////////
class iirButterLP : public iirfilter {
public:
	iirButterLP(int ord, double f, double dt) : iirfilter::iirfilter(IIR_LP, IIR_BU, ord, 0.0, 0.0, 0.0, f, dt)	{ }
	~iirButterLP() {}
};
class iirButterHP : public iirfilter {
public:
	iirButterHP(int ord, double f, double dt) : iirfilter::iirfilter(IIR_HP, IIR_BU, ord, 0.0, 0.0, f, 0.0, dt)	{ }
	~iirButterHP() {}
};
class iirButterBP : public iirfilter {
public:
	iirButterBP(int ord, double f1, double f2, double dt) : iirfilter::iirfilter(IIR_BP, IIR_BU, ord, 0.0, 0.0, f1, f2, dt) { }
	~iirButterBP() {}
};
class iirButterBR : public iirfilter {
public:
	iirButterBR(int ord, double f1, double f2, double dt) : iirfilter::iirfilter(IIR_BR, IIR_BU, ord, 0.0, 0.0, f1, f2, dt)	{ }
	~iirButterBR() {}
};
// BESSEL ////////////////////////////////////////////////////////////////
class iirBesselLP : public iirfilter {
public:
	iirBesselLP(int ord, double f, double dt) : iirfilter::iirfilter(IIR_LP, IIR_BE, ord, 0.0, 0.0, 0.0, f, dt) { }
	~iirBesselLP() {}
};
class iirBesselHP : public iirfilter {
public:
	iirBesselHP(int ord, double f, double dt) : iirfilter::iirfilter(IIR_HP, IIR_BE, ord, 0.0, 0.0, f, 0.0, dt) { }
	~iirBesselHP() {}
};
class iirBesselBP : public iirfilter {
public:
	iirBesselBP(int ord, double f1, double f2, double dt) : iirfilter::iirfilter(IIR_BP, IIR_BE, ord, 0.0, 0.0, f1, f2, dt) { }
	~iirBesselBP() {}
};
class iirBesselBR : public iirfilter {
public:
	iirBesselBR(int ord, double f1, double f2, double dt) : iirfilter::iirfilter(IIR_BR, IIR_BE, ord, 0.0, 0.0, f1, f2, dt) { }
	~iirBesselBR() {}
};
// C1 ////////////////////////////////////////////////////////////////
class iirCheb1_LP : public iirfilter {
public:
	iirCheb1_LP(int ord, double f, double dt, double a, double trbndw) : iirfilter::iirfilter(IIR_LP, IIR_C1, ord, a, trbndw, 0.0, f, dt) { }
	~iirCheb1_LP() {}
};
class iirCheb1_HP : public iirfilter {
public:
	iirCheb1_HP(int ord, double f, double dt, double a, double trbndw) : iirfilter::iirfilter(IIR_HP, IIR_C1, ord, a, trbndw, f, 0.0, dt) { }
	~iirCheb1_HP() {}
};
class iirCheb1_BP : public iirfilter {
public:
	iirCheb1_BP(int ord, double f1, double f2, double dt, double a, double trbndw) : iirfilter::iirfilter(IIR_BP, IIR_C1, ord, a, trbndw, f1, f2, dt) { }
	~iirCheb1_BP() {}
};
class iirCheb1_BR : public iirfilter {
public:
	iirCheb1_BR(int ord, double f1, double f2, double dt, double a, double trbndw) : iirfilter::iirfilter(IIR_BR, IIR_C1, ord, a, trbndw, f1, f2, dt) { }
	~iirCheb1_BR() {}
};
// C2 ////////////////////////////////////////////////////////////////
class iirCheb2_LP : public iirfilter {
public:
	iirCheb2_LP(int ord, double f, double dt, double a, double trbndw) : iirfilter::iirfilter(IIR_LP, IIR_C2, ord, a, trbndw, 0.0, f, dt) { }
	~iirCheb2_LP() {}
};
class iirCheb2_HP : public iirfilter {
public:
	iirCheb2_HP(int ord, double f, double dt, double a, double trbndw) : iirfilter::iirfilter(IIR_HP, IIR_C2, ord, a, trbndw, f, 0.0, dt) { }
	~iirCheb2_HP() {}
};
class iirCheb2_BP : public iirfilter {
public:
	iirCheb2_BP(int ord, double f1, double f2, double dt, double a, double trbndw) : iirfilter::iirfilter(IIR_BP, IIR_C2, ord, a, trbndw, f1, f2, dt) { }
	~iirCheb2_BP() {}
};
class iirCheb2_BR : public iirfilter {
public:
	iirCheb2_BR(int ord, double f1, double f2, double dt, double a, double trbndw) : iirfilter::iirfilter(IIR_BR, IIR_C2, ord, a, trbndw, f1, f2, dt) { }
	~iirCheb2_BR() {}
};
//////////////////////////////////////////////////////////////////
void xapiir(float *data, int nsamps, IIR_APROTO aproto, double trbndw, double a, int iord, IIR_TYPE type, double flo, double fhi, double ts, int passes);

#endif // !__IIR_FILTER__
