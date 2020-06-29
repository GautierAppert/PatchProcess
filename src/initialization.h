#include <Rcpp.h>

using namespace Rcpp;

class mergeParam { 
	public:
	mergeParam ( 
  	IntegerVector &A_v, 
		IntegerVector &Asizes_v,
  	IntegerVector &B_v, 
		IntegerVector &Bsizes_v,
		IntegerVector &Twx_v,
		IntegerVector &Jt_v,
  	List &paramList ) :
		A(&A_v[0]), 
		Asizes(&Asizes_v[0]), 
		B(&B_v[0]),
		Bsizes(&Bsizes_v[0]),
		Twx(&Twx_v[0]),
		Jt(&Jt_v[0]),
		n(as<int>(paramList["n"])), // number of images
		d(as<int>(paramList["d"])), // number of pixels
		betaCoeff(as<double>(paramList["betaCoeff"])),
		firstNewLabel(as<int>(paramList["firstNewLabel"])),
    lastNewLabel(as<int>(paramList["lastNewLabel"])),
    firstMerge(as<bool>(paramList["firstMerge"])),
		shift(0) {
			Crit = new double[lastNewLabel];
		}
	~mergeParam(void) {
		delete[] Crit;
	}
	void merge (int iter); 
	void output(void);
	void initialization(void); 
	double *Crit;
	int *A; 
	int *Asizes;
	int *B;
	int *Bsizes;
	int *Twx;
	int *Jt;
	int *Cont;
	int n;
	int d;
	double betaCoeff;
	int firstNewLabel;
	int lastNewLabel;
	bool firstMerge;
	int shift;
};


class Patches {
	public:
	Patches ( 
  NumericVector &m_v,
  NumericVector &v_v,
  NumericVector &w_v,
  NumericVector &C_v,
  IntegerVector &A_v,
  IntegerVector &B_v,
	IntegerVector &Wsx_v,
  List &paramList ) : 
		m(&m_v[0]),
		v(&v_v[0]),
		w(&w_v[0]),
		C(&C_v[0]),
		A(&A_v[0]), 
		B(&B_v[0]),
		Wsx(&Wsx_v[0]),
		betaCoeff(as<double>(paramList["betaCoeff"])),
		thresholdCoeff(as<double>(paramList["thresholdCoeff"])),
		thresholdValue(as<double>(paramList["thresholdValue"])),
		n(as<int>(paramList["n"])),
		d(as<int>(paramList["d"])),
		firstNewLabel(as<int>(paramList["firstNewLabel"])),
		lastNewLabel(as<int>(paramList["lastNewLabel"])),
		shift(0)
		{}	
	void computeBeta(void);
	void output(void);
	void split (int iter);
	double *m;
	double *v;
	double *w;
	double *C;
	int *A;
	int *B;
	int *Wsx;
	double betaCoeff;
	double thresholdValue;
	double thresholdCoeff;
	int n;
	int d;
	int firstNewLabel;
	int lastNewLabel;
	int shift;
};


