#define c 3e10
#define h1 1.0546e-34
#define h11 1.0546e-27//Plank's const for heaviside-lorentz units
#define eps0 8.85e-12
#define mu0 1.26e-6
#define E_from_SI_to_HL 1e-4/3
#define dab_from_SI_to_HL 1e2*3*1e9
#define na_from_SI_to_HL 1e-6

struct two_level_params
{
	double     p1;
	double    p2;
	double p3;

	double n1;
	double n2 ;
	double n3;


	double kp;
	double ke;
	double kn;
};


struct FdtdFields
{
	 double E;
	 double H;
	 double B;
	 double D;
	 double P;
	 double Pold;
	 double Pold2;
	 double N;	
};


struct SourceFields
{
	double Ein;
	double Hin;
};

struct PMLParams
{
	double sigma;
	double k;
	double a1;
	double b1;
	double c1;
	double d1;
};

struct MediaParamsIndex
{
	int two_level_par_index;
	int eps_r;
	int pump;
};

