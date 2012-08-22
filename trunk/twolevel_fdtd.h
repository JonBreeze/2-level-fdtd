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

