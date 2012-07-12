// twolevel_fdtd.cpp: определяет точку входа для консольного приложения.
//

#include <stdio.h>

#define M_PI       3.14159265358979323846
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <exception>


using namespace std;





//global vars

int time_points, space_points, media_start, media_end, media_number;
int time_samps;
int  ref_start, ref_length, ref_layers,ref_config;
static double omega,gama,na,delta,delta1,gama_nr,dab,N0,N11,Tp,w0,E0,step_lambda_rate;
double reflector_index;
int z_out_time;


void readln(FILE * f)
{
	while (fgetc(f) != '\n');
}

int  read_parameters()
{

	FILE *infile = NULL;
	infile = fopen("config.txt", "r");

	if (infile == NULL)
	{
		printf("file not found, exitting\n");
		return 0;
	}

	readln(infile);
	fscanf(infile,"%d\n", &time_points);
	readln(infile);
	fscanf(infile,"%d\n",&space_points);
	readln(infile);
	fscanf(infile,"%d\n",&media_start);
	readln(infile);
	fscanf(infile,"%d\n",&media_end);
	readln(infile);
	fscanf(infile,"%d\n",&media_number);
	readln(infile);
	fscanf(infile,"%d\n",&time_samps);
	readln(infile);
	fscanf(infile,"%d\n", &ref_start);
	readln(infile);
	fscanf(infile,"%d\n", &ref_length);
	readln(infile);
	fscanf(infile,"%d\n", &ref_layers);
	readln(infile);
	fscanf(infile,"%lf\n",&reflector_index);
	readln(infile);
	fscanf(infile,"%d\n",&ref_config);
	readln(infile);
	fscanf(infile,"%lf\n",&omega);
	readln(infile);
	fscanf(infile,"%lf\n",&gama);
	readln(infile);
	fscanf(infile,"%lf\n",&na);
	readln(infile);
	fscanf(infile,"%lf\n",&delta);
	readln(infile);
	fscanf(infile,"%lf\n",&delta1);
	readln(infile);
	fscanf(infile,"%lf\n",&gama_nr);
	readln(infile);
	fscanf(infile,"%lf\n",&dab);
	readln(infile);
	fscanf(infile,"%lf\n",&N0);
	readln(infile);
	fscanf(infile,"%lf\n",&N11);
	readln(infile);
	fscanf(infile,"%lf\n",&Tp);
	readln(infile);
	fscanf(infile,"%lf\n",&w0);
	readln(infile);
	fscanf(infile,"%lf\n",&E0);
	readln(infile);
	fscanf(infile,"%lf\n",&step_lambda_rate);
	readln(infile);
	fscanf(infile,"%d\n",&z_out_time);
	fclose(infile);
	return 1;

}


int main(int argc, char* argv[])
{
	time_t start = time(NULL);
	read_parameters();

	int time_i = 0;
	int Ntime = time_points / time_samps - 1 ;

	double c=3e10;
	double h1=1.0546e-34;
	double h11=1.0546e-27;//Plank's const for heaviside-lorentz units

	double lambda=2*M_PI*c/w0;
	double dz=step_lambda_rate*lambda;
	double dt=0.5*dz/c;
	double eps0=8.85e-12;
	double mu0=1.26e-6;
	double eta0 = sqrt(mu0/eps0);
	//Tp=2*PI*10/w0;
	double koren= sqrt(eps0);
	double E_from_SI_to_HL=(1/3)*1e-4;
	double dab_from_SI_to_HL=1e2*3*1e9;
	double na_from_SI_to_HL=1e-6;
	E0= E0*E_from_SI_to_HL;
	dab=dab*dab_from_SI_to_HL;
	na=na*na_from_SI_to_HL;
	double Hup=c*dt/dz;
	double Eup= c*dt/dz;
	double Dup=c*dt/dz;


	double omega0sqr = gama*gama+omega*omega;
	int Nabc = 50;      //Number of PML cells
	double Labc = Nabc*dz;
	int pulse_points = time_points;

	//allocate arrays
	int res = 0;
	vector<double> E(space_points+2*Nabc, 0);
	vector<double> H(space_points+2*Nabc, 0);
	vector<double> B(space_points+2*Nabc, 0);
	vector<double> D(space_points+2*Nabc, 0);
	vector<double> P(space_points+1+2*Nabc, 0);
	vector<double> Pold(space_points+1+2*Nabc,0);
	vector<double> Pold2(space_points+1+2*Nabc,0);
	vector<double> N(space_points+1+2*Nabc, 0);
	vector<double> Ein(pulse_points, 0);
	vector<double> Hin(pulse_points, 0);


	vector<double> sigma(Nabc, 0);
	vector<double> k(Nabc, 0);
	vector<double> a1(Nabc, 0);
	vector<double> b1(Nabc, 0);
	vector<double> c1(Nabc, 0);
	vector<double> d1(Nabc, 0);
	vector<int> eps_r(space_points+2*Nabc, 0);
	vector<int> pump(space_points+2*Nabc, 0);
	//store inverted epsilon data to speedup computations
	double media[3]= {1.0,1/reflector_index,1/(0.2*reflector_index)};
	double media_pump[3]= {0,delta,delta1};


	//wid=round(lambda/(4*dz*sqrt(media[1])));
	long wid2=(long)(lambda/(4*dz*sqrt(media[2])));
	long wid=wid2;
	//layers=20;

	switch(ref_config)
	{
	case  1 :
		{
			for (int j= 0; j<ref_layers;j++)
			{
				//left mirror--------
				//first material
				for (int i= ref_start+(wid+wid2)*j;i<=ref_start+(wid+wid2)*j+wid;i++)
					eps_r[i]=1;
				//second material
				for (int i= ref_start+(wid+wid2)*j+wid;i<=ref_start+(wid+wid2)*j+wid+wid2;i++)
					eps_r[i]=2;
				//right mirror-------
				for (int i= ref_start+2*wid*(ref_layers-1)+wid+ref_length+2*wid*j; i<= ref_start+2*wid*(ref_layers-1)+wid+ref_length+2*wid*j+wid; i++)
					eps_r[i]=1;
			}
			break;
		}

	case 2:
		{
			for (int j= 0 ; j<=ref_layers-1; j++)
			{
				//left mirror--------
				//first material
				for (int i= ref_start+(wid+wid2)*j;i<=ref_start+(wid+wid2)*j+wid;i++)
					eps_r[i]=1;
				//second material
				for (int i= ref_start+(wid+wid2)*j+wid;i<=ref_start+(wid+wid2)*j+wid+wid2;i++)
					eps_r[i]=2;
			}

		}
		break;
	default:
		{
			long media_length= (long)(media_end - media_start)*media_number;
			//bulk media with reflecting edges
			for (int i=0; i<=space_points-1+Nabc+Nabc; i++)
			{
				eps_r[i]= 0;
			}

			for (int j= media_start-media_length;j<=media_end+media_length-1;j++)
			{
				eps_r[j]= 1;
			}
		}


	}



	FILE *eps = fopen("eps.txt","wt");
	for(int i=0; i<space_points+2*Nabc; i++)
	{
		fprintf(eps,"%ld\n", eps_r[i]);
	}
	fclose(eps);

	//media with absorbing cell
	long delta_number=(long)((2/3)*(media_end-media_start));

	for (int j=media_start; j<=media_start+delta_number-1;j++)
	{
		pump[j]=1;
	}

	for (int j=media_start+delta_number; j<=media_end-1-delta_number;j++)
	{
		pump[j]=2;
	}

	for (int j=media_end-delta_number;j<=media_end-1;j++)
	{
		pump[j]=1;
	}

	/*
	for (int j= 0;j<=ref_layers-1; j++)
	{
	for (int i= ref_start+2*wid*j; ref_start+2*wid*j+wid;j++)
	eps_r[i]=1;
	for (int i= ref_start+2*wid*(ref_layers-1)+wid+ref_length+2*wid*j;i<=ref_start+2*wid*(ref_layers-1)+wid+ref_length+2*wid*j+wid;i++)
	eps_r[i]=1;
	}


	for j= 0 to layers-1 do
	for i= start+2*wid*(layers-1)+wid+length+2*wid*j to start+2*wid*(layers-1)+wid+length+2*wid*j+wid do
	eps_r[i]=1;
	*/

	/*
	for i=0 to space_points-1+Nabc+Nabc do
	writeln(eps,' ', 1/(media[eps_r[i]]));
	CloseFile(eps);
	*/
	

	for (int i = media_start; i<=media_end-1;i++)
	{
		N[i]=N11;
	}

	srand(  time(NULL));
	for (int i = media_start; i<= media_end-1; i++)
	{
		P[i]=sin(w0*i*dt/*+2*M_PI*((double)rand()/(double)RAND_MAX)*/);
	}

	FILE *pstart = fopen("po.txt","wt");
	for(int i=0; i< space_points+2*Nabc+1; i++)
	{
		fprintf(pstart,"%.15e\n",P[i]);
	}
	fclose(pstart);
	

	double Hdelay  = 0;//-dt/2;//-dz/2/c;


	//hyperbolic secant pulse
	/*
	{ for nq=0 to pulse_points-1 do
	{
	t=nq*dt;
	//G=0.5*(t-Tp)/(Tp);
	Ein[nq] = E0*sin( w0*t )/cosh(10*(t-0.75*Tp)/Tp);
	Hin[nq] = E0*sin( w0*(t+Hdelay) )/cosh(10*(t+Hdelay-0.75*Tp)/Tp);
	writeln(epulse,  Ein[nq], Hin[nq]);
	}
	closeFile(epulse);}
	*/

	//pml parameters calculation
	double abc_pow = 3;
	double kmax=2;
	double sigma_max = -(abc_pow+1)*log(0.0001)*c/(2*/*eta0**/Labc); //2e16;
	for (int i=0; i<= Nabc-1; i++)
	{
		sigma[i] =  sigma_max*pow(i/Nabc, abc_pow)/**mu0/eps0*/;
		k[i]=1+(kmax-1)*pow(i/Nabc, abc_pow);
		//sigma[i]=0;
		//k[i]=1;
		a1[i]= k[i]+0.5*sigma[i]*dt;
		b1[i]= 1/a1[i];//1/(k[i]+0.5*sigma[i]*dt);
		c1[i]= k[i]-0.5*sigma[i]*dt;
		d1[i]= c1[i]/a1[i];//(k[i]-0.5*sigma[i]*dt)/(k[i]+0.5*sigma[i]*dt);

		//writeln(polarization,' ', a1[i],'  ',b1[i],' ',c1[i],'  ',d1[i],' ',sigma[i],' ',k[i]);
	}


	//pumping rate cmdline input
	/*
	if ParamCount>0 then
	delta = StrToFloat(ParamStr(1));
	*/
	//precalculated media parameters
	const double     p1=(2-omega0sqr*dt*dt)/(1+gama*dt);
	const double    p2=(1-gama*dt)/(1+gama*dt);
	const double p3=omega0sqr*dt*dt/(h11*omega*(1+gama*dt));

	const double n1=(2-gama_nr*dt)/(2+gama_nr*dt);
	const double n2 = 2*dt*(-2*delta+gama_nr*N0)/(2+gama_nr*dt);
	const double n3=4*omega/(h11*omega0sqr*(2+gama_nr*dt));


	const double kp= p3*dab;
	const double ke= 2*na*dab;
	const double kn= n3*dab;

	//excitation pulse inserted at z = zin
	long zin = 100;
	
	FILE *etime = fopen("etime.txt","wt");
	FILE *inversion_time = fopen("inversion.txt","wt");
	FILE	*eout = fopen("eout.txt","wt");
	FILE	*inversion = fopen("inversion.txt","wt");

	
	printf("init completed in %d secs\n", start-time(NULL));
	//-------------------------------
	//time cycle {s here
	//-------------------------------
	try
	{
		for (int nq=0; nq<time_points; nq++)
		{
			double t=nq*dt;
			//#pragma omp parallel num_threads(2)
			{

				// left PML magnetic field
				for (int m= 0;m<=Nabc-1;m++)
				{
					int j=-m+Nabc-1;
					double Bold = B[m];
					B[m]= d1[j]*B[m]+(dt*c*b1[j]/dz)*(E[m+1]-E[m]); //E[m]-E[m-1]
					H[m]= d1[j]*H[m]+b1[j]*(a1[j]*B[m]-c1[j]*Bold);
				}

				//main grid magnetic field

				for (int m=Nabc; m<= Nabc+space_points-1;m++)
				{
					B[m]=B[m] + Hup*( E[m+1]-E[m] );
					H[m] = B[m];
				}

				// right PML magnetic field
				for (int m= space_points+Nabc; m<= space_points + Nabc+ Nabc-2;m++)
				{
					double Bold = B[m];
					B[m]= d1[m-space_points-Nabc]*B[m]+(dt*c*b1[m-space_points-Nabc]/dz)*(E[m+1]-E[m]); //E[m]-E[m-1]
					H[m]= d1[m-space_points-Nabc]*H[m]+b1[m-space_points-Nabc]*(a1[m-space_points-Nabc]*B[m]-c1[m-space_points-Nabc]*Bold);
				}



				//left PML electric field
				for (int m= 1; m<= Nabc-1; m++)
				{

					double Dold=D[m];
					D[m]= d1[-m+Nabc-1]*D[m]+(dt*c*b1[-m+Nabc-1]/dz)*(H[m]-H[m-1]);
					E[m]= d1[-m+Nabc-1]*E[m]+b1[-m+Nabc-1]*(a1[-m+Nabc-1]*D[m]-c1[-m+Nabc-1]*Dold);
				}


				//main grid electric field
				for (int m= Nabc ;m<=media_start-1;m++)
				{
					D[m]=D[m] + Dup*( H[m] - H[m-1] );
					E[m]=D[m]*media[eps_r[m]];
				}




				//#pragma omp parallel for 
				for (int m=media_start; m<= media_end-1;m++)
				{
					// n2 = 2*dt*(-2*media_pump[pump[m]]+gama_nr*N0)/(2+gama_nr*dt);

					D[m]=D[m] + Dup*( H[m] - H[m-1] );  //look at Taflove p. 247
					//if nq<>1 then
					Pold2[m]=Pold[m];
					//if nq<>0 then
					Pold[m]=P[m];
					//here product of sqrt(eps0) added

					P[m] = p1*P[m]-p2*Pold2[m]+kp*E[m]*N[m];
					//save electric field for previous time step
					double  Eprevt=E[m];
					E[m]=(D[m]-ke*P[m])*media[eps_r[m]];
					N[m]=n1*N[m]+n2-kn*(E[m]+Eprevt)*(P[m]-Pold[m]);

				}

				for (int m=media_end; m<= space_points+Nabc-1; m++)
				{
					D[m]=D[m] + Dup*( H[m] - H[m-1] );
					E[m]=D[m]*media[eps_r[m]];
				}



				// right PML electric field
				for (int m= space_points+Nabc; m<= space_points + Nabc+ Nabc - 2; m++)
				{
					double Dold=D[m];
					D[m]= d1[m-space_points-Nabc]*D[m]+(dt*c*b1[m-space_points-Nabc]/dz)*(H[m]-H[m-1]);
					E[m]= d1[m-space_points-Nabc]*E[m]+b1[m-space_points-Nabc]*(a1[m-space_points-Nabc]*D[m]-c1[m-space_points-Nabc]*Dold);
				}
			}

			//output of coordinate field, polarisation and inversion distributions
			//of count=Ntime equally spaced time intervals

			if (nq>Ntime*time_i)
			{

				for (int ic=0;  ic <= space_points-1+Nabc+Nabc; ic++)
				{
					fprintf(eout,"%.15e ", E[ic]);
					fprintf(inversion,"%.15e ", N[ic]);
					//write(polarization,' ', P[ic]);
				}
				fprintf(eout,"\n");
				fprintf(inversion,"\n");


				//calculate percentage of job completed
				printf("%ld %% completed %d sec\r",(long)(100*nq/time_points), time(NULL)-start);
				time_i++;
			}

			fprintf(etime, "%.15e\n", E[z_out_time]);

			//fprintf(inversion_time, "%ld\n", N[(long)(media_start+0.5*(media_end-media_start))]);
		}
	}
	catch (exception e)
	{
		printf( "%s\n", e.what());
	}



	/*CloseFile(eout);
	CloseFile(inversion);
	CloseFile(polarization);
	*/
	fclose(etime);   
	fclose(inversion_time);
	fclose(eout);
	fclose(inversion);

	printf("\ncompleted %d sec\n", time(NULL)-start);
	return 0;
}

