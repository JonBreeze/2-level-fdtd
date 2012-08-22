// twolevel_fdtd.cpp: определяет точку входа для консольного приложения.
//

#include <stdio.h>
#include <memory.h>

#define M_PI       3.14159265358979323846
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <exception>
#include "twolevel_fdtd.h"


#pragma warning(disable:4996)



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
	
	//store inverted epsilon data to speedup computations
	double media[3]= {1.0,1/reflector_index,1/(0.2*reflector_index)};
	double media_pump[3]= {0,delta,delta1};
	two_level_params active_medium_params[2];

	FdtdFields *fdtd_fields = (FdtdFields*)malloc((space_points+2*Nabc)*sizeof(FdtdFields));
	memset(fdtd_fields, 0, (space_points+2*Nabc)*sizeof(FdtdFields));
	PMLParams *pml = (PMLParams*)malloc(Nabc*sizeof(PMLParams));
	SourceFields * source = (SourceFields*)malloc(sizeof(SourceFields)*time_points);
	MediaParamsIndex * paramsIndex = (MediaParamsIndex *)malloc(sizeof(MediaParamsIndex)*(space_points+2*Nabc));
	memset(paramsIndex, 0, sizeof(MediaParamsIndex)*(space_points+2*Nabc));
	

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
					paramsIndex[i].eps_r=1;
				//second material
				for (int i= ref_start+(wid+wid2)*j+wid;i<=ref_start+(wid+wid2)*j+wid+wid2;i++)
					paramsIndex[i].eps_r=2;
				//right mirror-------
				for (int i= ref_start+2*wid*(ref_layers-1)+wid+ref_length+2*wid*j; i<= ref_start+2*wid*(ref_layers-1)+wid+ref_length+2*wid*j+wid; i++)
					paramsIndex[i].eps_r=1;
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
					paramsIndex[i].eps_r=1;
				//second material
				for (int i= ref_start+(wid+wid2)*j+wid;i<=ref_start+(wid+wid2)*j+wid+wid2;i++)
					paramsIndex[i].eps_r=2;
			}

		}
		break;
	default:
		{
			long media_length= (long)(media_end - media_start)*media_number;
			//bulk media with reflecting edges
			for (int i=0; i<=space_points-1+Nabc+Nabc; i++)
			{
				paramsIndex[i].eps_r= 0;
			}

			for (int j= media_start-media_length;j<=media_end+media_length-1;j++)
			{
				paramsIndex[j].eps_r= 1;
			}
		}


	}



	FILE *eps = fopen("eps.txt","wt");
	for(int i=0; i<space_points+2*Nabc; i++)
	{
		fprintf(eps,"%ld\n", paramsIndex[i].eps_r);
	}
	fclose(eps);

	//media with absorbing cell
	long delta_number=(long)((2/3)*(media_end-media_start));

	for (int j=media_start; j<=media_start+delta_number-1;j++)
	{
		paramsIndex[j].pump=1;
	}

	for (int j=media_start+delta_number; j<=media_end-1-delta_number;j++)
	{
		paramsIndex[j].pump=2;
	}

	for (int j=media_end-delta_number;j<=media_end-1;j++)
	{
		paramsIndex[j].pump=1;
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
		fdtd_fields[i].N=N11;
	}

	srand(  time(NULL));
	for (int i = media_start; i<= media_end-1; i++)
	{
		fdtd_fields[i].P=sin(w0*i*dt/*+2*M_PI*((double)rand()/(double)RAND_MAX)*/);
	}

	FILE *pstart = fopen("po.txt","wt");
	for(int i=0; i< space_points+2*Nabc+1; i++)
	{
		fprintf(pstart,"%.15e\n",fdtd_fields[i].P);
	}
	fclose(pstart);
	

	double Hdelay  = 0;//-dt/2;//-dz/2/c;


	//hyperbolic secant pulse
	
	for(int  nq=0; nq< pulse_points; nq++)
	{
		double t=nq*dt;
		//G=0.5*(t-Tp)/(Tp);
		source[nq].Ein = E0*sin( w0*t )/cosh(10*(t-0.75*Tp)/Tp);
		source[nq].Hin = E0*sin( w0*(t+Hdelay) )/cosh(10*(t+Hdelay-0.75*Tp)/Tp);

	}
	
	

	//pml parameters calculation
	double abc_pow = 3;
	double kmax=2;
	double sigma_max = -(abc_pow+1)*log(0.0001)*c/(2*/*eta0*/Labc); //2e16;
	for (int i=0; i<= Nabc-1; i++)
	{
		pml[i].sigma =  sigma_max*pow(i/Nabc, abc_pow)/**mu0/eps0*/;
		pml[i].k=1+(kmax-1)*pow(i/Nabc, abc_pow);
		//sigma[i]=0;
		//k[i]=1;
		pml[i].a1= pml[i].k+0.5*pml[i].sigma*dt;
		pml[i].b1= 1/pml[i].a1;//1/(k[i]+0.5*sigma[i]*dt);
		pml[i].c1= pml[i].k-0.5*pml[i].sigma*dt;
		pml[i].d1= pml[i].c1/pml[i].a1;//(k[i]-0.5*sigma[i]*dt)/(k[i]+0.5*sigma[i]*dt);

		//writeln(polarization,' ', a1[i],'  ',b1[i],' ',c1[i],'  ',d1[i],' ',sigma[i],' ',k[i]);
	}


	//pumping rate cmdline input
	/*
	if ParamCount>0 then
	delta = StrToFloat(ParamStr(1));
	*/
	//precalculated media parameters
/*
	const double     p1=(2-omega0sqr*dt*dt)/(1+gama*dt);
	const double    p2=(1-gama*dt)/(1+gama*dt);
	const double p3=omega0sqr*dt*dt/(h11*omega*(1+gama*dt));

	const double n1=(2-gama_nr*dt)/(2+gama_nr*dt);
	const double n2 = 2*dt*(-2*delta+gama_nr*N0)/(2+gama_nr*dt);
	const double n3=4*omega/(h11*omega0sqr*(2+gama_nr*dt));


	const double kp= p3*dab;
	const double ke= 2*na*dab;
	const double kn= n3*dab;
*/	

	memset(active_medium_params, 0, sizeof(active_medium_params));
	
	active_medium_params[1].p1=(2-omega0sqr*dt*dt)/(1+gama*dt);
	active_medium_params[1].p2=(1-gama*dt)/(1+gama*dt);
	active_medium_params[1].p3=omega0sqr*dt*dt/(h11*omega*(1+gama*dt));

	active_medium_params[1].n1=(2-gama_nr*dt)/(2+gama_nr*dt);
	active_medium_params[1].n2 = 2*dt*(-2*delta+gama_nr*N0)/(2+gama_nr*dt);
	active_medium_params[1].n3 = 4*omega/(h11*omega0sqr*(2+gama_nr*dt));


	active_medium_params[1].kp= active_medium_params[0].p3*dab;
	active_medium_params[1].ke= 2*na*dab;
	active_medium_params[1].kn= active_medium_params[0].n3*dab;

//init indices
	for (int m=media_start; m<= media_end-1;m++)
	{
		paramsIndex[m].two_level_par_index = 1;
	}

	
	

	//excitation pulse inserted at z = zin
	long zin = 150;
	
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
					double Bold = fdtd_fields[m].B;
					fdtd_fields[m].B= pml[j].d1*fdtd_fields[m].B + (dt*c*pml[j].b1/dz)*(fdtd_fields[m+1].E-fdtd_fields[m].E); //E[m]-E[m-1]
					fdtd_fields[m].H= pml[j].d1*fdtd_fields[m].H+pml[j].b1*(pml[j].a1*fdtd_fields[m].B-pml[j].c1*Bold);
				}

				//main grid magnetic field

				for (int m=Nabc; m<= Nabc+space_points-1;m++)
				{
					fdtd_fields[m].B = fdtd_fields[m].B + Hup*( fdtd_fields[m+1].E-fdtd_fields[m].E );
					fdtd_fields[m].H = fdtd_fields[m].B;
				}

				// right PML magnetic field
				for (int m= space_points+Nabc; m<= space_points + Nabc+ Nabc-2;m++)
				{
					double Bold = fdtd_fields[m].B;
					int j = m-space_points-Nabc;
					fdtd_fields[m].B =
						pml[j].d1*fdtd_fields[m].B+(dt*c*pml[j].b1/dz)*(fdtd_fields[m+1].E-fdtd_fields[m].E); //E[m]-E[m-1]
					fdtd_fields[m].H =
						pml[j].d1*fdtd_fields[m].H+pml[j].b1*(pml[j].a1*fdtd_fields[m].B-pml[j].c1*Bold);
				}


				//left PML electric field
				for (int m= 1; m<= Nabc-1; m++)
				{

					double Dold=fdtd_fields[m].D;
					int j = -m+Nabc-1;
					fdtd_fields[m].D= 
						pml[j].d1*fdtd_fields[m].D+(dt*c*pml[j].b1/dz)*(fdtd_fields[m].H-fdtd_fields[m-1].H);
					fdtd_fields[m].E=
						pml[j].d1*fdtd_fields[m].E+pml[j].b1*(pml[j].a1*fdtd_fields[m].D-pml[j].c1*Dold);
				}


				//main grid electric field
				
			 
				for (int m= Nabc; m<= space_points+Nabc-1;m++)
				{
					
					fdtd_fields[m].D=fdtd_fields[m].D + Dup*( fdtd_fields[m].H - fdtd_fields[m-1].H );  //look at Taflove p. 247
					//if nq<>1 then
					fdtd_fields[m].Pold2=fdtd_fields[m].Pold;
					//if nq<>0 then
					fdtd_fields[m].Pold=fdtd_fields[m].P;
					//here product of sqrt(eps0) added
					int iam =paramsIndex[m].two_level_par_index; 
					fdtd_fields[m].P = active_medium_params[iam].p1*fdtd_fields[m].P-
						active_medium_params[iam].p2*fdtd_fields[m].Pold2+
						active_medium_params[iam].kp*fdtd_fields[m].E*fdtd_fields[m].N;
					//save electric field for previous time step
					double  Eprevt=fdtd_fields[m].E;
					fdtd_fields[m].E=(fdtd_fields[m].D-active_medium_params[iam].ke*fdtd_fields[m].P)*media[paramsIndex[m].eps_r];
					fdtd_fields[m].N=active_medium_params[iam].n1*fdtd_fields[m].N+
						active_medium_params[iam].n2-
						active_medium_params[iam].kn*(fdtd_fields[m].E+Eprevt)*(fdtd_fields[m].P-fdtd_fields[m].Pold);

				}

				//insert hard source
				//fdtd_fields[zin].E =source[nq].Ein;



				// right PML electric field
				for (int m= space_points+Nabc; m<= space_points + Nabc+ Nabc - 2; m++)
				{
					double Dold=fdtd_fields[m].D;
					int j = m-space_points-Nabc;
					fdtd_fields[m].D= 
						pml[j].d1*fdtd_fields[m].D+(dt*c*pml[j].b1/dz)*(fdtd_fields[m].H-fdtd_fields[m].H);
					fdtd_fields[m].E= 
						pml[j].d1*fdtd_fields[m].E+pml[j].b1*(pml[j].a1*fdtd_fields[m].D-pml[j].c1*Dold);
				}
			}

			//output of coordinate field, polarisation and inversion distributions
			//of count=Ntime equally spaced time intervals

			if (nq>Ntime*time_i)
			{

				for (int ic=0;  ic <= space_points-1+Nabc+Nabc; ic++)
				{
					fprintf(eout,"%.15e ", fdtd_fields[ic].E);
					fprintf(inversion,"%.15e ", fdtd_fields[ic].N);
					//write(polarization,' ', P[ic]);
				}
				fprintf(eout,"\n");
				fprintf(inversion,"\n");


				//calculate percentage of job completed
				printf("%ld %% completed %d sec\r",(long)(100*nq/time_points), time(NULL)-start);
				time_i++;
			}

			fprintf(etime, "%.15e\n", fdtd_fields[z_out_time].E);

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
	free(fdtd_fields);
	fclose(etime);   
	fclose(inversion_time);
	fclose(eout);
	fclose(inversion);

	printf("\ncompleted %d sec\n", time(NULL)-start);
	return 0;
}

