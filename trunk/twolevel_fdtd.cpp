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
	clock_t start = clock();
	if (!read_parameters())
	{
		printf("can't load config.txt!");
		return 0;
	}

	int time_i = 0;
	int Ntime = time_points / time_samps ;

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
	double E_from_SI_to_HL=1e-4/3;
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


	//allocate arrays
	int res = 0;

	//store inverted epsilon data to speedup computations
	double media[3]= {1.0,1/reflector_index,1/(0.2*reflector_index)};
	double media_pump[3]= {0,delta,delta1};
	two_level_params active_medium_params[2];

	FdtdFields *fdtd_fields = (FdtdFields*)malloc((space_points)*sizeof(FdtdFields));
	memset(fdtd_fields, 0, (space_points)*sizeof(FdtdFields));
	PMLParams *pml = (PMLParams*)malloc(Nabc*sizeof(PMLParams));	
	MediaParamsIndex * paramsIndex = (MediaParamsIndex *)malloc(sizeof(MediaParamsIndex)*(space_points));
	memset(paramsIndex, 0, sizeof(MediaParamsIndex)*(space_points));


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
			for (int i=0; i<space_points; i++)
			{
				paramsIndex[i].eps_r= 0;
			}

			for (int j= media_start-media_length;j<=media_end+media_length-1;j++)
			{
				paramsIndex[j].eps_r= 1;
			}
		}


	}





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


	//initial conditions for population inversion and polarization rate
	for (int i = media_start; i<=media_end-1;i++)
	{
		fdtd_fields[i].N=N11;
	}

	srand( time(NULL));
	for (int i = media_start; i<= media_end-1; i++)
	{
		fdtd_fields[i].P=sin(w0*i*dt+2*M_PI*((double)rand()/(double)RAND_MAX));
	}

	FILE *pstart = fopen("po.txt","wt");
	for(int i=0; i< space_points; i++)
	{
		fprintf(pstart,"%.15e\n",fdtd_fields[i].P);
	}
	fclose(pstart);


	

	int pulse_points = time_points;//2*(int)(Tp/dt);
	SourceFields * source = (SourceFields*)malloc(sizeof(SourceFields)*pulse_points);
	
	double Hdelay  = dt/2;
	FILE * fsrc = fopen("esrc.txt", "wt");
	for(int  nq=0; nq< pulse_points; nq++)
	{
		double t=nq*dt;
		//G=0.5*(t-Tp)/(Tp);
		source[nq].Ein = E0*sin( w0*t );//cosh(10*(t-0.75*Tp)/Tp);
		source[nq].Hin = E0*sin( w0*(t+Hdelay) );//cosh(10*(t+Hdelay-0.75*Tp)/Tp);
		fprintf(fsrc, "%.15e\n", source[nq].Ein, source[nq].Hin);

	}
	fclose(fsrc);




	//pml parameters calculation
	double abc_pow = 3;
	double kmax=2;
	double sigma_max = -2e16;
	//-(abc_pow+1)*log(0.0001)*c/(2*/*eta0*/Labc); //2e16;
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




	memset(active_medium_params, 0, sizeof(active_medium_params));

	active_medium_params[1].p1=(2-omega0sqr*dt*dt)/(1+gama*dt);
	active_medium_params[1].p2=(1-gama*dt)/(1+gama*dt);
	active_medium_params[1].p3=omega0sqr*dt*dt/(h11*omega*(1+gama*dt));

	active_medium_params[1].n1=(2-gama_nr*dt)/(2+gama_nr*dt);
	active_medium_params[1].n2 = 2*dt*(-2*delta+gama_nr*N0)/(2+gama_nr*dt);
	active_medium_params[1].n3 = 4*omega/(h11*omega0sqr*(2+gama_nr*dt));


	active_medium_params[1].kp= active_medium_params[1].p3*dab;
	active_medium_params[1].ke= 2*na*dab;
	active_medium_params[1].kn= active_medium_params[1].n3*dab;

	//init indices
	for (int m=media_start; m<= media_end-1;m++)
	{
		paramsIndex[m].two_level_par_index = 0;
	}

	//read epsilon two-level params and pump rate indexes from a file
	FILE *eps = NULL;	
	eps =fopen("media_in.txt","rt");	
	if (eps!=NULL)
	{
		for(int i=0; i<space_points; i++)
		{
			fscanf(eps,"%d%d%d\n", &paramsIndex[i].eps_r, 
				&paramsIndex[i].two_level_par_index,
				&paramsIndex[i].pump);
		}
		fclose(eps);
	}

	FILE *fmedia = NULL;	
	fmedia =fopen("media_out.txt","wt");	
	for(int i=0; i<space_points; i++)
	{
		fprintf(fmedia,"%d %d %d\n", paramsIndex[i].eps_r, 
			paramsIndex[i].two_level_par_index,
			paramsIndex[i].pump);
	}
	fclose(fmedia);



	//write dt and dx  
	FILE *steps = NULL;
	steps =fopen("steps.txt","wt");
	if (steps!=NULL)
	{
		fprintf(steps, "%.15e\n%.15e\n", dt, dz);
		fclose(steps);
	}	


	//excitation pulse inserted at z = zin
	long zin = 200;
	long zin_right = 800;

	FILE *etime = fopen("e_time.txt","wt");
	FILE *inversion_time = fopen("inversion_time.txt","wt");		


	printf("init completed in %.0f msecs\n", ((double)clock()-(double)start)*1000/CLOCKS_PER_SEC);
/* TODO: implement auxiliary FDTD grid for source field calculation
	double *Esrc, *Hsrc;
	int srcLength = zin_right-zin+1;
	Esrc = (double*)malloc(srcLength*sizeof(double));
*/
	//-------------------------------
	//time cycle {s here
	//-------------------------------
	double mur_factor = (c*dt-dz)/(c*dt+dz);
	try
	{
		for (int nq=0; nq<time_points; nq++)
		{
			double t=nq*dt;			

			//main grid magnetic field
			double Hold = fdtd_fields[space_points-2].H;
			for (int m=0; m < space_points-1;m++)
			{
				fdtd_fields[m].H = fdtd_fields[m].H + Hup*( fdtd_fields[m+1].E-fdtd_fields[m].E );

			}
			//total/scattered field correction
			fdtd_fields[zin].H-=Hup*source[nq].Ein;
			//apply Mur boundary conditions
			fdtd_fields[space_points-1].H=
				Hold  +	mur_factor*(fdtd_fields[space_points-2].H - fdtd_fields[space_points-1].H);

			//main grid electric field

			double Eold = fdtd_fields[1].E;
			for (int m= 1; m<= space_points-1;m++)
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

			//total\scattered correction
			fdtd_fields[zin].E-=Eup*source[nq].Hin;

			//Mur boundary
			fdtd_fields[0].E=
				Eold  +	mur_factor*(fdtd_fields[1].E - fdtd_fields[0].E);

			//insert hard source
/*
			if (nq<pulse_points)
			{
				fdtd_fields[zin].E =source[nq].Ein;
			}
*/



			//}

			//output of coordinate field, polarisation and inversion distributions
			//of count=Ntime equally spaced time intervals

			if (nq>Ntime*time_i)
			{				
				char fname[256]="";;
				sprintf(fname,"eout%d.txt", time_i);
				FILE	*eout = fopen(fname,"wt");
				sprintf(fname,"inversion%d.txt", time_i);
				FILE	*inversion = fopen(fname,"wt");
				for (int ic=0;  ic < space_points; ic++)
				{					
					fprintf(eout,"%.15e \n", fdtd_fields[ic].E);
					fprintf(inversion,"%.15e \n", fdtd_fields[ic].N);																	
					
				}			
				fclose(inversion);
				fclose(eout);


				//calculate percentage of job completed
				printf("%ld %% completed %.0f msec\r",(long)(100*nq/time_points), ((double)clock()-(double)start)*1000/CLOCKS_PER_SEC);
				time_i++;
			}

			fprintf(etime, "%.15e\n", fdtd_fields[z_out_time].E);
			fprintf(inversion_time, "%.15e\n", fdtd_fields[(media_start+media_end)/2].N);

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
	free(pml);
	free(source);
	free(paramsIndex);
	fclose(etime);   
	fclose(inversion_time);

	printf("all completed in %.0f msecs\n", ((double)clock()-(double)start)*1000/CLOCKS_PER_SEC);
	return 0;
}

