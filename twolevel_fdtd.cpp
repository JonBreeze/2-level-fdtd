#include <stdio.h>
#include <memory.h>

#define M_PI       3.14159265358979323846
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <exception>
//#include <thread>
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
double Hup=c*dt/dz;
double Eup= c*dt/dz;
double Dup=c*dt/dz;


double omega0sqr = gama*gama+omega*omega;
int Nabc = 50;      //Number of PML cells
double Labc = Nabc*dz;
double mur_factor = 0;

//allocate arrays
int res = 0;

//store inverted epsilon data to speedup computations
double media[3]= {1.0,1/reflector_index,1/(0.2*reflector_index)};
double media_pump[3]= {0,delta,delta1};
two_level_params active_medium_params[2];

FdtdFields *fdtd_fields = 0;
MediaParamsIndex * paramsIndex = 0;
int source_xsize = space_points+2;
SourceFields * source = 0;
int pulse_points = time_points;//2*(int)(Tp/dt);

//excitation pulse inserted at z = zin
long zin = 200;
long zin_right = 700;

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

void mur_1d_source_grid_timestep(int xstart, int xend)
{
    //source time-stepping
    bool is_left_bound = (xstart == 0)?true:false;
    bool is_right_bound = (xend == space_points-1)?true:false;
    int e_xstart = is_left_bound ? 1 : xstart;
    int h_xend = is_right_bound ? xend-1 : xend;

    double Esrcold = source[1].Ein;
    for (int m= e_xstart; m < xend; m++)
    {
        source[m].Ein += Eup*(source[m].Hin - source[m-1].Hin);
    }
    //Mur boundary
    if (is_left_bound)
    {
        source[0].Ein=
                Esrcold  +	mur_factor*(source[1].Ein - source[0].Ein);
    }

    double Hsrcold = source[source_xsize-2].Hin;
    for (int m=0; m <= h_xend; m++)
    {
        source[m].Hin += Hup*(source[m+1].Ein-source[m].Ein);
    }
    //Mur boundary
    if (is_right_bound)
    {
        source[source_xsize-1].Hin=
                Hsrcold  +	mur_factor*(source[source_xsize-2].Hin - source[source_xsize-1].Hin);
    }
}

void mur_1d_main_grid_timestep(int xstart, int xend)
{
    //main grid electric field
    bool is_left_bound = (xstart == 0)?true:false;
    bool is_right_bound = (xend == space_points-1)?true:false;
    int e_xstart = is_left_bound ? 1 : xstart;
    int h_xend = is_right_bound ? xend-1 : xend;
    double Eold = fdtd_fields[1].E;
    for (int m = e_xstart; m <= xend; m++)
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

        fdtd_fields[m].E += Eup*( fdtd_fields[m].H - fdtd_fields[m-1].H );

    }

    //Mur boundary
    if (is_left_bound)
    {
        fdtd_fields[0].E=
                Eold  +	mur_factor*(fdtd_fields[1].E - fdtd_fields[0].E);
    }

    //total\scattered correction
    /*fdtd_fields[zin].E-=Eup*source[1].Hin;
    fdtd_fields[zin_right].E+=Eup*source[zin_right - zin + 2].Hin;
    */


    //main grid magnetic field
    double Hold = fdtd_fields[space_points-2].H;
    for (int m=0; m < h_xend; m++)
    {
        fdtd_fields[m].H = fdtd_fields[m].H + Hup*( fdtd_fields[m+1].E-fdtd_fields[m].E );
    }

    //apply Mur boundary conditions
    if (is_right_bound)
    {
        fdtd_fields[space_points-1].H=
                Hold  +	mur_factor*(fdtd_fields[space_points-2].H - fdtd_fields[space_points-1].H);
    }

    //total/scattered field correction
    /*
    fdtd_fields[zin-1].H -= Hup*source[2].Ein;
    fdtd_fields[zin_right].H += Eup*source[zin_right-zin+2].Ein;
    */

}


int main(int argc, char* argv[])
{
    clock_t start = clock();
    if (!read_parameters())
    {
        printf("can't load config.txt!");
        return 0;
    }
    //memory allocation
    fdtd_fields = (FdtdFields*)malloc((space_points)*sizeof(FdtdFields));
    memset(fdtd_fields, 0, (space_points)*sizeof(FdtdFields));
    paramsIndex = (MediaParamsIndex *)malloc(sizeof(MediaParamsIndex)*(space_points));
    memset(paramsIndex, 0, sizeof(MediaParamsIndex)*(space_points));

    //convert to heaviside-lorenz units
    E0= E0*E_from_SI_to_HL;
    dab=dab*dab_from_SI_to_HL;
    na=na*na_from_SI_to_HL;

    //initial conditions for population inversion and polarization rate
    srand( time(NULL));
    for (int i = media_start; i<= media_end-1; i++)
    {
        fdtd_fields[i].N=N11;
        fdtd_fields[i].P=0;//sin(w0*i*dt+2*M_PI*((double)rand()/(double)RAND_MAX));
    }

    FILE *pstart = fopen("po.txt","wt");
    for(int i=0; i< space_points; i++)
    {
        fprintf(pstart,"%.15e\n",fdtd_fields[i].P);
    }
    fclose(pstart);


    source = (SourceFields*)malloc(sizeof(SourceFields)*source_xsize);
    memset(source, 0, sizeof(SourceFields)*source_xsize);

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
        paramsIndex[m].two_level_par_index = 1;
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
        fprintf(steps, "time step=%.15e\nspace step=%.15e\nCourant factor=%.15e\n", dt, dz, c*dt/dz);
        fclose(steps);
    }


    FILE *etime = fopen("e_time.txt","wt");
    FILE *inversion_time = fopen("inversion_time.txt","wt");


    printf("init completed in %.0f msecs\n", ((double)clock()-(double)start)*1000/CLOCKS_PER_SEC);


    //-------------------------------
    //time cycle starts here
    //-------------------------------
    mur_factor = (c*dt-dz)/(c*dt+dz);
    try
    {
        for (int nq=0; nq<time_points; nq++)
        {
            //TODO: add source time step function call
            mur_1d_main_grid_timestep(0,space_points-1);
            //insert hard source into the main_grid
            //TODO: hard source should be on the source's grid and
            //total/scattered field correction on main grid must be used
            if (nq<pulse_points)
            {
                fdtd_fields[zin].E =source[nq].Ein;
            }
            //}

            //output of coordinate field, polarisation and inversion distributions
            //of count=Ntime equally spaced time intervals

            if (nq>Ntime*time_i)
            {
                char fname[256]="";
                sprintf(fname,"eout%d.txt", time_i);
                FILE	*eout = fopen(fname,"wt");
                sprintf(fname,"inversion%d.txt", time_i);
                FILE	*inversion = fopen(fname,"wt");
                sprintf(fname,"esrc%d.txt", time_i);
                FILE * fsrc = fopen(fname, "wt");
                sprintf(fname,"polarization%d.txt", time_i);
                FILE * fpol = fopen(fname, "wt");


                for (int ic=0;  ic < space_points; ic++)
                {
                    fprintf(eout,"%.15e \n", fdtd_fields[ic].E);
                    fprintf(inversion,"%.15e \n", fdtd_fields[ic].N);
                    fprintf(fsrc,"%.15e \n", source[ic].Ein);
                    fprintf(fpol,"%.15e \n", fdtd_fields[ic].P);

                }
                fclose(inversion);
                fclose(eout);
                fclose(fsrc);
                fclose(fpol);

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
    free(source);
    free(paramsIndex);
    fclose(etime);
    fclose(inversion_time);

    printf("all completed in %.0f msecs\n", ((double)clock()-(double)start)*1000/CLOCKS_PER_SEC);
    return 0;
}

