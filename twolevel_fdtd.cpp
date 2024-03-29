#include <stdio.h>
#include <memory.h>

#define M_PI       3.14159265358979323846
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <exception>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include "twolevel_fdtd.h"


#pragma warning(disable:4996)



using namespace std;


//vars to be read from config file
int time_points, space_points, media_start, media_end, media_number;
int time_samps;
int  ref_start, ref_length, ref_layers,ref_config;
double omega,gama,na,delta,delta1,gama_nr,dab,N0,N11,Tp,w0,E0,step_lambda_rate;
double reflector_index;
int z_out_time;
int time_i = 0;





//this values are calculated on the base of the above vars
double dz, dt;
double Hup, Eup, Dup;

double omega0sqr = 0.0;
double mur_factor = 0.0;


//store inverted epsilon data to speedup computations
double media[3]= {1.0,1.0,1.0};
double media_pump[3]= {0.0, 0.0, 0.0};
two_level_params active_medium_params[2];

FdtdFields *fdtd_fields = 0;
MediaParamsIndex * paramsIndex = 0;

SourceFields * source = 0;
int source_xsize = 0;

//excitation pulse inserted at z = zin and removed from zin_right
long zin = 50;
long zin_right = 700;


clock_t start;



mutex mx;
condition_variable cv;
condition_variable time_step_cv;
int total_thread_count = 1;
int working_thread_count(total_thread_count);
bool file_written = true;


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

void mur_1d_source_grid_timestep(int xstart, int xend, double atime)
{

    //source time-stepping
    bool is_left_bound = (xstart == 0)?true:false;
    bool is_right_bound = (xend == source_xsize-1)?true:false;
    int e_xstart = is_left_bound ? 1 : xstart;
    int h_xend = is_right_bound ? xend-1 : xend;


    double Esrcold_left = source[1].Ein;

    for (int m= e_xstart; m <= xend; m++)
    {
        source[m].Ein += Eup*(source[m].Hin - source[m-1].Hin);
    }
    //Mur boundary for E field
    if (is_left_bound)
    {
        source[0].Ein=
                Esrcold_left  +	mur_factor*(source[1].Ein - source[0].Ein);
    }

    source[zin].Ein = E0*exp(-(atime/Tp-1)*(atime/Tp-1))*sin(w0*atime);

    double Hsrcold_right = source[source_xsize-2].Hin;
    for (int m=0; m <= h_xend; m++)
    {
        source[m].Hin += Hup*(source[m+1].Ein-source[m].Ein);
    }
    //Mur boundary for H field
    if (is_right_bound)
    {
        source[source_xsize-1].Hin =
                Hsrcold_right - mur_factor*(source[source_xsize-1].Hin - source[source_xsize-2].Hin);
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
    for (int m=0; m <= h_xend; m++)
    {
        fdtd_fields[m].H = fdtd_fields[m].H + Hup*( fdtd_fields[m+1].E-fdtd_fields[m].E );
    }

    //apply Mur boundary conditions
    if (is_right_bound)
    {
        fdtd_fields[space_points-1].H=
                Hold  -	mur_factor*(fdtd_fields[space_points-2].H - fdtd_fields[space_points-1].H);
    }

    //total/scattered field correction
    /*
    fdtd_fields[zin-1].H -= Hup*source[2].Ein;
    fdtd_fields[zin_right].H += Eup*source[zin_right-zin+2].Ein;
    */

}


void thread_func(int xstart, int xend)
{
    unique_lock<std::mutex> ulock(mx);
    printf("fdtd thread: after ulock\n");

    for (int nq=0; nq<time_points; nq++)
    {
        printf("before lock: fdtd thread nq=%d\n", nq);
        while (!file_written)
        {
            cv.wait(ulock);
        }

        printf("fdtd thread nq=%d\n", nq);
        double time = nq*dt;
        // mur_1d_main_grid_timestep(0,space_points-1);

        //TODO: hard source should be on the source's grid and
        //total/scattered field correction on main grid must be used

        mur_1d_source_grid_timestep( xstart, xend, time);
        working_thread_count--;
        if (0==working_thread_count)
        {
            file_written = false;
            cv.notify_all();
            time_step_cv.notify_all();
        }
        printf("fdtd thread : working_thread_count=%d\n",static_cast<int>(working_thread_count));
        while (working_thread_count > 0)
        {
            time_step_cv.wait(ulock);
        }
    }
}

void file_thread_func()
{    
    unique_lock<std::mutex> ulock(mx);
    printf("file thread: after ulock\n");
    FILE *etime = fopen("e_time.txt","wt");
    FILE *inversion_time = fopen("inversion_time.txt","wt");
    int Ntime = time_points / time_samps;
    for (int nq=0; nq<time_points; nq++)
    {
        //output of coordinate field, polarisation and inversion distributions
        //of count=Ntime equally spaced time intervals
        printf("before lock: file thread nq=%d\n", nq);
        while (working_thread_count > 0)
        {
            cv.wait(ulock);
        }
        printf("file thread nq=%d\n", nq);

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
        //fprintf(inversion_time, "%.15e\n", fdtd_fields[(media_start+media_end)/2].N);
        file_written = true;
        working_thread_count = total_thread_count;
        cv.notify_all();
    }
    fclose(etime);
    fclose(inversion_time);
}


int main(int argc, char* argv[])
{


    start = clock();
    if (!read_parameters())
    {
        printf("can't load config.txt!");
        return 0;
    }


    source_xsize = space_points+2;

    //set medium non-resonant epsilon
    media[0]= 1.0;
    media[1] = 1/reflector_index;
    media[2] = 1/(0.2*reflector_index);
    //set pumping rate
    media_pump[0]= 0;
    media_pump[1] = delta;
    media_pump[2] = delta;
    //grid steps
    dz=step_lambda_rate*(2*M_PI*c/w0);
    dt=0.5*dz/c;
    //update coeffitients
    Hup=c*dt/dz;
    Eup= c*dt/dz;
    Dup=c*dt/dz;
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


    omega0sqr = gama*gama+omega*omega;
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

    //read epsilon two-level params and pump rate indices from a file
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




    printf("init completed in %.0f msecs\n", ((double)clock()-(double)start)*1000/CLOCKS_PER_SEC);



    //-------------------------------
    //time cycle starts here
    //-------------------------------
    mur_factor = (c*dt-dz)/(c*dt+dz);
    try
    {
        total_thread_count = 2;
        working_thread_count = total_thread_count;
        printf("running %d threads\n", total_thread_count);
        std::thread th( thread_func, 0, source_xsize/2);
        std::thread th1( thread_func, source_xsize/2+1, source_xsize-1);
        std::thread file_th(file_thread_func);



        th.join();
        th1.join();
        file_th.join();
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

    printf("all completed in %.0f msecs\n", ((double)clock()-(double)start)*1000/CLOCKS_PER_SEC);
    return 0;
}

