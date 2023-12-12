////20231024
//2021107 plate is not moving.   define is modified. not complete.
//   abs(vvv) < 1E-5	-> stopping
// no backward motion
// functionalized
// cleaned
// further modified
// measure the maginitude
//initial state is given by random numbers
//structured
//double xxxnew[NB+2],vvvnew[NB+2] in struct Matter are removed
//record the configuration just before the start of the slip
//skip the stopping interval
//kinetic frictional force is modified.  singma is introduced
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//
#define RAND_MAX 0x7fffffff
#define STEP_MAX 0x7fffffff
#define PI 3.141592653589793
//#define FILE (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/')	+ 1 : __FILE__)
//
#define DEBUG_INITIAL
#define DEBUG_ON
#define FLAGPRTOUT  1
#define FLAGPRTMAG  1
#define FLAGFPRXXX  0
#define FLAGF01OUT  1
#define FLAGF02MAG  1
#define FLAGF03XXX  1
//
//
#define FLAGINITRAND  0 // initial displacement is given by random numbers,	otherwise initial displacement is 0.0
#define NB  10	// # of Block
#define EVENTNO	 10   // max # of events
#define INITEVENTNO  10   // # of initial events
#define MAXTIME	200000.0
#define KPP 1.0
#define KXX 10.0
#define SIGMA 0.01
#define ALPHA 1.0
#define BETA 2.0
#define VPL 1.0E-2
#define FRICMAX	1.0
#define ISEED  728761879   // seed of random number
#define WIRAND 2*PI   // width of initial	displacement given	by random number
//#define DATAFILEDIR "/Users/hiro/Dropbox/Prog/BK210926/BK2109281D/Data/"っc
#define DATAFILEDIR "/Users/uzawataiki/data/data"
#define OUTPUTPARAMATERS "# %s	\n# %s \n# %s\n# NB= %d    VPL=%e	KPP=%e  KXX=%e \n# FRICMAX=%e  SIGMA=%e  ALPHA=%e  deltat=%e ISEED=%d \n", arg1,arg2,arg3,NB,VPL,KPP,KXX,FRICMAX,SIGMA,ALPHA,deltat,ISEED//
#define OUTPUTTITLE "#    Ev.N.    St_St     Ed_St  St_Time       End_Time	Plate Pos.    Seis.Mom.   Log10(Seis.Mom.)  \n"
#define OUTPUTMAG "%8d  %8d  %8d  %e  %e %e  %e	%e \n",no_of_earthquake,istep00,istep,istep00*deltat,istep*deltat,xppt,sumavxxx,magnitude
#define OUTPUTSTEP "istep = %6d  time = %e flagslip=%d xppt = %+e  avxxx	= %+e  avvvv = %+e maxvvv = %+e\n",istep,istep*deltat,flagslip,xppt,avxxx,avvvv,maxvvv
#define ACCURACY00 0.0E-3//  MUST	be greater than 1.0E-7.  Accuracy of velocity 0.
//#define ACCURACYTH 0.0//  MUST be greater than 1.0E-7.  Accuracy of static frictional force.
//#define FLAG00(xxx)   ( (int)( (xxx-ACCURACY00)*1.0E-7	+ 1.0 ) ) // FLAG of xxx >0 (->1) or <=0 (->0)
//#define FLAGTH(xxx,xth)   ( (int)( (fabs(xxx)-xth-ACCURACYTH)*1.0E-7 + 1.0 )	) // FLAG of abs(xxx) - xth  > 0 (->1) or	<=0 (->0)
//#define FLAGSLIP(forcel,maxstfric,vvv)  ( 1- ( 1-FLAGTH(forcel,maxstfric) )*( 1-FLAG00(fabs(vvv)) ) ) // FLAG of SLIPPING
#define ELFORCE(xx0,xx1,xx2,xxp)  (	KXX* (xx1 + xx2 - 2.0*xx0) + KPP*(xxp-xx0) ) // Elastic Force
//#define KINFRIC(vvv)  ( - ( FLAG00(vvv)-FLAG00(-vvv) ) * (1-SIGMA) / ( 1.0+ALPHA*fabs(vvv)/(1-SIGMA) ) ) // Kinetic	Friction Force
//#define FFXX(vvv)		( 0.5*( vvv + fabs(vvv) ) ) // Change of X
//#define ReLU(xxx)		( 0.5*( xxx + fabs(xxx) ) ) // ReLU
#define FFVV(forcel,forceip,vvv)    ( forcel +	forceip - vvv ) // Change of V,	Total Force
#define MAX(a, b) ( (a) > (b) ? (a) : (b) )
#define V0 1.0
//
struct Matter
{
    double xxx[NB+2],vvv[NB+2],qqx[NB+2];
};
//
//
void rk_rout(struct	Matter *M,double deltat,double xppt);
void initial(struct Matter *M);
void calaveq(struct	Matter *M, double *avxxx, double *avvvv, double *maxvvv);
void calmaxforce(struct Matter *M, double *maxforce, double xppt);
//
//
int main(void)
{
    int flagslip;
    int no_of_earthquake;
    double vplt;
    double deltat;
    double xppt;
    double avxxx, avvvv, maxvvv;
    double avxxx00,avvvv00;
    double xxs[NB+2],vvs[NB+2];
    double maxforce;
    struct Matter MMA, *MM;
    MM=&MMA;

    char arg1[200],arg2[200],arg3[200];
    FILE *file01,*file02,*file03;
    char dirname[]=DATAFILEDIR;
    time_t t1;
    struct tm *tp;

		/* 印字用文字列の格納のため */
		char date[200];

		/* 日付情報の取得 */
		time( &t1 );
		tp = localtime( &t1 );

		/* YYYY/MM/DD HH24:MI:SS 形式で印字 */
		strftime( date, sizeof(date), "%Y%m%d_%H%M%S", tp );

    //* ファイル名をargに書込む
    if(FLAGF01OUT)
        {
            //snprintf(arg1,sizeof(arg1),"%s%s_%s_%s",dirname, __FILE__, date,"OUT");
            snprintf(arg1,sizeof(arg1),"%s_%s_%s",dirname,date,"OUT");
            file01 = fopen(arg1,"w");
        }

    if(FLAGF02MAG)
        {
            //snprintf(arg2,sizeof(arg2),"%s%s_%s_%s",dirname, __FILE__, date,"MAG");
            snprintf(arg2,sizeof(arg2),"%s_%s_%s",dirname,date,"MAG");
            file02 = fopen(arg2,"w");
        }

    if(FLAGF03XXX)
        {
            //snprintf(arg3,sizeof(arg3),"%s%s_%s_%s",dirname, __FILE__, date,"XXX");
            snprintf(arg3,sizeof(arg3),"%s_%s_%s",dirname,date,"XXX");
            file03 = fopen(arg3,"w");
        }

//Initialization
    initial(MM);

    #ifdef DEBUG_INITIAL
        for(int ii = 1; ii	<= NB; ii++)
        {
            printf("ii = %5d  xxx= %10.5lf	vvv= %10.5lf \n",ii,MM->xxx[ii],MM->vvv[ii]);
        }
    #endif

    calaveq(MM, &avxxx, &avvvv, &maxvvv);

    // avxxx00 = avxxx;
    // avvvv00 = avvvv;

    deltat=1.0/pow(2.0,4);

    if( FLAGPRTOUT || FLAGPRTMAG )	printf(        OUTPUTPARAMATERS);

    #ifdef DEBUG_ON
        printf("X01 \n");
        printf("arg1 =	%s \n", arg1);
        printf("arg2 =	%s \n", arg2);
        printf("arg3 =	%s \n", arg3);
    #endif


    if( FLAGF01OUT )               fprintf(file01,OUTPUTPARAMATERS);
    if( FLAGF02MAG )               fprintf(file02,OUTPUTPARAMATERS);
    if( FLAGF03XXX	)              fprintf(file03,OUTPUTPARAMATERS);

    #ifdef DEBUG_ON
        printf("X02 \n");
    #endif


    if( FLAGPRTMAG ) printf (       OUTPUTTITLE);
    if( FLAGF02MAG ) fprintf(file02,OUTPUTTITLE);

    #ifdef DEBUG_ON
        printf("X03 \n");
    #endif

    if( FLAGPRTOUT )
    {
         int istep = 0;
         xppt = 0.0;
         printf (       OUTPUTSTEP);
     }

     if( FLAGF01OUT )
     {
          int istep = 0;
          xppt = 0.0;
          fprintf(file01,OUTPUTSTEP);
      }

    flagslip=0;
    no_of_earthquake=0;
    vplt = VPL;
    xppt = 0.0;

    for(int inoev =0; inoev <= ( INITEVENTNO + EVENTNO ); inoev++)
    {
        double maxforce;

        calmaxforce(MM, &maxforce, xppt);

        xppt += (FRICMAX - maxforce)/KPP;

        avxxx00 = avxxx;
        avvvv00 = avvvv;

        for(int ii = 1; ii	<= NB; ii++)
        {
            xxs[ii] = MM->xxx[ii];
        }

        for(int istep =	0; istep <= STEP_MAX; istep++)
        {
            int istep00;
            double delavxxx, delavvvv;
            double sumavxxx;
            double magnitude;
            double avxxxst;

            if( !flagslip)	avxxxst=avxxx00;

            if ( istep )  xppt += deltat*vplt;

            rk_rout(MM,deltat,xppt);

            calaveq(MM, &avxxx, &avvvv, &maxvvv);

            delavxxx =	avxxx - avxxx00;
            delavvvv =	avvvv - avvvv00;


            if( FLAGPRTOUT ) printf (       OUTPUTSTEP);
            if( FLAGF01OUT ) fprintf(file01,OUTPUTSTEP);

             //judgement of the beginning of	slipping
            if (flagslip == 0 && (maxvvv >=	ACCURACY00) )
            {
                #ifdef DEBUG_ON
                printf("A \n");
                #endif

                flagslip  =1;
                vplt=0.0;
                no_of_earthquake ++;
                sumavxxx  = 0.0;
                istep00=istep;
            }

             //judgement of the continuing of slipping
            if (flagslip == 1 &&  (maxvvv >=	ACCURACY00) )
            {
                #ifdef DEBUG_ON
                printf("B \n");
                #endif

                sumavxxx += delavxxx;
            }

            avxxx00 = avxxx;
            avvvv00 = avvvv;

            //judgement	of the end of slipping
            if (flagslip == 1 &&  (maxvvv < ACCURACY00) )
            {
                #ifdef DEBUG_ON
                printf("C \n");
                #endif

                flagslip	=0;
                vplt=VPL;
                sumavxxx  += delavxxx;
                magnitude = log10(sumavxxx);
                //magnitude = log(sumavxxx);

                if( FLAGPRTMAG ) printf (	  OUTPUTMAG);
                if( FLAGF02MAG ) fprintf(file02,OUTPUTMAG);

                //
                #ifdef DEBUG_ON
                    printf("avxxx-avxxxst = %e	\n", avxxx-avxxxst);
                #endif
                //

                if( FLAGFPRXXX && (no_of_earthquake > INITEVENTNO) )
                {
                    printf("%10.5lf", magnitude);
                    for(int ii = 1; ii <= NB; ii++)
                    {
                        printf(",%10.5lf", xxs[ii]);
                    }
                    printf("\n");
                }

                if( FLAGF03XXX && (no_of_earthquake > INITEVENTNO) )
                {
                    fprintf(file03,"%10.5lf", magnitude);
                    for(int ii = 1; ii <= NB; ii++)
                    {
                        fprintf(file03,",%10.5lf", xxs[ii]);
                    }
                    fprintf(file03,"\n");
                }

            break;

            }

            if( inoev >=		STEP_MAX )
            {
                printf(" step = %d exceeds max step \n", istep);
                return 0;
            }

        }

    }

    if(FLAGF01OUT) fclose(file01);
    if(FLAGF02MAG) fclose(file02);
    if(FLAGF03XXX) fclose(file03);

    return 0;

}
//
//
void initial(struct Matter *MM)
{
    if( FLAGINITRAND )
    {
        srand(ISEED);
        //printf("RAND_MAX: %d\n", RAND_MAX);
        for(int ii = 1; ii	<= NB; ii++)
        {
            double xrand;
            xrand= WIRAND*( rand()/(double)RAND_MAX-0.5 );
            MM->qqx[ii]	= xrand;
            MM->vvv[ii]	= 0.0;
        }
    }
    else
    {
        for(int ii = 1; ii	<= NB; ii++)
        {
            MM->xxx[ii]	= 0.0;
            MM->vvv[ii]	= 0.0;
        }
    }

    MM->xxx[0]     =	MM->xxx[1];
    MM->xxx[NB +1]	= MM->xxx[NB];
    MM->vvv[0]     =	MM->vvv[1];
    MM->vvv[NB +1]	= MM->vvv[NB];

}
//
//
void calaveq(struct	Matter *MM, double	*avxxx, double *avvvv, double *maxvvv)
{

    *avxxx  = 0.0;
    *avvvv  = 0.0;
    *maxvvv = 0.0;

    for(int ii = 1; ii <=	NB; ii++)
    {
        *avxxx += MM->xxx[ii];
        *avvvv += MM->vvv[ii];
        *maxvvv = MAX(*maxvvv, MM->vvv[ii]);
    }

    *avxxx = *avxxx/NB;
    *avvvv = *avvvv/NB;
}
//
//
void rk_rout(struct	Matter *MM,double deltat, double xppt)
{
    double xx1[NB+2],vv1[NB+2];
    double xx2[NB+2],vv2[NB+2];
    double xx3[NB+2],vv3[NB+2];
    double xxk[NB+2],vvk[NB+2];
    double deltatd2;
    double xpph,xppn;


    xpph = xppt +0.5*deltat*VPL;
    xppn = xppt +1.0*deltat*VPL;

    deltatd2=0.5*deltat;

    //RK 1st-step
    for(int ii = 1; ii <=	NB; ii++)
    {
        double kx1,kv1,forceel0,forceip0;

        kx1      = deltatd2*MM->vvv[ii];
        xxk[ii]  = kx1;
        xx1[ii]  = MM->xxx[ii] + kx1;
        forceel0 = ELFORCE(MM->xxx[ii],MM->xxx[ii-1],MM->xxx[ii+1],xppt) ;
        forceip0 = V0 *sin(MM->xxx[ii] + MM->qqx[ii]);
        kv1      = deltatd2*FFVV(forceel0,forceip0,MM->vvv[ii]);
        vvk[ii]  = kv1;
        vv1[ii]  = MM->vvv[ii] + kv1;
    }

    xx1[0]    = xx1[1];
    xx1[NB+1] = xx1[NB];

//RK 2nd-step
    for(int ii = 1; ii <=	NB; ii++)
    {
        double kx2,kv2,forceel1,forceip1;

        kx2      = deltatd2*vv1[ii];
        xxk[ii]  = xxk[ii] + 2.0*kx2;
        xx2[ii]  = xx1[ii] + kx2;
        forceel1 = ELFORCE(xx1[ii],xx1[ii-1],xx1[ii+1],xppt) ;
        forceip1 = V0 *sin(xx1[ii] + MM->qqx[ii]);
        kv2      = deltatd2*FFVV(forceel1,forceip1,vv1[ii]);
        vvk[ii]  = vvk[ii] + 2.0*kv2;
        vv2[ii]  = MM->vvv[ii] + kv2;
    }

    xx2[0]    = xx2[1];
    xx2[NB+1] = xx2[NB];


    //RK 3rd-step
    for(int ii = 1; ii <=	NB; ii++)
    {
        double kx3,kv3,forceel2,forceip2;

        kx3      = deltat*vv2[ii];
        xxk[ii]  = xxk[ii] + kx3;
        xx3[ii]  = MM->xxx[ii] + kx3;
        forceel2 = ELFORCE(xx2[ii],xx2[ii-1],xx2[ii+1],xppt) ;
        forceip2 = V0 *sin(xx2[ii] + MM->qqx[ii]);
        kv3      = deltat*FFVV(forceel2,forceip2,vv2[ii]);
        vvk[ii]  = vvk[ii] + kv3;
        vv3[ii]  = MM->vvv[ii] + kv3;
    }

    xx3[0]    = xx3[1];
    xx3[NB+1] = xx3[NB];

//RK 4th-step
    for(int ii = 1; ii <=	NB; ii++)
    {
        double forceel3,forceip3,vvvtmp;

        //Caution 20231024
        xxk[ii]     = 3.3333333333333333E-1*(xxk[ii] + deltatd2*vv3[ii] );
        MM->xxx[ii] =	MM->xxx[ii] + xxk[ii];
        forceel3 = ELFORCE(xx3[ii],xx3[ii-1],xx3[ii+1],xppt) ;
        forceip3 = V0 *sin(xx3[ii] + MM->qqx[ii]);
        vvk[ii]     = 3.3333333333333333E-1*(vvk[ii] + deltatd2*FFVV(forceel3,forceip3,vv3[ii]) );
        //vvvtmp  = vvv[ii] + vvk[ii];
        MM->vvv[ii] =	MM->vvv[ii] +	vvk[ii];
    }

    MM->xxx[0]    =	MM->xxx[1];
    MM->xxx[NB+1]	= MM->xxx[NB];

}
//
//
void calmaxforce(struct Matter *MM, double *maxforce, double xppt)
{

    *maxforce  = -1.0E10;

    for(int ii = 1; ii <=	NB; ii++)
    {
        *maxforce = MAX( *maxforce, ELFORCE(MM->xxx[ii],MM->xxx[ii-1],MM->xxx[ii+1],xppt) );
    }

}

//
