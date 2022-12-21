/** kinema.h ************/
/* 00/05/15 add new macro _KINE to avoid doublicate inclusion. */
/* 03/07/14 add new function calcp3lab. */
#include <math.h>

#ifndef _KINE
#define _KINE 1
#define AMU 931.494
#define MP 938.272
#define CV 299792458
#define HBARC 197.327053
#define D_TO_R (M_PI/180.0)
#define R_TO_D (180.0/M_PI)


struct skine {
  double m1,m2,m3,m4; /* rest mass for each particle *//********/
  /*double ex3,ex4;*/  /* excitation energy */
  double p1,p2,p3[2],p4[2]; /* momentum for each particle in lab */
  double p1c,p2c,p3c,p4c; /* momentum for each particle in cms */

  double b1c,b2c,b3c,b4c; /* velocity for each particle in cms */
  double g1c,g2c,g3c,g4c; /* gamma factor for each particle in cms */

  double E1,E2,E3,E4; /* energy of each particle in lab */
                      /* E1, E2 */
  double E_cm;        /* energy in cms */
  double E1c,E2c,E3c,E4c; /* energy of each particle in cms */
  double K1,K2,K3[2],K4[2];  /* kinematic energy of each particle in lab */
                          /* K1,K2 */
  double Q_val;           /* reaction Q-value */  /*****/
  double delta23;           /* delta-value delta23=beta_2/beta_3 */
  double delta24;           /* delta-value delta24=beta_2/beta_4 */
  double factor[2];          /* dOmega_lab/dOmega_cm */
  double qtra[2];            /* q-transfer */
  double curang;   /* current angle */  /*************/
  double curang_r;
  double labang;   /* current angle in lab*/
  double labang_r;

  double dexdk1[2],dexdk3[2],dexdt3[2];  /* Differential excitation energies */

  int numofang;   /* number of scattering angle (0 or 1 or 2)*/
  double ang_cm[2];   /* current angle in cms*/
  double ang_cm_r[2];   /* current angle in cms*/
  double ang_lb4_r[2];   /* Angle 4 in lab in radian. */
  double ang_lb4[2];     /* Angle 4 in lab. */
  double ang_cm4_r[2];   /* Angle 4 in cm in radian. */
  double ang_cm4[2];     /* Angle 4 in cm. */
  
  double W;  /* total cm energy */
  double gamma_cm; /* gamma of the cm */
  double beta_cm;  /* velocity of the cm */
};

int calccmang(double *cmang,double labang,double gamma,double delta);
double calclabk(double cmang_r,double gamma2,double p3c,double m3);
double calcfact(double cm_ang_r,double lab_ang_r,double gamma2,double delta);
double calcex4(double lab_ang_r,double m1,double m2,
                double m3,double m4,double K1,double K3);
double calcp3lab(double lab_ang_r,double m1,double m2,
				 double m3,double m4,double K1);
int calcscatt(struct skine *kine,int icmlab);
double calclabang(double *labang_r,double cmang_r,double gamma2,double delta);

#endif

