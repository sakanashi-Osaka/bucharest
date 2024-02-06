/** kinelib.c *************************/
/** calculate relation between q-value and momentum **/
/* 03/07/14 add new function calcp3lab. */

#include <stdio.h>
#include <math.h>
#include "kinema.h"

/**************************************************************************
int calccmang(double *cmang_r,double labang_r,double gamma2,double delta)

    Calculate scattering angle in the CM frame.
    double *cmang_r .. pointer to results.
                       calccmang return one or two angles.
                       results were expressed in radian.
    double labang_r .. scattering angle in the LAB frame in radian.
    double gamma2   .. gamma factor for center of mass.
                       (same with gamma factor for particle 2 in the CM frame)
    double delta    .. beta2/beta3 in the CM frame.
    return value    .. number of results.
    ***********************************************************************/

int calccmang(double *cmang_r,double labang_r,double gamma2,double delta){
  double firstterm,secondterm;
  int n=0,i;
  
  double xsq,tmpang;
  double cmcos[2];

  xsq=pow(gamma2*tan(labang_r),2.0);
  firstterm=-delta*xsq/(1.0+xsq);
  /*  secondterm=sqrt(pow(delta*xsq/(1.0+xsq),2.0)-(xsq-n1.0)/(xsq+1.0));*/
  secondterm=sqrt((1.0-delta*delta)*xsq+1)/(1.0+xsq);
  cmcos[0]=firstterm+secondterm;
  cmcos[1]=firstterm-secondterm;


  if(fabs(cmcos[0])< 1 || cmcos[0]==1) {
    tmpang=acos(cmcos[0]);
    if((cos(tmpang)+delta)*cos(labang_r)>0){
      *(cmang_r++)=tmpang;
      n++;
    }
  }

  if(((fabs(cmcos[1])<=1) && (fabs(secondterm)>10e-10)) || (cmcos[1]==-1)) {
    tmpang=acos(cmcos[1]);
    if((cos(tmpang)+delta)*cos(labang_r)>0){
      *(cmang_r++)=tmpang;
      n++;
    }
  }
  return(n);
}

/**************************************************************************
double calclabk(double cmang_r,double gamma2,double p3c,double m3)
    Calculate kinematic energy of particle 3 in the LAB frame.
     double cmang_r .. scattering angle of paticle 3 in the CM frame in radian.
                        (CM angle can be calculated by calccmang)
     double gamma2  .. gamma factor for center of mass.
                        (same with gamma factor for particle 2 in the CM frame)
     double p3c     .. momentum of particle 3 in the CM frame in MeV/c.
     double m3      .. mass of partilce 3 in MeV.
     return value   .. kinematic energy of paritcle 3in the LAB frame.
     *********************************************************************/

double calclabk(double cmang_r,double gamma2,double p3c,double m3){
  double E3,b2c,E3c;

  b2c=sqrt(1.0-1.0/(gamma2*gamma2));
  E3c=sqrt(m3*m3+p3c*p3c);
  E3=gamma2*(b2c*p3c*cos(cmang_r)+E3c);
  return(E3-m3);
}


/**************************************************************************
double calcfact(double cm_ang_r,double lab_ang_r,double gamma2,double delta)
   Caluculate conversion factor of cross section. (dOmega_lab/dOmega_cm)
     double cmang_r .. scattering angle of paticle 3 in the CM frame in radian.
                        (CM angle can be calculated by calccmang)
     double labang_r .. scattering angle in the LAB frame in radian.
     double gamma2   .. gamma factor for center of mass.
                        (same with gamma factor for particle 2 in the CM frame)
     double delta   ... beta2/beta3 in the CM frame.
     return value   ... dOmega_lab/dOmega_cm.
     ************************************************************************/

double calcfact(double cm_ang_r,double lab_ang_r,double gamma2,double delta){
  double fac1,fac2u,fac2d;
  /*  fac1=sin(lab_ang_r)/sin(cm_ang_r);
  fac2u=gamma2*sin(lab_ang_r)*sin(cm_ang_r)+cos(lab_ang_r)*cos(cm_ang_r);
  fac2d=gamma2*cos(lab_ang_r)*cos(cm_ang_r)+sin(lab_ang_r)*sin(cm_ang_r)
    +gamma2*delta;*/
  fac2u=gamma2*(1.0+delta*cos(cm_ang_r));
  fac1=pow(sin(cm_ang_r),2.0)+pow(gamma2*(delta+cos(cm_ang_r)),2.0);
  fac2d=pow(fac1,1.5);
  return(fac2u/fac2d);
}


/**************************************************************************
double calcex4(double lab_ang_r,double m1,double m2,
		double m3,double m4,double K1,double K3){
   Caluculate excitation energy of particle 4.
     double labang_r .. scattering angle in the LAB frame in radian.
     double m1       .. rest mass of projectile (particle 1)
     double m2       .. rest mass of target (particle 2)
     double m3       .. rest mass of scattered partile (particle 3)
     double m4       .. rest mass of ground state of
                                       recoil partile (particle 4)
     return value   ... reaction q-vaule
     m1, m2 and m4 should be given with intrinsic excitation energy.
     ************************************************************************/

double calcex4(double lab_ang_r,double m1,double m2,
		double m3,double m4,double K1,double K3){
  double s,t,u;
  double p1,p3;
  double E1,E3;
  double ex4,totalm4;
  double kernel;

  E1=m1+K1;
  E3=m3+K3;

  p1=sqrt(E1*E1-m1*m1);
  p3=sqrt(E3*E3-m3*m3);

  s=m1*m1+m2*m2+2*m2*E1;
  t=m1*m1+m3*m3+2*(p1*p3*cos(lab_ang_r)-E1*E3);
  u=m2*m2+m3*m3-2*m2*E3;

  kernel=s+t+u-m1*m1-m2*m2-m3*m3;
  if(kernel<0){
    ex4=-K1;
  }
  else{
    totalm4=sqrt(kernel);
    ex4=totalm4-m4;
  }
  return(ex4);
}


/**************************************************************************
double calcp3lab(double lab_ang_r,double m1,double m2,
		double m3,double m4,double K1){
   Caluculate momentum of particle 3.
     double labang_r .. scattering angle in the LAB frame in radian.
     double m1       .. rest mass of projectile (particle 1)
     double m2       .. rest mass of target (particle 2)
     double m3       .. rest mass of scattered partile (particle 3)
     double m4       .. rest mass of ground state of
                                       recoil partile (particle 4)
     return value   ... reaction q-vaule
     m1, m2 and m4 should be given with intrinsic excitation energy.
     ************************************************************************/

double calcp3lab(double lab_ang_r,double m1,double m2,
		double m3,double m4,double K1){
  double E1,E2,K3,W;
  double gamma_cm,beta_cm,g2c,delta;
  double p1,p2,p3;
  double p1c,p2c,p3c,p4c;
  double E1c,E2c,E3c,E4c;
  double b1c,b2c,b3c,b4c;
  double ang_cm_r[2];
  
  int numofang;

  E1=m1+K1;
  E2=m2;
  p1=sqrt(E1*E1-m1*m1);
  p2=0;
  W=sqrt(m1*m1+m2*m2+2*m2*E1);

  gamma_cm=(E1+m2)/W;
  beta_cm=p1/(E1+m2);

  p1c=1.0/(2*W)*sqrt((W*W-(m1+m2)*(m1+m2))*(W*W-(m1-m2)*(m1-m2)));
  p2c=p1c;
  p3c=1.0/(2*W)*sqrt((W*W-(m3+m4)*(m3+m4))*(W*W-(m3-m4)*(m3-m4)));
  p4c=p3c;

  E1c=sqrt(m1*m1+p1c*p1c);
  E2c=sqrt(m2*m2+p2c*p2c);
  E3c=sqrt(m3*m3+p3c*p3c);
  E4c=sqrt(m4*m4+p4c*p4c);
  g2c=E2c/m2;

  b1c=p1c/E1c;
  b2c=p2c/E2c;
  b3c=p3c/E3c;
  b4c=p4c/E4c;
  delta=b2c/b3c;

  numofang=calccmang(ang_cm_r,lab_ang_r,g2c,delta);
  switch(numofang){
  case 2:	
	ang_cm_r[0]=(ang_cm_r[0]<ang_cm_r[1]?ang_cm_r[0]:ang_cm_r[1]);
  case 1:
	K3=calclabk(ang_cm_r[0],g2c,p3c,m3);
	p3=sqrt(pow(K3+m3,2.0)-m3*m3);
	break;
  default:
	p3=0;
	break;
  }
  return(p3);
}


/**************************************************************************
double calcscatt(struct skine *kine,,int icmlab);
Input parameters
  icmlab ... Flag to define CM/LAB frame of the scattered angle 
  skine->m1  .... Mass of particle 1 including Ex (MeV).
  skine->m2  .... Mass of particle 2 including Ex (MeV).
  skine->m3  .... Mass of particle 3 including Ex (MeV).
  skine->m4  .... Mass of particle 4 including Ex (MeV).
  skine->K1  .... Energy of particle 1 in lab in MeV.
  skine->curang ... Scattering angle of particle 3 in lab (if icmlab=0)
                    or in CM (if icmlab=1).
Return value: 
   >=0 : No particle detected
    -1 : Incorrect icmlab
*/

int calcscatt(struct skine *kine,int icmlab){
  int i,j,k;

  kine->curang_r=kine->curang*D_TO_R;
  kine->K2=0;


/**** calculate some values********************************/
    kine->E1=kine->m1+kine->K1;
    kine->E2=kine->m2;
    kine->Q_val=kine->m3+kine->m4-kine->m1-kine->m2;
    kine->p1=sqrt(kine->E1*kine->E1-kine->m1*kine->m1);
    kine->p2=0;
    kine->W=sqrt(kine->m1*kine->m1+kine->m2*kine->m2+2*kine->m2*kine->E1);

    kine->gamma_cm=(kine->E1+kine->m2)/kine->W;
    kine->beta_cm=kine->p1/(kine->E1+kine->m2);
    kine->p1c=1.0/(2*kine->W)
      *sqrt((kine->W*kine->W-(kine->m1+kine->m2)
	     *(kine->m1+kine->m2))
	    *(kine->W*kine->W-(kine->m1-kine->m2)*(kine->m1-kine->m2)));
    kine->p2c=kine->p1c;
    kine->p3c=1.0/(2*kine->W)
      *sqrt((kine->W*kine->W-(kine->m3+kine->m4)
	     *(kine->m3+kine->m4))
	    *(kine->W*kine->W-(kine->m3-kine->m4)*(kine->m3-kine->m4)));
    kine->p4c=kine->p3c;
    kine->E1c=sqrt(kine->m1*kine->m1+kine->p1c*kine->p1c);
    kine->E2c=sqrt(kine->m2*kine->m2+kine->p2c*kine->p2c);
    kine->E3c=sqrt(kine->m3*kine->m3+kine->p3c*kine->p3c);
    kine->E4c=sqrt(kine->m4*kine->m4+kine->p4c*kine->p4c);
    kine->g2c=kine->E2c/kine->m2;

    kine->b1c=kine->p1c/kine->E1c;
    kine->b2c=kine->p2c/kine->E2c;
    kine->b3c=kine->p3c/kine->E3c;
    kine->b4c=kine->p4c/kine->E4c;
    kine->delta23=kine->b2c/kine->b3c;
    kine->delta24=kine->b2c/kine->b4c;

    switch(icmlab){
    case 1:  /* if cm angle is given */
      kine->ang_cm[0]=kine->curang;
      kine->ang_cm_r[0]=kine->curang_r;
      kine->labang=calclabang(&kine->labang_r,kine->curang_r,
			      kine->g2c,kine->delta23);
      kine->numofang=1;
      break;
    case 0:  /* if lab angle is given */
      kine->labang=kine->curang;
      kine->labang_r=kine->curang_r;
      kine->numofang=calccmang(kine->ang_cm_r,kine->curang_r,
			       kine->g2c,kine->delta23);
      if(kine->numofang<1) return(0); /* No particle detected */
      break;
    default:
      return(-1);
    }

    for(j=0;j<kine->numofang;j++){
      double qval;
      kine->ang_cm[j]=kine->ang_cm_r[j]*R_TO_D;
      kine->qtra[j]=sqrt(kine->p1c*kine->p1c+kine->p3c*kine->p3c
		      -2*kine->p1c*kine->p3c*cos(kine->ang_cm_r[j]))/HBARC;
      kine->K3[j]=calclabk(kine->ang_cm_r[j],kine->g2c,kine->p3c,kine->m3);
      kine->p3[j]=sqrt(2*kine->m3*kine->K3[j]+kine->K3[j]*kine->K3[j]);
      kine->factor[j]=calcfact(kine->ang_cm_r[j],kine->labang_r,
			    kine->g2c,kine->delta23);
      qval=calcex4(kine->labang_r,kine->m1,kine->m2,kine->m3,kine->m4,
		   kine->K1,kine->K3[j]);

      /* calculate recoil particle */
      kine->ang_cm4[j]=180-kine->ang_cm[j];
      kine->ang_cm4_r[j]=kine->ang_cm4[j]*D_TO_R;
      kine->ang_lb4[j]=calclabang(&kine->ang_lb4_r[j],kine->ang_cm4_r[j],
				  kine->g2c,kine->delta24);
      kine->K4[j]=kine->m1+kine->K1+kine->m2-(kine->m3+kine->K3[j]+kine->m4);
      /* differential factors */
      kine->dexdk1[j]=(kine->m2-kine->m3-kine->K3[j]+kine->p3[j]/kine->p1
		       *(kine->K1+kine->m1)*cos(kine->labang_r))/kine->m4;
      kine->dexdk3[j]=(-kine->m1-kine->m2-kine->K1+kine->p1/kine->p3[j]
		       *(kine->K3[j]+kine->m3)*cos(kine->labang_r))/kine->m4;
      kine->dexdt3[j]=-kine->p1*kine->p3[j]*sin(kine->labang_r)/kine->m4;
    }
    return(kine->numofang);
}

/**************************************************************************
double calclabang(double *labang_r,double cmang_r,double gamma2,double delta)

    Calculate scattering angle in the CM frame.
    double *labang_r .. pointer to results.
                       LAB angle in radian
    double cmang_r ..  scattering angle in the CM frame in radian.
    double gamma2   .. gamma factor for center of mass.
                       (same with gamma factor for particle 2 in the CM frame)
    double delta    .. beta2/beta3 in the CM frame.
    return value    .. LAB angle in degree
    ***********************************************************************/

double calclabang(double *labang_r,double cmang_r,double gamma2,double delta){

  double numera,denomi;
  double cocm,sicm;  
  double labcos;

  cocm=cos(cmang_r);sicm=sin(cmang_r);

  numera=gamma2*(cocm+delta);
  denomi=sqrt(sicm*sicm+pow(gamma2*(delta+cocm),2.0));

  labcos=numera/denomi;

  /*
  printf("\n cmang=%7.4f  delta=%7.4f  cos=%7.4f\n",
	 cmang_r*R_TO_D,delta,labcos);
	 printf(" cocm=%7.4f  sicm=%7.4f\n",cocm,sicm);*/

  *labang_r=acos(labcos);
  return(*labang_r*R_TO_D);
}
