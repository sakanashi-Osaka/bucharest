#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <strings.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdint>
#include <vector>
#include <signal.h>
#include <map>
#include <algorithm>
#include <regex>
#include <TCanvas.h>
#include <TEventList.h>
#include <TH1F.h>
#include <TH2.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include "lorlib.cpp"
#include "kinelib.cpp"


#define MAX_STAT 3 // number of state
#define N_ALPHA 10

using namespace std;

double mexes[N_ALPHA]={  /* Mass excess of 4N nuclei. */
 2.424911, 4.94166, 0.0000, -4.736998, -7.041929, -13.93340, -21.49283,  -26.01694, -30.2305, -34.8461
};
double mgs[N_ALPHA]; /* Mass of 4N nuclei in MeV */
double th_ene[N_ALPHA]; /* Threshold energy of alpha deay in 4N nuclei in MeV */


struct ssimval {
  /* General conditions **************************************/
  int nstat;  /* number of state */
  double enepa;
  double a1;    /* mass numbers */
  double mvec1[4];  /* Measured four vector of incident particle */
  double mvec2[4];  /* Measured four vector of target proton */
  double ex_ene[MAX_STAT]; /* Excitation energy */

/* Simulation condition **************************************/
  double dSi; /* distance between target and Si in mm */
  double diaSi; /* diameter of  Si in mm */
  double innSi; /* inner of  Si in mm */
  int nxSi; /* Number of horizontal strips in Si */
  int nySi; /* Number of vertical strips in Si */
  double thSi; /* Angle of Si in degree */
  double thrSi;   /* Angle of Si in radian */
  double tdeadSi;  /* Dead layer in Si in mg/cm2 */

  double sigK1;  /* Energy spread of the beam in MeV */
  double sigthr1; /* In-plane angular spread of the beam in radian */
  double sigx1;  /* Lateral spread of beam spot in mm */
  double offx1;  /* Lateral offset of beam spot in mm */
  double offy1;  /* Lateral offset of beam spot in mm */
  double tthick;   /* Target thickness in ug/cm^2 */
  double sigK3;  /* Energy resolution of recoiled 12C in MeV */
  double sigK4;  /* Energy resolution of recoiled 4He in MeV */
  double sigthxrC;  /* Angular resolution in theta of scattered 12C */
  double sigthyrC;  /* Angular resolution in phi of scattered 12C */
  /*************************************************************/
    
  /* Simulated variables **************************************/
  /* s: simulated, m: measured */
  /* Section 1: Up to the collision point */
  /* [0] ... before target, [1] ... before scattering */
  double sK1[2];  /* Energy of 4He Beam */
  double sthxr1[2]; /* Trajectory angle of 4He Beam in horizontal plane */
  double sthyr1[2]; /* Trajectory angle of 4He Beam in vertical plane */
  double sx1;    /* Lateral spread of 4He Beam in horizontal plane */
  double sy1;    /* Lateral spread of 4He Beam in vertical plane */
  double sdep1;  /* Depth where the reaction occurs */
  double sthr1[2];  /* Polar angle of 4He Beam */
  double sphr1[2];  /* Azimuthal angle of 4He Beam */
  double sth1[2];  /* Polar angle of 4He Beam in degree */
  double sph1[2];  /* Azimuthal angle of 4He Beam  in degree */

  
  /* Section 2: From the collision point to the exit of the target */
  double sthrcm; /* Polar scattering angle in CM */
  double sphrcm; /* Azimuthal scattering angle in CM */
  /* [0]  ... after the scattering, [1] ... after the target */
  double sK3[2]; /* Energy of particle 3 */
  double sK4[2]; /* Energy of particle 4 */
  double sthr3[2];  /* Polar angle of particle3 */
  double sphr3[2];  /* Azimuthal angle of particle3 */
  double sthr4[2];  /* Polar angle of 4 after the collision point */
  double sphr4[2];  /* Azimuthal angle of 4 after the collision point */
  double sth3[2];  /* Polar angle of particle3 in degree */
  double sph3[2];  /* Azimuthal angle of particle3 in degree */
  double sth4[2];  /* Polar angle of 4 after the collision point in deg. */
  double sph4[2];  /* Azimuthal angle of 4 after the collision point in deg. */
  double astra3;  /* In-plane ang. stragg. of particle 3 in target */
  double astra4;  /* In-plane ang. stragg. of particle 4 in target */
  double sthrstra3; /* Polar ang. stragg. of particle 3 in target */
  double sphrstra3; /* Azimuthal ang. stragg. of particle 3 in target */
  double sthrstra4; /* Polar ang. stragg. of particle 4 in target */
  double sphrstra4; /* Azimuthal ang. stragg. of particle 4 in target */

  
  /* Section 3: In-flight decay from 12C to 3a */
  /* Observables after collision point ********************/
  double savec[3][4]; /* 4 momentum of alpha particles */
                         /*             after collision point */
                         /* [0][] .. 1st alpha emitted into a + 8Be channel */
                         /* [1][] .. 2nd alpha emitted from 8Be */
                         /* [2][] .. 2nd alpha emitted from 8Be */
  double sbvec[4];    /* 4 momentum of 8Be */
  double sKa[3];    /* Energy of alpha particles */
  double sthra[3]; /* Polar angle of alpha particles in radian */
  double sphra[3]; /* Azimuthal angle of alpha particles in radian*/
  double stha[3];  /* Polar angle of alpha particles in degree */
  double spha[3];  /* Azimuthal ang. of alpha parcticles in degree */

  /* Observables after the target *******************/
  /* a ... alpha just after collision, af ... after the target */
  double safvec[3][4]; /* 4 momentum of alpha particles */
                         /*             after collision point */
                         /* [0][] .. 1st alpha emitted into a + 8Be channel */
                         /* [1][] .. 2nd alpha emitted from 8Be */
                         /* [2][] .. 2nd alpha emitted from 8Be */

  double sK3a[3];    /* Energy of alpha particles */
  double sthraf[3]; /* Polar angle of alpha particles in radian */
  double sphraf[3]; /* Azimuthal angle of alpha particles in radian*/
  double sthaf[3];  /* Polar angle of alpha particles in degree */
  double sphaf[3];  /* Azimuthal ang. of alpha parcticles in degree */

  double astraf[3];  /* In-plane ang. stragg. of alpha particle in target */
  double sthrstraf[3]; /* Polar ang. stragg. of alpha particle in target */
  double sphrstraf[3]; /* Azim. ang. stragg. of alpha particle in target */

  double stotea;  /* Total energy of 6 alpha particles in the target */
  double stoteaf;  /* Total energy of 6 alpha particles after the target */
  int sna; /* Number of alpha particles which escape from the target */

  
  /* Section 4: 4He detection */
  int idp4He; /* Detection flag for particle 4 0..Not detected  1..detected */
  double s4Hevec[2][3]; /* 3-dim momentum of 4He */
                      /* [0][] ... Detector frame */
                      /* [1][] ... Laboratory frame */
  double s4Hepos[3]; /* Position in the detector frame */
  double mK4He;  /* Measured energy of 4He */
  double m4Hepos[2][3]; /* Measured position of proton */
                      /* [0][] ... Detector frame */
                      /* [1][] ... Laboratory frame */
  double mthr4He;  /* Measured polar angle of 4He in radian */
  double mphr4He;  /* Measured azimuthal angle of 4He in radian */
  double mth4He;  /* Measured polar angle of 4He in degree */
  double mph4He;  /* Measured azimuthal angle of 4He in degree */
  double mthxr4He; /* Measured horizontal angle of 4He in radian */
  double mthyr4He; /* Measured vertical angle of 4He in radian */
  double mthx4He;  /* Measured horizontal angle of 4He in degree */
  double mthy4He;  /* Measured vertical angle of 4He in degree */
  double mexC;  /* Measured excitation energy of 12C */
  int m4Heseg;
  int isx4HeSi; /* ID of hit channel of Si (horizontal) */
  int isy4HeSi; /* ID of hit channel of Si (vertical) */
  double mvec4He[4];  /* Measured four vector of recoiled 4He */

  
  /* Section 5: 3a detection */
  int idp3a[3]; /* Detection flag for 3a 0..Not detected  1..detected */
  double s3avec[2][3][3]; /* 3-dim momentum of 3a */
  double s3apos[3][3]; /* Position in the detector frame */
  double mK3a[3];  /* Measured energy of 3a */
  double m3apos[2][3][3]; /* Measured position of 3a */
  double mthr3a[3];  /* Measured polar angle of 3a in radian */
  double mphr3a[3];  /* Measured azimuthal angle of 3a in radian */
  double mth3a[3];  /* Measured polar angle of 3a in degree */
  double mph3a[3];  /* Measured azimuthal angle of 3a in degree */
  double mthxr3a[3]; /* Measured horizontal angle of 3a in radian */
  double mthyr3a[3]; /* Measured vertical angle of 3a in radian */
  double mthx3a[3];  /* Measured horizontal angle of 3a in degree */
  double mthy3a[3];  /* Measured vertical angle of 3a in degree */
  double mex3a;  /* Measured excitation energy of 3a */
  double sex3a;  /* simutated excitation energy of 3a */
  int isx3aSi[3]; /* ID of hit channel of Si (horizontal) */
  int isy3aSi[3]; /* ID of hit channel of Si (vertical) */
  int m3aseg[3];
  double svec3a[3][4];  /* Simulated four vector of recoiled 3a */
  double mvec3a[3][4];  /* Measured four vector of recoiled 3a */

  /* Section 6: In-flight gamma decay of 12C */
  double scvec[MAX_STAT][4]; /* 4 momentum of 12C in the i th state */
  double sgvec[MAX_STAT][4]; 
  /* 4 mom. of gamma ray emitted from the i th state */
  double sKC[MAX_STAT]; /* Energy of 12C in the i th state */
  double spC[MAX_STAT]; /* Momentum of 12C in the i th state */
  double sthrC[MAX_STAT]; /* Polar angle of 12C in the i th state */
  double sphrC[MAX_STAT]; /* Azimuthal angle of 12C in the i th state */
  double sthC[MAX_STAT];  /* Polar angle of 12C in the i th state in degree */
  double sphC[MAX_STAT];  /* Azimuthal ang. of 12C in the i th state in deg */

  
    /* Section 7: 12C detection */
  int idp12C; /* Detection flag for particle 3 0..Not detected  1..detected */
  double s12Cvec[2][3]; /* 3-dim momentum of 12C */
  double s12Cpos[3]; /* Position in the detector frame */
  double mK12C;  /* Measured energy of 12C */
  double m12Cpos[2][3]; /* Measured position of 12C */
  double mthr12C;  /* Measured polar angle of 12C in radian */
  double mphr12C;  /* Measured azimuthal angle of 12C in radian */
  double mth12C;  /* Measured polar angle of 12C in degree */
  double mph12C;  /* Measured azimuthal angle of 12C in degree */
  double mthxr12C; /* Measured horizontal angle of 12C in radian */
  double mthyr12C; /* Measured vertical angle of 12C in radian */
  double mthx12C;  /* Measured horizontal angle of 12C in degree */
  double mthy12C;  /* Measured vertical angle of 12C in degree */
  double mex12C;  /* Measured excitation energy of 12C */
  int m12Cseg;
  int isx12CSi; /* ID of hit channel of Si (horizontal) */
  int isy12CSi; /* ID of hit channel of Si (vertical) */
  double mvec12C[4];  /* Measured four vector of recoiled 12C */  
};


int genalpha(int i, double *p0, double *p1, double *p2);
int gengamma(double eg,double *p0, double *p1, double *p2);
double deflect2(double *v1,double *v2, double th1,double ph1,double th2,double ph2);
double genvdir(double ang,double *v);
double genvthph(double th,double  ph,double *v);
double genini(double *p, double *v, double ex, double t, double m);
double gettheta4(double *p);
double gettheta(double *p);
double getphi4(double *p);
double getphi(double *p);
double genrndg(double *a, double *b);

/* simulation routine ***********************/
double create4He(struct skine *k, struct ssimval *s);
double scatt4He12c(struct skine *k, struct ssimval *s);
double decay3a(struct skine *k, struct ssimval *s);
double decay12C(struct skine *k, struct ssimval *s);
int detect4He(struct skine *k, struct ssimval *s);
int detect3a(struct skine *k, struct ssimval *s);
int detect12C(struct skine *k, struct ssimval *s);

/**multiple scattering************************************/
double eloss4He(double ene,double thick, double tthick);
double estra4He(double ene,double thick, double tthick);
double astra4He(double ene,double thick, double tthick);
double elossc(double ene,double thick, double tthick);
double estrac(double ene,double thick, double tthick);
double astrac(double ene,double thick, double tthick);

double tmp(double thetacm_i, double thetacm_f, int n_div, int i_div);



int main(int argc, char *argv[]){
 
  struct ssimval ssim;
  struct skine kine;

/* Initialization of simulation conditions ***********************/
  bool flag_ang = true; // true: apply the angular distribution, false: isotropic distribution 

  ssim.enepa=25; /* incident energy per */
  ssim.nstat=3; /* number of state */
  //  ssim.ex_ene[0]=9.64;
  ssim.ex_ene[0]=7.65;  /* default values */
  ssim.ex_ene[1]=4.44;
  ssim.ex_ene[2]=0;
  ssim.a1=4;    /* mass numbers of beam particle*/
  
  ssim.thSi=0; /* Angle of Si in degree */
  ssim.dSi=40.0; /* distance between target and Si in mm */
  ssim.diaSi=96; /* diameter of Si in mm */
  ssim.innSi=48; /* inner of Si in mm */
  ssim.nxSi=16*4; /* Number of horizontal strips in Si */
  ssim.nySi=16; /* Number of vertical strips in Si */
  ssim.thrSi=ssim.thSi*D_TO_R;   /* Angle of Si in radian */
  ssim.tdeadSi=0;  /* Dead layer in Si in mg/cm2 */

  ssim.sigK1=0.05;  /* Energy spread of the beam in MeV */
  ssim.sigthr1=0.2*D_TO_R; /* In-plane angular spread of the beam in radian */ 
  ssim.sigx1=0.5;  /* Lateral spread of beam spot in mm */
  ssim.offx1=0;  /* Lateral offset of beam spot in mm */
  ssim.offy1=2.4;  /* Lateral offset of beam spot in mm */
  ssim.tthick=50.0;   /* Target thickness in ug/cm^2 */
  ssim.sigK3=0.02;  /* Energy resolution of scatterd 12C in MeV */ //%%%
  ssim.sigK4=0.02;  /* Energy resolution of recoiled 4He in MeV */ //%%%
  ssim.sigthxrC=0.002;  /* Ang. resolution in theta of scattered 12C in rad.*/
  ssim.sigthyrC=0.002;  /* Ang. resolution in phi of scattered 12C in rad. */

  //#define _IDEAL
#ifdef _IDEAL    /* ideal conditions */
  ssim.tthick=0;
  ssim.sigx1=0;
  ssim.sigK1=0.00;
  ssim.sigthr1=0;
  ssim.sigK3=0.00;
  ssim.sigK4=0.00;
  ssim.tdeadSi=0;
  ssim.sigthxrC=0.0;
  ssim.sigthyrC=0.0;
#endif
  /*******************************************************************/

    
  /*sim. condition output*/  
  FILE *outputfile2;
  const char *txt_name2;
  if(argc>1){
    char str[100];
    sprintf(str, "condition/%s.txt", argv[1]);
    txt_name2=str;
  }else{
    printf("file name is not input\n");
  }
  outputfile2 = fopen(txt_name2, "w");  
  if (outputfile2 == NULL) {
    printf("cannot open %s\n",txt_name2);  
    exit(1);
  }
  
  time_t tnow;
  tnow=time(NULL);
  printf("\n/***** simulation condition *****/\n");
  printf("Time: %s",ctime(&tnow));
  printf("Output: %s\n",argv[1]);
  printf("Number of States: %d\n",ssim.nstat);
  printf("Mass Number: %d\n",(int)ssim.a1);
  printf("Incident energy: %7.1f MeV (%5.1f MeV/u)\n", ssim.enepa,ssim.enepa/ssim.a1);
  printf("  beam position: x=%6.3f mm, y=%6.3f mm\n",ssim.offx1, ssim.offy1);
  printf("  Energy spread: %6.3f MeV\n",ssim.sigK1);
  printf("  Angular spread: %6.3f rad\n",ssim.sigthr1);
  printf("  Lateral spread: %6.3f mm\n",ssim.sigx1);
  printf("  Target thickness: %5.1f mg/cm2\n",ssim.tthick);
  printf("Si\n");
  printf("  Angle: %6.2f deg\n",ssim.thSi);
  printf("  Distance: %6.1f mm\n",ssim.dSi);
  printf("  Resolution: %6.3f MeV\n",ssim.sigK4);
  for(int i=0;i<ssim.nstat;i++){
    printf("Excitation energy (%d): %7.3f MeV\n",i,ssim.ex_ene[i]);
  }
  if(flag_ang==1) printf("Angular distribution flag: true\n");
  if(flag_ang==0) printf("Angular distribution flag: false\n");
  
  
  std::string check;
  std::cout << "\ncontinue?(y/n) ";
  std::cin >> check;
  if(check!="Y" && check!="y") return 1;

  fprintf(outputfile2,"Time: %s",ctime(&tnow));
  fprintf(outputfile2,"Output: %s\n",argv[1]);
  fprintf(outputfile2,"Number of States: %d\n",ssim.nstat);
  fprintf(outputfile2,"Mass Number: %d\n",(int)ssim.a1);
  fprintf(outputfile2,"Incident energy: %7.1f MeV (%5.1f MeV/u)\n", ssim.enepa,ssim.enepa/ssim.a1);
  fprintf(outputfile2,"  beam position: x=%6.3f mm, y=%6.3f mm\n",ssim.offx1, ssim.offy1);
  fprintf(outputfile2,"  Energy spread: %6.3f MeV\n",ssim.sigK1);
  fprintf(outputfile2,"  Angular spread: %6.3f rad\n",ssim.sigthr1);
  fprintf(outputfile2,"  Lateral spread: %6.3f mm\n",ssim.sigx1);
  fprintf(outputfile2,"Target thickness: %5.1f mg/cm2\n",ssim.tthick);
  fprintf(outputfile2,"Si\n");
  fprintf(outputfile2,"  Angle: %6.2f deg\n",ssim.thSi);
  fprintf(outputfile2,"  Distance: %6.1f mm\n",ssim.dSi);
  fprintf(outputfile2,"  Resolution: %6.3f MeV\n",ssim.sigK4);
  if(flag_ang==1) printf("Angular distribution flag: true\n");
  if(flag_ang==0) printf("Angular distribution flag: false\n");
  for(int i=0;i<ssim.nstat;i++){
    fprintf(outputfile2,"Excitation energy (%d): %7.3f MeV\n",i,ssim.ex_ene[i]);
  }

  for(int i=0; i<N_ALPHA; i++){
    mgs[i]=(i+1.0)*4.0*AMU+mexes[i];
    th_ene[i]=mgs[0]*(i+1.0)-mgs[i];
  };

  /***** input kinematics ***********/
  kine.K1=ssim.enepa;
  kine.m1=mgs[0];
  kine.m2=mgs[2];
  kine.m3=mgs[0];
  kine.m4=mgs[2]+ssim.ex_ene[0];
  ssim.mvec1[0]=mgs[0]+ssim.enepa;
  ssim.mvec1[1]=0;
  ssim.mvec1[2]=0;
  ssim.mvec1[3]=sqrt(ssim.mvec1[0]*ssim.mvec1[0]-mgs[0]*mgs[0]);
  ssim.mvec2[0]=mgs[2];
  ssim.mvec2[1]=0;
  ssim.mvec2[2]=0;
  ssim.mvec2[3]=0;
  /***** input kinematics ***********/


  
  /*data output*/  
  TFile *fout =new TFile(Form("rootfile/%s.root",argv[1]),"recreate");
  TTree *tree = new TTree("tree","tree");
  
  //section 1
  tree->Branch("sK1",ssim.sK1,"1_sK1[2]/D");
  tree->Branch("sx1",&ssim.sx1,"1_sx1/D"); 
  tree->Branch("sy1",&ssim.sy1,"1_sy1/D"); 
  tree->Branch("sdep1",&ssim.sdep1,"1_sdep1/D"); 

  //section 2
  tree->Branch("sK3",ssim.sK3,"2_sK3[2]/D");
  tree->Branch("sth3",ssim.sth3,"2_sth3[2]/D");
  tree->Branch("sph3",ssim.sph3,"2_sph3[2]/D");
  tree->Branch("sK4",ssim.sK4,"2_sK4[2]/D");
  tree->Branch("sth4",ssim.sth4,"2_sth4[2]/D");
  tree->Branch("sph4",ssim.sph4,"2_sph4[2]/D");
  
  //section 4
  tree->Branch("idp4He",&ssim.idp4He,"4_idp4He/I");
  tree->Branch("mK4He",&ssim.mK4He,"4_mK4He/D");
  tree->Branch("mth4He",&ssim.mth4He,"4_mth4He/D");
  tree->Branch("mph4He",&ssim.mph4He,"4_mph4He/D");
  tree->Branch("s4Hepos",ssim.s4Hepos,"4_s4Hepos[3]/D");
  tree->Branch("m4Hepos",ssim.m4Hepos[1],"4_m4Hepos[3]/D");
  tree->Branch("m4Heseg",&ssim.m4Heseg,"4_m4Heseg/I");
  tree->Branch("m4HechF",&ssim.isx4HeSi,"4_m4HechF/I");
  tree->Branch("m4HechR",&ssim.isy4HeSi,"4_m4HechR/I");
  tree->Branch("mex12C",&ssim.mexC,"4_mex12C/D");

  // section 5
  tree->Branch("idp3a",ssim.idp3a,"5_idp3a[3]/I");
  tree->Branch("sK3a",ssim.sK3a,"5_sK3a[3]/D");
  tree->Branch("mK3a",ssim.mK3a,"5_mK3a[3]/D");
  tree->Branch("mth3a",ssim.mth3a,"5_mth3a[3]/D");
  tree->Branch("mph3a",ssim.mph3a,"5_mph3a[3]/D");
  tree->Branch("s3apos",ssim.s3apos,"5_s3apos[3][3]/D");
  tree->Branch("m3apos",ssim.m3apos[1],"5_m3apos[2][3][3]/D");
  tree->Branch("m3aseg",ssim.m3aseg,"5_m3aseg[3]/I");
  tree->Branch("m3achF",ssim.isx3aSi,"5_m3achF[3]/I");
  tree->Branch("m3achR",ssim.isy3aSi,"5_m3achR[3]/I");
  tree->Branch("sex3a",&ssim.sex3a,"5_sex3a/D");
  tree->Branch("mex3a",&ssim.mex3a,"5_mex3a/D");

  //section 6
  tree->Branch("sK12C",&ssim.sKC[MAX_STAT-1],"6_sK12C/D");
  tree->Branch("sth12C",&ssim.sthC[MAX_STAT-1],"6_sth12C/D");
  tree->Branch("sph12C",&ssim.sphC[MAX_STAT-1],"6_sph12C/D");
  
  //section 7
  tree->Branch("idp12C",&ssim.idp12C,"7_idp12C/I");
  tree->Branch("mK12C",&ssim.mK12C,"7_mK12C/D");
  tree->Branch("mth12C",&ssim.mth12C,"7_mth12C/D");
  tree->Branch("mph12C",&ssim.mph12C,"7_mph12C/D");
  tree->Branch("s12Cpos",ssim.s12Cpos,"7_s12Cpos[3]/D");
  tree->Branch("m12Cpos",ssim.m12Cpos[1],"7_m12Cpos[3]/D");
  tree->Branch("m12Cseg",&ssim.m12Cseg,"7_m12Cseg/I");
  tree->Branch("m12CchF",&ssim.isx12CSi,"7_m12CchF/I");
  tree->Branch("m12CchR",&ssim.isy12CSi,"7_m12CchR/I");


  
  int nev=100000;
  //  int nev=6e6;

  // angular distribution
  double thetacm_i = 40.; //*_i < thetacm < *_f
  double thetacm_f = 75.;
  const int n_div = 350; 
  double tmp_array[n_div];
  double theta_array[n_div];
  double total=0.;
  for(int i=0; i<n_div; i++){
    theta_array[i] = 0;
    tmp_array[i] = tmp(thetacm_i,thetacm_f,n_div,i);
    total += tmp_array[i]; 
  }
  int loop = (int)((double)nev*(thetacm_f-thetacm_i)/180/total);
  cout << loop * total << endl;

  
  srand48(time(NULL));  
  for(int iev=0;iev<nev*2;iev++){

    if(iev%10000==0) {
      std::cout << iev << " event finished" << endl;
      std::cout << std::flush; 
    }
    ssim.sx1=0; ssim.sy1=0;
    ssim.sdep1=0;
    ssim.sthrcm=0; ssim.sphrcm=0;
    ssim.sthrstra3=0; ssim.sphrstra4=0;
    ssim.astra3=0; ssim.astra4=0;
    for(int l=0;l<2;l++){
      ssim.sK1[l]=0; ssim.sthr1[l]=0; ssim.sphr1[l]=0; ssim.sth1[l]=0; ssim.sph1[l]=0;
      ssim.sthxr1[l]=0; ssim.sthyr1[l]=0;
      ssim.sK3[l]=0; ssim.sthr3[l]=0; ssim.sphr3[l]=0; ssim.sth3[l]=0; ssim.sph3[l]=0;
      ssim.sK4[l]=0; ssim.sthr4[l]=0; ssim.sphr4[l]=0; ssim.sth4[l]=0; ssim.sph4[l]=0;
    }
    for(int l=0;l<3;l++){
      ssim.sKa[l]=0; ssim.sthra[l]=0; ssim.sphra[l]=0; ssim.stha[l]=0; ssim.spha[l]=0;
      ssim.sKC[l]=0; ssim.sthrC[l]=0; ssim.sphrC[l]=0; ssim.sthC[l]=0; ssim.sphC[l]=0; ssim.spC[l]=0; 
      for(int m=0;m<4;m++){
	ssim.savec[l][m]=0;
	ssim.sbvec[m]=0;
	ssim.scvec[l][m]=0;
	ssim.sgvec[l][m]=0;
      }
    }
    ssim.idp4He=0;
    ssim.s4Hevec[0][0]=0; ssim.s4Hevec[0][1]=0; ssim.s4Hevec[0][2]=0; 
    ssim.s4Hevec[1][0]=0; ssim.s4Hevec[1][1]=0; ssim.s4Hevec[1][2]=0; 
    ssim.s4Hepos[0]=-100; ssim.s4Hepos[1]=-100; ssim.s4Hepos[2]=-100; 
    ssim.mK4He=0;
    ssim.m4Hepos[0][0]=-100; ssim.m4Hepos[0][1]=-100; ssim.m4Hepos[0][2]=-100; 
    ssim.m4Hepos[1][0]=-100; ssim.m4Hepos[1][1]=-100; ssim.m4Hepos[1][2]=-100; 
    ssim.mthr4He=0; ssim.mphr4He=0; ssim.mth4He=0; ssim.mph4He=0;
    ssim.mthxr4He=0; ssim.mthyr4He=0; ssim.mthx4He=0; ssim.mthy4He=0;
    ssim.mexC=0;
    ssim.m4Heseg=-1;
    ssim.isx4HeSi=-1; ssim.isy4HeSi=-1;
    ssim.mvec4He[0]=0; ssim.mvec4He[1]=0; ssim.mvec4He[2]=0; ssim.mvec4He[3]=0; 
    ssim.idp12C=0;
    ssim.s12Cvec[0][0]=0; ssim.s12Cvec[0][1]=0; ssim.s12Cvec[0][2]=0; 
    ssim.s12Cvec[1][0]=0; ssim.s12Cvec[1][1]=0; ssim.s12Cvec[1][2]=0; 
    ssim.s12Cpos[0]=-100; ssim.s12Cpos[1]=-100; ssim.s12Cpos[2]=-100; 
    ssim.mK12C=0;
    ssim.m12Cpos[0][0]=-100; ssim.m12Cpos[0][1]=-100; ssim.m12Cpos[0][2]=-100; 
    ssim.m12Cpos[1][0]=-100; ssim.m12Cpos[1][1]=-100; ssim.m12Cpos[1][2]=-100; 
    ssim.mthr12C=0; ssim.mphr12C=0; ssim.mth12C=0; ssim.mph12C=0;
    ssim.mthxr12C=0; ssim.mthyr12C=0; ssim.mthx12C=0; ssim.mthy12C=0;
    ssim.m12Cseg=-1;
    ssim.isx12CSi=-1; ssim.isy12CSi=-1;
    ssim.mvec12C[0]=0; ssim.mvec12C[1]=0; ssim.mvec12C[2]=0; ssim.mvec12C[3]=0; 
    for(int l=0;l<2;l++){
      for(int m=0;m<3;m++){
	ssim.idp3a[m]=0;
	ssim.mK3a[m]=0;	ssim.mthr3a[m]=0; ssim.mphr3a[m]=0; ssim.mth3a[m]=0; ssim.mph3a[m]=0;
	ssim.mthxr3a[m]=0; ssim.mthyr3a[m]=0; ssim.mthx3a[m]=0; ssim.mthy3a[m]=0;
	ssim.mex3a=0; ssim.sex3a=0;
	ssim.isx3aSi[m]=-1; ssim.isy3aSi[m]=-1;
	ssim.m3aseg[m]=-1;
	for(int n=0;n<3;n++){
	  ssim.s3avec[l][m][n]=0; ssim.s3apos[m][n]=-100; ssim.m3apos[l][m][n]=-100;
	}
	for(int n=0;n<4;n++){
	  ssim.svec3a[m][n]=0;	  
	  ssim.mvec3a[m][n]=0;	  
	}
      }
    }
    
    /** Create 4He Beam ***********************************************/
    create4He(&kine,&ssim); /* Up to the collision point */

    /** Scattering ***********************************************/
    scatt4He12c(&kine,&ssim); /* 4He + 12C Scattering */

    if(flag_ang==true){
      int array_n = (ssim.sthrcm*R_TO_D-thetacm_i)/((thetacm_f-thetacm_i)/n_div);
      //    cout << ssim.sthrcm*R_TO_D << " " << array_n << endl;
      if(array_n<0 || array_n-1>n_div) continue;
      theta_array[array_n]+=1;
      if(theta_array[array_n] > tmp_array[array_n]*loop) continue;
    }
    
    /** Decay 12C->3a ***********************************************/
    decay3a(&kine,&ssim);
    
    /** 4He detection *******************************************/
    detect4He(&kine,&ssim);

    /** 3a detection *******************************************/
    detect3a(&kine,&ssim);

    /** 12C gamma decay *******************************************/
    decay12C(&kine,&ssim);

    /** 12C detection *******************************************/
    detect12C(&kine,&ssim);

    tree->Fill();
  }
  tree->AutoSave();
  fout->Close();
  
  return(0);
}





/******************************************************/
/* double create4He(struct skine *k, struct ssimval *s)*/
/* Create initial 4He beam just before scattering in the target */
/* Input parameters:                                    */
/* s->tthick ... Target thickness in the unit of ug/cm^2  */
/* s->enepa ... Beam energy in MeV                            */
/* s->sigK1 ... Energy spread of Beam in MeV               */ 
/* s->sigthr1 ... In-plane angular spread of Beam in radian    */
/* s->sigx1 ... Lateral spread of beam spot in mm */
/* Output parameters:                                      */
/* s->sdep1 ... Depth where a collision occurs in the target */
/* s->sK1[0]  ... Beam energy before the target */
/* s->sK1[1]  ... Beam energy before the collision in the target */
/* s->thxr1[0] ... Trajectory angle of beam in hori. plane before tareget */
/* s->thyr1[0] ... Trajectory angle of beam in vert. plane before tareget */
/* s->thxr1[1] ... Trajectory angle of beam in hori. plane before collision */
/* s->thyr1[1] ... Trajectory angle of beam in vert. plane before collision */
/* s->sx1[1] ... Postion of beam in horizontal plane before collision */
/* s->sy1[1] ... Position of beam in vertical plane before collision */
/* s->thr1[0] ... Polar angle of beam before target */
/* s->thr1[1] ... Polar angle of beam before collision */
/* s->phr1[0] ... Azimuthal angle of beam before target */
/* s->phr1[1] ... Azimuthal angle of beam before collision */
/* Return value: kinetic energy of 4He beam just before the scattering */

double create4He(struct skine *k, struct ssimval *s){
  int i;
  double tmpr1,tmpr2;

  s->sdep1=drand48()*s->tthick;  /* Collision depth */

  genrndg(&tmpr1,&tmpr2);
  s->sK1[0]=s->enepa+tmpr1*s->sigK1;  /* Beam energy */
  s->sK1[1]=s->sK1[0]-eloss4He(s->sK1[0],s->sdep1,s->tthick)+tmpr2*estra4He(s->sK1[0],s->sdep1,s->tthick) ;
  //  s->sK1[1]=s->sK1[0]-eloss4He(s->sK1[0],s->sdep1,s->tthick)+tmpr2*estra4He(s->sK1[0]-eloss4He(s->sK1[0],s->sdep1,s->tthick),s->sdep1,s->tthick) ;
  /* energy loss in target before scattering */

  genrndg(&tmpr1,&tmpr2);  /* Angular spread of beam */
  s->sthxr1[0]=tmpr1*s->sigthr1;
  s->sthyr1[0]=tmpr2*s->sigthr1;
  genrndg(&tmpr1,&tmpr2);  /* Angular straggling in the taret */
  s->sthxr1[1]=s->sthxr1[0]+tmpr1*astrac(s->sK1[0],s->sdep1,s->tthick);
  s->sthyr1[1]=s->sthyr1[0]+tmpr2*astrac(s->sK1[0],s->sdep1,s->tthick);

  genrndg(&tmpr1,&tmpr2);  /* Lateral spread of beam */
  s->sx1=tmpr1*s->sigx1+s->offx1;
  s->sy1=tmpr2*s->sigx1+s->offy1;

  /* Calculate polar and azimuthal angles */
  for(i=0;i<2;i++){
    double vtra,tanx,tany;
    tanx=tan(s->sthxr1[i]);
    tany=tan(s->sthyr1[i]);
    s->sphr1[i]=atan2(tany,tanx);
    s->sph1[i]=s->sphr1[i]*R_TO_D;
    vtra=sqrt(tanx*tanx+tany*tany+1.0);
    s->sthr1[i]=acos(1.0/vtra);
    s->sth1[i]=s->sthr1[i]*R_TO_D;
  }
  return(s->sK1[1]);
}



/*******************************************************/
/* double scatt4He12c(struct skine *k, struct ssimval *s) */
/* Scattering 4He + 12C and create 4He and 12C vector       */
/* Input parameters:                                    */
/* s->tthick ... Target thickness in the unit of mg/cm^2  */
/* s->sdep1 ... Depth where a collision occurs in the target */
/* s->sK1[1]  ... Beam energy before the collision in the target */
/* s->sthr1[1] ... Polar angle of beam before collision */
/* s->sphr1[1] ... Azimuthal angle of beam before collision */
/* Output parameters:                                  */
/* s->sthrcm ... Polar scattering angle in CM */
/* s->sphrcm ... Azimuthal scattering angle in CM */
/* s->sK3[0]  ... Energy of particle 3 after the collision point */
/* s->sK4[0]  ... Energy of particle 4 after the collision point */
/* s->sthr3[0] ... Polar angle of 3 after the collision point */
/* s->sphr3[0] ... Azimuthal angle of 3 after the collision point */
/* s->sthr4[0] ... Polar angle of 4 after the collision point */
/* s->sphr4[0] ... Azimuthal angle of 4 after the collision point */
/* s->astra3 ... In-plane ang. stragg. of particle 3 in target */
/* s->astra4 ... In-plane ang. stragg. of particle 4 in target */
/* s->sthrstra3 ... Polar ang. stragg. of particle 3 in target */
/* s->sphrstra3 ... Azimuthal ang. stragg. of particle 3 in target */
/* s->sthrstra4 ... Polar ang. stragg. of particle 4 in target */
/* s->sphrstra4 ... Azimuthal ang. stragg. of particle 4 in target */
/* s->sK3[1]  ... Energy of particle 3 after the target */
/* s->sK4[1]  ... Energy of particle 4 after the target */
/* s->sthr3[1] ... Polar scatt. ang. of 3 after the target */
/* s->sphr3[1] ... Azimuthal scatt. ang. of 3 after the target */
/* s->sth3[1] ... Polar scatt. ang. of 3 after the target in deg. */
/* s->sph3[1] ... Azimuthal scatt. ang. of 3 after the target in deg. */
/* s->sthr4[1] ... Polar scatt. ang. of 4 after the target */
/* s->sphr4[1] ... Azimuthal scatt. ang. of 4 after the target */
/* s->sth4[1] ... Polar scatt. ang. of 4 after the target in deg. */
/* s->sph4[1] ... Azimuthal scatt. ang. of 4 after the target in deg. */
/* Return value: kinetic energy of 4He beam after the target */

double scatt4He12c(struct skine *k, struct ssimval *s){
  int i;
  double v1[3],v2[3];
  double tmpr1,tmpr2,sdep2;
    
  s->sthrcm=acos(2.0*drand48()-1.0);//%% 0-180 deg
  //  s->sthrcm=1.6;
  s->sphrcm=drand48()*2.0*M_PI; //0-360 deg
  //  s->sphrcm=1.5; //0-360 deg
  k->K1=s->sK1[1];
  k->curang=s->sthrcm*R_TO_D;
  i=calcscatt(k,1);
  if(i!=1) printf("Number of particle(s) %d\n",i); /* calculation in CM */
  s->sK3[0]=k->K3[0];
  s->sK4[0]=k->K4[0];

  /* scattering of particle 3 */
  deflect2(v1,v2,s->sthr1[1],s->sphr1[1],k->labang_r,s->sphrcm);
  s->sth3[0]=gettheta(v2);
  s->sph3[0]=getphi(v2);
  s->sthr3[0]=s->sth3[0]*D_TO_R;
  s->sphr3[0]=s->sph3[0]*D_TO_R;

  /* scattering of particle 4 */
  //  deflect2(v1,v2,s->sthr1[1],s->sphr1[1],k->ang_lb4[0]*D_TO_R,s->sphrcm);
  deflect2(v1,v2,s->sthr1[1],s->sphr1[1],k->ang_lb4[0]*D_TO_R,s->sphrcm+180*D_TO_R);//%%%
  s->sth4[0]=gettheta(v2);
  s->sph4[0]=getphi(v2);
  s->sthr4[0]=s->sth4[0]*D_TO_R;
  s->sphr4[0]=s->sph4[0]*D_TO_R;

  //  printf("%f,%f,%f,%f\n",s->sth3[0],s->sph3[0],s->sth4[0],s->sph4[0]);//%%%%
  
  sdep2=s->tthick-s->sdep1;
  genrndg(&tmpr1,&tmpr2);  /* Energy loss in the target */
  //  s->sK3[1]=s->sK3[0]-elossc(s->sK3[0],sdep2,s->tthick)+tmpr1*estrac(s->sK3[0],sdep2,s->tthick) ;
  //  s->sK4[1]=s->sK4[0]-eloss4He(s->sK4[0],sdep2,s->tthick)+tmpr2*estra4He(s->sK4[0],sdep2,s->tthick) ;
  s->sK3[1]=s->sK3[0];
  s->sK4[1]=s->sK4[0];

  //if(s->sK3[1]<0)  printf("%f,%f,%f\n",s->sK3[0],elossc(s->sK3[0],sdep2,s->tthick),estrac(s->sK3[0],sdep2,s->tthick));//%%%%
  //  if(s->sK4[1]<1000)  printf("%f,%f,%f\n",s->sK4[0],eloss4He(s->sK4[0],sdep2,s->tthick),estra4He(s->sK4[0],sdep2,s->tthick));//%%%%
    
  genrndg(&tmpr1,&tmpr2);  /* Angular spread of beam */

  
  s->astra3=tmpr1*astrac(s->sK3[0],sdep2,s->tthick);
  s->sthrstra3=s->astra3*sqrt(2.0);
  s->sphrstra3=drand48()*2.0*M_PI; //%
  deflect2(v1,v2,s->sthr3[0],s->sphr3[0],s->sthrstra3,s->sphrstra3);
  s->sthr3[1]=s->sthr3[0];
  s->sphr3[1]=s->sphr3[0];
  s->sth3[1]= s->sthr3[1]*R_TO_D;
  s->sph3[1]= s->sphr3[1]*R_TO_D;
  
  /*
  s->astra3=tmpr1*astrac(s->sK3[0],sdep2,s->tthick);
  s->sthrstra3=s->astra3*sqrt(2.0);
  s->sphrstra3=drand48()*2.0*M_PI; //%
  deflect2(v1,v2,s->sthr3[0],s->sphr3[0],s->sthrstra3,s->sphrstra3);
  s->sth3[1]=gettheta(v2);
  s->sph3[1]=getphi(v2);
  s->sthr3[1]=s->sth3[1]*D_TO_R;
  s->sphr3[1]=s->sph3[1]*D_TO_R;
  */

  
  if(s->sK4[0]<0||s->sK4[0]>50){
    s->sK4[1]=-1000;
    s->sth4[1]=-1000;
    s->sph4[1]=-10000;
    s->sthr4[1]=-10000;
    s->sphr4[1]=-10000;
  }
  else{
    s->astra4=tmpr2*astra4He(s->sK4[0],sdep2,s->tthick);
    s->sthrstra4=s->astra4*sqrt(2.0);
    s->sphrstra4=drand48()*2.0*M_PI;
    deflect2(v1,v2,s->sthr4[0],s->sphr4[0],s->sthrstra4,s->sphrstra4);
    s->sth4[1]=gettheta(v2);
    s->sph4[1]=getphi(v2);
    s->sthr4[1]=s->sth4[1]*D_TO_R;
    s->sphr4[1]=s->sph4[1]*D_TO_R;
  }

  return(s->sK3[1]);
}



/*******************************************************/
/* double decay3a(struct skine *k, struct ssimval *s) */
/* In-flight decay of 12C into 3 alpha particles via 8Be    */
/* Input parameters:                                        */
/* k->m1 ... Rest mass of 12C */
/* s->sdep1 ... Depth where a collision occurs in the target */
/* s->tthick ... Target thickness in the unit of mg/cm^2  */
/* s->sK3vec[0] ... Four vector of particle 3 after the collision point */
/* s->sK4vec[0] ... Four vector of particle 4 after the collision point */
/* s->sthr3[0] ... Polar angle of 3 after the collision point */
/* s->sphr3[0] ... Azimuthal angle of 3 after the collision point */
/* s->sthr4[0] ... Polar angle of 4 after the collision point */
/* s->sphr4[0] ... Azimuthal angle of 4 after the collision point */
/* Output parameters:                                  */
/* Just after the collision point ******************/
/* s->savec[3][4] ... 4 momentum of alpha particles */
                /* [][0][] .. 1st alpha emitted into a + 8Be channel */
                /* [][1][] .. 2nd alpha emitted from 8Be */
                /* [][2][] .. 2nd alpha emitted from 8Be */
/* s->sKa[2][3] ... Energy of alpha particle */
/* s->sbvec[2][4] ... 4 momentum of 8Be      */
/* s->sthra[2][3] ... Polar angle of alpha particles in radian */
/* s->sphra[2][3] ... Azimuthal angle of alpha particles in radian */
/* s->stha[2][3] ... Polar angle of alpha particles in degree */
/* s->spha[2][3] ... Azimuthal ang. of alpha parcticles in degree */
/* Just after the target ******************/
/* s->safvec[3][4] ... 4 momentum of alpha particles */
                /* [][0][] .. 1st alpha emitted into a + 8Be channel */
                /* [][1][] .. 2nd alpha emitted from 8Be */
                /* [][2][] .. 2nd alpha emitted from 8Be */
/* s->sK3a[2][3] ... Energy of alpha particle */
/* s->sthraf[2][3] ... Polar angle of alpha particles in radian */
/* s->sphraf[2][3] ... Azimuthal angle of alpha particles in radian */
/* s->sthaf[2][3] ... Polar angle of alpha particles in degree */
/* s->sphaf[2][3] ... Azimuthal ang. of alpha parcticles in degree */
/* s->astraf[2][3] ... In-plane ang. stragg. of alpha particle in target */
/* s->sthrstraf[2][3] ... Polar ang. stragg. of alpha particle in target */
/* s->sphrstraf[2][3] ... Azim. ang. stragg. of alpha particle in target */
/* s->stotea   ... Sum of alpha energies in the target */
/* s->stoteaf  ... Sum of alpha energies after the target */
/* s->sna  ... Number of alpha particles which escape from the target */
/* Return value: Sum of kinetic energies of 6 alpha particles */
/*               after the target */


double decay3a(struct skine *k, struct ssimval *s){
  int i,j,l;
  double betag[3];
  double m0p;
  double m8be=mgs[1];
  double m4he=mgs[0];
  double c12vec[4];
  double p8bec[4],p4hec[4],p4he1c[4],p4he2c[4];
  double sdep2,spath;
  double tmpr1,tmpr2;
  double tmpmom,tmpeloss;
  double v1[3],v2[3];
  

  c12vec[0]=s->sK4[1]+mgs[2]+s->ex_ene[0];
  tmpmom=sqrt(c12vec[0]*c12vec[0]-(mgs[2]+s->ex_ene[0])*(mgs[2]+s->ex_ene[0]));
  
  deflect2(v1,v2,s->sthr1[1],s->sphr1[1],k->ang_lb4[0]*D_TO_R,s->sphrcm+180*D_TO_R);
  for(i=0;i<3;i++) c12vec[i+1]=v2[i]*tmpmom; 
  
  sdep2=s->tthick-s->sdep1;

  s->stotea=0;
  s->stoteaf=0;
  s->sna=0;

  for(i=0;i<4;i++) {
    p8bec[i]=0.; 
    p4hec[i]=0.; 
  }
  
  /* 8Be + 4He  */
  getbeta(c12vec,betag);
  vecadd(betag,betag,betag,0.0,-1.0);
  m0p=sqrt(scapro4(c12vec,c12vec));

  
  if(decay2body(p8bec,p4hec,m0p,m8be,m4he)<0) {
    printf("Cannot decay to 8Be + 4He channel !!\n");
    return(0);
  }
  lortra(p8bec,s->sbvec,betag);
  lortra(p4hec,s->savec[0],betag);
  
  /* 4He + 4He  */
  getbeta(s->sbvec,betag);
  vecadd(betag,betag,betag,0.0,-1.0);
  if(decay2body(p4he1c,p4he2c,m8be,m4he,m4he)<0) {
    printf("Cannot decay to 4He + 4He channel !!\n");
    return(0);
  }
  lortra(p4he1c,s->savec[1],betag);
  lortra(p4he2c,s->savec[2],betag);
  
  
  // Energy loss in the target
  for(j=0;j<3;j++){
    s->sKa[j]=s->savec[j][0]-m4he;
    s->stotea+=s->sKa[j];
    s->stha[j]=gettheta4(s->savec[j]);
    s->spha[j]=getphi4(s->savec[j]);
    s->sthra[j]=s->stha[j]*D_TO_R;
    s->sphra[j]=s->spha[j]*D_TO_R;

    genrndg(&tmpr1,&tmpr2);
    if(s->savec[j][3]>0) spath=sdep2/cos(s->sthra[j]);
    else {
      //      printf("Alpha emitted in backward direction !!!\n");
      spath=fabs(s->sdep1/cos(s->sthra[j]));
    }
    // Anguar straggling 
    s->safvec[j][0]=-1000;
    s->safvec[j][1]=-1000;
    s->safvec[j][2]=-1000;
    s->safvec[j][3]=-1000;
    
    s->astraf[j]=tmpr1*astra4He(s->sKa[j],spath,s->tthick);
    //    s->astraf[j]=0;
    s->sthrstraf[j]=s->astraf[j]*sqrt(2.0);
    s->sphrstraf[j]=drand48()*2.0*M_PI;
    deflect2(v1,v2,s->sthra[j],s->sphra[j],
	     s->sthrstraf[j],s->sphrstraf[j]);
    s->sthaf[j]=gettheta(v2);
    s->sphaf[j]=getphi(v2);
    s->sthraf[j]=s->sthaf[j]*D_TO_R;
    s->sphraf[j]=s->sphaf[j]*D_TO_R;
    
    tmpeloss=eloss4He(s->sKa[j],spath,s->tthick);
    
    // Energy straggling 
    s->sK3a[j]=s->sKa[j]-tmpeloss+tmpr2*estra4He(s->sKa[j]-tmpeloss,spath,s->tthick);
    //    s->sK3a[j]=s->sKa[j];
    if(s->sK3a[j]<0){
      s->sK3a[j]=-10000;
      s->sthaf[j]=-10000;
      s->sphaf[j]=-10000;
      s->sthraf[j]=-10000;
      s->sphraf[j]=-10000;
    }else{
      s->sna++;
      s->stoteaf+=s->sK3a[j];
      
      // Reconstruct 4 momentum 
      s->safvec[j][0]=s->sK3a[j]+m4he;
      tmpmom=sqrt(s->safvec[j][0]*s->safvec[j][0]-m4he*m4he);
      for(l=0;l<3;l++) s->safvec[j][l+1]=v2[l]*tmpmom; 

    }
  }
  
  return(s->stoteaf);
}




/*******************************************************/
/* int detect4He(struct skine *k, struct ssimval *s) */
/* Proton detection by Si-CsI                        */
/* Input parameters:                                    */
/* k->m3 ... Rest mass of particle 3 (proton) */
/* s->nstat ... Number of states in 12C */
/* k->m1 ... Rest mass of 12C in the ground state */
/* s->sK3[1]  ... Energy of particle 3 after the target */
/* s->sthr3[1] ... Polar scatt. ang. of 3 after the target */
/* s->sphr3[1] ... Azimuthal scatt. ang. of 3 after the target */
/* s->sth3[1] ... Polar scatt. ang. of 3 after the target in deg. */
/* s->sph3[1] ... Azimuthal scatt. ang. of 3 after the target in deg. */
/* s->sx1 ... Lateral spread of 12C Beam in horizontal plane */
/* s->sy1 ... Lateral spread of 12C Beam in vertical plane */

/* s->dSi ... distance between target and Si in mm */
/* s->vSi ... Vertical size of Si in mm */
/* s->hSi ... Horizontal size of Si in mm */
/* s->nxSi ... Number of horizontal strips in Si */
/* s->nySi ...  Number of vertical strips in Si */
/* s->thSi ...  Angle of Si in degree */
/* s->thrSi ... Angle of Si in radian */
/* s->sigK3 ... Energy resolution of recoiled proton in MeV */

/* Output parameters:                                  */
/* s->idp4He ... Detection flag for particle 3 ..Not detected  1..detected */
/* s->s4Hevec[2][3] ... 3-dim direction of 4He trajectory */
                      /* [0][] ... Detector frame */
                      /* [1][] ... Laboratory frame */
/* s->s4Hepos[3] ... Position in the detector frame */
/* s->mK4He ... Measured energy of proton */
/* s->mppos[2][3] ... Measured position of proton */
                      /* [0][] ... Detector frame */
                      /* [1][] ... Laboratory frame */
/* s->mthrp ... Measured polar angle of proton in radian */
/* s->mphrp ... Measured azimuthal angle of proton in radian */
/* s->mthp ... Measured polar angle of proton in degree */
/* s->mphp ... Measured azimuthal angle of proton in degree */
/* s->mthxrp ... Measured horizontal angle of proton in radian */
/* s->mthyrp ... Measured vertical angle of proton in radian */
/* s->mthxp ... Measured horizontal angle of proton in degree */
/* s->mthyp ... Measured vertical angle of proton in degree */
/* s->mexC ... Measured excitation energy of 12C */
/* s->mvec4He[4] ... Measured four vector of recoiled proton */
/* Return value: Detection flag for particle 3 (s->idp4He) */

int detect4He(struct skine *k, struct ssimval *s){
  double slatex[2][3];
  double null[3]={0,0,0};
  double tmpr1=0,tmpr2=0,tmpx=0,tmpy=0,tmpz=0;
  double tmpr=0;
  double tmptheta=0;
  int quadrant=-1;
  int i,j;
  genvthph(s->sth3[1],s->sph3[1],s->s4Hevec[1]);
  rotvec(s->s4Hevec[1],s->s4Hevec[0],1,-1.0*s->thrSi);

  
  genrndg(&tmpr1,&tmpr2);
  //  s->mK4He=s->sK3[1]+tmpr1*s->sigK3-eloss4He(s->sK3[1]+tmpr1*s->sigK3,s->tdeadSi,s->tthick);
  s->mK4He=s->sK3[1]+tmpr1*s->sigK3-eloss4He(s->sK3[1],s->tdeadSi,s->tthick);
  if(s->s4Hevec[0][2]>0){
    slatex[1][0]=s->sx1;
    slatex[1][1]=s->sy1;
    slatex[1][2]=0;
    rotvec(slatex[1],slatex[0],1,-1.0*s->thrSi);
    vecadd(s->s4Hevec[0],slatex[0],s->s4Hepos,
           (s->dSi-slatex[0][2])/s->s4Hevec[0][2],1.0);    
    tmpr = pow((pow(s->s4Hepos[0],2)+pow(s->s4Hepos[1],2)),0.5);
    if(s->s4Hepos[0]>0) tmptheta=atan(s->s4Hepos[1]/s->s4Hepos[0])*R_TO_D;
    if(s->s4Hepos[0]<0) tmptheta=atan(s->s4Hepos[1]/s->s4Hepos[0])*R_TO_D+180;
    //    printf("4He:x,y,z,r:%f,%f,%f,%f\n",s->s4Hepos[0],s->s4Hepos[1],s->s4Hepos[2],tmpr);

    if((tmpr < s->diaSi/2) && (tmpr > s->innSi/2)){
      s->idp4He=1;
      if(s->s4Hepos[0]>0 && s->s4Hepos[1]>0){
        quadrant=0;
      }else if(s->s4Hepos[0]<0 && s->s4Hepos[1]>0){
        quadrant=1;
      }else if(s->s4Hepos[0]<0 && s->s4Hepos[1]<0){
        quadrant=2;
      }else if(s->s4Hepos[0]>0 && s->s4Hepos[1]<0){
        quadrant=3;
      }else{
        quadrant=-100;
        printf("detection process error");
      }
      s->m4Heseg = quadrant;
      //      printf("quadrant:%d\n",quadrant);

      s->isx4HeSi=(int)((tmpr-s->innSi/2)/((s->diaSi-s->innSi)/2)*(s->nxSi/4)) + 16*quadrant;//%% 
      if(quadrant==0) s->isy4HeSi = (int)(((atan(s->s4Hepos[1]/s->s4Hepos[0])*R_TO_D)/22.5)) +4*0;
      if(quadrant==1) s->isy4HeSi = (int)(((atan(s->s4Hepos[1]/s->s4Hepos[0])*R_TO_D)/22.5)-1)+4 +4*1;
      if(quadrant==2) s->isy4HeSi = (int)(((atan(s->s4Hepos[1]/s->s4Hepos[0])*R_TO_D)/22.5)) +4*2;
      if(quadrant==3) s->isy4HeSi = (int)(((atan(s->s4Hepos[1]/s->s4Hepos[0])*R_TO_D)/22.5)-1)+4 +4*3;
      //printf("x,y:%d,%d\n",s->isx4HeSi,s->isy4HeSi);
      
      s->m4Hepos[0][2]=s->s4Hepos[2];
      s->m4Hepos[0][0]=(((s->isx4HeSi%16+0.5)*((s->diaSi-s->innSi)/2)/(s->nxSi/4))+s->innSi/2) * cos(((s->isy4HeSi+0.5)*22.5)*D_TO_R);
      s->m4Hepos[0][1]=(((s->isx4HeSi%16+0.5)*((s->diaSi-s->innSi)/2)/(s->nxSi/4))+s->innSi/2) * sin(((s->isy4HeSi+0.5)*22.5)*D_TO_R);
      rotvec(s->m4Hepos[0],s->m4Hepos[1],1,s->thrSi);
      //      printf("posx, posy, posz : %f, %f %f\n",s->m4Hepos[0][0],s->m4Hepos[0][1],s->m4Hepos[0][2]);

      double tmppos[3]; 
      tmppos[0]=s->m4Hepos[1][0]-s->offx1;
      tmppos[1]=s->m4Hepos[1][1]-s->offy1;
      tmppos[2]=s->m4Hepos[1][2]-0;
     
      s->mth4He=gettheta(tmppos);
      s->mph4He=getphi(tmppos);
      s->mthr4He=s->mth4He*D_TO_R;
      s->mphr4He=s->mph4He*D_TO_R;
      //      printf("pos,theta:%f,%f\n",s->m4Hepos[1],s->mth4He);
      
      tmpx=sin(s->mthr4He)*cos(s->mphr4He);
      tmpy=sin(s->mthr4He)*sin(s->mphr4He);
      tmpz=cos(s->mthr4He);
      s->mthxr4He=atan(tmpx/tmpz);
      s->mthyr4He=atan(tmpy/tmpz);
      s->mthx4He=s->mthxr4He*R_TO_D;
      s->mthy4He=s->mthyr4He*R_TO_D;
      s->mexC=calcex4(s->mthr4He,k->m1,k->m2,k->m1,k->m2,s->enepa,s->mK4He);
      //      s->mexC=calcex4(s->mthr4He,mgs[0],mgs[2],mgs[0],mgs[2],s->enepa,s->mK4He);
      //      printf("%f %f\n",s->mK4He,s->sK3[1]);
      
      s->mvec4He[0]=k->m3+s->mK4He;
      //  /*Opposite direction from 12C*/
      genvthph(-1.0*s->mth4He,s->mph4He,&s->mvec4He[1]);
      vecadd(&s->mvec4He[1],null,&s->mvec4He[1],
	     sqrt(s->mvec4He[0]*s->mvec4He[0]-k->m3*k->m3),1);
      
    }else{
      s->idp4He=0;
    }
  }
  //  printf("idp4He: %f\n",s->idp4He);
  return(s->idp4He);
};



/*******************************************************/
/* double decay12C(struct skine *k, struct ssimval *s) */
/* In-flight decay of 12C                               */
/* Input parameters:                                    */
/* k->m2 ... Rest mass of 12C */
/* k->m4 ... Rest mass of 12C in the ex_ene[0]  state */
/* s->sK4[1]  ... Energy of particle 4 after the target */
/* s->sth4[1] ... Polar scatt. ang. of 4 after the target */	 
/* s->sph4[1] ... Azimuthal scatt. ang. of 4 after the target */
/* s->sthr4[1] ... Polar scatt. ang. of 4 after the target */	 
/* s->sphr4[1] ... Azimuthal scatt. ang. of 4 after the target */
/* s->nstat ... Number of states in 12C */
/* s->ex_ene[i] ... Excitation energy of the i th state */

/* Output parameters:                                  */
/* s->scvec[i][4] ... Four momentum of 12C in the i th state */
/* s->sgvec[i][4] ... Four momentum of gamma ray emitted from the i th state */
/* s->sKC[i] ... Energy of 12C in the i th state */
/* s->spC[i] ... Momentum of 12C in the i th state */
/* s->sthrC[i] ... Polar angle of 12C in the i th state */
/* s->sphrC[i] ... Azimuthal angle of 12C in the i th state */
/* s->sthC[i] ... Polar angle of 12C in the i th state in degree */
/* s->sphC[i] ... Azimuthal angle of 12C in the i th state in degree */
/* Return value: kinetic energy of 12C in the ground state */
double decay12C(struct skine *k, struct ssimval *s){
  int i;
  double vdummy[3];
  double nullv[4]={0,0,0,0};
  s->sKC[0]=s->sK4[1];
  s->sthrC[0]=s->sthr4[1];
  s->sphrC[0]=s->sphr4[1];
  s->sthC[0]=s->sth4[1];
  s->sphC[0]=s->sph4[1];


  /* Create four memntum of the scattered particle */
  s->scvec[0][0]=k->m4+s->sKC[0];
  s->spC[0]=sqrt(s->scvec[0][0]*s->scvec[0][0]-k->m4*k->m4);
  genvthph(s->sthC[0],s->sphC[0],vdummy);
  vecadd(vdummy,nullv,&s->scvec[0][1],s->spC[0],1.0);

  for(i=1;i<s->nstat;i++){
    double segc; /* energy of emitted gamma ray in CM */
    segc=s->ex_ene[i-1]-s->ex_ene[i];
    gengamma(segc,s->scvec[i-1],s->sgvec[i-1],s->scvec[i]);
    s->sthC[i]=gettheta4(s->scvec[i]);
    s->sphC[i]=getphi4(s->scvec[i]);
    s->sthrC[i]=s->sthC[i]*D_TO_R;
    s->sphrC[i]=s->sphC[i]*D_TO_R;
    s->sKC[i]=s->scvec[i][0]-sqrt(scapro4(s->scvec[i],s->scvec[i]));
    s->spC[i]=sqrt(scapro(&s->scvec[i][1],&s->scvec[i][1]));
  }

  return(s->sKC[i]);
}



int detect12C(struct skine *k, struct ssimval *s){
  double slatex[2][3];
  double null[3]={0,0,0};
  double tmpr1=0,tmpr2=0,tmpx=0,tmpy=0,tmpz=0;
  double tmpr=0;
  double tmptheta=0;
  int quadrant=-1;
  int i,j;
  genvthph(s->sthC[2],s->sphC[2],s->s12Cvec[1]);
  rotvec(s->s12Cvec[1],s->s12Cvec[0],1,-1.0*s->thrSi);

  double spath;
  spath=(s->tthick-s->sdep1)/cos(s->sthrC[2]);
  
  genrndg(&tmpr1,&tmpr2);
  s->mK12C=s->sKC[2]+tmpr1*s->sigK4-elossc(s->sKC[2]+tmpr1*s->sigK4,spath,s->tthick);
  if(s->s12Cvec[0][2]>0){
    slatex[1][0]=s->sx1;
    slatex[1][1]=s->sy1;
    slatex[1][2]=0;
    rotvec(slatex[1],slatex[0],1,-1.0*s->thrSi);
    vecadd(s->s12Cvec[0],slatex[0],s->s12Cpos,
           (s->dSi-slatex[0][2])/s->s12Cvec[0][2],1.0);    
    tmpr = pow((pow(s->s12Cpos[0],2)+pow(s->s12Cpos[1],2)),0.5);
    if(s->s12Cpos[0]>0) tmptheta=atan(s->s12Cpos[1]/s->s12Cpos[0])*R_TO_D;
    if(s->s12Cpos[0]<0) tmptheta=atan(s->s12Cpos[1]/s->s12Cpos[0])*R_TO_D+180;
    //    printf("12C:x,y,z,r:%f,%f,%f,%f\n",s->s12Cpos[0],s->s12Cpos[1],s->s12Cpos[2],tmpr);

    if((tmpr < s->diaSi/2) && (tmpr > s->innSi/2)){
      s->idp12C=1;
      if(s->s12Cpos[0]>0 && s->s12Cpos[1]>0){
        quadrant=0;
      }else if(s->s12Cpos[0]<0 && s->s12Cpos[1]>0){
        quadrant=1;
      }else if(s->s12Cpos[0]<0 && s->s12Cpos[1]<0){
        quadrant=2;
      }else if(s->s12Cpos[0]>0 && s->s12Cpos[1]<0){
        quadrant=3;
      }else{
        quadrant=-100;
        printf("detection process error");
      }
      s->m12Cseg = quadrant;

      //      printf("quadrant:%d\n",quadrant);
      
      s->isx12CSi=(int)((tmpr-s->innSi/2)/((s->diaSi-s->innSi)/2)*(s->nxSi/4)) + 16*quadrant;//%% 
      if(quadrant==0) s->isy12CSi = (int)(((atan(s->s12Cpos[1]/s->s12Cpos[0])*R_TO_D)/22.5)) +4*0;
      if(quadrant==1) s->isy12CSi = (int)(((atan(s->s12Cpos[1]/s->s12Cpos[0])*R_TO_D)/22.5)-1)+4 +4*1;
      if(quadrant==2) s->isy12CSi = (int)(((atan(s->s12Cpos[1]/s->s12Cpos[0])*R_TO_D)/22.5)) +4*2;
      if(quadrant==3) s->isy12CSi = (int)(((atan(s->s12Cpos[1]/s->s12Cpos[0])*R_TO_D)/22.5)-1)+4 +4*3;
      //	printf("x,y:%d,%d\n",s->isx12CSi,s->isy12CSi);
     
	  
      s->m12Cpos[0][2]=s->s12Cpos[2];
      s->m12Cpos[0][0]=(((s->isx12CSi%16+0.5)*((s->diaSi-s->innSi)/2)/(s->nxSi/4))+s->innSi/2) * cos(((s->isy12CSi+0.5)*22.5)*D_TO_R);
      s->m12Cpos[0][1]=(((s->isx12CSi%16+0.5)*((s->diaSi-s->innSi)/2)/(s->nxSi/4))+s->innSi/2) * sin(((s->isy12CSi+0.5)*22.5)*D_TO_R);
      rotvec(s->m12Cpos[0],s->m12Cpos[1],1,s->thrSi);
      //      printf("posx, posy, posz : %f, %f %f\n",s->m12Cpos[0][0],s->m12Cpos[0][1],s->m12Cpos[0][2]);
      
      double tmppos[3]; 
      tmppos[0]=s->m12Cpos[1][0]-s->offx1;
      tmppos[1]=s->m12Cpos[1][1]-s->offy1;
      tmppos[2]=s->m12Cpos[1][2]-0;

      s->mth12C=gettheta(tmppos);
      s->mph12C=getphi(tmppos);
      s->mthr12C=s->mth12C*D_TO_R;
      s->mphr12C=s->mph12C*D_TO_R;
      //      printf("pos,theta:%f,%f\n",s->m4Hepos[1],s->mth4He);
      
      tmpx=sin(s->mthr12C)*cos(s->mphr12C);
      tmpy=sin(s->mthr12C)*sin(s->mphr12C);
      tmpz=cos(s->mthr12C);
      s->mthxr12C=atan(tmpx/tmpz);
      s->mthyr12C=atan(tmpy/tmpz);
      s->mthx12C=s->mthxr12C*R_TO_D;
      s->mthy12C=s->mthyr12C*R_TO_D;
      //      s->mex12C=calcex4(s->mthr12C,k->m1,k->m2,k->m1,k->m2,s->enepa,s->mK12C); //need correction
      s->mex12C=s->enepa-s->mK12C-s->mK4He;

      s->mvec12C[0]=k->m4+s->mK12C;
      //  /*Opposite direction from 12C*/
      genvthph(-1.0*s->mth12C,s->mph12C,&s->mvec12C[1]);
      vecadd(&s->mvec12C[1],null,&s->mvec12C[1],
	     sqrt(s->mvec12C[0]*s->mvec12C[0]-k->m4*k->m4),1);
    }
    else s->idp12C=0;
  }
  return(s->idp12C);
};




/*******************************************************/
/* int detect3a(struct skine *k, struct ssimval *s) */
/* 12C detection by Si-CsI                        */
/* Input parameters:                                    */
/* k->m4 ... Rest mass of particle 4 (12C) */
/* s->nstat ... Number of states in 12C */
/* k->m2 ... Rest mass of 12C in the ground state */
/* s->sK4[1]  ... Energy of particle 4 after the target */
/* s->sthr4[1] ... Polar scatt. ang. of 4 after the target */
/* s->sphr4[1] ... Azimuthal scatt. ang. of 4 after the target */
/* s->sth4[1] ... Polar scatt. ang. of 4 after the target in deg. */
/* s->sph4[1] ... Azimuthal scatt. ang. of 4 after the target in deg. */
/* s->sx1 ... Lateral spread of 12C Beam in horizontal plane */
/* s->sy1 ... Lateral spread of 12C Beam in vertical plane */

/* s->dSi ... distance between target and Si in mm */
/* s->vSi ... Vertical size of Si in mm */
/* s->hSi ... Horizontal size of Si in mm */
/* s->nxSi ... Number of horizontal strips in Si */
/* s->nySi ...  Number of vertical strips in Si */
/* s->thSi ...  Angle of Si in degree */
/* s->thrSi ... Angle of Si in radian */
/* s->sigK4 ... Energy resolution of recoiled proton in MeV */

/* Output parameters:                                  */
/* s->idp3a ... Detection flag for particle 4 ..Not detected  1..detected */
/* s->s12Cvec[2][3] ... 3-dim direction of 4He trajectory */
                      /* [0][] ... Detector frame */
                      /* [1][] ... Laboratory frame */
/* s->s12Cpos[3] ... Position in the detector frame */
/* s->mK12C ... Measured energy of proton */
/* s->m12Cpos[2][3] ... Measured position of proton */
                      /* [0][] ... Detector frame */
                      /* [1][] ... Laboratory frame */
/* s->mthr12C ... Measured polar angle of proton in radian */
/* s->mphr12C ... Measured azimuthal angle of proton in radian */
/* s->mth12C ... Measured polar angle of proton in degree */
/* s->mph12C ... Measured azimuthal angle of proton in degree */
/* s->mthxr12C ... Measured horizontal angle of proton in radian */
/* s->mthyr12C ... Measured vertical angle of proton in radian */
/* s->mthx12C ... Measured horizontal angle of proton in degree */
/* s->mthy12C ... Measured vertical angle of proton in degree */
/* s->mex12C ... Measured excitation energy of 12C */
/* s->mvec4He[4] ... Measured four vector of recoiled proton */
/* Return value: Detection flag for particle 4 (s->idp3a) */



int detect3a(struct skine *k, struct ssimval *s){
  double m4he=mgs[0];
  double slatex[2][3][3];
  double null[3]={0,0,0};
  double tmpr1=0,tmpr2=0,tmpx=0,tmpy=0,tmpz=0;
  double tmpr[3]={0,0,0};
  double tmptheta=0;
  int quadrant=-1;
  int i,j;

  for(j=0;j<3;j++){
    genvthph(s->sthaf[j],s->sphaf[j],s->s3avec[1][j]);
    rotvec(s->s3avec[1][j],s->s3avec[0][j],1,-1.0*s->thrSi);

    //    cout << s->sthaf[j] << " " << s->sphaf[j] << endl;
   
    genrndg(&tmpr1,&tmpr2);
    //  s->mK12C=s->sK4[1]+tmpr1*s->sigK4-elossc(s->sK4[1]+tmpr1*s->sigK4,s->tdeadSi,s->tthick);
    s->mK3a[j]=s->sK3a[j]+tmpr1*s->sigK4;
    if(s->s3avec[0][j][2]>0){
      slatex[1][j][0]=s->sx1;
      slatex[1][j][1]=s->sy1;
      slatex[1][j][2]=0;
      rotvec(slatex[1][j],slatex[0][j],1,-1.0*s->thrSi);
      vecadd(s->s3avec[0][j],slatex[0][j],s->s3apos[j],
	     (s->dSi-slatex[0][j][2])/s->s3avec[0][j][2],1.0);    
      tmpr[j] = pow((pow(s->s3apos[j][0],2)+pow(s->s3apos[j][1],2)),0.5);
      if(s->s3apos[j][0]>0) tmptheta=atan(s->s3apos[j][1]/s->s3apos[j][0])*R_TO_D;
      if(s->s3apos[j][0]<0) tmptheta=atan(s->s3apos[j][1]/s->s3apos[j][0])*R_TO_D+180;

      
      if((tmpr[j] < s->diaSi/2) && (tmpr[j] > s->innSi/2)){
	s->idp3a[j]=1;
	if(s->s3apos[j][0]>0 && s->s3apos[j][1]>0){
	  quadrant=0;
	}else if(s->s3apos[j][0]<0 && s->s3apos[j][1]>0){
	  quadrant=1;
	}else if(s->s3apos[j][0]<0 && s->s3apos[j][1]<0){
	  quadrant=2;
	}else if(s->s3apos[j][0]>0 && s->s3apos[j][1]<0){
	  quadrant=3;
	}else{
	  quadrant=-100;
	  printf("detection process error");
	}
	s->m3aseg[j] = quadrant;
	//      printf("quadrant:%d\n",quadrant);
      
	s->isx3aSi[j]=(int)((tmpr[j]-s->innSi/2)/((s->diaSi-s->innSi)/2)*(s->nxSi/4)) + 16*quadrant;//%% 
	if(quadrant==0) s->isy3aSi[j] = (int)(((int)(atan(s->s3apos[j][1]/s->s3apos[j][0])*R_TO_D)+0)/22.5) +4*0;
	if(quadrant==1) s->isy3aSi[j] = (int)(((int)(atan(s->s3apos[j][1]/s->s3apos[j][0])*R_TO_D)+90)/22.5) +4*1;
	if(quadrant==2) s->isy3aSi[j] = (int)(((int)(atan(s->s3apos[j][1]/s->s3apos[j][0])*R_TO_D)+0)/22.5) +4*2;
	if(quadrant==3) s->isy3aSi[j] = (int)(((int)(atan(s->s3apos[j][1]/s->s3apos[j][0])*R_TO_D)+90)/22.5) +4*3;
	//      printf("x,y:%d,%d\n",s->isx12CSi,s->isy12CSi);
      
	s->m3apos[0][j][2]=s->s3apos[j][2];
	s->m3apos[0][j][0]=(((s->isx3aSi[j]%16+0.5)*((s->diaSi-s->innSi)/2)/(s->nxSi/4))+s->innSi/2) * cos(((s->isy3aSi[j]+0.5)*22.5)*D_TO_R);
	s->m3apos[0][j][1]=(((s->isx3aSi[j]%16+0.5)*((s->diaSi-s->innSi)/2)/(s->nxSi/4))+s->innSi/2) * sin(((s->isy3aSi[j]+0.5)*22.5)*D_TO_R);
	rotvec(s->m3apos[0][j],s->m3apos[1][j],1,s->thrSi);
	//	printf("posx, posy, posz : %f, %f %f\n",s->m3apos[1][j][0],s->m3apos[1][j][1],s->m3apos[1][j][2]);
      
	double tmppos[3][3]; 
	tmppos[j][0]=s->m3apos[1][j][0]-s->offx1;
	tmppos[j][1]=s->m3apos[1][j][1]-s->offy1;
	tmppos[j][2]=s->m3apos[1][j][2]-0;


	s->mth3a[j]=gettheta(tmppos[j]);
	s->mph3a[j]=getphi(tmppos[j]);
	s->mthr3a[j]=s->mth3a[j]*D_TO_R;
	s->mphr3a[j]=s->mph3a[j]*D_TO_R;
	//	printf("%f %f\n",s->mth3a[j], s->mph3a[j]);
      
	tmpx=sin(s->mthr3a[j])*cos(s->mphr3a[j]);
	tmpy=sin(s->mthr3a[j])*sin(s->mphr3a[j]);
	tmpz=cos(s->mthr3a[j]);
	s->mthxr3a[j]=atan(tmpx/tmpz);
	s->mthyr3a[j]=atan(tmpy/tmpz);
	s->mthx3a[j]=s->mthxr3a[j]*R_TO_D;
	s->mthy3a[j]=s->mthyr3a[j]*R_TO_D;

	s->mvec3a[j][0]=m4he+s->mK3a[j];
	genvthph(s->mth3a[j],s->mph3a[j],&s->mvec3a[j][1]);
	vecadd(&s->mvec3a[j][1],null,&s->mvec3a[j][1],sqrt(s->mvec3a[j][0]*s->mvec3a[j][0]-m4he*m4he),1);

      }
    }else{
      s->idp3a[j]=0;
    }  
    s->svec3a[j][0]=m4he+s->sK3a[j];
    genvthph(s->sthaf[j],s->sphaf[j],&s->svec3a[j][1]);
    vecadd(&s->svec3a[j][1],null,&s->svec3a[j][1],sqrt(s->svec3a[j][0]*s->svec3a[j][0]-m4he*m4he),1);
        
  }
  if(s->idp3a[0]==1 && s->idp3a[1]==1 && s->idp3a[2]==1){
    double tmp_num = pow(s->mvec3a[0][1]+s->mvec3a[1][1]+s->mvec3a[2][1],2)+pow(s->mvec3a[0][2]+s->mvec3a[1][2]+s->mvec3a[2][2],2)+pow(s->mvec3a[0][3]+s->mvec3a[1][3]+s->mvec3a[2][3],2)+pow(mgs[2],2);
    if(tmp_num<0) printf("Error");
    double total3a = pow(tmp_num,0.5);
    s->mex3a = s->mK3a[0]+s->mK3a[1]+s->mK3a[2]+mgs[0]*3-total3a;      
  }
  
  double tmp_num1 = pow(s->svec3a[0][1]+s->svec3a[1][1]+s->svec3a[2][1],2)+pow(s->svec3a[0][2]+s->svec3a[1][2]+s->svec3a[2][2],2)+pow(s->svec3a[0][3]+s->svec3a[1][3]+s->svec3a[2][3],2)+pow(mgs[2],2);
  double stotal3a = pow(tmp_num1,0.5);
  if(tmp_num1<0) printf("Error");
  s->sex3a = s->sK3a[0]+s->sK3a[1]+s->sK3a[2]+mgs[0]*3-stotal3a;

  
  return(s->idp3a[0]);
};






double culc_3a(double *enea, int *tha, int *pha){

  double m4He=931.494*4+2.4249;
  double m12C=931.494*12;
  double ene3a[3]={};
  double th3a[3]={};
  double ph3a[3]={};
  double moma[3][3]={};
  double mom12C[3]={0,0,0};
  double total12C=-1;
  double ex12C=-1;
  
  for(int n=0; n<3; n++){
    ene3a[n] = enea[n];
    th3a[n] = tha[n];
    ph3a[n] = pha[n];
    //    cout << tha[n] << " " << pha[n] << endl;
    //    cout << tha[n]*180/3.14 << " " << pha[n]*180/3.14 << endl;
    
    moma[n][0] = pow(2*m4He*ene3a[n]+ene3a[n]*ene3a[n],0.5) * sin(tha[n]) * cos(pha[n]);
    moma[n][1] = pow(2*m4He*ene3a[n]+ene3a[n]*ene3a[n],0.5) * sin(tha[n]) * sin(pha[n]);
    moma[n][2] = pow(2*m4He*ene3a[n]+ene3a[n]*ene3a[n],0.5) * cos(tha[n]);
    //    cout << ene3a[n] << " " << moma[n][0]  << " " << moma[n][1] << " "  << moma[n][2] << endl;
    //    cout << m4He+ene3a[n] << " " << moma[n][0]  << " " << moma[n][1] << " "  << moma[n][2] << endl;

    mom12C[0] += moma[n][0];
    mom12C[1] += moma[n][1];
    mom12C[2] += moma[n][2];
  }
  //  cout << mom12C[0]  << " " << mom12C[1] << " "  << mom12C[2] << endl;

  total12C = pow(mom12C[0]*mom12C[0]+mom12C[1]*mom12C[1]+mom12C[2]*mom12C[2]+m12C*m12C,0.5);
  ex12C = ene3a[0]+ene3a[1]+ene3a[2]+m4He*3 - total12C;
  
  return ex12C;
}








/******************************************************/
/*  double genrndg(double *a, double *b)              */
/*  Generate random standard distribution             */
/*  double a, b; random numbers (out)                 */
/*  Retrun value: one of two random numbers          */
double genrndg(double *a, double *b){
  double a1,a2,c;
  do{
    a1=2.0*drand48()-1.0;
    a2=2.0*drand48()-1.0;
    c=a1*a1+a2*a2;
  }while(c >= 1.0);
  c = sqrt((-2.0*log(c))/c);
  *a=a1*c;
  *b=a2*c;
  return(*a);
}

/******************************************************/
/*
double eloss4He(double ene,double thick);
double estra4He(double ene,double thick);
double astr4Hep(double ene,double thick);
double elossc(double ene,double thick);
double estrac(double ene,double thick);
double astrac(double ene,double thick);
  Get energy loss, energy straggring, angular straggring for
  proton with 5 MeV < E < 100 MeV
  3a with 100 MeV < E < 1000 MeV

  double ene ... energy of the particle in MeV
  double thick ... mass thickness of the material in the unit of 1 mg/cm^2
  
  All the values are estimated based on LISE++ for C 1 mg/cm^2 foil.
 
  Return value: Energy loss in MeV
                Energy straggling in MeV
                Angular straggling in plane in radian  
                -1000 if ene is out of range.
*/

double eloss4He(double ene,double thick, double tthick){
  int i;
  static int np=5;
  static double par[5]={8.588e-2, -1.368e-2, 1.126e-3, -4.383e-5, 6.4e-7};
  double eloss=0;
  //  if(ene<0.5 || ene>100) return(ene*thick/tthick/2);
  for(i=0;i<np;i++) eloss+=par[i]*pow(ene,(double)i);
  if(eloss>ene) eloss=ene;
  //  eloss=0;
  return(eloss*thick/tthick);
}


double estra4He(double ene,double thick, double tthick){
  int i;
  double eloss=0;
  //  if(ene<0.5 || ene>100) return(0);
  if(0.012<ene && ene<10){
    eloss=0.012;
  }else if(ene>10){
    eloss=0.0108;
  }else{
    eloss=ene;
  }    
  //  eloss=0;
  return(eloss*thick/tthick);
}

double astra4He(double ene,double thick, double tthick){
  int i;
  static int np=2;
  static double par[2]={-4.38e-2, 1.127e+1};
  double eloss=0,eneinv;
  if(ene<0.2 || ene>100) return(0);
  eneinv=1.0/ene;
  for(i=0;i<np;i++) eloss+=par[i]*pow(eneinv,(double)i);
  return(eloss*thick/tthick/1000.);
}

double elossc(double ene,double thick, double tthick){
  int i;
  static int np=3;
  static double par[3]={0.394, -9.628e-3, 9.2e-5};
  double eloss=0;
  for(i=0;i<np;i++) eloss+=par[i]*pow(ene,(double)i);
  if(eloss>ene) eloss=ene;
  //  eloss=0;
  return(eloss*thick/tthick);
}

double estrac(double ene,double thick, double tthick){
  int i;
  double eloss=0;
  //  if(ene<0.5 || ene>100) return(0);
  if(0.012<ene && ene<10){
    eloss=0.012;
  }else if(ene>10){
    eloss=0.0108;
  }else{
    eloss=ene;
  }    
  //  eloss=0;
  return(eloss*thick/tthick);
}

double astrac(double ene,double thick, double tthick){
  int i;
  static int np=2;
  static double par[2]={-1.187e-2, 3.444e+1};
  double eloss=0,eneinv;
  if(ene<0.2 || ene>100) return(0);
  eneinv=1.0/ene;
  for(i=0;i<np;i++) eloss+=par[i]*pow(eneinv,(double)i);
  return(eloss*thick/tthick/1000);
}


/******************************************************/
double gettheta4(double *p){
  /* Get scattering angle of particle */
  /* double *p ... four momentum of A (in) */
  /* return value ... theta (deg, angle from z-axis) */
  double vtra,theta;
  vtra=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
  theta=acos(p[3]/vtra)*R_TO_D;
  return(theta);
}

/******************************************************/
double gettheta(double *p){
  /* Get scattering angle of particle */
  /* double *p ... three momentum of A (in) */
  /* return value ... theta (deg, angle from z-axis) */
  double vtra,theta;
  vtra=sqrt(p[1]*p[1]+p[2]*p[2]+p[0]*p[0]);
  theta=acos(p[2]/vtra)*R_TO_D;
  return(theta);
}


/******************************************************/
double deflect2(double *v1,double *v2,
		double th1,double ph1,double th2,double ph2){
  /* Deflect vector twice from z-axis*/
  /* double *v1 ... vector after 1st deflection (out)
     double *v2 ... vector after 2nd deflection (out)
     double th1 ... theta in 1st deflection in radian (in)
     double ph1 ... phi in 1st deflection in radian (in)
     double th2 ... theta in 2nd deflection in radian (in)
     double ph2 ... phi in 2nd deflection in radian (in) */
  double v0[3]={0,0,1};
  double vtmp[4][3];
  rotvec(v0,vtmp[0],1,th2);
  rotvec(vtmp[0],vtmp[1],2,ph2);
  rotvec(vtmp[1],vtmp[2],1,th1);
  rotvec(vtmp[2],v2,2,ph1);

  rotvec(v0,vtmp[0],1,th1);
  rotvec(vtmp[0],v1,2,ph1);
  return(gettheta(v2));
}



/******************************************************/
double getphi4(double *p){
  /* Get azumithal scattering angle of particle */
  /* double *p ... four momentum of A (in) */
  /* return value ... phi (deg, angle from z-axis) */
  double phi;
  phi=atan2(p[2],p[1])*R_TO_D;
  return(phi);
}

/******************************************************/
double getphi(double *p){
  /* Get azumithal scattering angle of particle */
  /* double *p ... three momentum of A (in) */
  /* return value ... phi (deg, angle from z-axis) */
  double phi;
  phi=atan2(p[1],p[0])*R_TO_D;
  return(phi);
}

/******************************************************/
 double genvdir(double ang,double *v){
  /* Generate scattered particle (A)
     double ang ... scattering angle theta (in)
     double *v ... direction of A (out)
     return value ... azumuthal scattering phi */
   double phir,angr;
   phir=drand48()*2.0*M_PI;
   angr=ang*D_TO_R;
   v[0]=sin(angr)*cos(phir);
   v[1]=sin(angr)*sin(phir);
   v[2]=cos(angr);
   return(phir*R_TO_D);
 }

/******************************************************/
double genvthph(double th,double  ph,double *v){
  /* Generate scattered particle (A) acconding to given th and ph
     double th ... scattering angle theta (in)
     double ph ... scattering angle phi (in)
     double *v ... direction of A (out)
     return value ... azumuthal scattering phi */
   double phir,angr;
   phir=ph*D_TO_R;
   angr=th*D_TO_R;
   v[0]=sin(angr)*cos(phir);
   v[1]=sin(angr)*sin(phir);
   v[2]=cos(angr);
   return(phir*R_TO_D);
 }

/******************************************************/
double genini(double *p, double *v, double ex, double t, double m){
  /* Generate initial particle
     double *p ... four momentum of A (out)
     double *v ... direction of A (in)
     double ex ... Excitation energy of A (in)
     double t ... Kinetic energy of A (in)
     double m ... Mass of Ground State of A (in)
     return value ... Invariant mass of A */
  int i;
  double pmag,univ[3],massp;
  unitvec(v,univ);
  massp=m+ex;
  p[0]=massp+t;
  pmag=sqrt(p[0]*p[0]-massp*massp);
  for(i=1;i<4;i++) p[i]=pmag*univ[i-1];
  return(sqrt(scapro4(p,p)));
}

/******************************************************/
int gengamma(double eg,double *p0, double *p1, double *p2){
  /* Gammna emission from mother nucleus.
     double eg ... Energy of gamma (in)
     double *p0 ... four momentum of mother nucleus (in)
     double *p1 ... four momentum of gamma (out)
     double *p2 ... four momentum of daughter nucleus (out)
     return value ... 1 successfully generated.
                      0 fail
*/

  double betag[3];
  double m0p,m1p,m2p;  /* mass of paticles */
  double p1c[4],p2c[4];
  int j;

  getbeta(p0,betag);
  vecadd(betag,betag,betag,0.0,-1.0);

  m0p=sqrt(scapro4(p0,p0));
  m1p=0.0;
  m2p=m0p-eg;


  if(decay2body(p1c,p2c,m0p,m1p,m2p)<0){
    return(0); //%
  }
  lortra(p1c,p1,betag);
  lortra(p2c,p2,betag);

  return(1);
}


/******************************************************/
int genalpha(int i, double *p0, double *p1, double *p2){
  /* Alpha emission from mother nucleus.
     int i ... id of mother nuclus (must be >= 1)
     double *p0 ... four momentum of mother nucleus
     double *p1 ... four momentum of alpha (emitted particle)
     double *p2 ... four momentum of daughter nucleus
     return value ... 1 successfully generated.
                      0 fail
*/

  double betag[3];
  double m0p,m1p,m2p;  /* mass of paticles */
  double p1c[4],p2c[4];
  int j;

  getbeta(p0,betag);
  vecadd(betag,betag,betag,0.0,-1.0);
  m0p=sqrt(scapro4(p0,p0));
  m1p=mgs[0];
  /*  m2p=mgs[i-1]+ex_ene[i-1];*/

  if(decay2body(p1c,p2c,m0p,m1p,m2p)<0) return(0);

  lortra(p1c,p1,betag);
  lortra(p2c,p2,betag);

  return(1);
}





double tmp(double thetacm_i=40., double thetacm_f=75., int n_div=350, int i_div=0){
  double p0 = -241.512; //fitted using eventbuild_cor/hit/count/hist.C                                                                                                                                                               
  double p1 = 3.3118;
  double p2 = 112.729;
  double p3 = 54.2054;
  double p4 = 8.41897;
  double p5 = 130.893;
  double p6 = 36.5291;
  double p7 = 5.08719;

  double x = thetacm_i + (thetacm_f-thetacm_i)/n_div * (double)i_div;
  double num;
  num = p0 + p1*x + p2*exp(-pow(x-p3,2)/(2*pow(p4,2))) + p5*exp(-pow(x-p6,2)/(2*pow(p7,2)));

  return num;
}
