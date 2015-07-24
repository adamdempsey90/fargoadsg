#include "mp.h"

//#define POWERPROF
//#define GAUSSPROF
//#define TAPERPROF
#define MKLINPROF
//#define LAUGHLINPROF

extern real ScalingFactor;


#ifdef GAUSSPROF
static const real rmax1 = .45;
static const real pwidth1 = .05;
#endif

#ifdef TAPERPROF
static const real inner_slope = 3;
static const real rmax1 = 1;
#endif

#ifdef MKLINPROF
const static real h_p = .05;
const static real eps_p = .1;
const static real delta_R = 5.0;
const static real R1 = 1;
const static real R2 = 2;

#define s_p  (1 + FLARINGINDEX)
#define delta1 (delta_R * h_p * pow(R1,FLARINGINDEX + 1))
#define  delta2  (delta_R * h_p * pow(R2,FLARINGINDEX + 1))
real sech ();
real bump_function ();
real f1_func ();
real drf1_func ();
real d2rf1_func ();
real f2_func ();
real drf2_func ();
real d2rf2_func ();
#endif

#ifdef LAUGHLINPROF
const real dens_width2 = .1;
const real dens_peak_rad = .3;
#endif

/* Surface density */
real Sigma(r)
     real r;
{
  real res;
  real cavity = 1.0;
  if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO;
  /* This is *not* a steady state */
  /* profile, if a cavity is defined. It first needs */
  /* to relax towards steady state, on a viscous time scale */

#ifdef POWERPROF
  res=cavity*ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE);
#endif

#ifdef GAUSSPROF
  res =  SIGMA0*exp(-(r-rmax1)*(r-rmax1)/pwidth1);
#endif

#ifdef TAPERPROF
  real outer_slope = -SIGMASLOPE;
  res = SIGMA0 /( pow(r/rmax1,-inner_slope) + pow(r/rmax1,-outer_slope));
#endif

#ifdef MKLINPROF
  real sigfac = h_p /(2 * M_PI * bump_function(2.0)); // Qout = 2.0
  res = sigfac * bump_function(r) * pow(r,-s_p);
#endif

#ifdef LAUGHLINPROF
   res = SIGMA0 * exp(-(r-dens_peak_rad)*(r-dens_peak_rad)/(dens_width2));
#endif

//
  return res;
}

real dlogsigma_func(r)
  real r;
{
  real res;
#ifdef POWERPROF
  res = -SIGMASLOPE;
#endif

#ifdef GAUSSPROF
  res =  2*(rmax1 - r)*r/pwidth1;
#endif

#ifdef TAPERPROF
  real outer_slope = -SIGMASLOPE;
  real denom = pow(r/rmax1,inner_slope) + pow(r/rmax1,outer_slope);
  res = outer_slope + (inner_slope - outer_slope)*pow(r/rmax1,outer_slope)/denom;
#endif

#ifdef MKLINPROF
  res = -s_p + drf1_func(r)/f1_func(r) + drf2_func(r)/f2_func(r);
#endif

#ifdef LAUGHLINPROF
  res = 2*(dens_peak_rad - r)*r/dens_width2;
#endif

	return res;
}
real d2logsigma_func(r)
  real r;
{
  real res;
#ifdef POWERPROF
  res = 0;
#endif

#ifdef GAUSSPROF
  res =   2*(rmax1 - 2*r)*r/pwidth1;
#endif

#ifdef TAPERPROF
  real outer_slope = -SIGMASLOPE;
  real denom = pow(r/rmax1,inner_slope) + pow(r/rmax1,outer_slope);
  denom *= denom;
  res =  -(outer_slope-inner_slope)*(outer_slope-inner_slope)*pow(r/rmax1,inner_slope+outer_slope)/denom;
#endif
#ifdef MKLINPROF
  real d2lf1 = d2rf1_func(r) - drf1_func(r)*drf1_func(r)/f1_func(r);
  real d2lf2 = d2rf2_func(r) - drf2_func(r)*drf2_func(r)/f2_func(r);

  res = d2lf1/f1_func(r) + d2lf2/f2_func(r);
#endif

#ifdef LAUGHLINPROF
  res = 2*(dens_peak_rad - 2*r)*r/dens_width2;
#endif
  return res;
}


void FillSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    SigmaMed[i] = Sigma(Rmed[i]);
    SigmaInf[i] = Sigma(Rinf[i]);
  }
}

void RefillSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    SigmaMed[i] = moy;
  }
  SigmaInf[0] = SigmaMed[0];
  for (i = 1; i < nr; i++) {
    SigmaInf[i] = (SigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   SigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}

real scaleH_func(r)
  real r;
{
    return sqrt(temp_func(r))*pow(r,1.5);
}

real temp_func(r)
  real r;
{
  real res = ASPECTRATIO*ASPECTRATIO*pow(r, 2.0*FLARINGINDEX - 1.0);
#ifdef LAUGHLINPROF
   res =  0.74528*SIGMA0*FLARINGINDEX*pow(Sigma(r),FLARINGINDEX-1);
#endif
}

real dlogtemp_func(r)
  real r;
{
  real res = 2*FLARINGINDEX - 1;
#ifdef LAUGHLINPROF
  res = (FLARINGINDEX - 1)*dlogsigma_func(r);
#endif
  return res;
}
real d2logtemp_func(r)
  real r;
{
  real res = 0;
#ifdef LAUGHLINPROF
  res = (FLARINGINDEX - 1)*d2logsigma_func(r);
#endif
  return res;
}


/* Thermal energy */
real Energy(r)
     real r;
{
  real energy0;
  if (ADIABATICINDEX == 1.0) {
    fprintf (stderr, "The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.\n");
    prs_exit (1);
  }
  else
    energy0 = R/MU/(ADIABATICINDEX-1.0)*Sigma(r)*temp_func(r); //*pow(ASPECTRATIO,2.0)*pow(r,-1.0+2.0*FLARINGINDEX);
  return energy0;
}

void FillEnergy() {
  int i;
  for (i = 0; i < NRAD; i++)
    EnergyMed[i] = Energy(Rmed[i]);
}


void RefillEnergy (energy)
     PolarGrid *energy;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = energy->Nrad;
  ns = energy->Nsec;
  field = energy->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    EnergyMed[i] = moy;
  }
}

/* Cooling time */
real CoolingTime(r)
     real r;
{
  real ct0;
  ct0 = COOLINGTIME0*pow(r,1.5); // t_cool = 1/(\beta \Omega) = COOLINGTIME0/\Omega
//  ct0 = COOLINGTIME0*pow(r,2.0+2.0*FLARINGINDEX);
  return ct0;
}

void FillCoolingTime() {
  int i;
  for (i = 0; i < NRAD; i++)
    CoolingTimeMed[i] = CoolingTime(Rmed[i]);
}

/* Heating source term */
real Qplusinit(r)
     real r;
{
  real qp0, viscosity;
  viscosity = FViscosity(r);
  qp0 = 2.25*viscosity*Sigma(r)*pow(r,-3.0);
  return qp0;
}

void FillQplus() {
  int i;
  for (i = 0; i < NRAD; i++)
    QplusMed[i] = Qplusinit(Rmed[i]);
}

/* AMD 7/7/15 */
void AddUserIC (Rho, Vr, Vt, Energy)
	PolarGrid *Rho, *Vr, *Vt, *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *vr, *vt, *energy;
  real lambdaWave, initAmp, infsig;
  dens = Rho->Field;
  vr = Vr->Field;
  vt = Vt->Field;
  energy = Energy->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  lambdaWave = GlobalRmed[nr-1] - GlobalRmed[0];
  srand((unsigned int)time(NULL));
  initAmp = .0001;
 // infsig = Sigma(GlobalRmed[GLOBALNRAD-1]);
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
  //    dens[l] += infsig*initAmp*cos(2*M_PI*j/ns); // m=1 perturbation
      dens[l] = dens[l]*(1 + initAmp * (2*(real)rand()/(real)(RAND_MAX) - 1)); // random noise
      if (dens[l] < 0)	printf("Error Sigma < 0, CPU: %d,\t(%g,%g,%g)\n",CPU_Rank,Rmed[i],2*M_PI*j/ns,dens[l]);
    }
  }
  return;
}

#ifdef MKLINPROF

real sech(r)
  real r;
{
    return 1/cosh(r);
}
real f1_func(r)
  real r;
{
	return .5*(1-eps_p)*(1+tanh((r-R1)/delta1))+eps_p;

}

real f2_func(r)
  real r;
{
	return  .5*(1-eps_p)*(1-tanh((r-R2)/delta2))+eps_p;
}

real drf1_func(r)
  real r;
{
  real arg = (r - R1)/delta1;
	return  -r * (eps_p - 1) * sech(arg)*sech(arg)/(2*delta1);
}
real drf2_func(r)
  real r;
{
	real arg = (r - R2)/delta2;

	return  r * (eps_p - 1) * sech(arg)*sech(arg)/(2*delta2);

}

real d2rf1_func(r)
  real r;
{
	real arg = (r - R1)/delta1;
	return r*(eps_p - 1)*sech(arg)*sech(arg)*(2*r*tanh(arg)-delta1)/(2*delta1*delta1);
}

real d2rf2_func(r)
  real r;
{
	real arg = (r - R2)/delta2;
	return r*(eps_p - 1)*sech(arg)*sech(arg)*(-2*r*tanh(arg)+delta2)/(2*delta2*delta2);
}

real bump_function(r)
  real r;
{
	return f1_func(r) * f2_func(r);
}
#endif
