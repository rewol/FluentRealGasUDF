/**************************************************************/
/* UDRGM : Volume Corrected Peng Robinson Equation of State   */
/* Wan der Waals Mixture rules are utilized to calculate      */
/* density and thermodynamic properties. As for transport     */
/* properties, CHUNG method is used to calculate viscosity    */
/* and thermal conductivity. For mass diffusivity, Takahasi   */
/* method is used to account for high-pressure affects.       */
/* Volume correction is achieved with volume-translation      */
/* method, suggested by Abudour et al.                        */
/**************************************************************/

/**************************************************************/
/*															  */
/* This code includes 9 species and was created to utilized   */
/* the modified Jones-Lindstedt (JL-R) kinetics model. One    */
/* needs to remember that addition/extraction of the species  */
/* requires species data to be modified down below. A manual  */
/* will be provided in the future models.                     */
/* Species are as follows :									  */
/* CH4, O2, CO, H2, H20, O, H, OH, CO2                        */
/**************************************************************/

/**************************************************************/
/*															  */
/* AUTHOR : Refik Alper Tuncer                                */
/* Date: February, 2018                                       */
/* Version: 1.2                                               */
/*                                                            */
/**************************************************************/

#include "udf.h"
#include "stdio.h"
#include "ctype.h"
#include "stdarg.h"

#define N_SPECIES 9
#define N_SPECIES_NAME 80
#define R 8314.34
#define PI 3.141592654

static char gas[N_SPECIES][N_SPECIES_NAME];
static double ref_p;
static double ref_T;
static int (*usersMessage)(char *,...);
static void (*usersError)(char *,...);

/* STATIC PROPERTY PARAMETERS */
static double mw[N_SPECIES]; 
static double hf[N_SPECIES];
static double t_crit[N_SPECIES];
static double v_crit[N_SPECIES];
static double p_crit[N_SPECIES];
static double z_crit[N_SPECIES];
static double w[N_SPECIES];

/* STATIC PARAMETERS ASSOCIATED WITH REDLICH-KWONG REAL GAS MODEL */
static double a0[N_SPECIES];
static double b[N_SPECIES];
static double m[N_SPECIES];
static double rgas[N_SPECIES];

/* STRUCTURE DEFINITION FOR MIXTURE FUNCTIONS */
struct thermo {
	double temperature;
	double den;
	double press;
	double species[N_SPECIES];
	double molefrac[N_SPECIES];
	double parameter1;
	double parameter2;
	double parameter3;
	int parameter4;
	double parameter5;
	double parameter6; /* dpdv for volume correction */
};

/* DECLARATION FOR FUNCTIONS THAT ARE NEEDED FOR SPECIES CALCULATIONS */
double a_mixture_wanderwaals(struct thermo *argument);
double b_mixture_wanderwaals(struct thermo *argument);
double dvdt_mixture(struct thermo *argument);
double dvdp_mixture(struct thermo *argument);
double a_mixture_firder(struct thermo *argument);
double a_mixture_secder(struct thermo *argument);
double dpdt_mixture(struct thermo *argument);
double dpdv_mixture(struct thermo *argument);
double RK_Ideal_cp(struct thermo *argument);
double RK_Ideal_Enthalpy(struct thermo *argument);
double RK_Ideal_Entropy(struct thermo *argument);
double collision_integral(double temp, double e_k_mixture);
double mixture_sigma(double sig[N_SPECIES][N_SPECIES], double yi[]);
double mixture_critical_volume(double mixture_sigma);
double mixture_critical_temp(double mix_e_k);
double mixture_e_k(double sig[N_SPECIES][N_SPECIES], double yi[], double e_k[N_SPECIES][N_SPECIES]);
double mixture_mol_weight(double sigma[N_SPECIES][N_SPECIES], double yi[], double e_k[N_SPECIES][N_SPECIES], double MW[N_SPECIES][N_SPECIES]);
double mixture_acc_fac(double sigma[N_SPECIES][N_SPECIES], double yi[], double ww[N_SPECIES][N_SPECIES]);
double mixture_dip_mom(double sigma[N_SPECIES][N_SPECIES], double yi[], double dipole_moment[], double e_k[N_SPECIES][N_SPECIES]);
double mixture_corr_fac(double yi[], double k[N_SPECIES][N_SPECIES]);
double low_pressure_vis(double temp, double sigma[N_SPECIES][N_SPECIES], double yi[], double e_k[N_SPECIES][N_SPECIES], double MW[N_SPECIES][N_SPECIES], double dipole_moment[], double ww[N_SPECIES][N_SPECIES], double k[N_SPECIES][N_SPECIES]);
double ideal_specific_volume(double temp, double yi[]);
double species_sum(struct thermo *argument);
double species_sum_two(int i, double temp, double yi[]);
double volume_corr(struct thermo *argument);


DEFINE_ON_DEMAND(I_do_nothing)
{
}
void Mixture_error(int err, char *f, char *msg)
{
if (err)
usersError("Mixture_error (%d) from function: %s\n%s\n",err,f,msg);
}

/*******************************************************************/
/* Transport Properties Combination values in CHUNG et al.         */
/* These are the only functions called from ANSYS FLUENT Code      */
/*******************************************************************/

/* CH4 = 0, O2 = 1, CO = 2, H2 = 3, H20 = 4, O = 5, H = 6, OH = 7, CO2 = 8 */

double mw[N_SPECIES] = {16.04303, 31.99880, 28.01, 2.016, 18.01534, 15.994, 1.00797, 17.007, 44.00995};
double hf[N_SPECIES] = {-74895176, 0, -3946090, 0, -2.418379e+08, 249197000, 217977000, 38985000, -3.9353235e+08};
double t_crit[N_SPECIES] = {190.56, 154.58, 132.85, 32.98, 647.14, 44.5, 33.2, 400, 304.12};
double v_crit[N_SPECIES] = {0.006146, 0.002294, 0.003223, 0.0318, 0.003106, 0.0022, 0.0185, 0.0027, 0.002136};
double p_crit[N_SPECIES] = {45.99*1e5, 50.43*1e5, 34.94*1e5, 12.93*1e5, 220.64*1e5, 26.90*1e5, 1360000, 8200000, 73.74*1e5};
double w[N_SPECIES] = {0.0114, 0.021, 0.045, -0.217, 0.344, 0, 0, 0.2, 0.225};
double z_crit[N_SPECIES] = {0.286, 0.288, 0.292, 0.303, 0.229, 0.255, 0.091, 0.113, 0.273};
 

/*******************************************************************/
/* GLOBILIZED EoS parameters to be used in MIXTURE functions       */
/*                                                                 */
/*******************************************************************/

double a_mixture = 0;
double b_mixture = 0;
double dadt = 0;

/*******************************************************************/
/* Mixture Functions                                               */
/* These are the only functions called from ANSYS FLUENT Code      */
/*******************************************************************/

void MIXTURE_Setup(Domain *domain, cxboolean vapor_phase, char *specielist, int (*messagefunc)(char *format,...), void (*errorfunc)(char *format,...))
{
unsigned int i;
usersMessage = messagefunc;
usersError = errorfunc;
ref_p = ABS_P(RP_Get_Real("reference-pressure"),op_pres);
ref_T = 298.15;
Message0("\n MIXTURE_Setup: Peng-Robinson equation of State with Wan der Waals mixing rules \n");
Message0("\n MIXTURE_Setup: reference-temperature is %f \n", ref_T);

if (ref_p == 0.0)
{
Message0("\n MIXTURE_Setup: reference-pressure was not set by user \n");
Message0("\n MIXTURE_Setup: setting reference-pressure to 101325 Pa \n");
ref_p = 101325.0;
}

(void)strcpy(gas[0],"CH4");
(void)strcpy(gas[1],"O2");
(void)strcpy(gas[2],"CO2") ;
(void)strcpy(gas[3],"H2O");
(void)strcpy(gas[4],"N2") ;

Message0("\n MIXTURE_Setup: RealGas mixture initialization \n");
Message0("\n MIXTURE_Setup: Number of Species = %d \n",N_SPECIES);

for (i=0; i<N_SPECIES; ++i)
{
Message0("\n MIXTURE_Setup: Specie[%d] = %s \n",i,gas[i]);
}

strcat(specielist,gas[0]);

for (i=1; i<N_SPECIES; ++i)
{
strcat(specielist," ");
strcat(specielist,gas[i]);
}


for (i=0; i<N_SPECIES; i++)
{
	rgas[i] = R/mw[i];
	a0[i] = 0.457235*rgas[i]*rgas[i]*t_crit[i]*t_crit[i]/p_crit[i];
	b[i] = 0.07779*rgas[i]*t_crit[i]/p_crit[i];
	m[i] = 0.37464+1.54226*w[i]-0.26992*w[i]*w[i];
}

}

/*********************************************/
/************ MIXTURE FUNCTIONS **************/
/*********************************************/

/* MIXTURE MOLECULAR WEIGHT */

double MIXTURE_Molecular_Weight(double yi[]) /*ok */
{
	double sum = 0;
	double mixture_molecular_weight;
	unsigned int i;
	
	for(i=0;i<N_SPECIES;i++)
	{
		sum = sum + yi[i]/mw[i];
	}
	
	mixture_molecular_weight = 1/sum;
	return mixture_molecular_weight; /* kg/kmol */
}

/* MIXTURE DENSITY */

double MIXTURE_Density(double temp, double pressure, double yi[]) /*ok */
{
	unsigned int i,j;
	struct thermo argument;
	argument.temperature = temp;
	argument.press = pressure;
	for(i=0;i<N_SPECIES;i++){
	argument.species[i] = yi[i];
	}
	double mixture_mw = MIXTURE_Molecular_Weight(yi);
		for(i=0;i<N_SPECIES;i++)
	{
		argument.molefrac[i] = yi[i]*mixture_mw/mw[i];
	}		
	argument.parameter5 = mixture_mw;
	double mixture_density;
	a_mixture = a_mixture_wanderwaals(&argument);
	b_mixture = b_mixture_wanderwaals(&argument);
	dadt = a_mixture_firder(&argument);
	double a1,a2,a3;
	double vv,vv1,vv2,vv3;
	double qq,qq3,sqq,rr,tt,dd;
	double rrgas = (R/MIXTURE_Molecular_Weight(yi));
	a1 = b_mixture - rrgas*temp/pressure;
	a2 = -3*pow(b_mixture,2) - (2*temp*b_mixture*rrgas/pressure) + a_mixture/pressure;
	a3 = pow(b_mixture,3) + (rrgas*temp*pow(b_mixture,2)/pressure) - a_mixture*b_mixture/pressure;
	
	qq = (a1*a1-3*a2)/9;
	rr = (2*a1*a1*a1-9*a1*a2+27*a3)/54;
	qq3 = qq*qq*qq;
	dd = qq3-rr*rr;
	
	if (dd < 0) 
	{
	tt = -SIGN(rr)*(pow(pow(-dd,0.5)+fabs(rr),0.333333));
	vv = (tt+qq/tt)-a1/3;
	} 
	else 
	{
	if (rr/pow(qq3,0.5)<-1) 
	{
	tt = PI;
	} 
	else if (rr/pow(qq3,0.5)>1) 
	{
	tt = 0;
	}
	else 
	{
	tt = acos(rr/pow(qq3,0.5));
	}
	sqq = pow(qq,0.5);
	vv1 = -2*sqq*cos(tt/3)-a1/3;
	vv2 = -2*sqq*cos((tt+2*PI)/3)-a1/3;
	vv3 = -2*sqq*cos((tt+4*PI)/3)-a1/3;
	vv = (vv1 > vv2) ? vv1 : vv2;
	vv = (vv > vv3) ? vv : vv3;
}
	
	argument.den = 1/vv;
	argument.parameter6 = dpdv_mixture(&argument);
	argument.den = vv;
	double correction = volume_corr(&argument);
	
	mixture_density = 1/(vv+correction);
	return mixture_density;	
}

/* MIXTURE SPECIFIC HEAT */

double MIXTURE_Specific_Heat(double temp, double density, double pressure, double yi[]) 
{
		
	unsigned int i,j;
	
	struct thermo argument1;
	argument1.temperature = temp;
	argument1.den = density;
	argument1.press = pressure;
	for(i=0;i<N_SPECIES;i++){
	argument1.species[i] = yi[i];
	}
	double mixture_mw = MIXTURE_Molecular_Weight(yi);
		for(i=0;i<N_SPECIES;i++)
	{
		argument1.molefrac[i] = yi[i]*mixture_mw/mw[i];
	}	
	argument1.parameter5 = mixture_mw;
	
	/*double a_mixture = a_mixture_wanderwaals(&argument1);
	double b_mixture = b_mixture_wanderwaals(&argument1);
	double dadt = a_mixture_firder(&argument1);*/
	double da2dt2 = a_mixture_secder(&argument1);
	argument1.parameter1 = a_mixture;
	argument1.parameter2 = b_mixture;
	argument1.parameter3 = dadt;
	
	double mixture_ideal_gas_cp = 0;
	double mixture_departure_cp;
	double v = 1/density;
	double dpdt = dpdt_mixture(&argument1);
	double dpdv = dpdv_mixture(&argument1);
	double rrgas = (R/MIXTURE_Molecular_Weight(yi));
	
	
	
	for(i=0;i<N_SPECIES;i++)
	{
		argument1.parameter4 = i;
		mixture_ideal_gas_cp = mixture_ideal_gas_cp + yi[i]*RK_Ideal_cp(&argument1);             /* BURADA R gibi SAYILARIN BİRİMİ KONTROL EDİLECEK */
	}
	/*Message("ideal cp is %g \n",mixture_ideal_gas_cp);*/
	mixture_departure_cp = -temp*pow(dpdt,2)/dpdv - rrgas - da2dt2*temp/(2.82842*b_mixture)*log((v-0.4142*b_mixture)/(v+2.4142*b_mixture));
	
	return mixture_ideal_gas_cp + mixture_departure_cp;
}

/* MIXTURE SENSIBLE ENTHALPY */

double MIXTURE_Enthalpy(double temp, double density, double pressure, double yi[])
{
	unsigned int i,j;
	
	struct thermo argument1;
	argument1.temperature = temp;
	argument1.den = density;
	argument1.press = pressure;
	for(i=0;i<N_SPECIES;i++){
	argument1.species[i] = yi[i];
	}
	double mixture_mw = MIXTURE_Molecular_Weight(yi);
		for(i=0;i<N_SPECIES;i++)
	{
		argument1.molefrac[i] = yi[i]*mixture_mw/mw[i];
	}		
	argument1.parameter5 = mixture_mw;
	/*double a_mixture = a_mixture_wanderwaals(&argument1);
	double b_mixture = b_mixture_wanderwaals(&argument1);
	double dadt = a_mixture_firder(&argument1);*/
	
	double mixture_ideal_gas_h = 0;
	double mixture_departure_h;
	double v = 1/density; 
	double rrgas = (R/MIXTURE_Molecular_Weight(yi));
	for(i=0;i<N_SPECIES;i++)
	{
		argument1.parameter4 = i;
		mixture_ideal_gas_h = mixture_ideal_gas_h + yi[i]*RK_Ideal_Enthalpy(&argument1);         /* BURADA R gibi SAYILARIN BİRİMİ KONTROL EDİLECEK */
	}
	mixture_departure_h = pressure*v - rrgas*temp + (a_mixture - temp*dadt)*log((v-0.4142*b_mixture)/(v+2.4142*b_mixture))/(2.82842*b_mixture);
	
	return mixture_ideal_gas_h + mixture_departure_h;	
}

/* MIXTURE ENTHALPY */

double MIXTURE_enthalpy_prime(double temp, double density, double pressure, double yi[], double hi[])
{
	unsigned int i,j;

	struct thermo argument1;
	argument1.temperature = temp;
	argument1.den = density;
	argument1.press = pressure;
	for(i=0;i<N_SPECIES;i++){
	argument1.species[i] = yi[i];
	}
	double mixture_mw = MIXTURE_Molecular_Weight(yi);
		for(i=0;i<N_SPECIES;i++)
	{
		argument1.molefrac[i] = yi[i]*mixture_mw/mw[i];
	}
	argument1.parameter5 = mixture_mw;
	/*double a_mixture = a_mixture_wanderwaals(&argument1);
	double b_mixture = b_mixture_wanderwaals(&argument1);
	double dadt = a_mixture_firder(&argument1);*/
	argument1.parameter1 = a_mixture;
	argument1.parameter2 = b_mixture;
	argument1.parameter3 = dadt;
	
		
	double mixture_ideal_gas_h = 0;
	double mixture_formation_h = 0;
	double mixture_departure_h;
	double v = 1/density;
	double vi[N_SPECIES];
	double dpdv = dpdv_mixture(&argument1);
	double rrgas = (R/MIXTURE_Molecular_Weight(yi));
	double K1 = log((v-0.4142*b_mixture)/(v+2.4142*b_mixture));
	double s = v*v + 2*v*b_mixture - b_mixture*b_mixture;
	double ss;
	for(i=0;i<N_SPECIES;i++) 
	{
		/*vi[i] = (-1/dpdv)*(rgas[i]*temp/(v-b_mixture) + rgas[i]*temp*b[i]/pow(v-b_mixture,2) - 2*species_sum(i,temp,yi)/(v*v + 2*v*b_mixture - pow(b_mixture,2)) + 2*a_mixture*(v-b_mixture)*b[i]/pow(v*v + 2*v*b_mixture - b_mixture*b_mixture,2)); /* işaretler değişebilir */
		/*hi[i] = hf[i]/mw[i] + ((species_sum(i,temp,yi) - temp*species_sum_two(i,temp,yi))/(2.8284*b_mixture) - (a_mixture - temp*dadt)*b[i]/(2.8284*pow(b_mixture,2)))*K1 + pressure*vi[i] - rgas[i]*temp + ((a_mixture - temp*dadt)/(b_mixture))*((b_mixture*vi[i] - b[i]*v)/((v+2.4142*b_mixture)*(v-0.4142*b_mixture))) + RK_Ideal_Enthalpy(temp,i);*/
		argument1.parameter4 = i;
		ss=species_sum(&argument1);
		vi[i] = (-1/dpdv)*(rrgas*temp/(v-b_mixture) + rrgas*temp*b[i]/pow(v-b_mixture,2)) - (1/dpdv)*(2*a_mixture*(v-b_mixture)*b[i]/pow(s,2) - 2*ss/s);
		hi[i] = hf[i]/mw[i] + RK_Ideal_Enthalpy(&argument1) + pressure*vi[i] - rrgas*temp + (a_mixture - temp*dadt)*(vi[i] - v*b[i]/b_mixture)/s - K1*(0.3535/b_mixture)*(ss - temp*species_sum_two(i,temp,yi) - (a_mixture - temp*dadt)*b[i]/b_mixture);
		mixture_ideal_gas_h = mixture_ideal_gas_h + yi[i]*RK_Ideal_Enthalpy(&argument1);  /* BURADA R gibi SAYILARIN BİRİMİ KONTROL EDİLECEK */
		mixture_formation_h = mixture_formation_h + hf[i]*yi[i]/mw[i];
	}
	

	/*for(i=0;i<N_SPECIES;i++)
	{
		mixture_ideal_gas_h = mixture_ideal_gas_h + yi[i]*RK_Ideal_Enthalpy(temp,i);  
		mixture_formation_h = mixture_formation_h + hf[i]*yi[i]/mw[i];
	}*/
	
	mixture_departure_h = pressure*v - rrgas*temp + (a_mixture - temp*dadt)*K1/(2.8284*b_mixture);
	return mixture_formation_h + mixture_ideal_gas_h + mixture_departure_h;
}

/* MIXTURE ENTROPY */

double MIXTURE_Entropy(double temp, double density, double pressure, double yi[])
{
	unsigned int i,j;

	struct thermo argument1;
	argument1.temperature = temp;
	argument1.den = density;
	argument1.press = pressure;
	for(i=0;i<N_SPECIES;i++){
	argument1.species[i] = yi[i];
	}
	double mixture_mw = MIXTURE_Molecular_Weight(yi);
		for(i=0;i<N_SPECIES;i++)
	{
		argument1.molefrac[i] = yi[i]*mixture_mw/mw[i];
	}	
	argument1.parameter5 = mixture_mw;
	/*double a_mixture = a_mixture_wanderwaals(&argument1);
	double b_mixture = b_mixture_wanderwaals(&argument1);
	double dadt = a_mixture_firder(&argument1);*/
	double mixture_ideal_gas_s = 0;
	double mixture_departure_s;
	double v = 1/density;
	double v0 = (R/MIXTURE_Molecular_Weight(yi))*temp/ref_p;
	double rrgas = (R/MIXTURE_Molecular_Weight(yi));
	for(i=0;i<N_SPECIES;i++)
	{
		argument1.parameter4 = i;
		mixture_ideal_gas_s = mixture_ideal_gas_s + yi[i]*RK_Ideal_Entropy(&argument1);  /* BURADA R gibi SAYILARIN BİRİMİ KONTROL EDİLECEK */
	}
	
	mixture_departure_s = rrgas*log((v-b_mixture)/v0) - dadt/(2.8284*b_mixture)*log((v-0.4142*b_mixture)/(v+2.4142*b_mixture)); 
	return mixture_ideal_gas_s + mixture_departure_s;
}

/* MIXTURE SPEED OF SOUND */

double MIXTURE_Speed_of_sound(double temp, double density, double pressure, double yi[]) /* ok */
{
	int i;
	struct thermo argument1;
	argument1.temperature = temp;
	argument1.den = density;
	argument1.press = pressure;
	for(i=0;i<N_SPECIES;i++){
	argument1.species[i] = yi[i];
	}
	double mixture_mw = MIXTURE_Molecular_Weight(yi);
		for(i=0;i<N_SPECIES;i++)
	{
		argument1.molefrac[i] = yi[i]*mixture_mw/mw[i];
	}
	argument1.parameter5 = mixture_mw;
	/*double a_mixture = a_mixture_wanderwaals(&argument1);
	double b_mixture = b_mixture_wanderwaals(&argument1);
	double dadt = a_mixture_firder(&argument1);*/
	argument1.parameter1 = a_mixture;
	argument1.parameter2 = b_mixture;
	argument1.parameter3 = dadt;
	double v = 1/density;
	double speed_of_sound;
	double cp = MIXTURE_Specific_Heat(temp, density, pressure, yi);
	double dvdt = dvdt_mixture(&argument1);
	double dvdp = dvdp_mixture(&argument1);
	double dpdv = dpdv_mixture(&argument1);
	speed_of_sound = v*pow((-1*cp/(cp + temp*dvdt*dvdt*dpdv))*(1/dvdp), 0.5);
	
	return speed_of_sound;
}

/* MIXTURE d(rho)/d(T) */

double MIXTURE_rho_t(double temp, double density, double pressure, double yi[])
{
	int i;
	struct thermo argument1;
	argument1.temperature = temp;
	argument1.den = density;
	argument1.press = pressure;
	for(i=0;i<N_SPECIES;i++){
	argument1.species[i] = yi[i];
	}
	double mixture_mw = MIXTURE_Molecular_Weight(yi);
		for(i=0;i<N_SPECIES;i++)
	{
		argument1.molefrac[i] = yi[i]*mixture_mw/mw[i];
	}
	argument1.parameter5 = mixture_mw;
	/*double a_mixture = a_mixture_wanderwaals(&argument1);
	double b_mixture = b_mixture_wanderwaals(&argument1);
	double dadt = a_mixture_firder(&argument1);*/
	argument1.parameter1 = a_mixture;
	argument1.parameter2 = b_mixture;
	argument1.parameter3 = dadt;
	
	double dvdt = dvdt_mixture(&argument1);
	double mixture_rho_t;
	mixture_rho_t = -density*density*dvdt;
	return mixture_rho_t;
}

/* MIXTURE d(rho)/d(P) */

double MIXTURE_rho_p(double temp, double density, double pressure, double yi[])
{
	int i;
	struct thermo argument1;
	argument1.temperature = temp;
	argument1.den = density;
	argument1.press = pressure;
	for(i=0;i<N_SPECIES;i++){
	argument1.species[i] = yi[i];
	}
	double mixture_mw = MIXTURE_Molecular_Weight(yi);
		for(i=0;i<N_SPECIES;i++)
	{
		argument1.molefrac[i] = yi[i]*mixture_mw/mw[i];
	}
	argument1.parameter5 = mixture_mw;
	/*double a_mixture = a_mixture_wanderwaals(&argument1);
	double b_mixture = b_mixture_wanderwaals(&argument1);
	double dadt = a_mixture_firder(&argument1);*/
	argument1.parameter1 = a_mixture;
	argument1.parameter2 = b_mixture;
	argument1.parameter3 = dadt;
	
	double mixture_rho_p;
	double dvdp = dvdp_mixture(&argument1);
	mixture_rho_p = -density*density*dvdp;
	return mixture_rho_p;
}

/* MIXTURE d(H)/d(T) */

double MIXTURE_enthalpy_t(double temp, double density, double pressure, double yi[])
{
	return MIXTURE_Specific_Heat(temp, density, pressure, yi);
}

/* MIXTURE d(H)/d(P) */

double MIXTURE_enthalpy_p(double temp, double density, double pressure, double yi[])
{
	int i;
	double v = 1/density;	
	struct thermo argument1;
	argument1.temperature = temp;
	argument1.den = density;
	argument1.press = pressure;
	for(i=0;i<N_SPECIES;i++){
	argument1.species[i] = yi[i];
	}
	double mixture_mw = MIXTURE_Molecular_Weight(yi);
		for(i=0;i<N_SPECIES;i++)
	{
		argument1.molefrac[i] = yi[i]*mixture_mw/mw[i];
	}
	argument1.parameter5 = mixture_mw;
	/*double a_mixture = a_mixture_wanderwaals(&argument1);
	double b_mixture = b_mixture_wanderwaals(&argument1);
	double dadt = a_mixture_firder(&argument1);*/
	argument1.parameter1 = a_mixture;
	argument1.parameter2 = b_mixture;
	argument1.parameter3 = dadt;
	double dvdt = dvdt_mixture(&argument1);
	return v - temp*dvdt;
}

/* MIXTURE VISCOSITY */
double MIXTURE_Viscosity(double temp, double density, double pressure, double yi[]) /* ok */
{
	double e_k[N_SPECIES][N_SPECIES] = {{151.3222,	136.2898,	126.3477,	62.9524,	278.8597,	73.1251,	63.1620,	219.2384,	191.1653	}, {136.2898,	122.7507,	113.7963,	56.6987,	251.1577,	65.8609,	56.8875,	197.4592,	172.1749	}, {126.3477,	113.7963,	105.4951,	52.5626,	232.8363,	61.0565,	52.7376,	183.0550,	159.6151	}, {62.9524,	56.6987,	52.5626,	26.1892,	116.0100,	30.4212,	26.2764,	91.2066,	79.5277	}, {278.8597,	251.1577,	232.8363,	116.0100,	513.8887,	134.7566,	116.3962,	404.0173,	352.2836	}, {73.1251,	65.8609,	61.0565,	30.4212,	134.7566,	35.3371,	30.5225,	105.9451,	92.3790	}, {63.1620,	56.8875,	52.7376,	26.2764,	116.3962,	30.5225,	26.3639,	91.5103,	79.7925	}, {219.2384,	197.4592,	183.0550,	91.2066,	404.0173,	105.9451,	91.5103,	317.6368,	276.9640	}, {191.1653,	172.1749,	159.6151,	79.5277,	352.2836,	92.3790,	79.7925,	276.9640,	241.4992	}}; 
	double k[N_SPECIES][N_SPECIES] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
	double MW[N_SPECIES][N_SPECIES] = {{16.0430,	21.3713,	20.4011,	3.5819,	16.9721,	16.0185,	1.8968,	16.5110,	23.5143	}, {21.3713,	31.9988,	29.8718,	3.7930,	23.0523,	21.3277,	1.9544,	22.2098,	37.0554	}, {20.4011,	29.8718,	28.0100,	3.7613,	21.9275,	20.3614,	1.9459,	21.1638,	34.2327	}, {3.5819,	3.7930,	3.7613,	2.0160,	3.6262,	3.5807,	1.3440,	3.6047,	3.8554	}, {16.9721,	23.0523,	21.9275,	3.6262,	18.0153,	16.9446,	1.9091,	17.4967,	25.5655	}, {16.0185,	21.3277,	20.3614,	3.5807,	16.9446,	15.9940,	1.8964,	16.4850,	23.4616	}, {1.8968,	1.9544,	1.9459,	1.3440,	1.9091,	1.8964,	1.0080,	1.9031,	1.9708	}, {16.5110,	22.2098,	21.1638,	3.6047,	17.4967,	16.4850,	1.9031,	17.0070,	24.5334	}, {23.5143,	37.0554,	34.2327,	3.8554,	25.5655,	23.4616,	1.9708,	24.5334,	44.0100	}}; 
	double sigma[N_SPECIES][N_SPECIES] = {{3.7374,	3.5581,	3.6829,	3.4787,	3.4007,	3.1477,	2.8316,	3.2905,	3.7078	}, {3.5581,	3.3873,	3.5062,	3.3117,	3.2375,	2.9966,	2.6957,	3.1326,	3.5299	}, {3.6829,	3.5062,	3.6292,	3.4279,	3.3511,	3.1018,	2.7903,	3.2425,	3.6537	}, {3.4787,	3.3117,	3.4279,	3.2378,	3.1653,	2.9298,	2.6356,	3.0627,	3.4511	}, {3.4007,	3.2375,	3.3511,	3.1653,	3.0943,	2.8641,	2.5765,	2.9940,	3.3738	}, {3.1477,	2.9966,	3.1018,	2.9298,	2.8641,	2.6510,	2.3848,	2.7713,	3.1227	}, {2.8316,	2.6957,	2.7903,	2.6356,	2.5765,	2.3848,	2.1453,	2.4930,	2.8092	}, {3.2905,	3.1326,	3.2425,	3.0627,	2.9940,	2.7713,	2.4930,	2.8970,	3.2644	}, {3.7078,	3.5299,	3.6537,	3.4511,	3.3738,	3.1227,	2.8092,	3.2644,	3.6785	}};   
	double ww[N_SPECIES][N_SPECIES] = {{0.0114,	0.0162,	0.0282,	-0.1028,	0.1777,	0.0057,	0.0057,	0.1057,	0.1182	}, {0.0162,	0.0210,	0.0330,	-0.0980,	0.1825,	0.0105,	0.0105,	0.1105,	0.1230	}, {0.0282,	0.0330,	0.0450,	-0.0860,	0.1945,	0.0225,	0.0225,	0.1225,	0.1350	}, {-0.1028,	-0.0980,	-0.0860,	-0.2170,	0.0635,	-0.1085,	-0.1085,	-0.0085,	0.0040	}, {0.1777,	0.1825,	0.1945,	0.0635,	0.3440,	0.1720,	0.1720,	0.2720,	0.2845	}, {0.0057,	0.0105,	0.0225,	-0.1085,	0.1720,	0.0000,	0.0000,	0.1000,	0.1125	}, {0.0057,	0.0105,	0.0225,	-0.1085,	0.1720,	0.0000,	0.0000,	0.1000,	0.1125	}, {0.1057,	0.1105,	0.1225,	-0.0085,	0.2720,	0.1000,	0.1000,	0.2000,	0.2125	}, {0.1182,	0.1230,	0.1350,	0.0040,	0.2845,	0.1125,	0.1125,	0.2125,	0.2250	}}; 
	double dipole_moment[5] = {0, 0, 0, 1.8, 0, 0, 0, 0, 0};
	
	double a_v[10] = {6.324, 0.00121, 5.283, 6.623, 19.745, -1.9, 24.275, 0.7972, -0.2382, 0.06863};
	double b_v[10] = {50.412, -0.001154, 254.209, 38.096, 7.630, -12.537, 3.45, 1.117, 0.0677, 0.3479};
	double c_v[10] = {-51.68, -0.006257, -168.48, -8.464, -14.354, 4.985, -11.291, 0.01235, -0.8163, 0.5926};
	double d_v[10] = {1189, 0.03728, 3989, 31.42, 31.53, -18.15, 69.35, -4.117, 4.025, -0.727};
	
	double sigma_mix = mixture_sigma(sigma, yi);
	double crit_volume_mix = mixture_critical_volume(sigma_mix);
	double e_k_mix = mixture_e_k(sigma, yi, e_k);
	double crit_temperature_mix = mixture_critical_temp(e_k_mix);
	double lp_viscosity = low_pressure_vis(temp, sigma, yi, e_k, MW, dipole_moment, ww, k);
	double dipole_moment_mix = mixture_dip_mom(sigma, yi, dipole_moment, e_k);               
	double accentric_factor_mix = mixture_acc_fac(sigma, yi, ww);                   
	double correction_factor_mix = mixture_corr_fac(yi, k);
	double dimensionless_moment_mix = (131.3*dipole_moment_mix)/(pow(crit_volume_mix*crit_temperature_mix,0.5));
	double E[10];
	unsigned int i;

	for(i=0; i<10; i++)
		{
			E[i] = a_v[i] + b_v[i]*accentric_factor_mix + c_v[i]*pow(dimensionless_moment_mix,4) + d_v[i]*correction_factor_mix;
		}
	
	double molecular_weight_mix = mixture_mol_weight(sigma, yi, e_k, MW);
	double y = (density/(molecular_weight_mix*1000)*crit_volume_mix)/6;
	double G_one = (1-0.5*y)/pow((1-y),3);
	double G_two = (E[0]*((1-exp(-E[3]*y))/y) + E[1]*G_one*exp(E[4]*y) + E[2]*G_one)/((E[0]*E[3] + E[1] + E[2]));
	double hp_corr_1 = ((1/G_two) + E[5]*y)*lp_viscosity;
	double hp_corr_2 = (36.344*1e-6*pow(molecular_weight_mix*crit_temperature_mix,0.5))/pow(crit_volume_mix,0.666667)*E[6]*pow(y,2)*G_two*exp(E[7] + E[8]*pow(temp/e_k_mix,-1) + E[9]*pow(temp/e_k_mix,-2));
	double hp_viscosity = (hp_corr_1 + hp_corr_2)*1e-1; 

	return hp_viscosity;	
}
/* CH4 = 0, O2 = 1, CO = 2, H2 = 3, H20 = 4, O = 5, H = 6, OH = 7, CO2 = 8 */
/*MIXTURE THERMAL CONDUCTIVITY */
double MIXTURE_ThermalConductivity(double temp, double density, double pressure, double yi[])  /* ok */
{
	double e_k[N_SPECIES][N_SPECIES] = {{151.3222,	136.2898,	126.3477,	62.9524,	278.8597,	73.1251,	63.1620,	219.2384,	191.1653	}, {136.2898,	122.7507,	113.7963,	56.6987,	251.1577,	65.8609,	56.8875,	197.4592,	172.1749	}, {126.3477,	113.7963,	105.4951,	52.5626,	232.8363,	61.0565,	52.7376,	183.0550,	159.6151	}, {62.9524,	56.6987,	52.5626,	26.1892,	116.0100,	30.4212,	26.2764,	91.2066,	79.5277	}, {278.8597,	251.1577,	232.8363,	116.0100,	513.8887,	134.7566,	116.3962,	404.0173,	352.2836	}, {73.1251,	65.8609,	61.0565,	30.4212,	134.7566,	35.3371,	30.5225,	105.9451,	92.3790	}, {63.1620,	56.8875,	52.7376,	26.2764,	116.3962,	30.5225,	26.3639,	91.5103,	79.7925	}, {219.2384,	197.4592,	183.0550,	91.2066,	404.0173,	105.9451,	91.5103,	317.6368,	276.9640	}, {191.1653,	172.1749,	159.6151,	79.5277,	352.2836,	92.3790,	79.7925,	276.9640,	241.4992	}}; 
	double k[N_SPECIES][N_SPECIES] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
	double MW[N_SPECIES][N_SPECIES] = {{16.0430,	21.3713,	20.4011,	3.5819,	16.9721,	16.0185,	1.8968,	16.5110,	23.5143	}, {21.3713,	31.9988,	29.8718,	3.7930,	23.0523,	21.3277,	1.9544,	22.2098,	37.0554	}, {20.4011,	29.8718,	28.0100,	3.7613,	21.9275,	20.3614,	1.9459,	21.1638,	34.2327	}, {3.5819,	3.7930,	3.7613,	2.0160,	3.6262,	3.5807,	1.3440,	3.6047,	3.8554	}, {16.9721,	23.0523,	21.9275,	3.6262,	18.0153,	16.9446,	1.9091,	17.4967,	25.5655	}, {16.0185,	21.3277,	20.3614,	3.5807,	16.9446,	15.9940,	1.8964,	16.4850,	23.4616	}, {1.8968,	1.9544,	1.9459,	1.3440,	1.9091,	1.8964,	1.0080,	1.9031,	1.9708	}, {16.5110,	22.2098,	21.1638,	3.6047,	17.4967,	16.4850,	1.9031,	17.0070,	24.5334	}, {23.5143,	37.0554,	34.2327,	3.8554,	25.5655,	23.4616,	1.9708,	24.5334,	44.0100	}}; 
	double sigma[N_SPECIES][N_SPECIES] = {{3.7374,	3.5581,	3.6829,	3.4787,	3.4007,	3.1477,	2.8316,	3.2905,	3.7078	}, {3.5581,	3.3873,	3.5062,	3.3117,	3.2375,	2.9966,	2.6957,	3.1326,	3.5299	}, {3.6829,	3.5062,	3.6292,	3.4279,	3.3511,	3.1018,	2.7903,	3.2425,	3.6537	}, {3.4787,	3.3117,	3.4279,	3.2378,	3.1653,	2.9298,	2.6356,	3.0627,	3.4511	}, {3.4007,	3.2375,	3.3511,	3.1653,	3.0943,	2.8641,	2.5765,	2.9940,	3.3738	}, {3.1477,	2.9966,	3.1018,	2.9298,	2.8641,	2.6510,	2.3848,	2.7713,	3.1227	}, {2.8316,	2.6957,	2.7903,	2.6356,	2.5765,	2.3848,	2.1453,	2.4930,	2.8092	}, {3.2905,	3.1326,	3.2425,	3.0627,	2.9940,	2.7713,	2.4930,	2.8970,	3.2644	}, {3.7078,	3.5299,	3.6537,	3.4511,	3.3738,	3.1227,	2.8092,	3.2644,	3.6785	}};   
	double ww[N_SPECIES][N_SPECIES] = {{0.0114,	0.0162,	0.0282,	-0.1028, 0.1777,	0.0057,	0.0057,	0.1057,	0.1182	}, {0.0162,	0.0210,	0.0330,	-0.0980, 0.1825, 0.0105, 0.0105, 0.1105, 0.1230	}, {0.0282,	0.0330,	0.0450,	-0.0860,	0.1945,	0.0225,	0.0225,	0.1225,	0.1350	}, {-0.1028,	-0.0980,	-0.0860,	-0.2170,	0.0635,	-0.1085,	-0.1085,	-0.0085,	0.0040	}, {0.1777,	0.1825,	0.1945,	0.0635,	0.3440,	0.1720,	0.1720,	0.2720,	0.2845	}, {0.0057,	0.0105,	0.0225,	-0.1085,	0.1720,	0.0000,	0.0000,	0.1000,	0.1125	}, {0.0057,	0.0105,	0.0225,	-0.1085,	0.1720,	0.0000,	0.0000,	0.1000,	0.1125	}, {0.1057,	0.1105,	0.1225,	-0.0085,	0.2720,	0.1000,	0.1000,	0.2000,	0.2125	}, {0.1182,	0.1230,	0.1350,	0.0040,	0.2845,	0.1125,	0.1125,	0.2125,	0.2250	}}; 
	double dipole_moment[5] = {0, 0, 0, 1.8, 0, 0, 0, 0, 0};
	
	double a_k[7] = {2.4166, -0.50924, 6.6107, 14.543, 0.79274, -5.8634, 91.089};
	double b_k[7] = {0.74824, -1.5094, 5.6207, -8.9139, 0.82019, 12.801, 128.11};
	double c_k[7] = {-0.91858, -49.991, 64.760, -5.6379, -0.69369, 9.5893, -60.841};
	double d_k[7] = {121.72, 69.983, 27.039, 74.344, 6.3173, 65.529, 523.81};
	double dipole_moment[5] = {0, 0, 0, 1.8, 0};
	double sigma_mix = mixture_sigma(sigma, yi);
	double crit_volume_mix = mixture_critical_volume(sigma_mix);
	double e_k_mix = mixture_e_k(sigma, yi, e_k);
	double crit_temperature_mix = mixture_critical_temp(e_k_mix);
	double molecular_weight_mix = mixture_mol_weight(sigma, yi, e_k, MW);
	double accentric_factor_mix = mixture_acc_fac(sigma, yi, ww);
	double correction_factor_mix = mixture_corr_fac(yi, k);                              
	double dipole_moment_mix = mixture_dip_mom(sigma, yi, dipole_moment, e_k);               
	double cv = ideal_specific_volume(temp, yi);
	double Z = 2 + 10.5*pow((temp/crit_temperature_mix),2);
	double alpha = (cv/8.314) - 1.5;
	double beta = 0.7862 - 0.7109*accentric_factor_mix + 1.3168*accentric_factor_mix*accentric_factor_mix;
	double eta = 1 + alpha*((0.215 + 0.028288*alpha - 1.061*beta + 0.26665*Z)/(0.6366 + beta*Z + 1.061*alpha*beta));
	double y = (density/(molecular_weight_mix*1000)*crit_volume_mix)/6;
	double G_one = (1-0.5*y)/pow((1-y),3);
	double dimensionless_moment_mix = (131.3*dipole_moment_mix)/(pow(crit_volume_mix*crit_temperature_mix,0.5));
	double B[7];
	unsigned int i;

	for(i=0; i<7; i++)
	{
		B[i] = a_k[i] + b_k[i]*accentric_factor_mix + c_k[i]*pow(dimensionless_moment_mix,4) + d_k[i]*correction_factor_mix;
	}
	
	double lp_viscosity = low_pressure_vis(temp, sigma, yi, e_k, MW, dipole_moment, ww, k);
	double GG_two = (  (B[0]/y)*(1-exp(-1*B[3]*y)) + B[1]*G_one*exp(B[4]*y) + B[2]*G_one)/((B[0]*B[3] + B[1] + B[2]));
	double q = 3.586*1e-3*pow(crit_temperature_mix*1000/molecular_weight_mix,0.5)/pow(crit_volume_mix,0.6667);
	double thermal_conductivity = (31.2*1e2*lp_viscosity*eta/molecular_weight_mix)*(1/GG_two + B[5]*y) + q*B[6]*y*y*pow(temp/crit_temperature_mix,0.5)*GG_two;
		
	return thermal_conductivity;
}

/*******************************************************************/
/* Real Gas Associated Coefficient Functions :                     */
/* 0 = CH4, 1 = O2, 2 = CO2, 3 = H2O, 4 = N2                       */
/*******************************************************************/
/* a_mixture */
double a_mixture_wanderwaals(struct thermo *argument )
{
	double xi[N_SPECIES]; 
	double sum = 0;
	double mixture_mw;
	double a_mix = 0;
	unsigned int i,j;
	double sss;
	
	
	
	for(i=0;i<N_SPECIES;i++)
		for(j=0;j<N_SPECIES;j++)
	{
		sss = argument->molefrac[i]*argument->molefrac[j]*pow(a0[i]*a0[j],0.5)*pow(pow(1+m[i]*(1-pow(argument->temperature/t_crit[i],0.5)),2)*pow(1+m[j]*(1-pow(argument->temperature/t_crit[j],0.5)),2),0.5);
		a_mix = a_mix + sss;
	}
	
	return a_mix;
}

double b_mixture_wanderwaals(struct thermo *argument)
{
	unsigned int i;
	double b_mix=0;
	double sum = 0;
	double mixture_mw;
	double xi[5];
		
	for(i=0;i<N_SPECIES;i++)
	{
		/*xi[i] = argument->species[i]*mixture_mw/mw[i];*/
		b_mix = b_mix + argument->molefrac[i]*b[i];
	}
	
	
	return b_mix;
}
/* dvdt mixture */
double dvdt_mixture(struct thermo *argument)
{
	double dvdt;
	double v;
	v = 1/argument->den;
	
	unsigned int i,j;
	double da1,da2,da3;
	double a1,a2,a3; 
	double rrgas = (R/argument->parameter5);
	
	a1 = argument->parameter2 - rrgas*argument->temperature/argument->press;
	a2 = -3*pow(argument->parameter2,2) - (2*argument->temperature*argument->parameter2*rrgas/argument->press) + argument->parameter1/argument->press;
	a3 = pow(argument->parameter2,3) + (2*rrgas*argument->parameter2*argument->parameter2/argument->press) - argument->parameter1*argument->parameter2/argument->press;
	
	da1 = -rrgas/argument->press;
	da2 = -2*argument->parameter2*rrgas/argument->press + argument->parameter3/argument->press;
	da3 = rrgas*pow(argument->parameter2,2)/argument->press - (argument->parameter2/argument->press)*argument->parameter3;
	dvdt = -1*(da1*pow(v,2) + da2*v + da3)/(3*pow(v,2) + 2*v*a1 + a2);
	
	return dvdt;	
}
/* dvdp mixture */
double dvdp_mixture(struct thermo *argument)
{
	double dvdp;
	double v;
	v = 1/argument->den;
	
	unsigned int i,j;
	double da1,da2,da3;
	double a1,a2,a3;
	double rrgas = (R/argument->parameter5);
	
	a1 = argument->parameter2 - rrgas*argument->temperature/argument->press;
	a2 = -3*pow(argument->parameter2,2) - (2*argument->temperature*argument->parameter2*rrgas/argument->press) + argument->parameter1/argument->press;
	a3 = pow(argument->parameter2,3) + (2*rrgas*argument->parameter2*argument->parameter2/argument->press) - argument->parameter1*argument->parameter2/argument->press;
	
	da1 = rrgas*argument->temperature/pow(argument->press,2);
	da2 = (2*argument->temperature*argument->parameter2*rrgas)/pow(argument->press,2) - argument->parameter1/pow(argument->press,2);
	da3 = -1*rrgas*argument->temperature*pow(argument->parameter2,2)/pow(argument->press,2) + argument->parameter1*argument->parameter2/pow(argument->press,2);
	dvdp = -1*(da1*pow(v,2) + da2*v + da3)/(3*pow(v,2) + 2*v*a1 + a2);
	
	return dvdp;
}

/* dadt mixture */
double a_mixture_firder(struct thermo *argument) 
{
	double dadtt = 0;
	unsigned int i,j;
	
	double tcrit_ij,vcrit_ij,pcrit_ij,zcrit_ij,c_ij,aij,w_ij;
	double rrgas = (R/argument->parameter5);
	
	for(i=0;i<N_SPECIES;i++)
	{
		for(j=0;j<N_SPECIES;j++)
		{
			tcrit_ij = pow(t_crit[i]*t_crit[j],0.5);
			vcrit_ij = 0.125*(pow(pow(v_crit[i],0.333333) + pow(v_crit[j],0.333333),3));
			zcrit_ij = (z_crit[i] + z_crit[j])*0.5;
			pcrit_ij = zcrit_ij*rrgas*tcrit_ij/vcrit_ij;
			w_ij = (w[i]+w[j])*0.5;
			c_ij = 0.37464 + 1.52226*w_ij - 0.26992*w_ij*w_ij;
			aij = 0.457235*rrgas*rrgas*tcrit_ij*tcrit_ij*pow(1+c_ij*(1-pow(argument->temperature/tcrit_ij,0.5)),2)/pcrit_ij;
			dadtt = dadtt + (-1/argument->temperature)*argument->molefrac[i]*argument->molefrac[j]*aij*(c_ij*pow(argument->temperature/tcrit_ij,0.5))/(1+c_ij*(1-pow(argument->temperature/tcrit_ij,0.5)));
			
		}
	}
	return dadtt;
}

/* da2dat2 mixture */
double a_mixture_secder(struct thermo *argument) /* ok */
{
	double da2dt2 = 0;
	unsigned int i,j;

	double tcrit_ij,vcrit_ij,pcrit_ij,zcrit_ij,c_ij,aij,w_ij;
	double rrgas = (R/argument->parameter5);
	
	for(i=0;i<N_SPECIES;i++)
	{
		for(j=0;j<N_SPECIES;j++)
		{
			tcrit_ij = pow(t_crit[i]*t_crit[j],0.5);
			vcrit_ij = 0.125*(pow(pow(v_crit[i],0.333333) + pow(v_crit[j],0.333333),3));
			zcrit_ij = (z_crit[i] + z_crit[j])*0.5;
			pcrit_ij = zcrit_ij*rrgas*tcrit_ij/vcrit_ij;
			w_ij = (w[i]+w[j])*0.5;
			c_ij = 0.37464 + 1.52226*w_ij - 0.26992*w_ij*w_ij;
			da2dt2 = da2dt2 + (0.457235*rrgas*rrgas/(2*argument->temperature))*argument->molefrac[i]*argument->molefrac[j]*(1-c_ij)*(tcrit_ij/pcrit_ij)*pow(tcrit_ij/argument->temperature,0.5); /* 1-c_ij yerine c_ij*(1+c_ij) yazılabilir */
		}
	}	
	
	return da2dt2;
}

/*dpdt mixture */
double dpdt_mixture(struct thermo *argument) /* ok */
{
	double dpdt;
	double v;
	v = 1/argument->den;
	double rrgas = (R/argument->parameter5);
	dpdt = rrgas/(v-argument->parameter2) - argument->parameter3/(v*v + 2*argument->parameter2*v - argument->parameter2*argument->parameter2);
	return dpdt;
}

/* dpdv mixture */
double dpdv_mixture(struct thermo *argument) /* ok */
{
	double dpdv;
	double v;
	v = 1/argument->den;
	double rrgas = (R/argument->parameter5);
	dpdv = -1*rrgas*argument->temperature/pow(v-argument->parameter2,2) + 2*argument->parameter1*(v + argument->parameter2)/pow(v*v + 2*argument->parameter2*v - argument->parameter2*argument->parameter2,2);
	return dpdv;
}

/******************************************************************************/
/* Species Ideal Thermodynamic Property functions                             */
/* CH4 = 0, O2 = 1, CO = 2, H2 = 3, H20 = 4, O = 5, H = 6, OH = 7, CO2 = 8    */ 
/******************************************************************************/

/******************* IDEAL GAS SPECIFIC HEAT COEFS ******************/
double RK_Ideal_cp(struct thermo *argument)
{
	double cp;
	double t = argument->temperature*0.001;
	
	if (argument->parameter4==0) { /* CH4 */
		if(argument->temperature<=300){
			cp=(8.314*(4.568 - 0.008976*argument->temperature + 0.00003631*pow(argument->temperature,2) - 0.00000003407*pow(argument->temperature,3) + 0.00000000001091*pow(argument->temperature,4)));
		}
		/*Message("cp is %g",cp*1000/mw[i]);*/
		else if (argument->temperature>300 && argument->temperature<1300) {
			cp = (-0.703029 + 108.4773*t - 42.52157*pow(t,2) + 5.862788*pow(t,3) + 0.678565*pow(t,-2));
		}
		else {
			cp = (85.81217 + 11.26467*t - 2.114146*pow(t,2) + 0.138190*pow(t,3) - 26.42221*pow(t,-2));
		}
	}
	if (argument->parameter4==1) { /* O2 */
		if(argument->temperature<=300){
			cp=(8.314*(3.630 - 0.001794*argument->temperature + 0.00000658*pow(argument->temperature,2) - 0.00000000601*pow(argument->temperature,3) + 0.00000000000179*pow(argument->temperature,4)));
		}
		/*Message("cp is %g",cp*1000/mw[i]);*/
		else if (argument->temperature>300 &&argument->temperature<=700) {
		    cp = (31.32234 - 20.23531*t + 57.86644*pow(t,2) - 36.50624*pow(t,3) - 0.007374*pow(t,-2));
		}
		else if(argument->temperature>700 && argument->temperature <=2000){
			cp = (30.03235 + 8.772972*t - 3.988133*pow(t,2) + 0.788313*pow(t,3) - 0.741599*pow(t,-2));
		}
		else {
			cp = (20.91111 + 10.72071*t - 2.020498*pow(t,2) + 0.146449*pow(t,3) + 9.245722*pow(t,-2));
		}
	}
	if (argument->parameter4==2) { /* CO */
		if(argument->temperature<=300){
			cp=(8.314*(3.912 - 3.913*argument->temperature + 0.00001182*pow(argument->temperature,2) - 0.00000001302*pow(argument->temperature,3) + 0.515*1e-11*pow(argument->temperature,4)));
		}
		else if (argument->temperature>300 &&argument->temperature<=1300) {
			cp = (25.56 + 6.09613*t + 4.054656*pow(t,2) - 2.671301*pow(t,3) + 0.131021*pow(t,-2));
		}
		else {
			cp = (35.1507 + 1.300095*t - 0.205921*pow(t,2) 0.01355*pow(t,3) - 3.28278*pow(t,-2));
		}
	}
	if (argument->parameter4==3) { /* H2 */
		if(argument->temperature<=300){
			cp=(8.314*(2.883 - 3.681*argument->temperature - 0.772*1e-5*pow(argument->temperature,2) + 0.692*1e-8*pow(argument->temperature,3) - 0.213*1e-11*pow(argument->temperature,4)));
		}
		else if (argument->temperature>300 &&argument->temperature<=1000) {
			cp = (33.066 - 11.3634*t + 11.4328*pow(t,2) - 2.77287*pow(t,3) - 0.15855*pow(t,-2));
		}
		else if(argument->temperature>1000 &&argument->temperature<=2500){
			cp = (18.563083 + 12.25735*t - 2.85978*pow(t,2) + 0.26823*pow(t,3) + 1.978*pow(t,-2));
		}
		else {
			cp = (43.41356 -  4.293079*t + 1.272428*pow(t,2) - 0.096876*pow(t,3) - 20.53386*pow(t,-2));
		}
	}
	if(argument->parameter4==4) {	/* H2O */
		if(argument->temperature<=300){
			cp=(8.314*(4.395 - 0.004186*argument->temperature + 0.00001405*pow(argument->temperature,2) - 0.00000001564*pow(argument->temperature,3) + 0.00000000000632*pow(argument->temperature,4)));
		}
		/*Message("cp is %g",cp*1000/mw[i]);*/
		else if(argument->temperature>300 && argument->temperature<1700) {
			cp = (30.09200 + 6.832514*t + 6.793435*pow(t,2) - 2.534480*pow(t,3) + 0.082139*pow(t,-2));
		}
		else {
			cp = (41.96426 + 8.622053*t - 1.499780*pow(t,2) + 0.098119*pow(t,3) - 11.15764*pow(t,-2));
		}
	}
	if(argument->parameter4==5) { /* O */
		if(argument->temperature<=1000){
			cp = 1000*8.314*(2.542059 -  0.02755*1e-3*argument->temperature - 0.031028*1e-7*pow(argument->temperature,2) + 0.04551067*1e-10*pow(argument->temperature,3) - 0.043680*1e-14*pow(argument->temperature,4));
		}
		else {
			cp = 1000*8.314*(2.9464 - 0.1638166*1e-2*argument->temperature + 0.0242103*1e-4*pow(argument->temperature,2) - 0.0.1602843*1e-8*pow(argument->temperature,3) + 0.0389069*1e-11*pow(argument->temperature,4));
		}
	}	
	if(argument->parameter4==6) { /* H */
			cp = 1000*2.5*8.314;
	}
	if(argument->parameter4==7) { /* OH */
		if(argument->temperature<=1300){
			cp = (32.277 - 11.3629*t + 13.60545*pow(t,2) - 3.846486*pow(t,3) - 0.001335*pow(t,-2));
		}
		else {
			cp = (28.7470 + 4.714489*t - 0.814725*pow(t,2) + 0.054748*pow(t,3) - 2.747829*pow(t,-2));
		}
	}
	if(argument->parameter4==8) { /* CO2 */
		if(argument->temperature<=300){
			cp=(8.314*(3.259 + 0.001356*argument->temperature + 0.00001502*pow(argument->temperature,2) - 0.00000002374*pow(argument->temperature,3) + 0.00000000001056*pow(argument->temperature,4)));
		}
		/*Message("cp is %g",cp*1000/mw[i]);*/
		else if(argument->temperature>300 && argument->temperature<1200) {
			cp = (24.99735 + 55.18696*t - 33.69137*pow(t,2) + 7.948387*pow(t,3) - 0.136638*pow(t,-2));
		}
		else {
			cp = (58.16639 + 2.720074*t - 0.492289*pow(t,2) + 0.038844*pow(t,3) - 6.447293*pow(t,-2));
		}
	}
	return cp*1000/mw[argument->parameter4];
}

/******************* IDEAL GAS ENTHALPY COEFS ******************/

double RK_Ideal_Enthalpy(struct thermo *argument)
{
	double h;
	double t = argument->temperature/1000;
	
	if (argument->parameter4==0) { /* CH4 */
		if (argument->temperature<1300) {
			h = (-0.703029*t + 108.4773*t*t*0.5 - 42.52157*pow(t,3)*0.3333 + 5.862788*pow(t,4)*0.25 - 0.678565*pow(t,-1) - 76.84376 + 74.87310);
		}
		else {
			h = (85.81217*t + 11.26467*t*t*0.5  - 2.114146*pow(t,3)*0.3333 + 0.138190*pow(t,4)*0.25 + 26.42221*pow(t,-1) - 153.5327 + 74.87310);
		}
	}
	if (argument->parameter4==1) {
		if (argument->temperature<=700) { /* O2 */
		    h = (31.32234*t - 20.23531*t*t*0.5  + 57.86644*pow(t,3)*0.3333 - 36.50624*pow(t,4)*0.25 + 0.007374*pow(t,-1) - 8.903471);
		}
		else if(argument->temperature>700 && argument->temperature <=2000){
			h = (30.03235*t + 8.772972*t*t*0.5  - 3.988133*pow(t,3)*0.3333 + 0.788313*pow(t,4)*0.25 + 0.741599*pow(t,-1) - 11.32468);
		}
		else {
			h = (20.91111*t + 10.72071*t*t*0.5  - 2.020498*pow(t,3)*0.3333 + 0.146449*pow(t,4)*0.25 - 9.245722*pow(t,-1) + 5.337651);
		}
	}
	if (argument->parameter4==2) { /* CO */
		if (argument->temperature>300 &&argument->temperature<=1300) {
			h = (25.56*t + 6.09613*t*t*0.5 + 4.054656*pow(t,3)*0.3333 - 2.671301*pow(t,4)*0.25 - 0.131021*pow(t,-1) - 118.0089 + 110.5271) ;
		}
		else {
			h = (35.1507*t + 1.300095*t*t*0.5 - 0.205921*pow(t,3)*0.3333 0.01355*pow(t,4)*0.25 + 3.28278*pow(t,-1) - 127.8375 + 110.5271);
		}
	}
	if (argument->parameter4==3) { /* H2 */
		if (argument->temperature>300 &&argument->temperature<=1000) {
			h = (33.066*t - 11.3634*t*t*0.5 + 11.4328*pow(t,3)*0.3333 - 2.77287*pow(t,4)*0.25 + 0.15855*pow(t,-1) - 9.980797);
		}
		else if(argument->temperature>1000 &&argument->temperature<=2500){
			h = (18.563083*t + 12.25735*t*t*0.5 - 2.85978*pow(t,3)*0.3333 + 0.26823*pow(t,4)*0.25 - 1.978*pow(t,-1) - 1.147438);
		}
		else {
			h = (43.41356*t -  4.293079*t*t*0.5 + 1.272428*pow(t,3)*0.3333 - 0.096876*pow(t,4)*0.25 + 20.53386*pow(t,-1) - 38.515158);
		}
	}
	if(argument->parameter4==4) { /* H2O */		
		if(argument->temperature<1700) {
			h = (30.09200*t + 6.832514*t*t*0.5  + 6.793435*pow(t,3)*0.3333 - 2.534480*pow(t,4)*0.25 - 0.082139*pow(t,-1) - 250.8810 + 241.8264);
		}
		else {
			h = (41.96426*t + 8.622053*t*t*0.5  - 1.499780*pow(t,3)*0.3333 + 0.098119*pow(t,4)*0.25 + 11.15764*pow(t,-1) - 272.1797 + 241.8264);
		}
	}
	if(argument->parameter4==5) { /* O */
		if(argument->temperature<=1000){
			h = 1000*8.314*argument->temperature*(2.542059 -  0.02755*1e-3*argument->temperature*0.5 - 0.031028*1e-7*pow(argument->temperature,2)*0.3333 + 0.04551067*1e-10*pow(argument->temperature,3)*0.25 - 0.043680*1e-14*pow(argument->temperature,4)*0.2 + 0.02923*1e6/argument->temperature);
		}
		else {
			h = 1000*8.314*argument->temperature*(2.9464 - 0.1638166*1e-2*argument->temperature*0.5 + 0.0242103*1e-4*pow(argument->temperature,2)*0.3333 - 0.0.1602843*1e-8*pow(argument->temperature,3)*0.25 + 0.0389069*1e-11*pow(argument->temperature,4)*0.2 + 0.02914*1e6/argument->temperature);
		}
	}
	if(argument->parameter4==6) { /* H */
			h = 1000*8.314*argument->temperature(2.5 + 0.025471*1e6/argument->temperature);
	}
	if(argument->parameter4==7) { /* OH */
		if(argument->temperature<=1300){
			h = (32.277*t - 11.3629*t*t*0.5 + 13.60545*pow(t,3)*0.3333 - 3.846486*pow(t,4)*0.25 + 0.001335*pow(t,-1) + 29.75113 - 38.98);
		}
		else {
			h = (28.7470*t + 4.714489*t*t*0.5 - 0.814725*pow(t,3)*0.3333 + 0.054748*pow(t,4)*0.25 + 2.747829*pow(t,-1) + 26.41439 - 38.98);
		}
	}
	if(argument->parameter4==8) { /* C02 */
		if(argument->temperature<1200) {
			h = (24.99735*t + 55.18696*t*t*0.5  - 33.69137*pow(t,3)*0.3333 + 7.948387*pow(t,4)*0.25 + 0.136638*pow(t,-1) - 403.6075 + 393.5224);
		}
		else {
			h = (58.16639*t + 2.720074*t*t*0.5  - 0.492289*pow(t,3)*0.3333 + 0.038844*pow(t,4)*0.25 + 6.447293*pow(t,-1) - 425.9186 + 393.5224);
		}
	}
	
	return h*1e6/(mw[argument->parameter4]);
	
}

/******************* IDEAL GAS ENTROPY COEFS ******************/

double RK_Ideal_Entropy(struct thermo *argument)
{
	double s;
	double t = argument->temperature/1000;
	
	if (argument->parameter4==0) { /* CH4 */
		if (argument->temperature<1300) {
			s = (-0.703029*log(t) + 108.4773*t - 42.52157*pow(t,2)*0.5  + 5.862788*pow(t,3)*0.3333 - 0.678565*pow(2*t*t,-1) + 158.7163);
		}
		else {
			s = (85.81217*log(t) + 11.26467*t - 2.114146*pow(t,2)*0.5  + 0.138190*pow(t,3)*0.3333 + 26.42221*pow(2*t*t,-1) + 224.4143);
		}
	}
	if (argument->parameter4==1) { /* O2 */
		if (argument->temperature<=700) {
		    s = (31.32234*log(t) - 20.23531*t + 57.86644*pow(t,2)*0.5  - 36.50624*pow(t,3)*0.3333 + 0.007374*pow(2*t*t,-1) + 246.7945);
		}
		else if(argument->temperature>700 && argument->temperature <=2000){
			s = (30.03235*log(t) + 8.772972*t - 3.988133*pow(t,2)*0.5  + 0.788313*pow(t,3)*0.3333 + 0.741599*pow(2*t*t,-1) + 236.1663);
		}
		else {
			s = (20.91111*log(t) + 10.72071*t - 2.020498*pow(t,2)*0.5  + 0.146449*pow(t,3)*0.3333 - 9.245722*pow(2*t*t,-1) + 237.6185);
		}
	}
	if (argument->parameter4==2) { /* CO */
		if (argument->temperature>300 &&argument->temperature<=1300) {
			s = (25.56*log(t) + 6.09613*t + 4.054656*pow(t,2)*0.5 - 2.671301*pow(t,3)*0.3333 - 0.131021*pow(2*t*t,-1) + 227.3665);
		}
		else {
			s = (35.1507*log(t) + 1.300095*t - 0.205921*pow(t,2)*0.5 0.01355*pow(t,3)*0.3333 + 3.28278*0.5*pow(2*t*t,-1) + 231.7120);
		}
	}
	if (argument->parameter4==3) { /* H2 */
		if (argument->temperature>300 &&argument->temperature<=1000) {
			s = (33.066*log(t) - 11.3634*t + 11.4328*pow(t,2)*0.5 - 2.77287*pow(t,3)*0.3333 + 0.15855*pow(2*t*t,-1) + 172.707974);
		}
		else if(argument->temperature>1000 &&argument->temperature<=2500){
			s = (18.563083*log(t) + 12.25735*t - 2.85978*pow(t,2)*0.5 + 0.26823*pow(t,3)*0.3333 - 1.978*pow(2*t*t,-1) + 156.288133);
		}
		else {
			s = (43.41356*log(t) -  4.293079*t + 1.272428*pow(t,2)*0.5 - 0.096876*pow(t,3)*0.3333 - 20.53386*pow(2*t*t,-1) + 162.081354);
		}
	}
	if(argument->parameter4==4) { /* h20 */
		if(argument->temperature<1700) {
			s = (30.09200*log(t) + 6.832514*t + 6.793435*pow(t,2)*0.5  - 2.534480*pow(t,3)*0.3333 - 0.082139*pow(2*t*t,-1) + 223.3967);
		}
		else {
			s = (41.96426*log(t) + 8.622053*t - 1.499780*pow(t,2)*0.5  + 0.098119*pow(t,3)*0.3333 + 11.15764*pow(2*t*t,-1) + 219.7809);
		}
	}
	if(argument->parameter4==5) { /* O */
		if(argument->temperature<=1000){
			s = 1000*8.314*(2.542059*log(argument->temperature) -  0.02755*1e-3*argument->temperature - 0.031028*1e-7*pow(argument->temperature,2)*0.5 + 0.04551067*1e-10*pow(argument->temperature,3)*0.3333 - 0.043680*1e-14*pow(argument->temperature,4)*0.25 + 4.9203);
		}
		else {
			s = 1000*8.314*(2.9464*log(argument->temperature) - 0.1638166*1e-2*argument->temperature + 0.0242103*1e-4*pow(argument->temperature,2)*0.5 - 0.1602843*1e-8*pow(argument->temperature,3)*0.3333 + 0.0389069*1e-11*pow(argument->temperature,4)*0.25 + 2.964);
		}
	}
	if(argument->parameter4==6) { /* H */
			s = 1000*8.314*(2.5*log(argument->temperature) - 0.460);
	}
	if(argument->parameter4==7) { /* OH */
		if(argument->temperature<=1300){
			s = (32.277*log(t) - 11.3629*t + 13.60545*pow(t,2)*0.5 - 3.846486*pow(t,3)*0.3333 + 0.001335*pow(2*t*t,-1) + 225.5783);
		}
		else {
			s = (28.7470*log(t) + 4.714489*t - 0.814725*pow(t,2)*0.5 + 0.054748*pow(t,3)*0.3333 + 2.747829*pow(2*t*t,-1) + 214.1166);
		}
	}
	if(argument->parameter4==2) { /* co2 */
		if(argument->temperature<1200) {
			s = (24.99735*log(t) + 55.18696*t - 33.69137*pow(t,2)*0.5  + 7.948387*pow(t,3)*0.3333 + 0.136638*pow(2*t*t,-1) + 228.2431);
		}
		else {
			s = (58.16639*log(t) + 2.720074*t - 0.492289*pow(t,2)*0.5  + 0.038844*pow(t,3)*0.3333 + 6.447293*pow(2*t*t,-1) + 263.6125);
		}
	}	

	return s*1000/mw[argument->parameter4];	
}


/*******************************************************************/
/* VISCOSITY AND THERMAL CONDUCTIVITY FUNCTIONS                    */
/* 0 = CH4, 1 = O2, 2 = CO2, 3 = H2O, 4 = N2                       */
/*******************************************************************/

/* COLLISION INTEGRAL */
double collision_integral(double temp, double e_k_mixture)  /* G,W, vs gibi değerlerde gelebilir ama ok diyelim */
{
double A, B, C, D, E, F, T_STAR, CI ;
A = 1.16145;
B = 0.14874;
C = 0.52487;
D = 0.77320;
E = 2.16178;
F = 2.43787;	
T_STAR = temp/e_k_mixture;	
CI = (A*pow(T_STAR,-B)) + C*exp(-D*T_STAR) + E*exp(-F*T_STAR);
return CI;	
}

/* 	MIXTURE SIGMA FUNCTION */

double mixture_sigma(double sig[N_SPECIES][N_SPECIES], double yi[])  /* ok */
{
double mix_sum = 0;
unsigned int i,j;
for(i=0; i<N_SPECIES; i++)
{
	for(j=0; j<N_SPECIES; j++)
	{
		mix_sum = mix_sum + yi[i]*yi[j]*pow(sig[i][j],3);
	}
}
return pow(mix_sum, 0.33333333);
}

/* MIXTURE CRITICAL VOLUME FUNCTION*/
double mixture_critical_volume(double mixture_sigma)  /* ok, cm3/mol olcak */
{
double mix_crit_volume;
mix_crit_volume = pow((mixture_sigma*1.2361),3);
return mix_crit_volume;
}

/* MIXTURE CRITICAL TEMPERATURE FUNCTION*/  /*ok */

double mixture_critical_temp(double mix_e_k)
{
double mix_crit_temp;
mix_crit_temp = 1.2593*mix_e_k;
return mix_crit_temp;
}

/* MIXTURE E/K FUNCTION */  
double mixture_e_k(double sig[N_SPECIES][N_SPECIES], double yi[], double e_k[N_SPECIES][N_SPECIES]) /* ok */
{
unsigned int i,j;
double mix_e_k = 0;
double sigma_cube = mixture_sigma(sig, yi);
sigma_cube = pow(sigma_cube,3);

for(i=0; i<N_SPECIES; i++)
{
	for(j=0; j<N_SPECIES; j++)
	{
		mix_e_k = mix_e_k + yi[i]*yi[j]*e_k[i][j]*pow(sig[i][j],3);
	}
}
return mix_e_k/sigma_cube;
}

/* MIXTURE MOLECULAR WEIGHT FUNCTION */

double mixture_mol_weight(double sigma[N_SPECIES][N_SPECIES], double yi[], double e_k[N_SPECIES][N_SPECIES], double MW[N_SPECIES][N_SPECIES])  /* ok */
{
unsigned int i,j;
double mix_mol_weight = 0;
double sigma_cube = mixture_sigma(sigma, yi);
double sigma_square = pow(sigma_cube,2);
double mix_e_k = mixture_e_k(sigma, yi, e_k);

for(i=0; i<N_SPECIES; i++)
{
	for(j=0; j<N_SPECIES; j++)
	{
		mix_mol_weight = mix_mol_weight + yi[i]*yi[j]*e_k[i][j]*pow(sigma[i][j],2)*pow(MW[i][j],0.5);
	}
}
return pow((mix_mol_weight/(sigma_square*mix_e_k)),2); 
}

/* MIXTURE ACCENTRIC FACTOR FUNCTION */  

double mixture_acc_fac(double sigma[N_SPECIES][N_SPECIES], double yi[], double ww[N_SPECIES][N_SPECIES])  /*ok */
{
unsigned int i,j;
double mix_acc_fac = 0;
double sigma_cube = mixture_sigma(sigma, yi);
sigma_cube = pow(sigma_cube,3);
for(i=0; i<N_SPECIES; i++)
{
	for(j=0; j<N_SPECIES; j++)
	{
		mix_acc_fac = mix_acc_fac + yi[i]*yi[j]*ww[i][j]*pow(sigma[i][j],3);
	}
}
return mix_acc_fac/sigma_cube;
}

/* MIXTURE DIPOLE MOMENT FUNCTION */

double mixture_dip_mom(double sigma[N_SPECIES][N_SPECIES], double yi[], double dipole_moment[], double e_k[N_SPECIES][N_SPECIES])  /* ok */
{
unsigned int i,j;
double mix_dipole_moment_fourth = 0;
double sigma_cube = mixture_sigma(sigma, yi);
sigma_cube = pow(sigma_cube,3);
double mix_e_k = mixture_e_k(sigma, yi, e_k);
for(i=0; i<N_SPECIES; i++)
{
	for(j=0; j<N_SPECIES; j++)
	{
		mix_dipole_moment_fourth = mix_dipole_moment_fourth + (yi[i]*yi[j]*pow(dipole_moment[i],2)*pow(dipole_moment[j],2))/(pow(sigma[i][j],3)*e_k[i][j]);
	}
}
mix_dipole_moment_fourth = mix_dipole_moment_fourth*sigma_cube*mix_e_k;
return pow(mix_dipole_moment_fourth,0.25);
}

/* MIXTURE SPECIES CORRECTION FACTOR FUNCTION */
	
double mixture_corr_fac(double yi[], double k[N_SPECIES][N_SPECIES])  /* ok */
{
unsigned int i,j;
double mix_corr_fac = 0;

for(i=0; i<N_SPECIES; i++)
{
	for(j=0; j<N_SPECIES; j++)
	{
		mix_corr_fac = mix_corr_fac + yi[i]*yi[j]*k[i][j];
	}
}
return mix_corr_fac;
}

double low_pressure_vis(double temp, double sigma[N_SPECIES][N_SPECIES], double yi[], double e_k[N_SPECIES][N_SPECIES], double MW[N_SPECIES][N_SPECIES], double dipole_moment[], double ww[N_SPECIES][N_SPECIES], double k[N_SPECIES][N_SPECIES])
{
double sigma_value = mixture_sigma(sigma, yi);                                        /* mixture sigma */
double e_k_mix = mixture_e_k(sigma, yi, e_k);                                   /* mixture e_k */
double crit_volume_mix = mixture_critical_volume(sigma_value);                                /* mixture critical volume*/
double crit_temperature_mix = mixture_critical_temp(e_k_mix);  /* mixture critical temperature */
double dipole_moment_mix = mixture_dip_mom(sigma, yi, dipole_moment, e_k);               /* mixture dipole moment */
double molecular_weight_mix = mixture_mol_weight(sigma, yi, e_k, MW);      /* mixture molecular weight*/
double accentric_factor_mix = mixture_acc_fac(sigma, yi, ww);                   /* mixture accentric factor*/
double correction_factor_mix = mixture_corr_fac(yi, k);                               /* mixture correction factor */
double dimensionless_moment_mix = (131.3*dipole_moment_mix)/(pow((crit_volume_mix*crit_temperature_mix),0.5));
double CI = collision_integral(temp, e_k_mix);                                                /* collision integral*/ 
double F = 1 - 0.2756*accentric_factor_mix + 0.059035*pow(dimensionless_moment_mix,4) + correction_factor_mix;

return 4.0785*1e-5*F*pow(molecular_weight_mix*temp,0.5)/(pow(crit_volume_mix,0.666667)*CI);
}

double ideal_specific_volume(double temp, double yi[]) /* ok */
{
	unsigned int i;
	double specific_volume = 0;
	struct thermo argument1;
	argument1.temperature = temp;
	for(i=0;i<N_SPECIES;i++)
	{
		argument1.parameter4 = i;
		specific_volume = specific_volume + yi[i]*((RK_Ideal_cp(&argument1)*mw[i]*0.001) - 8.3145);
	}
	
	return specific_volume;
}
/* summation term in partial specific volume and partial enthalpy */
double species_sum(struct thermo *argument)
{
	unsigned int j;
	/*double sum = 0;*/
	double species_summation = 0;
	/*double mixture_mw;
	double xi[5];
	
	
	for(j=0;j<N_SPECIES;j++)
	{
		sum = sum + yi[j]/mw[j];
	}
	
	mixture_mw = 1/sum;
	
	for(j=0;j<N_SPECIES;j++)
	{
		xi[j] = yi[j]*mixture_mw/mw[j];
	}*/
	
	for(j=0;j<N_SPECIES;j++)
		
		{
			species_summation = species_summation + argument->molefrac[j]*pow(a0[argument->parameter4]*a0[j],0.5)*pow(pow(1+m[argument->parameter4]*(1-pow(argument->temperature/t_crit[argument->parameter4],0.5)),2)*pow(1+m[j]*(1-pow(argument->temperature/t_crit[j],0.5)),2),0.5);
		}
		return species_summation;
}
/* summation term in species enthaly that contains derivative of aij */
double species_sum_two(int i, double temp, double yi[])
{
	unsigned int j;
	double sum = 0;
	double mixture_mw;
	double xi[5];
	double species_summation = 0;
	double tcrit_ij,vcrit_ij,pcrit_ij,zcrit_ij,c_ij,aij,w_ij;
	double rrgas = (R/MIXTURE_Molecular_Weight(yi));
	for(j=0;j<N_SPECIES;j++)
	{
		sum = sum + yi[j]/mw[j];
	}
	
	mixture_mw = 1/sum;
	
	for(j=0;j<N_SPECIES;j++)
	{
		xi[j] = yi[j]*mixture_mw/mw[j];
	}
	
	for(j=0;j<N_SPECIES;j++)
		{
			tcrit_ij = pow(t_crit[i]*t_crit[j],0.5);
			vcrit_ij = 0.125*(pow(pow(v_crit[i],0.333) + pow(v_crit[j],0.333),3));
			zcrit_ij = (z_crit[i] + z_crit[j])*0.5;
			pcrit_ij = zcrit_ij*rrgas*tcrit_ij/vcrit_ij;
			w_ij = (w[i]+w[j])*0.5;
			c_ij = 0.37464 + 1.52226*w_ij - 0.26992*w_ij*w_ij;
			aij = 0.457235*rrgas*rrgas*tcrit_ij*tcrit_ij*pow(1+c_ij*(1-pow(temp/tcrit_ij,0.5)),2)/pcrit_ij;
			species_summation = species_summation + (-1/temp)*xi[j]*aij*(c_ij*pow(temp/tcrit_ij,0.5))/(1+c_ij*(1-pow(temp/tcrit_ij,0.5)));
		}
		return species_summation;
	
}

double volume_corr(struct thermo *argument)
{
	unsigned int i;
	double belowsum = 0;
	double c1m = 0;
	double wm = 0;
	double rrgas = (R/argument->parameter5);
	double tcm = 0;
	double vcm = 0;
	double s;
	for(i=0;i<N_SPECIES;i++)
	{
		s = argument->molefrac[i]*pow(v_crit[i],0.6667);
		wm = argument->molefrac[i]*w[i];
		c1m = c1m + 0.4266*z_crit[i] - 0.1101;
		belowsum = belowsum + s;
		tcm = tcm + s*t_crit[i];
		vcm = vcm + s*v_crit[i];
	}
	
	tcm = tcm/belowsum;
	vcm = vcm/belowsum;
	double r = rrgas*tcm;
	double dm = -1*pow(argument->den,2)*argument->parameter6/r;
	double pcm = (0.2905-0.085*wm)*r/vcm;
	double vcm_pr = 0.3074*r/pcm;
	double sig = vcm_pr - vcm;
	double cm = r*(c1m - (0.004 + c1m)*exp(-2*dm))/pcm;
	
	return cm + sig*(0.35/(0.35 + dm));	
}

/*******************************************************************/
/* Mixture Functions Structure                                     */
/*******************************************************************/

UDF_EXPORT RGAS_Functions RealGasFunctionList =
{
	MIXTURE_Setup,
	MIXTURE_Density,
	MIXTURE_Enthalpy,
	MIXTURE_Entropy,
	MIXTURE_Specific_Heat,
	MIXTURE_Molecular_Weight,
	MIXTURE_Speed_of_sound,
	MIXTURE_Viscosity,
	MIXTURE_ThermalConductivity,
	MIXTURE_rho_t,
	MIXTURE_rho_p,
	MIXTURE_enthalpy_t,
	MIXTURE_enthalpy_p,
	MIXTURE_enthalpy_prime
};

/* CH4 = 0, O2 = 1, CO = 2, H2 = 3, H20 = 4, O = 5, H = 6, OH = 7, CO2 = 8 */
real scale_parameter[N_SPECIES] = {25.14, 16.3, 18, 4.62, 13.1, 6.11, 2.31, 8.42, 26.8};

/* TAKAHASHI MASS DIFFUSIVITY */
DEFINE_DIFFUSIVITY(takahashi, cell, thread, i)
{
real takahashi;
real DPlow = 1.01;
real crit_pressure[N_SPECIES] = {45.99*1e5, 50.43*1e5, 34.94*1e5, 12.93*1e5, 220.64*1e5, 26.90*1e5, 1360000, 8200000, 73.74*1e5}; /* IN PASCALS */
real crit_temperature[N_SPECIES]={190.56, 154.58, 132.85, 32.98, 647.14, 44.5, 33.2, 400, 304.12};
real mol_weight[N_SPECIES] = {16.04303, 31.99880, 28.01, 2.016, 18.01534, 15.994, 1.00797, 17.007, 44.00995};
real T_reduced;
real P_reduced;
real pressure = RP_Get_Real("operating-pressure");
temp = C_T(cell, thread);
real above_sum =0;
real below_sum = 0;
real X_I[N_SPECIES];
real Y_I[N_SPECIES];
real sum = 0;
real fuller;
real AB_mw;
real a,b,c,f;
int j;

for (j=0 ; j<N_SPECIES ; j++)
{
Y_I[j] = C_YI(cell,thread,j);
}

for(j=0 ; j<N_SPECIES ; j++)
{
sum = sum + Y_I[j]/mol_weight[j];
}

sum = 1/sum;

/* Mole Fractions of the Species in mixture */
for(j=0 ; j<N_SPECIES ; j++)
{
	X_I[j] = (Y_I[j]*sum/mol_weight[j]);
}

for(j=0;j<N_SPECIES;j++)
{
	if(j != i)

	{
		above_sum = above_sum + (X_I[j] + 1e-12)*mol_weight[j];
	}
}

for(j=0;j<N_SPECIES;j++)
{
	if(j != i)

	{
		AB_mw = 2*mol_weight[i]*mol_weight[j]/(mol_weight[i] + mol_weight[j]);
		P_reduced = pressure/(crit_pressure[i]*X_I[i] + crit_pressure[j]*X_I[j]);
		T_reduced = temp/(crit_temperature[i]*X_I[i] + crit_temperature[j]*X_I[j]);
		a = (T_reduced-2.4)/1.5;
		b = 6.293*T_reduced*T_reduced - 9.0433*T_reduced + 2.9334;
		c = 0.015*T_reduced - 0.036;
		if(T_reduced <= 2.4)
		{
			
			f = (exp(a*P_reduced)+b)/(1+b);
		}
		else
		{
			f = 1+c*P_reduced;
		}
		fuller = (0.00143*pow(temp,1.75))/(1.01325*pow(AB_mw,0.5)*pow(pow(scale_parameter[i],0.333) + pow(scale_parameter[j],0.333),2));
		below_sum = below_sum + (X_I[j] + 1e-12) / (f*(101325/pressure)*fuller);
	}
}
below_sum = below_sum*sum;
takahashi = above_sum/below_sum;
return takahashi*0.0001;
}
