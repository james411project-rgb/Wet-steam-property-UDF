#include "udf.h"
#include "stdio.h"
#include "ctype.h"
#include "stdarg.h"
#include "math.h"

#define R_const_UDRGM 0.46151     /* KJ/(kg¡¤K) specific gas constant */
#define M_PI_UDRGM 3.14159265358979323846

// Property calculation parameters
const double alpha_B_UDRGM = 10000.0;
const double a1_UDRGM = 0.0015;
const double a2_UDRGM = -0.000942;
const double a3_UDRGM = -0.0004882;
const double tau0_UDRGM = 0.8978;
const double alpha_C_UDRGM = 11.16;
const double aC_UDRGM = 1.772;
const double Tg_ref_UDRGM = 647.286;
const double R_gas_UDRGM = 0.46151;

static int (*usersMessage)(const char *, ...);
static void (*usersError)(const char *, ...);

// B calculation
double B_UDRGM(double Tg){
    double tau = 1500.0/Tg;
    double result = a1_UDRGM * pow((1+Tg/alpha_B_UDRGM),-1.0) + a2_UDRGM*exp(tau)*pow((1-exp(-tau)),2.5)*pow((tau),-0.5)+a3_UDRGM*tau;
    return result;
} 

// C calculation
double C_UDRGM(double Tg){
    double tau = Tg/647.286;
    double result = aC_UDRGM*(tau-tau0_UDRGM)*exp(-alpha_C_UDRGM*tau) + 1.5*pow(10.0,-6.0);
    return result;
} 

// First derivative of B
double dB_dTg_UDRGM(double Tg) {
    double tau = 1500.0 / Tg;
    double term1 = -a1_UDRGM / alpha_B_UDRGM * pow(1.0 + Tg / alpha_B_UDRGM, -2.0);

    double exp_tau = exp(tau);
    double exp_neg_tau = exp(-tau);

    double term2 = a2_UDRGM * exp_tau * (pow(tau, -0.5) - 0.5 * pow(tau, -1.5)) * pow(1 - exp_neg_tau, 2.5);
    double term3 = 2.5 * a2_UDRGM * pow(tau, -0.5) * pow(1 - exp_neg_tau, 1.5);
    double term4 = a3_UDRGM;

    double result = term1 - 1500.0 * pow(Tg, -2.0) * (term2 + term3 + term4);
    return result;
}

// Second derivative of B
double d2B_dTg2_UDRGM(double Tg) {
    double tau = 1500.0 / Tg;
    double exp_tau = exp(tau);
    double exp_neg_tau = exp(-tau);

    double term1 = 2.0 * a1_UDRGM / (alpha_B_UDRGM * alpha_B_UDRGM) * pow(1.0 + Tg / alpha_B_UDRGM, -3.0);

    double part1 = a2_UDRGM * exp_tau * (pow(tau, -0.5) - 0.5 * pow(tau, -1.5)) * pow(1 - exp_neg_tau, 2.5)
                 + 2.5 * a2_UDRGM * pow(tau, -0.5) * pow(1 - exp_neg_tau, 1.5)
                 + a3_UDRGM;

    double part2 = a2_UDRGM * exp_tau * (pow(tau, -0.5) - pow(tau, -1.5) + 0.75 * pow(tau, -2.5)) * pow(1 - exp_neg_tau, 2.5)
                 + 2.5 * a2_UDRGM * (pow(tau, -0.5) - 0.5 * pow(tau, -1.5)) * pow(1 - exp_neg_tau, 1.5)
                 - 1.25 * a2_UDRGM * pow(tau, -1.5) * pow(1 - exp_neg_tau, 1.5)
                 + 3.75 * a2_UDRGM * pow(tau, -0.5) * exp_neg_tau * pow(1 - exp_neg_tau, 0.5);

    double result = term1 + 3000.0 * pow(Tg, -3.0) * part1 + pow(1500.0,2) * pow(Tg, -4.0) * part2;
    return result;
}

// First derivative of C (dC/d¦Ó)
double dC_dtau_UDRGM(double Tg) {
    double tau = Tg / Tg_ref_UDRGM;
    double exp_neg = exp(-alpha_C_UDRGM * tau);
    double result = (aC_UDRGM * exp_neg - aC_UDRGM * alpha_C_UDRGM * (tau - tau0_UDRGM) * exp_neg)/647.286;
    return result;
}

// Second derivative of C (d2C/d¦Ó2)
double d2C_dtau2_UDRGM(double Tg) {
    double tau = Tg / Tg_ref_UDRGM;
    double exp_neg = exp(-alpha_C_UDRGM * tau);
    double result = -(aC_UDRGM * alpha_C_UDRGM / (Tg_ref_UDRGM * Tg_ref_UDRGM)) * exp_neg * (2.0 - alpha_C_UDRGM * (tau - tau0_UDRGM));
    return result;
}

// B1 calculation
double B1_UDRGM(double Tg){
    double result = Tg*dB_dTg_UDRGM(Tg);
    return result;
}

// B2 calculation
double B2_UDRGM(double Tg){
    double result = Tg*Tg*d2B_dTg2_UDRGM(Tg);
    return result;
} 

// C1 calculation
double C1_UDRGM(double Tg){
    double result = Tg*dC_dtau_UDRGM(Tg);
    return result;
}

// C2 calculation
double C2_UDRGM(double Tg){
    double result = Tg*Tg*d2C_dtau2_UDRGM(Tg);
    return result;
}

// Density calculation 
double RHO_UDRGM(double Tg, double P){
    double term1 = 2.0*P/(R_gas_UDRGM*Tg);
    double term2 = 1.0 + pow((1.0 +(4*P*B_UDRGM(Tg)/(R_gas_UDRGM*Tg))),0.5);
    return term1 / term2;
}

// alpha_T calculation
double alpha_T_UDRGM(double Tg,double P){
    return (1+(B_UDRGM(Tg)+B1_UDRGM(Tg))*RHO_UDRGM(Tg,P)+(C_UDRGM(Tg)+C1_UDRGM(Tg))*RHO_UDRGM(Tg,P)*RHO_UDRGM(Tg,P))/(1+2.0*B_UDRGM(Tg)*RHO_UDRGM(Tg,P)+3.0*C_UDRGM(Tg)*RHO_UDRGM(Tg,P)*RHO_UDRGM(Tg,P));
} 

// beta_P calculation
double beta_P_UDRGM(double Tg,double P){
    return(1.0+B_UDRGM(Tg)*RHO_UDRGM(Tg,P)+C_UDRGM(Tg)*RHO_UDRGM(Tg,P)*RHO_UDRGM(Tg,P))/(1.0+2.0*B_UDRGM(Tg)*RHO_UDRGM(Tg,P)+3.0*C_UDRGM(Tg)*RHO_UDRGM(Tg,P)*RHO_UDRGM(Tg,P));
} 

// c_p0 calculation
double c_p0_UDRGM(double Tg){
    return 46.0*pow(Tg,-1.0)+1.47276
    +8.3893*pow(10.0,-4.0)*Tg
    -2.19989*pow(10.0,-7.0)*pow(Tg,2.0)
    +2.46619*pow(10.0,-10.0)*pow(Tg,3.0)
    -9.70466*pow(10.0,-14.0)*pow(Tg,4.0);
} 

// Isobaric specific heat c_pg calculation
double c_pg_UDRGM(double Tg,double P){
    double term1 = c_p0_UDRGM(Tg);
    double term2 = ((1-alpha_T_UDRGM(Tg,P))*(B_UDRGM(Tg)-B1_UDRGM(Tg))-B2_UDRGM(Tg))*R_gas_UDRGM*RHO_UDRGM(Tg,P);
    double term3 = ((1-2.0*alpha_T_UDRGM(Tg,P))*C_UDRGM(Tg)+alpha_T_UDRGM(Tg,P)*C1_UDRGM(Tg)-0.5*C2_UDRGM(Tg))*R_gas_UDRGM*RHO_UDRGM(Tg,P)*RHO_UDRGM(Tg,P);
    return term1 + term2 + term3;
} 

// Isochoric specific heat c_vg calculation
double c_vg_UDRGM(double Tg,double P){
    double term1 = c_p0_UDRGM(Tg);
    double term2 = -(1.0 + (2.0*B1_UDRGM(Tg)+B2_UDRGM(Tg))*RHO_UDRGM(Tg,P))*R_gas_UDRGM;
    double term3 = -(C1_UDRGM(Tg)+0.5*C2_UDRGM(Tg))*R_gas_UDRGM*RHO_UDRGM(Tg,P)*RHO_UDRGM(Tg,P);
    return term1+term2+term3;
}

// Specific heat ratio Gama
double Gama_UDRGM(double Tg,double P){
    double result = c_pg_UDRGM(Tg,P)/(c_vg_UDRGM(Tg,P)*beta_P_UDRGM(Tg,P));
    return result;
} 

// Steam dynamic viscosity mu_g calculation
double mu_g_UDRGM(double Tg){
    double Tr=Tg/647.286;
    double result = pow(10.0,-6.0)*(-15.371+99.871*Tr-133.933*Tr*Tr+75.8226*Tr*Tr*Tr);
    return result;
} 

// Steam thermal conductivity lam_g calculation 
double lam_g_UDRGM(double Tg){
    double Tr=Tg/647.286;
    double result = -0.2045+1.137*Tr-1.939*Tr*Tr+1.1418*Tr*Tr*Tr;
    return result;
} 

// Steam enthalpy h_g calculation
double h_g_UDRGM(double Tg,double P){
    double result = 46.0*log(Tg) + 1.47276*Tg 
                    +4.19465*pow(10.0,-4.0)*pow(Tg,2.0)
                    -7.3329667*pow(10.0,-8.0)*pow(Tg,3.0)
                    +6.165475*pow(10.0,-11.0)*pow(Tg,4.0)
                    -1.940932*pow(10.0,-14.0)*pow(Tg,5.0)
                    +1811.6
                    +R_gas_UDRGM*Tg*((B_UDRGM(Tg)-B1_UDRGM(Tg))*RHO_UDRGM(Tg,P)+(C_UDRGM(Tg)-0.5*C1_UDRGM(Tg))*RHO_UDRGM(Tg,P)*RHO_UDRGM(Tg,P));
    return result;              
}                    

// Steam entropy s_g calculation
double s_g_UDRGM(double Tg,double P){
    double result = -46.0/Tg + 1.47276*log(Tg)
                    +8.3893*pow(10.0,-4.0)*Tg
                    -1.099945*pow(10.0,-7.0)*pow(Tg,2.0)
                    +8.22063333*pow(10.0,-11.0)*pow(Tg,3.0)
                    -2.426165*pow(10.0,-14.0)*pow(Tg,4.0)
                    +0.97012
                    -R_gas_UDRGM*(log(RHO_UDRGM(Tg,P))+(B_UDRGM(Tg)+B1_UDRGM(Tg))*RHO_UDRGM(Tg,P)+0.5*(C_UDRGM(Tg)+C1_UDRGM(Tg))*RHO_UDRGM(Tg,P)*RHO_UDRGM(Tg,P));
    return result;              
}

DEFINE_ON_DEMAND(I_do_nothing_UDRGM)
{
    /* This is a dummy function to allow us to use */
    /* the Compiled UDFs utility                   */
}

void IDEAL_error(int err, char *f, char *msg)
{
    if (err)
        usersError("IDEAL_error (%d) from function: %s\n%s\n", err, f, msg);
}

void IDEAL_Setup(Domain *domain, cxboolean vapor_phase, char *filename,
                 int (*messagefunc)(const char *format, ...),
                 void (*errorfunc)(const char *format, ...))
{
    usersMessage = messagefunc;
    usersError  = errorfunc;
    usersMessage("\nLoading Real Gas Model for Water Steam: %s\n", filename);
}

double IDEAL_density(cell_t cell, Thread *thread,
                     cxboolean vapor_phase, double Temp, double press, double yi[])
{
    double P_kPa = press / 1000.0;
    double rho = RHO_UDRGM(Temp, P_kPa);
    return rho;
}

double IDEAL_specific_heat(cell_t cell, Thread *thread,
                           double Temp, double density, double P, double yi[])
{
    double P_kPa = P / 1000.0;
    double cp = c_pg_UDRGM(Temp, P_kPa) * 1000.0;
    return cp;
}

double IDEAL_enthalpy(cell_t cell, Thread *thread,
                      double Temp, double density, double P, double yi[])
{
    double P_kPa = P / 1000.0;
    double h = h_g_UDRGM(Temp, P_kPa) * 1000.0;
    return h;
}

#define TDatum 288.15
#define PDatum 1.01325e5

double IDEAL_entropy(cell_t cell, Thread *thread,
                     double Temp, double density, double P, double yi[])
{
    double P_kPa = P / 1000.0;
    double s = s_g_UDRGM(Temp, P_kPa) * 1000.0;
    return s;
}

double IDEAL_mw(double yi[])
{
    return 18.015;
}

double IDEAL_speed_of_sound(cell_t cell, Thread *thread,
                            double Temp, double density, double P, double yi[])
{
    double P_kPa = P / 1000.0;
    double cp = c_pg_UDRGM(Temp, P_kPa) * 1000.0;
    double cv = c_vg_UDRGM(Temp, P_kPa) * 1000.0;
    double gamma = cp / cv;
    return sqrt(gamma*R_gas_UDRGM*1000.0*Temp);
}

double IDEAL_viscosity(cell_t cell, Thread *thread,
                       double Temp, double density, double P, double yi[])
{
    double mu = mu_g_UDRGM(Temp);
    return mu;
}

double IDEAL_thermal_conductivity(cell_t cell, Thread *thread,
                                  double Temp, double density, double P, double yi[])
{
    double ktc = lam_g_UDRGM(Temp);
    return ktc;
}

double IDEAL_rho_t(cell_t cell, Thread *thread,
                   double Temp, double density, double P, double yi[])
{
    double delta_T = 0.1;
    double P_kPa = P / 1000.0;
    double rho1 = RHO_UDRGM(Temp + delta_T/2.0, P_kPa);
    double rho2 = RHO_UDRGM(Temp - delta_T/2.0, P_kPa);
    double rho_t = (rho1 - rho2) / delta_T;
    return rho_t;
}

double IDEAL_rho_p(cell_t cell, Thread *thread,
                   double Temp, double density, double P, double yi[])
{
    double delta_P = 10.0;
    double P_kPa1 = (P + delta_P/2.0) / 1000.0;
    double P_kPa2 = (P - delta_P/2.0) / 1000.0;
    double rho1 = RHO_UDRGM(Temp, P_kPa1);
    double rho2 = RHO_UDRGM(Temp, P_kPa2);
    double rho_p = (rho1 - rho2) / delta_P;
    return rho_p;
}

double IDEAL_enthalpy_t(cell_t cell, Thread *thread,
                        double Temp, double density, double P, double yi[])
{
    double delta_T = 0.1;
    double P_kPa = P / 1000.0;
    double h1 = h_g_UDRGM(Temp + delta_T/2.0, P_kPa) * 1000.0;
    double h2 = h_g_UDRGM(Temp - delta_T/2.0, P_kPa) * 1000.0;
    double h_t = (h1 - h2) / delta_T;
    return h_t;
}

double IDEAL_enthalpy_p(cell_t cell, Thread *thread,
                        double Temp, double density, double P, double yi[])
{
    double delta_P = 10.0;
    double P_kPa1 = (P + delta_P/2.0) / 1000.0;
    double P_kPa2 = (P - delta_P/2.0) / 1000.0;
    double h1 = h_g_UDRGM(Temp, P_kPa1) * 1000.0;
    double h2 = h_g_UDRGM(Temp, P_kPa2) * 1000.0;
    double h_p = (h1 - h2) / delta_P;
    return h_p;
}

UDF_EXPORT RGAS_Functions RealGasFunctionList =
{
    IDEAL_Setup,
    IDEAL_density,
    IDEAL_enthalpy,
    IDEAL_entropy,
    IDEAL_specific_heat,
    IDEAL_mw,
    IDEAL_speed_of_sound,
    IDEAL_viscosity,
    IDEAL_thermal_conductivity,
    IDEAL_rho_t,
    IDEAL_rho_p,
    IDEAL_enthalpy_t,
    IDEAL_enthalpy_p
};
