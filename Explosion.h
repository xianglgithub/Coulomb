#ifndef EXPLOSION_H
#define EXPLOSION_H

#include "coulomb.h"
#include <math.h>
#include <vector>

#define PI 3.14159265
#define AU2AG 0.52917721092
#define AMU2AU 1836.15
#define HALFC 180
#define AU2FS 0.02418884326505


class Explosion
{
  
  const double theta120 = PI*120/HALFC;
  
  const double theta240 = PI*240/HALFC;  
    
  double fcx, fcy, fcz, fix, fiy, fiz;
  double fh1x, fh1y, fh1z, fh2x, fh2y, fh2z, fh3x, fh3y, fh3z;
  
  int numion=5;
  int dims=30;  
  double y[30];
  const double tempmass[5]={AMU2AU*12, AMU2AU*127, AMU2AU*1, AMU2AU*1, AMU2AU*1};//{1836.15*6,1836.15*6};
  double tempcharge[5];
  
  double hstart,epsabs,epsrel,iternum;
  double tao,rate,tmax;
  double t0,t1,dt1,dt2,po;
  
  std::string sep="\t";

  double fcgC,fcgI,fcgI_s,fcgHa,fcgHb,fcgHc,fcg;

  double time;
  bool charged_up;
  bool exploded;
  std::vector<double> y_vec,cg_vec;
  
  public:
  Explosion();
  ~Explosion();
   
void init_charge(double cgC,double cgI,double cgH1,double cgH2,double cgH3);
void init_fcharge(double fcgC1,double fcgI1,double fcgI_s1,double fcgH1,double fcgH2,double fcgH3);
  void init_geometry(double r_ci,double r_ch,double theta_hci);
  
  void init_coulomb(double hstart,double epsabs,double epsrel,int iternum);  
  
  void init_model(double tao,double rate,double tmax); 
  
  void init_time(double t0,double t1,double dt1,double dt2);  
  
  bool GetCharges(double tao, double ctrsf, double cg, double Ccg, double Icg, double Icg_s, double Hcga,double Hcgb,double Hcgc,double dt, double time, double *charges);
  void Expl();
  std::vector<double> getY();
  std::vector<double> getCG();  
  
};


#endif