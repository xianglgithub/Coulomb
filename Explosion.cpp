#include "Explosion.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Explosion::Explosion()
{
}

Explosion::~Explosion()
{
}


void Explosion::init_charge(double cgC,double cgI,double cgH1,double cgH2,double cgH3)
{
tempcharge[0] = cgC;tempcharge[1] = cgI;
tempcharge[2] = cgH1;tempcharge[3] = cgH2;tempcharge[4] = cgH3;
}

void Explosion::init_fcharge(double fcgC1,double fcgI1,double fcgI_s1,double fcgH1,double fcgH2,double fcgH3)
{
fcgC=fcgC1; fcgI=fcgI1; fcgHa=fcgH1; fcgI_s=fcgI_s1; fcgHb=fcgH2; fcgHc=fcgH3;
fcg = fcgC+fcgI+fcgHa+fcgHb+fcgHc;
} 
  
void Explosion::init_geometry(double r_ci,double r_ch,double theta_hci)
{
    

for (int i=0; i<dims; ++i)
{
y[i] = 0;
}

  theta_hci = PI*theta_hci/HALFC;
    
  r_ci = r_ci/AU2AG;
  r_ch = r_ch/AU2AG;

  y[10] = r_ci;
  y[12] = r_ch*sin(theta_hci); y[14] = 0; y[16] = r_ch*cos(theta_hci);
  fh1x = r_ch*sin(theta_hci); fh1y = 0; fh1z = r_ch*cos(theta_hci);    
  y[18] = cos(theta120)*fh1x-sin(theta120)*fh1y; y[20] = sin(theta120)*fh1x+cos(theta120)*fh1y; y[22] = fh1z;
  y[24] = cos(theta240)*fh1x-sin(theta240)*fh1y; y[26]= sin(theta240)*fh1x+cos(theta240)*fh1y; y[28] = fh1z; 

  //  for (int i=0; i<dims; ++i)
//{
//std::cout<<y[i]<<sep;
//}
// std::cout<<std::endl;     
}

  
void Explosion::init_coulomb(double hstart1,double epsabs1,double epsrel1,int iternum1) 
{
  hstart=hstart1; epsabs=epsabs1; epsrel=epsrel1; iternum=iternum1;
  charged_up = false;
  exploded = false;
  
}
  
void Explosion::init_model(double tao1,double rate1,double tmax1)
{
tao=tao1/AU2FS;rate=rate1*AU2FS;tmax=tmax1/AU2FS;
}
  
void Explosion::init_time(double t01,double t11,double dt11,double dt21)
{
t0=t01;t1=t11;dt1=dt11;dt2=dt21;
time=0;
}
  

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Explosion::GetCharges(double tao, double ctrsf, double cg, double Ccg, double Icg, double Icg_s, double Hcga,double Hcgb,double Hcgc,double dt, double time, double *charges)
{
double Tcg,Mcg;

if (charges[0]>=Ccg && (time>=tmax || charges[1]>=Icg_s))
{
return true;
}
else
{
 if (charges[0]>=Ccg)
 {
     Tcg = cg*(1-exp(-time/tao));
Mcg = charges[0]+charges[2]+charges[3]+charges[4];

charges[1] = Tcg-Mcg;     
 }
     else
     {
Tcg = cg*(1-exp(-time/tao));
Mcg = charges[0]+charges[2]+charges[3]+charges[4];

charges[1] = Tcg-Mcg;
//Mcg += ctrsf*exp(-kappa*d)*charges[1]*dt;
Mcg += ctrsf*charges[1]*dt;
charges[0] = Mcg*(Ccg/(Ccg+Hcga+Hcgb+Hcgc));
charges[2] = Mcg*(Hcga/(Ccg+Hcga+Hcgb+Hcgc));
charges[3] = Mcg*(Hcgb/(Ccg+Hcga+Hcgb+Hcgc));
charges[4] = Mcg*(Hcgc/(Ccg+Hcga+Hcgb+Hcgc));
     }
return false;
}
    
}


std::vector<double> Explosion::getY()
{
    y_vec.clear();
    if (exploded)
    {
        for (int i=0; i<dims; ++i)
        {
            y_vec.push_back(y[i]);
        }
               
    }
    else
    {
                for (int i=0; i<dims; ++i)
        {
            y_vec.push_back(0);
        }

    }
    

   return y_vec;
    
}


std::vector<double> Explosion::getCG()
{
    cg_vec.clear();

        
    for (int j=0; j<numion; ++j)
    {
        cg_vec.push_back(tempcharge[j]);
    }        
    
   return cg_vec;
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Explosion::Expl()
{

    
  while (!charged_up && time<t1)
 {
      
  charged_up = GetCharges(tao, rate, fcg, fcgC, fcgI, fcgI_s, fcgHa, fcgHb, fcgHc, dt1, time,tempcharge);

  my_params para;
  para.num=numion;
  para.mass.insert (para.mass.begin(), tempmass, tempmass+numion);

  
  para.charge.insert (para.charge.begin(), tempcharge, tempcharge+numion);

  po=potential(&para, y);

  rk8(dims,&para,t0,dt1,iternum,hstart,epsabs,epsrel,y);



  time += dt1;
  
 }
    
  while (time<t1)
 {

  my_params para;
  para.num=numion;
  para.mass.insert (para.mass.begin(), tempmass, tempmass+numion);

  
  para.charge.insert (para.charge.begin(), tempcharge, tempcharge+numion);

  po=potential(&para, y);

  rk8(dims,&para,t0,dt2,iternum,hstart,epsabs,epsrel,y);


  time += dt2;
  
 }
        
    
  exploded = true;
}

