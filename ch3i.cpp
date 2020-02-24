#include "Explosion.h"

int main(void)
{
     std::string sep=",";
Explosion ep;
ep.init_charge(0,0,0,0,0);
ep.init_fcharge(2,4,3.94,1,1,1);
ep.init_geometry(2.136,1.084,107.47);

ep.init_coulomb(1,1e-6,1e-8,1);
ep.init_model(9,0.37,100);    
ep.init_time(0.0,1e8,3,1e6);    
ep.Expl();
vector<double> y_vec,cg_vec;
y_vec = ep.getY();
cg_vec = ep.getCG();
    
   std::cout<<"Charges:";
   for (int i=0; i<cg_vec.size();++i)
   {
       std::cout<<cg_vec[i]<<sep;
   }
    std::cout<<std::endl;    

   std::cout<<"Y:";
    
    
    
return 0;
}

