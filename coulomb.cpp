#include "coulomb.h"


int
func (double t, const double y[], double f[],
      void *params)
{
  struct my_params * pa 
    =(struct my_params *)params;
  int n=(pa->num);
  const std::vector<double> massarr=pa->mass;
  const std::vector<double> chargearr=pa->charge;
  double *r=new double[n*n];
  for (int i=0;i<n;i++)
{
  for (int j=0;j<n;j++)
{
  if(i!=j)
{
  r[i*n+j]=sqrt(pow(y[i*6]-y[j*6],2)+pow(y[i*6+2]-y[j*6+2],2)+pow(y[i*6+4]-y[j*6+4],2));
}
}
}
  for (int i=0;i<n;i++)
{
  f[i*6]=y[i*6+1];
  f[i*6+2]=y[i*6+3];
  f[i*6+4]=y[i*6+5];
  f[i*6+1]=0;f[i*6+3]=0;f[i*6+5]=0;
///////////////////////////////////////////////////////////////////////////////////////////////
  for(int j=0;j<n;j++)
{
  if(i!=j)
{
  f[i*6+1]=f[i*6+1]+(chargearr[i]*chargearr[j]/massarr[i])*(y[i*6]-y[j*6])/pow(r[i*n+j],3);
  f[i*6+3]=f[i*6+3]+(chargearr[i]*chargearr[j]/massarr[i])*(y[i*6+2]-y[j*6+2])/pow(r[i*n+j],3);
  f[i*6+5]=f[i*6+5]+(chargearr[i]*chargearr[j]/massarr[i])*(y[i*6+4]-y[j*6+4])/pow(r[i*n+j],3);
}
}
}
  delete[] r;
  return GSL_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{
  struct my_params * pa 
    =(struct my_params *)params;
  int n=(pa->num);
  const std::vector<double> massarr=pa->mass;
  const std::vector<double> chargearr=pa->charge;
  double *r=new double[n*n];
  double *dfdyarr=new double[6*n*6*n];//=

  std::fill_n(dfdyarr,6*n*6*n,0);//=
  //fill_n(dfdyarr,30*30,0);//=
  dfdy=dfdyarr;

//////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i=0;i<n;i++)
{
  for (int j=0;j<n;j++)
{
  if(i!=j)
{
  r[i*n+j]=sqrt(pow(y[i*6]-y[j*6],2)+pow(y[i*6+2]-y[j*6+2],2)+pow(y[i*6+4]-y[j*6+4],2));
}
}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy,6*n, 6*n);//add 30//=
  gsl_matrix * m = &dfdy_mat.matrix; 

  double tempjac;
  for (int i=0;i<n;i++)
{
  gsl_matrix_set(m,i*6,i*6+1,1);
  gsl_matrix_set(m,i*6+2,i*6+3,1);
  gsl_matrix_set(m,i*6+4,i*6+5,1);

  for (int j=0;j<n;j++)
{
  if (i==j)
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////for the row Vx
  tempjac=0;
  for (int k=0;k<n;k++)
{
  if (i!=k)
{
  tempjac=tempjac+(chargearr[i]*chargearr[k]/massarr[i])*(pow(r[i*n+k],3)-3*(y[i*6]-y[k*6])*(y[i*6]-y[k*6])*r[i*n+k])/pow(r[i*n+k],6);
}
}
  gsl_matrix_set(m,i*6+1,j*6,tempjac);  //x

  tempjac=0;
  for (int k=0;k<n;k++)
{
  if (i!=k)
{
  tempjac=tempjac-(chargearr[i]*chargearr[k]/massarr[i])*(3*(y[i*6]-y[k*6])*(y[i*6+2]-y[k*6+2])*r[i*n+k])/pow(r[i*n+k],6);
}
}
  gsl_matrix_set(m,i*6+1,j*6+2,tempjac);  //y

  tempjac=0;
  for (int k=0;k<n;k++)
{
  if (i!=k)  
{
  tempjac=tempjac-(chargearr[i]*chargearr[k]/massarr[i])*(3*(y[i*6]-y[k*6])*(y[i*6+4]-y[k*6+4])*r[i*n+k])/pow(r[i*n+k],6);
}
}
  gsl_matrix_set(m,i*6+1,j*6+4,tempjac);  //z
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////for the row Vy
  tempjac=0;
  for (int k=0;k<n;k++)
{
  if (i!=k)
{
  tempjac=tempjac+(chargearr[i]*chargearr[k]/massarr[i])*(pow(r[i*n+k],3)-3*(y[i*6+2]-y[k*6+2])*(y[i*6+2]-y[k*6+2])*r[i*n+k])/pow(r[i*n+k],6);
}
}
  gsl_matrix_set(m,i*6+3,j*6+2,tempjac);//Vy/y

  tempjac=0;
  for (int k=0;k<n;k++)
{
  if (i!=k)
{
  tempjac=tempjac-(chargearr[i]*chargearr[k]/massarr[i])*(3*(y[i*6+2]-y[k*6+2])*(y[i*6]-y[k*6])*r[i*n+k])/pow(r[i*n+k],6);
}
}
  gsl_matrix_set(m,i*6+3,j*6,tempjac);//Vy/x

  tempjac=0;
  for (int k=0;k<n;k++)
{
  if (i!=k)  
{
  tempjac=tempjac-(chargearr[i]*chargearr[k]/massarr[i])*(3*(y[i*6+2]-y[k*6+2])*(y[i*6+4]-y[k*6+4])*r[i*n+k])/pow(r[i*n+k],6);
}
}
  gsl_matrix_set(m,i*6+3,j*6+4,tempjac);//Vy/z

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////for the row Vz
  tempjac=0;
  for (int k=0;k<n;k++)
{
  if (i!=k)
{
  tempjac=tempjac+(chargearr[i]*chargearr[k]/massarr[i])*(pow(r[i*n+k],3)-3*(y[i*6+4]-y[k*6+4])*(y[i*6+4]-y[k*6+4])*r[i*n+k])/pow(r[i*n+k],6);
}
}
  gsl_matrix_set(m,i*6+5,j*6+4,tempjac);//Vz/z

  tempjac=0;
  for (int k=0;k<n;k++)
{
  if (i!=k)
{
  tempjac=tempjac-(chargearr[i]*chargearr[k]/massarr[i])*(3*(y[i*6+4]-y[k*6+4])*(y[i*6]-y[k*6])*r[i*n+k])/pow(r[i*n+k],6);
}
}
  gsl_matrix_set(m,i*6+5,j*6,tempjac);//Vz/x

  tempjac=0;
  for (int k=0;k<n;k++)
{
  if (i!=k)  
{
  tempjac=tempjac-(chargearr[i]*chargearr[k]/massarr[i])*(3*(y[i*6+4]-y[k*6+4])*(y[i*6+2]-y[k*6+2])*r[i*n+k])/pow(r[i*n+k],6);
}
}
  gsl_matrix_set(m,i*6+5,j*6+2,tempjac);//Vz/y

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
else
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////for the row Vx
  tempjac=-(chargearr[i]*chargearr[j]/massarr[i])*(pow(r[i*n+j],3)+3*(y[i*6]-y[j*6])*(y[i*6]-y[j*6])*r[i*n+j])/pow(r[i*n+j],6);//x
  gsl_matrix_set(m,i*6+1,j*6,tempjac);
  tempjac=(chargearr[i]*chargearr[j]/massarr[i])*(3*(y[i*6]-y[j*6])*(y[i*6+2]-y[j*6+2])*r[i*n+j])/pow(r[i*n+j],6);//y
  gsl_matrix_set(m,i*6+1,j*6+2,tempjac);
  tempjac=(chargearr[i]*chargearr[j]/massarr[i])*(3*(y[i*6]-y[j*6])*(y[i*6+4]-y[j*6+4])*r[i*n+j])/pow(r[i*n+j],6);//z
  gsl_matrix_set(m,i*6+1,j*6+4,tempjac);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////for the row Vy
  tempjac=-(chargearr[i]*chargearr[j]/massarr[i])*(pow(r[i*n+j],3)+3*(y[i*6+2]-y[j*6+2])*(y[i*6+2]-y[j*6+2])*r[i*n+j])/pow(r[i*n+j],6);//y
  gsl_matrix_set(m,i*6+3,j*6+2,tempjac);
  tempjac=(chargearr[i]*chargearr[j]/massarr[i])*(3*(y[i*6+2]-y[j*6+2])*(y[i*6]-y[j*6])*r[i*n+j])/pow(r[i*n+j],6);//x
  gsl_matrix_set(m,i*6+3,j*6,tempjac);
  tempjac=(chargearr[i]*chargearr[j]/massarr[i])*(3*(y[i*6+2]-y[j*6+2])*(y[i*6+4]-y[j*6+4])*r[i*n+j])/pow(r[i*n+j],6);//z
  gsl_matrix_set(m,i*6+3,j*6+4,tempjac);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////for the row Vz
  tempjac=-(chargearr[i]*chargearr[j]/massarr[i])*(pow(r[i*n+j],3)+3*(y[i*6+4]-y[j*6+4])*(y[i*6+4]-y[j*6+4])*r[i*n+j])/pow(r[i*n+j],6);//z
  gsl_matrix_set(m,i*6+5,j*6+4,tempjac);
  tempjac=(chargearr[i]*chargearr[j]/massarr[i])*(3*(y[i*6+4]-y[j*6+4])*(y[i*6]-y[j*6])*r[i*n+j])/pow(r[i*n+j],6);//x
  gsl_matrix_set(m,i*6+5,j*6,tempjac);
  tempjac=(chargearr[i]*chargearr[j]/massarr[i])*(3*(y[i*6+4]-y[j*6+4])*(y[i*6+2]-y[j*6+2])*r[i*n+j])/pow(r[i*n+j],6);//y
  gsl_matrix_set(m,i*6+5,j*6+2,tempjac);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
  std::fill_n(dfdt,2*3*n,0);
  delete[] r;
  delete[] dfdyarr;
  return GSL_SUCCESS;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double potential(my_params *para, const double d[])
{
  int n = para->num;
  std::vector<double> charge = para->charge;
  double poten=0;
  double *r=new double[n*n];
//////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i=0;i<n;i++)
{
  for (int j=0;j<n;j++)
{
  if(i!=j)
{
  r[i*n+j]=sqrt(pow(d[i*6]-d[j*6],2)+pow(d[i*6+2]-d[j*6+2],2)+pow(d[i*6+4]-d[j*6+4],2));
}
}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i=0;i<n;i++)
{
  for (int j=0;j<n;j++)
{
  if (i!=j)
{
  poten=poten+(charge[i]*charge[j])/r[i*n+j];
}
}
}
  delete[] r;
  return 27.211*poten/2;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int rk8(size_t dims, my_params *para, double t, double t1, int iternum, double hstart, double epsabs, double epsrel, double y[])
{

//  string sep="\t";

  gsl_odeiv2_system sys = {func, jac, dims, para}; //add 30//=
  
//  gsl_odeiv2_driver * d = 
 //   gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
//				  hstart, epsabs, epsrel);
				  
    
  gsl_odeiv2_driver * d = 
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
				  hstart, epsabs, epsrel);
				      
 // double enc,enia,enib,en,po;
  

  int i;
 // ofstream myfile("coulomb.dat");
 // if (!myfile.is_open())
//{
 // cout << "Unable to open file <coulomb.dat>"<<endl;
//}

  for (i = 1; i <= iternum; i++)
    {
      double ti = i * t1 / iternum;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
	{
	  printf ("error, return value=%d\n", status);
	  break;
	}
 //     enc=27.211*0.5*1836.15*14*(pow(y[1],2)+pow(y[3],2)+pow(y[5],2));
 //     enia=27.211*0.5*1836.15*127*(pow(y[7],2)+pow(y[9],2)+pow(y[11],2));
 //     enib=27.211*0.5*1836.15*127*(pow(y[13],2)+pow(y[15],2)+pow(y[17],2));
 //     en=enc+enia+enib;
 //     po=potential(&para, y);
      //myfile<<t<<"\t"<<po<<"\t"<<en<<"\t"<<enc<<"\t"<<3*enh1<<"\t"<<eni<<endl;
 //     myfile<<0.02418884326505*t<<sep<<po<<sep<<en<<sep<<enc<<sep<<enia<<sep<<enib<<endl;
      //myfile<<0.02418884326505*t<<"|"<<po<<"|"<<en<<"|"<<enc<<"|"<<enh2<<"|"<<eni<<endl;
      //printf ("%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", t, en, enc, enh1, enh2, enh3, eni);
    }
 // myfile.close();
//  double energyev=0.5*1836.15*6*(pow(y[1],2)+pow(y[3],2)+pow(y[5],2)+pow(y[7],2)+pow(y[9],2)+pow(y[11],2));
 // double energyev=27.211*0.5*1836.15*6*(pow(y[1],2)+pow(y[3],2)+pow(y[5],2)+pow(y[7],2)+pow(y[9],2)+pow(y[11],2)); 
  gsl_odeiv2_driver_free (d);
  return 0;
}

