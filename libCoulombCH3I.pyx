# distutils: language = c++


# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector
# c++ interface to cython
cdef extern from "Explosion.h":
  cdef cppclass Explosion:
        Explosion() except +      
        void init_charge(double,double,double,double,double)
        void init_fcharge(double,double,double,double,double,double)
        void init_geometry(double,double,double) 
        void init_coulomb(double,double,double,int)    
        void init_model(double,double,double)   
        void init_time(double,double,double,double)  
        void Expl()
        vector[double] getY()
        vector[double] getCG()  
	    
# creating a cython wrapper class
cdef class PyCoulombCH3I:
    cdef Explosion *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new Explosion()
    def __dealloc__(self):
        del self.thisptr
        
    def init_charge(self,cgC,cgI,cgH1,cgH2,cgH3):
        self.thisptr.init_charge(cgC,cgI,cgH1,cgH2,cgH3)        
        
    def init_fcharge(self,fcgC1,fcgI1,fcgI_s1,fcgH1,fcgH2,fcgH3):
        self.thisptr.init_fcharge(fcgC1,fcgI1,fcgI_s1,fcgH1,fcgH2,fcgH3)
        
    def init_geometry(self,r_ci,r_ch,theta_hci):
        self.thisptr.init_geometry(r_ci,r_ch,theta_hci)
        
    def init_coulomb(self,hstart,epsabs,epsrel,iternum):
        self.thisptr.init_coulomb(hstart,epsabs,epsrel,iternum)

    def init_model(self,tao,rate,tmax):
        self.thisptr.init_model(tao,rate,tmax)   
        
    def init_time(self,t0,t1,dt1,dt2):
        self.thisptr.init_time(t0,t1,dt1,dt2)

    def Expl(self):
        self.thisptr.Expl()          
        
    def getY(self):
        return self.thisptr.getY()                       
        
    def getCG(self):
        return self.thisptr.getCG()       

        
        
        
        
        
        
        	         
        
