#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>

using namespace std;

//masses in reduced units
#define m 1.0

//lattice constant in reduced units
#define lc 3.5E+00

#define p 17.0
#define q 1.19
#define e 1.0
#define ro 1.0
#define A 0.03514






double Fr (double r, int N, int i, int j){
    
 double F = 0.0;
    
    
    double rr = r/ro - 1.0;
    for (int i = 0; i < N; i++){
        
        double Fel1 = 0.0;
        double Fel2 = 0.0;
        
        double Frep = 0.0;
        
        for (int j = 0; j < N; j++){
        
            if (j!=i){
                double v = -2.0*q* (rr);
                double u = -p*rr;
                
                double dv = -2.0*q/ro;
                double du = -p/ro;
                
                Frep += A*exp(u)*du;
                
                
                Fel1 += 0.5*exp(v)*dv;
                Fel2 += exp(v);
                
            }
        }
        
        
           F = Frep - (Fel1 / sqrt (Fel2));
           F = -F;
           
          
    }
        
        
        
        return F;
        
    
    
    
}


double Ur(double r,int N, int i, int j){
    
    double U = 0.0;
    
 
    
    double rr = r/ro - 1.0;
    for (int i = 0; i < N; i++){
        
        double Uel = 0.0;
        double Urep = 0.0;
        
        for (int j = 0; j < N; j++){
        
            if (j!=i){
                
                Uel += exp(-2.0*q*rr);
                Urep += A*exp(-p*rr);
            }
        }  
           
           U += e*0.5*(Urep - sqrt(Uel));
    }
        
        
        
        return U;
        
}

int main(){

    ofstream file; file.open("force-potential.txt");
    file.precision(8);
    
    file<< "r" <<"\t" << "U(Ni-Ni)" << "\t" << "F(Ni-Ni)" << "\t" << "-dU(Ni-Ni)/dr" << "\n";
    
    double N = 2; 
    double r = 0.0;
    
    double UNi= 0.0;
    double dUNi= 0.0;
    
    for (; r < 7; r+=0.001){
    
        UNi = Ur(r, N,1,1);
       
        
        
        dUNi = Fr(r, N,1,1);
        
        double h = 0.01;
     
        double fNi = -(Ur(r+h, N, 1,1)-Ur(r-h, N,1,1))/ (2*h) ;
       
        file<< r <<"\t" << UNi/N << "\t" << fNi << "\t" << dUNi << "\n";
        
        
    }
    
    
   return 0; 
    
}
