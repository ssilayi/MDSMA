/**
//  sma.cpp
//  md_modeling
//
//  Created by Swabir Silayi on 2/19/18.
//  Copyright © 2018 Swabir Silayi. All rights reserved.
*/


#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <random>


using namespace std;

random_device rd;
mt19937 engine(rd());
uniform_real_distribution<double> dist(0,1);


bool read_in_parameters(string &el, double &m, double &lc,
                        double &e, double &ro, double &A, double &p, double &q,
                        int &N, double &T, double &dt, int &nSteps){
    
    string hold;
    ifstream infile; infile.open("mdin");
    if(infile){
        infile>>hold>>hold>>hold>>hold>>el;
        infile>>hold>>hold>>hold>>hold>>m;
        infile>>hold>>hold>>hold>>hold>>lc;
        infile>>hold>>hold>>hold>>e>>ro>>A>>p>>q;
        infile>>hold>>hold>>hold>>hold>>N;
        infile>>hold>>hold>>hold>>T;
        infile>>hold>>hold>>dt;
        infile>>hold>>hold>>hold>>nSteps;
        
    }
    else{
        cout << "INPUT FILE 'mdin' NOT FOUND! EXITING. \n";
        return false;
    }
    
    
    return true;
}

void initialize_Positions(int N, double lc, double &L, double ** r) {
    
    // find M large enough to fit N atoms on an fcc lattice
    int M = 0.0;
    while (4*M*M*M < N)
        M = M+1;
    
    L = M*lc;
    
    // positions in fcc unit cell
    double dx[4] = {0.0, 0.5, 0.5, 0.0};
    double dy[4] = {0.0, 0.5, 0.0, 0.5};
    double dz[4] = {0.0, 0.0, 0.5, 0.5};
    
    int n = 0;
    
    for (int x = 0; x < M; x++){
        for (int y = 0; y < M; y++){
            for (int z = 0; z < M; z++)
                for (int i = 0; i < 4; i++){
                    if (n < N){
                        r[n][0] = (x + dx[i]) * lc;
                        r[n][1] = (y + dy[i]) * lc;
                        r[n][2] = (z + dz[i]) * lc;
                        ++n;
                    }
                }
        }
    }
    
    
    double temp[3] ={ 0.0 };
    
    for (int j = 0; j < n; j++){
        for (int i = 0; i < n-1; i++){
            if(r[i][0] > r[i+1][0] ){
                
                for (int k = 0; k < 3; k++){
                    temp[k] = r[i][k];
                    r[i][k] = r[i+1][k];
                    r[i+1][k] = temp[k];
                }
            }
        }
    }
    
    
    // Adjust  so center-of-mass is zero
    double rCM[3] = {0, 0, 0};
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            rCM[i] += r[n][i];
    for (int i = 0; i < 3; i++)
        rCM[i] /= N;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            r[n][i] -= rCM[i];
    
    
    return;
    
}


void initialize_Velocities(int N, double ** v, double m, double T, double &K) {
    
    //Uniform random
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] = dist(engine);
    // Adjust velocities so center-of-mass velocity is zero
    double vCM[3] = {0, 0, 0};
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vCM[i] += v[n][i];
    for (int i = 0; i < 3; i++)
        vCM[i] /= N;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] -= vCM[i];
    
    // Rescale velocities to get the desired instantaneous temperature
    K = 0.0;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            K += m*v[n][i] * v[n][i];
    
    K *=0.5;
    
    double Ts = 2.0 * K / (3.0 * (N - 1));
    
    double lambda = sqrt( T / Ts );
    
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] *= lambda;
    
    //initial KE
    K = 0.0;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            K += m*v[n][i] * v[n][i];
    
    K *=0.5;
    
    
    return;
}

double instantaneous_temperature(int N, double K) {
    
    return 2.0 * K / ( 3.0 * (N - 1));
    
}
void rescaleVelocities(int N, double K, double T, double **v) {
    
    
    double lambda = sqrt( T / instantaneous_temperature(N, K) );
    
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] *=  lambda;
    
}

void initialize_configuration(int N, double m, double lc, double &L,
                              double **r, double **v, double &K, double T){
    
    initialize_Positions( N, lc, L, r) ;
    initialize_Velocities( N, v, m, T, K);
    
    return;
}

void minimum_image(double &dr, double L){
    
    if (dr > 0.5*L)
        dr = dr- L;
    
    if (dr < -0.5*L)
        dr= dr+L;
    
    return;
}

void computeAccelerations(int N, double ** r, double ** a, double **f,
                          double L, double &U, double m, double A, double p, double q ) {
    
    int i, j;
    U = 0.0;
    double dx=0.0, dy=0.0, dz=0.0;
    
    double rij=0.0;
    double rr=0.0;
    double rCutOff = 0.49*L;
    
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    
    
    //setting to zero
    for (i = 0; i < N; i++)
    {
        
        f[i][0] = 0.0; f[i][1] = 0.0; f[i][2] = 0.0;
        
        a[i][0] = 0.0; a[i][1] = 0.0; a[i][2] = 0.0;
    }
    
    
    for (i = 0; i < N; i++)
    {
        sum1 = 0.0;
        sum2 = 0.0;
        sum3 = 0.0;
        
        for (j = 0; j < i; j++)
        {
            
            if (j != i){
                
                dx = r[i][0] - r[j][0];
                dy = r[i][1] - r[j][1];
                dz = r[i][2] - r[j][2];
                
                
                /*minimum image convention*/
                minimum_image(dx, L);
                minimum_image(dy, L);
                minimum_image(dz, L);
                
                
                rij = sqrt (dx * dx + dy * dy + dz * dz);
                
                if (rij <= rCutOff)
                {
                    rr = rij  - 1.0 ;
                    
                    sum1 += A  * exp (-p * rr );
                    sum2 += exp (-2.0* q * rr );
                    
                    double f1 = -(A * p * exp(-p*rr) )/rij;
                    double f2 = ( q*exp(-2.0*q*rr ) )/rij ;
                    double f3 =  sqrt( exp(-2.0*q*rr  ) );
                    
                    double ff = ( f1 + f2/ f3 ) ;
                    ff = ff;
                    
                    f[i][0] = f[i][0] - (ff*dx );
                    f[j][0] = f[j][0] + (ff*dx );
                    
                    f[i][1] = f[i][1] - (ff* dy) ;
                    f[j][1] = f[j][1] + (ff* dy) ;
                    
                    f[i][2] = f[i][2] - (ff* dz );
                    f[j][2] = f[j][2] + (ff* dz) ;
                    
                }//end if
                
            }//end if j .ne. i
            
        }//end for j
        
        sum3 = sqrt (sum2);
        
        U += 0.5 * ( sum1 - sum3 );
        
        
        a[i][0] = f[i][0]/m;
        a[i][1] = f[i][1]/m;
        a[i][2] = f[i][2]/m;
        
    }//end loop over i
    
    
    return;
}

void initialize_gr(double * gr, int n_bins, double &bin_size, double L){
    
    bin_size = (L/(double)(2.0*n_bins) );
    
    for (int i = 0; i < n_bins; i++)
        gr[i] = 0.0;
    
    return;
    
}

void update_gr (int N, double **r, double L, double * gr, double bin_size ) {
    
    double rCutOff = 0.49*L;
    
    double dr[3] = {0.0}; double rSqd = 0.0;
    
    for (int i = 0; i < N-1; i++)
    {
        for (int j = i+1; j < N; j++)
        {
            
            dr[0] = r[i][0] - r[j][0];
            dr[1] = r[i][1] - r[j][1];
            dr[2] = r[i][2] - r[j][2];
            
            
            for (int k = 0; k < 3; k++){
                if (dr[k] > 0.5*L)
                    dr[k] = dr[k]- L;
                
                if (dr[k] < -0.5*L)
                    dr[k] = dr[k]+L;
            }
            
            rSqd = dr[0]*dr[0]+ dr[1]*dr[1]+ dr[2]*dr[2];
            
            if (rSqd <= rCutOff*rCutOff)
            {
                int bin=(int)(sqrt(rSqd)/bin_size);
                gr[bin]+=2;
            }
        }
    }
    
    return;
}


void normalize_gr(int N, int ngr, double * gr, int nbins, double dbin, double L, string filename){
    
    ofstream gfile; gfile.open (filename.c_str());
    double rho = (double) N/ (L*L*L);
    
    /* Normalize radial distribution g(r) and save  to file*/
    for (int i=0; i< nbins;i++) {
        double rr = dbin*(i+0.5);
        double vb=((i+1)*(i+1)*(i+1)-i*i*i)*dbin*dbin*dbin;
        double nid=(4./3.)*M_PI*vb*rho;
        
        gfile <<rr << "\t"<< i*dbin<<"\t" <<(double)(gr[i])/(ngr*N*nid) <<"\n";
    }
    
    
}


void mean_square_displacement(int N, double ** r, double ** po, double * msd){
    
    msd[0] = 0.0;
    msd[1] = 0.0;
    msd[2] = 0.0;
    
    double dx, dy, dz;
    for (int i = 0; i < N; i++){
        
        dx = r[i][0] - po[i][0];
        dy = r[i][1] - po[i][1];
        dz = r[i][2] - po[i][2];
        
        msd[0] += dx*dx;
        msd[1] += dy*dy;
        msd[2] += dz*dz;
    }
    
    msd[0] = msd[0]/ (double) N;
    msd[1] = msd[1]/ (double) N;
    msd[2] = msd[2]/ (double) N;
    
    return;
}


void ave_var(int cc, double &Avg, double &Var, double U)
{
    
    double old_Avg = Avg;
    
    Avg = cc == 0 ? U : Avg + ((U - Avg) / cc);
    
    Var = cc == 0 ? 0 : Var * cc + (U - old_Avg) * (U - Avg);
    Var /= (cc + 1);
}



void save_configuration(int N, string el, double **r, double ** v, string filename){
    
    ofstream file; file.open(filename.c_str());
    
    file << N<<"\n";
    file << "Lattice=0 0 0 0.5 0.5 0.5 Properties=species:S:1:pos:R:3 Time=0 \n";
    //L <<"\t" << L << "\t" << L << "\n";;
    for (int i = 0; i < N; i++)
    {
        file << el <<"\t";
        for (int j = 0; j < 3; j++)
            file << r[i][j] <<setw(12)<< "\t";
        for (int j = 0; j < 3; j++)
            file << v[i][j] << setw(12)<<"\t";
        file << "\n";
        
    }
    
    file.close();
    
    return;
    
}

int main(int argc, const char * argv[]){
    
    string el;
    double m;
    double lc;
    double e,ro,A,p,q;
    int N;
    double Tt = 900.0;
    double dt;
    int nSteps;
    
    bool check = read_in_parameters(el, m, lc, e, ro, A, p, q, N, Tt, dt, nSteps);
    if(!check) return 0;
    
    
    double ** r;
    double ** v;
    double ** f;
    double ** a;
    double L;
    
    r = new double * [N];
    v = new double * [N];
    f = new double * [N];
    a = new double * [N];
    
    for (int i = 0; i < N; i++){
        r[i] = new double [3];
        v[i] = new double [3];
        f[i] = new double [3];
        a[i] = new double [3];
    }
    
    
    
    /*open file for writing*/
    system("mkdir -p Results/");
    ofstream hfile;
    std::ostringstream hf;
    hf <<"Results/data.txt" ;
    std::string hmf = hf.str();
    hfile.open(hmf.c_str(), ios_base::app);
    
    const double kb = 8.6173303e-5;  // eV K-1
    
    double T = kb*Tt;
    
    
    double U = 0.0; double UAvg = 0.0, UVar=0.0;
    double K = 0.0; double KAvg = 0.0, KVar=0.0;
    double E = 0.0; double EAvg = 0.0, EVar=0.0;
    double iT = 0.0; double iTAvg = 0.0, iTVar=0.0;
    
    
    initialize_configuration(N, m, lc, L, r, v, K, T);
    
    
    computeAccelerations( N, r, a, f, L, U, m, A, p, q);
    
    
    save_configuration(N, el, r, v, "Results/init.xyz");
    
    const int n_bins = 80;
    double bin_size;
    double *gr;
    gr = new double [n_bins];
    initialize_gr( gr, n_bins, bin_size, L);
    
    
    
    //mean square displacement
    //initial positions
    double ** po;
    po = new double * [N];
    for (int i = 0; i < N; i++){
        po[i] = new double [3];
    }
    double * msd = new double [3];
    
    ofstream mfile;
    std::ostringstream mf;
    mf <<"Results/msd_" << Tt << ".txt" ;
    std::string msf = mf.str();
    mfile.open(msf.c_str());
    for (int i = 0; i < N; i++){
        for (int k = 0; k < 3; k++){
            po[i][k] = r[i][k];
        }
    }
    
    
    
    int cc = 0; int ngr = 0;
    for (int n = 0; n < nSteps; n++)
    {
        
        //loop over atoms
        for (int i = 0; i < N; i++){
            
            for (int k = 0; k < 3; k++){
                
                r[i][k] = r[i][k] + v[i][k]*dt + 0.5* a[i][k]* dt*dt;
                
                
                //periodic boundary conditions
                while (r[i][k] > 0.5*L)
                    r[i][k] = r[i][k] - L;
                
                while (r[i][k] < -0.5*L)
                    r[i][k] = r[i][k] + L;
                
                
                v[i][k] = v[i][k] + 0.5 * a[i][k]*dt;
                
            }
            
        }//end loop over atoms
        
        computeAccelerations( N, r, a, f, L, U, m, A, p, q);
        
        K = 0.0;
        for (int i = 0; i < N; i++)
        {
            for (int k = 0; k < 3; k++)
            {
                v[i][k] =  v[i][k] + 0.5 * a[i][k] * dt;
                
                K = K + 0.5 * m * v[i][k] * v[i][k];
            }
        }
        
        E = K+U;
        
        iT = instantaneous_temperature(N, K)/kb;
        
        
        
        //production steps
        if (n > int(0.5*nSteps))
        {
            
            //calculate averages
            ave_var(cc, UAvg, UVar, U);
            ave_var(cc, KAvg, KVar, K);
            ave_var(cc, EAvg, EVar, E);
            ave_var(cc, iTAvg,iTVar, iT);
            
            cc++;
            
            ngr +=1;
            update_gr (N, r, L, gr, bin_size );
            
            //if (n > int(0.75*nSteps))
            //dt  = dt/10.0;
        }
        
        
       
        if ( n % 150 == 0  )
            rescaleVelocities(N, K, T, v);
        
        cout << n << "\t" << iTAvg << "\t" << KAvg << "\t" << UAvg << "\t" << EAvg << "\n";
        
        //mean square dispalcement
        mean_square_displacement( N, r, po, msd );
        mfile << n <<"\t" << msd[0] << "\t" << msd[1] << "\t" << msd[2] << "\t" << msd[0]+msd[1]+msd[2] << "\n";
        
        
    }//end loop over nSteps
    
    //write to file
    hfile << Tt << "\t"
    << iTAvg << "\t" << iTVar << "\t" << KAvg << "\t" << KVar << "\t"
    << UAvg << "\t" << UVar << "\t" << EAvg << "\t" << EVar  <<   "\n" ;
    
    
    //save temperature configuration
    std::ostringstream finf;
    finf <<"Results/fin_"<<Tt<<".xyz" ;
    std::string fnf = finf.str();
    save_configuration(N, el, r, v, fnf.c_str());
    
    
    cout << Tt << "\t" << iTAvg << "\t" << KAvg << "\t" << UAvg << "\t" << EAvg << "\n";
    
    
    std::ostringstream gn;
    gn <<"Results/rd_"<<Tt<<".txt";
    std::string gn1 = gn.str();
    normalize_gr( N, ngr, gr, n_bins, bin_size, L, gn1.c_str());
    
    
    //  }//end loop over temperature
    
    
    
    
    free (r);
    free (v);
    free (f);
    free (a);
    free (po);
    
    hfile.close();
    
   
   
    return 0;
    
}
