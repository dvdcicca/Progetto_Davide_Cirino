#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <io.h>
#include <cubegnr.h>


#define PI 3.14159265358979323846
#define N 1000

//Dichiarazione delle varie funzioni
double dist(double R[][3], int , int , double );
void dist_2(double R[][3], double Rij[3], double , int , int );
double Vol(double, double );
double V(double , double , double , double );
double W(double , double , double , double);

int main(){

    flush("CC2.dat");

    int n = 10;
    double rho = 0.8;
    double L = pow(N/rho, 1./3);
    double R[N][3];
    double dr=5;
    double M=L/(2*dr);
    double r=1;
    double D=0.18;
    double sigma = 1;
    double eps=1;
    double T1=1.1;
    double P=T1*rho;
    double c = 0;
    r_initiator(n, rho, R);
    
    return 0;
}



double dist(double R[][3], int p, int q, double L){    
    return sqrt((R[p][0]-R[q][0]-L*rint((R[p][0] - R[q][0])/L))*(R[p][0]-R[q][0]-L*rint((R[p][0] - R[q][0])/L))+(R[p][1]-R[q][1]-L*rint((R[p][1] - R[q][1])/L))*(R[p][1]-R[q][1]-L*rint((R[p][1] - R[q][1])/L))+(R[p][2]-R[q][2]-L*rint((R[p][2] - R[q][2])/L))*(R[p][2]-R[q][2]-L*rint((R[p][2] - R[q][2])/L)));
} 


void dist_2(double R[][3], double Rpq[3], double L, int p, int q){
    Rpq[0] = R[p][0] - R[q][0]-L*rint((R[p][0] - R[q][0])/L);
    Rpq[1] = R[p][1] - R[q][1]-L*rint((R[p][1] - R[q][1])/L);
    Rpq[2] = R[p][2] - R[q][2]-L*rint((R[p][2] - R[q][2])/L);
}


double Vol(double r, double dr){
    return (4./3)*PI*(r+dr/2)*(r+dr/2)*(r+dr/2)-(4./3)*PI*(r-dr/2)*(r-dr/2)*(r-dr/2);
}


//Definizione del potenziale di coppia
double V(double r, double L, double sigma, double epsilon){
    double r6 = (sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r);
    double L6 = (2*sigma/L)*(2*sigma/L)*(2*sigma/L)*(2*sigma/L)*(2*sigma/L)*(2*sigma/L);
    if(r <= L/2){
        return 4*epsilon*r6*(r6 - 1) - 4*epsilon*L6*(L6 - 1);
    }
    else{
        return 0;
    }
}


double W(double rij, double sigma, double L, double epsilon){
    double r6 = (sigma/rij)*(sigma/rij)*(sigma/rij)*(sigma/rij)*(sigma/rij)*(sigma/rij);
    double dV;
    if(rij <= L/2){
        dV = -6*(V(rij, L, sigma, epsilon) + 4*epsilon*r6*r6)/rij;
    }
    else{
        dV = 0;
    }
    return -dV*rij/N;

}