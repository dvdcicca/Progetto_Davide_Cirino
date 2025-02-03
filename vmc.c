#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "headers/cubegnr.h" //funzioni per creare il mio reticolo cubico e che mi include anche le funzioni per stampare in un file


#define PI 3.14159265358979323846
#define N 64

//Dichiarazione delle varie funzioni
double dist(double R[][3], int , int , double );
void dist_2(double R[][3], double Rij[3], double , int , int );
double Vol(double, double );
double V(double , double , double , double );
void montecarlo(double R[][3], double , double *, double, double, double, double);

int main(){

    flush("CC.dat");
    int n = 4;
    double rho = 0.0218;
    double L = pow(N/rho, 1./3);
    double R[N][3];
    double r=1;
    double D=0.18;
    double eps = 10.22;
    double sigma = 2.556;
    double h = 6.0596;
    double alpha = 2.5;
    r_initiator(n, rho, R, "CC.dat");
    
    return 0;
}

void montecarlo(double R[][3], double D, double *cumulative_E, double L, double alpha, double sigma, double epsilon){
    int accept, reject;
    double E, wf_old, wf_new;
    accept = reject = 0;
    E = 0;
    wf_old = wf(R, alpha, L);
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            double r = dist(R, i, j, L);
            E+= V(r, L, sigma, epsilon);
        }
    }
    for (int i = 1; i < 1000000; i++){

        double R_new[N][3];

        for (int j = 0; j < N; j++){
            for (int k = 0; k < 3; k++){
                R_new[j][k] = R[j][k] + D * (rand() / (RAND_MAX + 1.) - 0.5);
            }
        }
        wf_new = wf(R_new, alpha, L);

        if ((wf_new / wf_old) > rand() / (RAND_MAX + 1.)){
            for (int j = 0; j < N; j++){
                for (int k = 0; k < 3; k++){
                    R[j][k] = R_new[j][k];
                }
            }
            wf_old = wf_new;
            accept++;
        }
        else{
            reject++;
        }
    }
}

double dist(double R[][3], int p, int q, double L){    
    return sqrt((R[p][0]-R[q][0]-L*rint((R[p][0] - R[q][0])/L))*(R[p][0]-R[q][0]-L*rint((R[p][0] - R[q][0])/L))+(R[p][1]-R[q][1]-L*rint((R[p][1] - R[q][1])/L))*(R[p][1]-R[q][1]-L*rint((R[p][1] - R[q][1])/L))+(R[p][2]-R[q][2]-L*rint((R[p][2] - R[q][2])/L))*(R[p][2]-R[q][2]-L*rint((R[p][2] - R[q][2])/L)));
} 


void dist_2(double R[][3], double Rpq[3], double L, int p, int q){
    Rpq[0] = R[p][0] - R[q][0]-L*rint((R[p][0] - R[q][0])/L);
    Rpq[1] = R[p][1] - R[q][1]-L*rint((R[p][1] - R[q][1])/L);
    Rpq[2] = R[p][2] - R[q][2]-L*rint((R[p][2] - R[q][2])/L);
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

//Definizione della funzione u
double u(double r_ij, double alpha1){
    return (alpha1/r_ij)*(alpha1/r_ij)*(alpha1/r_ij)*(alpha1/r_ij)*(alpha1/r_ij);
}

//Definizione della funzione u con correzione per le condizioni al contorno
double u_2(double r_ij, double alpha1, double L){
    if(r_ij <= L/2){
        return u(r_ij, alpha1) + u(L-r_ij, alpha1) - 2*u(L/2, alpha1);
    }
    else{
        return 0;
    }
}

//Definizione derivata prima della funzione u
double u_prime(double r_ij, double alpha1, double L){
    return -5*u_2(r_ij, alpha1, L)/r_ij;
}


//Definizione derivata seconda della funzione u
double u_second(double r_ij, double alpha1, double L){
    return 30*u_2(r_ij, alpha1, L)/(r_ij*r_ij);
}

//Definizione della funzione d'onda
double psi(double R[][3], double alpha1, double L){
    double u_f = 0;
    for(int i = 0; i < N-1; i++){
        for(int j = i+1; j<N; j++){
            double r_ij = dist(R, i, j, L);
            u_f += u_2(r_ij, alpha1, L);
        }
    }
    return exp((-1./2)*u_f);
}

//Calcolo del modulo quadro della funzione d'onda
double wf(double R[][3], double alpha1, double L){
    double psii = psi(R, alpha1, L); 
    return psii*psii;
}


