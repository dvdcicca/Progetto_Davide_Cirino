#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include "headers/cubegnr.h" //funzioni per creare il mio reticolo cubico e che mi include anche le funzioni per stampare in un file

//Definizione delle costanti
#define PI 3.14159265358979323846
#define N 64
#define h 6.0596
#define rho 0.02186
#define eps 10.22
#define sigma 2.556
#define a2 5
#define stepn 20000
#define stepb 5000


//Dichiarazione delle varie funzioni
double dist(double R[][3], int , int , double );
void dist_2(double R[][3], double Rij[3], double , int , int );
double V(double , double );
void montecarlo(double R[][3], double, double, double);
double u(double , double);
double u_prime(double, double);
double u2_prime(double, double, double);
double u_second(double, double);
double u2_second(double, double, double);
double log_psi(double R[][3], double , double );
double u_2(double , double , double );

int main(){

    flush("CC.dat"); //pulisce il file se esistente
    flush("qmc.dat");
    int n = 4;
    double L = pow(N/rho, 1./3);
    double R[N][3];
    double r=1;
    double D=0.27;
    double alpha = 2.5;
    r_initiator(n, rho, R, "CC.dat");
    montecarlo(R, D, L, alpha);
    return 0;
}

//Definizione della funzione contenente il metodo M(RT)²
void montecarlo(double R[][3], double D, double L, double alpha){
    int accept, reject;
    double E, wf_old, wf_new;
    accept = reject = 0;
    E = 0;
    wf_old = log_psi(R, alpha, L);
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            double r = dist(R, i, j, L);
            E+= V(r, L);
        }
    }
    for(int i = 0; i < N; i++){
        double K_i = 0;
        for(int j = 0; j < N; j++){
            double r_ij = dist(R, i, j, L);
            if(j != i){
                K_i -= u2_prime(r_ij, alpha, L) + u2_second(r_ij, alpha, L)/2;
                for(int l = 0; l < N; l++){
                    double R_ij[3];
                    double R_il[3];
                    dist_2(R, R_ij, L, i, j);
                    dist_2(R, R_il, L, i, l);
                    double r_il = dist(R, i, l, L);
                    if(l != i){
                        K_i += (R_ij[0]*R_il[0]+R_ij[1]*R_il[1]+R_ij[2]*R_il[2])*u2_prime(r_ij, alpha, L)*u2_prime(r_il, alpha, L)/4;
                    }
                }
            }
        }
        E -= h*K_i;
    }

    output(0, E, 0.0, "qmc.dat");
    for (int i = 1; i <stepn; i++){
        double E_new = 0;
        double R_new[N][3];

        for (int j = 0; j < N; j++){
            for (int k = 0; k < 3; k++){
                R_new[j][k] = R[j][k] + D * (rand() / (RAND_MAX + 1.) - 0.5);
            }
        }
        wf_new = log_psi(R_new, alpha, L);

        for(int i = 0; i < N; i++){
            for(int j = i+1; j < N; j++){
                double r = dist(R_new, i, j, L);
                E_new+= V(r, L);
            }
        }
        for(int i = 0; i < N; i++){
            double K_i = 0;
            for(int j = 0; j < N; j++){
                double r_ij = dist(R, i, j, L);
                if(j != i){
                    K_i -= u2_prime(r_ij, alpha, L) + u2_second(r_ij, alpha, L)/2;
                    for(int l = 0; l < N; l++){
                        double R_ij[3];
                        double R_il[3];
                        dist_2(R_new, R_ij, L, i, j);
                        dist_2(R_new, R_il, L, i, l);
                        double r_il = dist(R_new, i, l, L);
                        if(l != i){
                            K_i += (R_ij[0]*R_il[0]+R_ij[1]*R_il[1]+R_ij[2]*R_il[2])*u2_prime(r_ij, alpha, L)*u2_prime(r_il, alpha, L)/4;
                        }
                    }
                }
            }
            E_new -= h*K_i;
        }


        if (2.0*(wf_new - wf_old) > log(rand() / (RAND_MAX + 1.))){
            for (int j = 0; j < N; j++){
                for (int k = 0; k < 3; k++){
                    R[j][k] = R_new[j][k];
                }
            }
            E = E_new;
            wf_old = wf_new;
            accept++;

        }
        else{
            reject++;
        }
        output(i, E, accept*1./i, "qmc.dat");
    }
}

//Funzione che calcola la distanza tra due particelle
double dist(double R[][3], int p, int q, double L){    
    return sqrt((R[p][0]-R[q][0]-L*rint((R[p][0] - R[q][0])/L))*(R[p][0]-R[q][0]-L*rint((R[p][0] - R[q][0])/L))+(R[p][1]-R[q][1]-L*rint((R[p][1] - R[q][1])/L))*(R[p][1]-R[q][1]-L*rint((R[p][1] - R[q][1])/L))+(R[p][2]-R[q][2]-L*rint((R[p][2] - R[q][2])/L))*(R[p][2]-R[q][2]-L*rint((R[p][2] - R[q][2])/L)));
} 

//Funzione che calcola la distanza sui tre assi tra due particelle
void dist_2(double R[][3], double Rpq[3], double L, int p, int q){
    Rpq[0] = R[p][0] - R[q][0]-L*rint((R[p][0] - R[q][0])/L);
    Rpq[1] = R[p][1] - R[q][1]-L*rint((R[p][1] - R[q][1])/L);
    Rpq[2] = R[p][2] - R[q][2]-L*rint((R[p][2] - R[q][2])/L);
}

//Definizione del potenziale di coppia
double V(double r, double L){
    double r6 = (sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r);
    double L6 = (2*sigma/L)*(2*sigma/L)*(2*sigma/L)*(2*sigma/L)*(2*sigma/L)*(2*sigma/L);
    if(r <= L/2){
        return 4*eps*r6*(r6 - 1) - 4*eps*L6*(L6 - 1);
    }
    else{
        return 0;
    }
}

//Definizione della funzione u
double u(double r, double a1){
    return (a1/r)*(a1/r)*(a1/r)*(a1/r)*(a1/r);
}

//Definizione della funzione u con correzione per le condizioni al contorno
double u_2(double r, double a1, double L){
    if(r <= L/2){
        return u(r, a1) + u(L-r, a1) - 2*u(L/2, a1);
    }
    else{
        return 0;
    }
}

//Definizione derivata prima della funzione u
double u_prime(double r, double a1){
    return -5*(a1*a1*a1*a1*a1)*1./(r*r*r*r*r*r);
}

double u2_prime(double r, double a1, double L){
    if (r<=L/2){
        return u_prime(r, a1) - u_prime(L-r, a1);
    }
    else{
        return 0;
    }
}

//Definizione derivata seconda della funzione u
double u_second(double r, double a1){
    return 30*(a1*a1*a1*a1*a1)*1./(r*r*r*r*r*r*r);
}

double u2_second(double r, double a1, double L){
    if (r<=L/2){
        return u_second(r, a1) + u_second(L-r, a1);
    }
    else{
        return 0;
    }
}

//Definizione del logaritmo della funzione d'onda
double log_psi(double R[][3], double a1, double L){
    double u_f = 0;
    for(int i = 0; i < N-1; i++){
        for(int j = i+1; j<N; j++){
            double r_ij = dist(R, i, j, L);
            u_f += u_2(r_ij, a1, L);
        }
    }
    return (-1./2)*u_f;
}

//Definizione della funzione che calcola l'energia locale


