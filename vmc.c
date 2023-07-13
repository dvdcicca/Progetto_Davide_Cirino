#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>

#define PI 3.14159265358979323846
#define N 1000

//Dichiarazione delle varie funzioni
void r_initiator(int , double , double R[][3] );
void r_initiator_BCC(int , double , double R[][3] );
void outputCC(double , double , double );
double dist(double R[][3], int , int , double );
void dist_2(double R[][3], double Rij[3], double , int , int );
double Vol(double, double );
double V(double , double , double , double );
double W(double , double , double , double);
void output(double , double , double , double , char *);
void flush(char *);

int main(){

    flush("CC2.dat");
    flush("monte_lq.dat");

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


//Definizione dell'inizializzatore delle posizioni'
void r_initiator(int n, double rho, double R[][3]){
    int i = 0;
    double a = pow(1./rho, 1./3);
    for(int k = 0; k<=n-1;k++){
        for(int l = 0; l<=n-1; l++){
            for(int m = 0; m<=n-1; m++){
                R[i][0]=k*a;
                R[i][1]=l*a;
                R[i][2]=m*a;
                outputCC( R[i][0],  R[i][1],  R[i][2]);
                i++;  
            }
        }
    }
}


void r_initiator_BCC(int n, double rho, double R[][3]){
    int i = 0;
    double a = pow(2./rho, 1./3);
    for(int k = 0; k<=n-1;k++){
        for(int l = 0; l<=n-1; l++){
            for(int m = 0; m<=n-1; m++){
                R[i][0]=k*a;
                R[i][1]=l*a;
                R[i][2]=m*a;
                outputCC( R[i][0],  R[i][1],  R[i][2]);
                i++;  
                R[i][0]=k*a+1./2;
                R[i][1]=l*a+1./2;
                R[i][2]=m*a+1./2;
                outputCC( R[i][0],  R[i][1],  R[i][2]);
                i++;  
            }
        }
    }
}


//Definizione della funzione per stampare in un file le posizioni iniziali
void outputCC(double x1, double x2, double x3){
    FILE *fp;
    fp = fopen("CC2.dat", "a");
    fprintf(fp, "%lf %lf %lf\n", x1, x2, x3);
    fclose(fp);
    return;
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


void output(double t, double P, double E, double acc, char *title){
    FILE *fp;
    fp = fopen(title, "a");
    fprintf(fp, "%lf %lf %lf %lf \n", t, P, E, acc);
    fclose(fp);
    return;
}

void flush(char *title){
FILE *fp;
fp = fopen(title, "w");
fclose(fp);
}
