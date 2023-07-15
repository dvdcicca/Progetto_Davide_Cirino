#include <math.h>
#include "io.h"

//Definizione dell'inizializzatore delle posizioni'
void r_initiator(int n, double rho, double R[][3], char *title){
    int i = 0;
    double a = pow(1./rho, 1./3);
    for(int k = 0; k<=n-1;k++){
        for(int l = 0; l<=n-1; l++){
            for(int m = 0; m<=n-1; m++){
                R[i][0]=k*a;
                R[i][1]=l*a;
                R[i][2]=m*a;
                outputCC( R[i][0],  R[i][1],  R[i][2], title);
                i++;  
            }
        }
    }
}


void r_initiator_BCC(int n, double rho, double R[][3], char *title){
    int i = 0;
    double a = pow(2./rho, 1./3);
    for(int k = 0; k<=n-1;k++){
        for(int l = 0; l<=n-1; l++){
            for(int m = 0; m<=n-1; m++){
                R[i][0]=k*a;
                R[i][1]=l*a;
                R[i][2]=m*a;
                outputCC( R[i][0],  R[i][1],  R[i][2], title);
                i++;  
                R[i][0]=k*a+1./2;
                R[i][1]=l*a+1./2;
                R[i][2]=m*a+1./2;
                outputCC( R[i][0],  R[i][1],  R[i][2], title);
                i++;  
            }
        }
    }
}