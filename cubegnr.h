#include <math.h>
#include <stdio.h>
void r_initiator(int , double , double R[][3] );
void r_initiator_BCC(int , double , double R[][3] );
void outputCC(double , double , double );


void outputCC(double x1, double x2, double x3){
    FILE *fp;
    fp = fopen("CC2.dat", "a");
    fprintf(fp, "%lf %lf %lf\n", x1, x2, x3);
    fclose(fp);
    return;
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