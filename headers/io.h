#include <stdio.h>

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


void outputCC(double x1, double x2, double x3, char *title){
    FILE *fp;
    fp = fopen(title, "a");
    fprintf(fp, "%lf %lf %lf\n", x1, x2, x3);
    fclose(fp);
    return;
}