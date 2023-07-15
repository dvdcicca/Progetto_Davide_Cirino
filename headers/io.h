#include <stdio.h>

void output(double , double , double , double , char *);
void flush(char *);



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