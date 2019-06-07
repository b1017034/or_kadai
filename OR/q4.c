#include <stdio.h>
#include <math.h>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define ABS(x)   ((x)<0.0 ? -(x): (x))


double nablaLagrangian_x(double x, double a, double b, double c, double d);
double func(double x, double a, double b, double c, double d);

int main() {
    int k, maxloop, id, prmod;
    double dx, x, alpha, eps, a, b, c, d;
    FILE *fpw;

    fpw = fopen("output.text", "w");

    maxloop = 40000;

    alpha = 0.01;
    eps = 1.0e-3;
    prmod = 50;

    id=4;

    fprintf(stdout,"------- id = %d -------\n", id);
    a = (double)id + 1.0;
    b = (double)id + 2.0;
    c = (double)id + 3.0;
    d = (double)id + 5.0;

    x = a-1;
    for (k = 0; k < maxloop; k++){
        dx = -nablaLagrangian_x(x, a, b, c, d);

        x = x + alpha * dx;

        if( k % prmod == 0){
            fprintf(stdout, "k=%d, dx_k=%9.2le,",k, dx);
            fprintf(stdout, "k_{k+1}=%9.2le, func(x_{k+1})=%9.2le\n",x, func(x, a, b, c, d));
        }
        if(ABS(dx)<eps){
            break;
        }
    }
    fprintf(stdout, "k=%d, dx_k=%9.2le, ",k, dx);
    fprintf(stdout, "k_{k+1}=%9.2le, func(x_{k+1})=%9.2le\n",x, func(x, a, b, c, d));
    fprintf(fpw, "id=%d, x=%9.2le, F(x)=%9.2le, k=%d\n",id, x, func(x, a, b, c, d), k);

    fclose(fpw);
    return 0;
}

double nablaLagrangian_x(double x, double a, double b, double c, double d){
    return (-a + x)*(-b + x)*(-c + x) + (-a + x)*(-b + x)*(-d + x) + (-a + x)*(-c + x)*(-d + x) + (-b + x)*(-c + x)*(-d + x);
}

double func(double x, double a, double b, double c, double d){
    return (x - a) * (x - b) * (x - c) * (x - d);
}