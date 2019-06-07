#include <stdio.h>
#include <math.h>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define ABS(x)   ((x)<0.0 ? -(x): (x))


double nablaLagrangian_x(double x, double lambda, double a, double b, double c);
double nablaLagrangian_lambda(double x, double lambda, double a, double b, double c);
double func(double x, double a, double b, double c);

int main() {
    int k, maxloop, id, prmod;
    double dx, dlambda, x, lambda, alpha, eps, a, b, c;
    FILE *fpw;

    fpw = fopen("output.text", "w");

    maxloop = 40000;

    alpha = 0.02;
    eps = 1.0e-3;
    prmod = 50;

    id=34;

    fprintf(stdout,"------- id = %d -------\n", id);
    a = (double)id + 1.0;
    b = (double)id / 1000.0;
    c = (double)id / 10.0;

    x = 0.0;
    lambda = 1.0;
    for (k = 0; k < maxloop; k++){
        dx = -nablaLagrangian_x(x, lambda, a, b, c);
        dlambda = -nablaLagrangian_lambda(x, lambda, a, b, c);

        x = x + alpha * dx;
        lambda = MAX(lambda - alpha * dlambda, 0.0);

        if( k % prmod == 0){
            fprintf(stdout, "k=%d, dx_k=%9.2le,",k, dx);
            fprintf(stdout, "dlambda_k=%9.2le, lambda_{k+1}=%9.2le, ",dlambda, lambda);
            fprintf(stdout, "k_{k+1}=%9.2le, func(x_{k+1})=%9.2le\n",x, func(x, a, b, c));
        }
        if((ABS(dx)<eps)&&(ABS(dlambda)<eps)){
            break;
        }
    }
    fprintf(stdout, "k=%d, dx_k=%9.2le, dlambda_k=%9.2le, lambda_{k+1}=%9.2le, ",k, dx, dlambda, lambda);
    fprintf(stdout, "k_{k+1}=%9.2le, func(x_{k+1})=%9.2le\n",x, func(x, a, b, c));
    fprintf(fpw, "id=%d, x=%9.2le,F(x)=%9.2le, k=%d\n",id, x, func(x, a, b, c), k);

    fclose(fpw);
    return 0;
}

double nablaLagrangian_x(double x, double lambda, double a, double b, double c){
    return (2.0 * (x - a) + 3.0 * x * x + lambda);
}

double nablaLagrangian_lambda(double x, double lambda, double a, double b, double c){
    return x - c;
}

double func(double x, double a, double b, double c){
    return (x - a) * (x - a) + b * pow(x, 3.0);
}