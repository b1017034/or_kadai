#include <stdio.h>
#include <math.h>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define ABS(x)   ((x)<0.0? -(x): (x))

double Lagrange_x(double x, double lamda);
double Lagrange_lamda(double x, double lamda);

int main(){
    int max_loop = 1000;
    double x,lamda,eps;
    double dx,dlamda;

    x=1.5;
    lamda=1.0;
    eps = 1.0e-3;
    for(int i = 0; i < max_loop; i++){
        dx = Lagrange_x(x,lamda);
        dlamda = Lagrange_lamda(x,lamda);

        x=x+dx;
        lamda=MAX(lamda-dlamda, 0.0);

        printf("i=%d, x=%9.2le, lamda=%9.2le, dx=%9.2le, dlamda=%9.2le\n",i,x,lamda,dx,dlamda);

        if((ABS(dx) < eps)&&(ABS(dlamda) < eps)){
            break;
        }
    }

    return 0;
}


double Lagrange_x(double x, double lamda){
    return -(2.0*(x-2.0) + lamda);
}
double Lagrange_lamda(double x, double lamda){
    return -(x-1.0);
}