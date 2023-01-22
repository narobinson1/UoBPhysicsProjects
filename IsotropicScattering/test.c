#include <stdint.h>
#include <stdio.h>
#include <math.h>

const float dx = 0.01;
const double g = 1.4;
const int N = 102;
// first order

int main()
{                               
                double Q[3][N];
                double F[3][N];
                double F_plus_half[3][N];
                double F_minus_half[3][N];

              
                double dt;
                double end_t = 0.2;
                double t;

                double x[N];
                

                double r[N];
                double v[N];
                double p[N];
                double E[N];

               

                double flux_1;
                double flux_2;
                double flux_3;

                double max_v;
                double c_s;

                double max_a;
        

                for(int j=0; j < N; j++)
                {

                                x[j]=j*0.01;

                                if (x[j] <= 0.3)
                                {

                                r[j]=1;
                                p[j]=1;
                                v[j]=0.75;
                                // E_0 is total energy density
                                E[j]= p[j]/((g-1)*r[j]) + 0.5*pow(v[j],2);

                                Q[0][j] = r[j];
                                Q[1][j] = r[j]*v[j];
                                Q[2][j] = r[j]*E[j];

                                }
                                

                                if (x[j] > 0.3 & x[j]<=1)
                                {
                                r[j]=0.125;
                                p[j]=0.1;
                                v[j]=0;
                                E[j]= p[j]/((g-1)*r[j]) + 0.5*pow(v[j],2);

                                // density, momentum density, total energy for each x (j)
                                Q[0][j] = r[j];
                                Q[1][j] = r[j]*v[j];

                                // total energy is total energy density multiplied by r
                                Q[2][j] = r[j]*E[j];
                                }               
                }
                

//    for(int d = 0; d < N; d++) {
//            printf("%f, %f, %f\n", Q[0][d], Q[1][d], Q[2][d]);
//    }


                c_s = sqrt(g*p[0]/r[0]);
                max_a = fabs(v[0]) + fabs(c_s);
        
                for(int d = 1; d < N-1; d++) {
                    if(max_a < fabs(v[d])+fabs(sqrt(g*p[d]/r[d]))){
                        max_a = fabs(v[d])+fabs(sqrt(g*p[d]/r[d]));
                    }
                }  

                dt = 0.5*dx/max_a;

                t=0;

             

                while(t<end_t){             

                    // out flow boundaries 

                    

                    
                    

                    for(int d = 1; d < N-1; d++) {
                        flux_1= r[d]*v[d];
                        flux_2= r[d]*v[d]*v[d] + p[d];
                        flux_3= v[d]*(E[d]*r[d] + p[d]);

                        F[0][d] = flux_1;
                        F[1][d] = flux_2;
                        F[2][d] = flux_3;
                    }

                    for(int d = 1; d < N-1; d++) {
                        F_plus_half[0][d] = 0.5*(F[0][d]+F[0][d+1]) + 0.5*(dx/dt)*(Q[0][d]-Q[0][d+1]);
                        F_minus_half[0][d] = 0.5*(F[0][d-1]+F[0][d]) + 0.5*(dx/dt)*(Q[0][d-1]-Q[0][d]);

                        F_plus_half[1][d] = 0.5*(F[1][d]+F[1][d+1]) + 0.5*(dx/dt)*(Q[1][d]-Q[1][d+1]);
                        F_minus_half[1][d] = 0.5*(F[1][d-1]+F[1][d]) + 0.5*(dx/dt)*(Q[1][d-1]-Q[1][d]);

                        F_plus_half[2][d] = 0.5*(F[2][d]+F[2][d+1]) + 0.5*(dx/dt)*(Q[2][d]-Q[2][d+1]);
                        F_minus_half[2][d] = 0.5*(F[2][d-1]+F[2][d]) + 0.5*(dx/dt)*(Q[2][d-1]-Q[2][d]);
                        

                        Q[0][d] = Q[0][d] - (dt/dx)*(F_plus_half[0][d]-F_minus_half[0][d]);
                        Q[1][d] = Q[1][d] - (dt/dx)*(F_plus_half[1][d]-F_minus_half[1][d]);
                        Q[2][d] = Q[2][d] - (dt/dx)*(F_plus_half[2][d]-F_minus_half[2][d]);

                  

                        r[d]=Q[0][d];
                        v[d]=Q[1][d]/r[d];
                        E[d]=Q[2][d]/r[d];
                        p[d]=(g-1)*r[d]*(E[d]-0.5*pow(v[d],2));
                    
                    }

                        // SUM OF V AND C_MAX is HIGHEST TAKE ABS  and filter for the sums!!

                        
                        c_s = sqrt(g*p[0]/r[0]);
                        max_a = fabs(v[0]) + fabs(c_s);

                        for(int d = 1; d < N-1; d++) {
                            if(max_a < fabs(v[d])+fabs(sqrt(g*p[d]/r[d]))){
                                max_a = fabs(v[d])+fabs(sqrt(g*p[d]/r[d]));
                            }
                        }  

                        dt = 0.5*dx/max_a;

                        t=t+dt;
}


        for(int i = 1; i < N-1; i++)
            printf("%f, %f, %f, %f, %f\n", x[i], r[i], v[i], E[i], p[i]);

        return 0;
}