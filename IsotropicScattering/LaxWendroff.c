#include <stdint.h>
#include <stdio.h>
#include <math.h>

const float dx = 0.01;
const double g = 1.4;
const int N = 100;
// first order

int main()
{                    
                // declaration of state vector array, flux array, average rightwards and average leftward flux array           
                double Q[3][N];
                double F[3][N];
                double F_plus_half[3][N];
                double F_minus_half[3][N];
                double Q_intermediate_plus_half[3][N];
                double Q_intermediate_minus_half[3][N];
              
                double dt;
                double end_t = 0.2;
                double t;

                double x[N];
                
                // arrays for each variable 

                double r[N];    // density
                double v[N];    // velocity
                double p[N];    // pressure
                double E[N];    // total energy density 

                double rho;
                double total_energy;
                double velocity;        
                double pressure;       
            
                double flux_1;
                double flux_2;
                double flux_3;

                double max_v;
                double c_s;

                double max_a;
        
                // initialization of variable arrays and state vectors using initial conditions
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
                

                // algorithm to determine the maximum value for a (v + c_s), and determine dt using CFL condition.
                c_s = sqrt(g*p[0]/r[0]);
                max_a = fabs(v[0]) + fabs(c_s);
        
                for(int d = 0; d < N; d++) {
                    if(max_a < fabs(v[d])+fabs(sqrt(g*p[d]/r[d]))){
                        max_a = fabs(v[d])+fabs(sqrt(g*p[d]/r[d]));
                    }
                }  

                dt = 0.5*dx/max_a;  // constant of 0.5 used as CFL condition constant

                t=0;

                while(t<end_t){             

                    // boundary conditions

                    r[0]=1;
                    p[0]=1;
                    v[0]=0.75;
                    E[0]= p[0]/((g-1)*r[0]) + 0.5*pow(v[0],2);

                    r[N-1]=0.125;
                    p[N-1]=0.1;
                    v[N-1]=0;
                    E[N-1]= p[N-1]/((g-1)*r[N-1]) + 0.5*pow(v[N-1],2);

                    
                    // flux of current state vectors determined 

                    for(int d = 0; d < N; d++) {
                        flux_1= r[d]*v[d];
                        flux_2= r[d]*v[d]*v[d] + p[d];
                        flux_3= v[d]*(E[d]*r[d] + p[d]);

                        F[0][d] = flux_1;
                        F[1][d] = flux_2;
                        F[2][d] = flux_3;
                    }

                    // reconstruction of flux and update of new state vector 

                    for(int d = 1; d < N-1; d++) {
                        Q_intermediate_plus_half[0][d] = 0.5*(Q[0][d]+Q[0][d+1]) + 0.5*(dt/dx)*(F[0][d]-F[0][d+1]);
                        Q_intermediate_minus_half[0][d] = 0.5*(Q[0][d-1]+Q[0][d]) + 0.5*(dt/dx)*(F[0][d-1]-F[0][d]);

                        Q_intermediate_plus_half[1][d] = 0.5*(Q[1][d]+Q[1][d+1]) + 0.5*(dt/dx)*(F[1][d]-F[1][d+1]);
                        Q_intermediate_minus_half[1][d] = 0.5*(Q[1][d-1]+Q[1][d]) + 0.5*(dt/dx)*(F[1][d-1]-F[1][d]);

                        Q_intermediate_plus_half[2][d] = 0.5*(Q[2][d]+Q[2][d+1]) + 0.5*(dt/dx)*(F[2][d]-F[2][d+1]);
                        Q_intermediate_minus_half[2][d] = 0.5*(Q[2][d-1]+Q[2][d]) + 0.5*(dt/dx)*(F[2][d-1]-F[2][d]);

                        

                        // determining average fluxes of Q intermediates 

                        F_plus_half[0][d] = Q_intermediate_plus_half[1][d];
                        F_minus_half[0][d] = Q_intermediate_minus_half[1][d];

                        // (used to calculate pressure in the flux terms)

                        rho = Q_intermediate_plus_half[0][d];
                        total_energy = Q_intermediate_plus_half[2][d]/rho;
                        velocity = Q_intermediate_plus_half[1][d];
                        pressure = (g-1)*rho*(total_energy-0.5*pow(velocity,2));

                        F_plus_half[1][d] = (pow(Q_intermediate_plus_half[1][d],2)/Q_intermediate_plus_half[0][d]) + pressure;
                        F_minus_half[1][d] = (pow(Q_intermediate_minus_half[1][d],2)/Q_intermediate_minus_half[0][d]) + pressure;

                        rho = Q_intermediate_plus_half[0][d];
                        total_energy = Q_intermediate_plus_half[2][d]/rho;
                        velocity = Q_intermediate_plus_half[1][d];
                        pressure = (g-1)*rho*(total_energy-0.5*pow(velocity,2));

                        F_plus_half[2][d] = Q_intermediate_plus_half[1][d]*(Q_intermediate_plus_half[2][d]+pressure)/Q_intermediate_plus_half[0][d];
                        F_minus_half[2][d] = Q_intermediate_minus_half[1][d]*(Q_intermediate_minus_half[2][d]+pressure)/Q_intermediate_minus_half[0][d];
                        

                        Q[0][d] = Q[0][d] - (dt/dx)*(F_plus_half[0][d]-F_minus_half[0][d]);
                        Q[1][d] = Q[1][d] - (dt/dx)*(F_plus_half[1][d]-F_minus_half[1][d]);
                        Q[2][d] = Q[2][d] - (dt/dx)*(F_plus_half[2][d]-F_minus_half[2][d]);

                        // setting variables for time t using updated state vectors 

                        r[d]=Q[0][d];
                        v[d]=Q[1][d]/r[d];
                        E[d]=Q[2][d]/r[d];
                        p[d]=(g-1)*r[d]*(E[d]-0.5*pow(v[d],2));
                    
                    }

                        // dt determined by finding CFL

                        c_s = sqrt(g*p[0]/r[0]);
                        max_a = fabs(v[0]) + fabs(c_s);

                        for(int d = 0; d < N; d++) {
                            if(max_a < fabs(v[d])+fabs(sqrt(g*p[d]/r[d]))){
                                max_a = fabs(v[d])+fabs(sqrt(g*p[d]/r[d]));
                            }
                        }  

                        dt = 0.5*dx/max_a;

                        // time updated by dt (smallest time interval possible)
                        t=t+dt;
}

        // variables printed for t = t
        for(int i = 0; i < N; i++)
            printf("%f, %f, %f, %f, %f\n", x[i], r[i], v[i], p[i], E[i]);

        return 0;
}