#pragma once

namespace mrcz
{
    // Evaluate Legendre P_m(x) and P_{m-1}(x) via recurrence
    static inline void legendre_Pm_and_Pm1(int m,double x,double &Pm,double &Pm1)
    {
        double P0=1;
        double P1=x;
        if (m==0)
        {
            Pm=P0;
            Pm1=0;
            return;
        }
        if (m==1)
        {
            Pm=P1;
            Pm1=P0;
            return;
        }

        for (int k=2;k<=m;k++)
        {
            double Pk=((2*k-1)*x*P1-(k-1)*P0)/k;
            P0=P1;
            P1=Pk;
        }
        Pm=P1; // P_m
        Pm1=P0; // P_{m-1}
    }

    // First derivative P'_m(x) using: P'_m = m/(x^2-1) * (x P_m - P_{m-1})
    static inline double legendre_dPm(int m,double x,double Pm,double Pm1)
    {
        return (m*(x*Pm-Pm1))/(x*x-1);
    }

    // Second derivative from Legendre ODE: (1-x^2)P'' - 2x P' + m(m+1)P = 0
    static inline double legendre_ddPm(int m,double x,double Pm,double dPm)
    {
        return (2*x*dPm-m*(m+1)*Pm)/(1-x*x);
    }

    template<typename T>
    vector<T> compute_GLL_nodes(uint32_t n)
    {
        int max_iter=100;
        double tol=1e-14;

        if (n<2)
            return {0};
        if (n==2) 
            return {-1,1};
        int m=n-1; // we need roots of P'_m
        vector<T> x(n);
        x[0]=-1;
        x[n-1]=1;

        // Initial guesses: Chebyshev-Lobatto nodes are a good start
        // x_k ~ -cos(pi*k/(n-1)), k=0..n-1 (already includes endpoints)
        // We'll refine only interior k = 1..n-2
        for (uint32_t k=1;k<=n-2;k++)
        {
            double xi=-cos(M_PI*k/m); // initial guess
            for (int it=0;it<max_iter;it++)
            {
                double Pm,Pm1;
                legendre_Pm_and_Pm1(m,xi,Pm,Pm1);
                double dPm=legendre_dPm(m,xi,Pm,Pm1);
                double ddPm=legendre_ddPm(m,xi,Pm,dPm);
                double dx=dPm/ddPm;
                xi-=dx;
                if (std::abs(dx)<tol)
                    break;
            }
            x[k]=xi;
        }
        
        x.front()=-1;
        x.back()=1;
        return x;
    }
}