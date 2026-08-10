#pragma once

using namespace std;

namespace mrcz
{
    // 通过迭代，计算坐标为x时的m阶和m-1阶勒让德多项式的值 // P_m(x)和P_{m-1}(x)
    static inline void legendre_Pm_and_Pm1(int m,double x,double &Pm,double &P_prev)
    {
        if (m==0)
        {
            Pm=1;
            P_prev=0;
            return;
        }
        Pm=x;
        P_prev=1;
        for (int k=2;k<=m;k++)
        {
            double P_next=((2*k-1)*x*Pm-(k-1)*P_prev)/k;
            P_prev=Pm;
            Pm=P_next;
        }
    }

    // 计算坐标为x时m阶勒让德多项式的一阶导数 // P'_m(x) = m*(x*P_m-P_{m-1})/(x^2-1) 
    static inline double legendre_dPm(int m,double x,double Pm,double P_prev)
    {
        return (m*(x*Pm-P_prev))/(x*x-1);
    }

    // 计算坐标为x时m阶勒让德多项式的二阶导数 // P''_m(x) = (2*x*P'_m-m*(m+1)*P_m)/(1-x^2)
    static inline double legendre_ddPm(int m,double x,double Pm,double dPm)
    {
        return (2*x*dPm-m*(m+1)*Pm)/(1-x*x);
    }

    template<typename T>
    vector<T> compute_GLL_nodes(uint32_t n)
    {
        int max_iter=100;
        double tol=1e-14;

        if (n==0)
            return {};
        if (n==1)
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
                double Pm,P_prev;
                legendre_Pm_and_Pm1(m,xi,Pm,P_prev);
                double dPm=legendre_dPm(m,xi,Pm,P_prev);
                double ddPm=legendre_ddPm(m,xi,Pm,dPm);
                double dx=dPm/ddPm;
                xi-=dx;
                if (abs(dx)<tol)
                    break;
            }
            x[k]=xi;
        }
        return x;
    }
}