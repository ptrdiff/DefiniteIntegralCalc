using System;
using MathNet.Numerics;

namespace DefiniteIntegralCalc
{
    class CalcDefInteg
    {
        private static double Function(double x) => (2 * Math.Cos(3.5 * x) * Math.Exp(5.0 * x / 3.0) + 3.0 * Math.Sin(1.5 * x) * Math.Exp(-4.0 * x) + 3.0);
        private static double Integral(double x) => (2 * Math.Cos(3.5 * x) * Math.Exp(5.0 * x / 3.0) + 3.0 * Math.Sin(1.5 * x) * Math.Exp(-4.0 * x) + 3.0) / (Math.Pow(x - 1.5, 1.0 / 5.0) * Math.Pow(2.3 - x, 0.0));
        private const double B = 2.3;
        private const double A = 1.5;

        public static double RealResult()
        {
            return Integrate.OnClosedInterval(Integral, 1.5, 2.3, 1e-15);
        }

        public static double RiemannSum(int n)
        {
            double h = (B - A) / n;
            double intSum = 0;
            for (int i = 1; i <= n; ++i)
            {
                intSum += Integral(A + (i - 1 / 2) * h);
            }
            return h * intSum;
        }

        public static double NewtonCotesIQF(double a = A, double b = B)
        {
            double x1 = a;
            double x2 = (a + b) / 2.0;
            double x3 = b;

            double nu0 = (5.0 / 4.0) * (Math.Pow(b - 1.5, 4.0 / 5.0) - Math.Pow(a - 1.5, 4.0 / 5.0));
            double nu1 = (5.0 / 72.0) * (Math.Pow(b - 1.5, 4.0 / 5.0) * (8 * b + 15) - Math.Pow(a - 1.5, 4.0 / 5.0) * (8 * a + 15));
            double nu2 = (5.0 / 336.0) * (Math.Pow(b - 1.5, 4.0 / 5.0) * (24.0 * Math.Pow(b, 2.0) + 40.0 * b + 75) - Math.Pow(a - 1.5, 4.0 / 5.0) * (24.0 * Math.Pow(a, 2.0) + 40.0 * a + 75));

            double a1 = (nu2 - nu1 * (x2 + x3) + nu0 * x2 * x3) / ((x2 - x1) * (x3 - x1));
            double a2 = -(nu2 - nu1 * (x1 + x3) + nu0 * x1 * x3) / ((x2 - x1) * (x3 - x2));
            double a3 = (nu2 - nu1 * (x2 + x1) + nu0 * x2 * x1) / ((x3 - x2) * (x3 - x1));

            return a1 * Function(x1) + a2 * Function(x2) + a3 * Function(x3);
        }

        public static double NewtonCotesCQF(int n)
        {
            double h = (B - A) / n;
            double Sum = 0;
            for (int i = 0; i < n; ++i)
            {
                double a = A + i * h;
                double b = A + (i + 1) * h;
                Sum += NewtonCotesIQF(a,b);
            }

            return Sum;
        }

        public static double NewtonCotesCQFHalf(int n , double eps)
        {
            int L = 2;
            double Sum1 = NewtonCotesCQF(n);
            double Sum2 = NewtonCotesCQF(n * L);
            double Sum3 = NewtonCotesCQF(n * L * L);

            double m = -Math.Log((Sum3 - Sum2) / (Sum2 - Sum1)) / Math.Log(L);
            double richardson = (Sum3 - Sum2) / (Math.Pow(L, m) - 1);
            int l = n * L * L;

            while (Math.Abs(richardson) > eps)
            {
                Sum1 = Sum2;
                Sum2 = Sum3;
                l *= L;
                Sum3 = NewtonCotesCQF(l);

                m = -Math.Log((Sum3 - Sum2) / (Sum2 - Sum1)) / Math.Log(L);
                richardson = (Sum3 - Sum2) / (Math.Pow(L, m) - 1);
            }

            return Sum3 + richardson;
        }

        public static double NewtonCotesCQFOptimal(double eps)
        {
            int L = 2;
            double Sum1 = NewtonCotesCQF(1);
            double Sum2 = NewtonCotesCQF(2);
            double Sum3 = NewtonCotesCQF(4);
            double m = 3;
            double richardson = (Sum3 - Sum2) / (Math.Pow(L, m) - 1);

            while (Math.Abs(richardson) > eps)
            {
                double hopt = Math.Pow((B - A) / L * (eps * (1 - Math.Pow(L, -m))) / (Math.Abs(Sum2 - Sum1)), 1.0 / m);
                int n = (int)Math.Ceiling((B - A) / (hopt));

                Sum1 = Sum2;
                Sum2 = Sum3;
                Sum3 = NewtonCotesCQF(n); 

                richardson = (Sum3 - Sum2) / (Math.Pow(L, m) - 1);
            }

            return Sum3 + richardson;
        }



    }
}
