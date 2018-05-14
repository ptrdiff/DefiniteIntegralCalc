using System;
using System.Collections.Generic;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace DefiniteIntegralCalc
{
    class CalcDefInteg
    {
        private static double Function(double x) => (2 * Math.Cos(3.5 * x) * Math.Exp(5.0 * x / 3.0) + 3.0 * Math.Sin(1.5 * x) * Math.Exp(-4.0 * x) + 3.0);
        private static double Integral(double x) => (2 * Math.Cos(3.5 * x) * Math.Exp(5.0 * x / 3.0) + 3.0 * Math.Sin(1.5 * x) * Math.Exp(-4.0 * x) + 3.0) / (Math.Pow(x - 1.5, 1.0 / 5.0) * Math.Pow(2.3 - x, 0.0));
        private const double B = 2.3;
        private const double A = 1.5;
        private const double Al = 0.2;

        private static List<double> SolvePoly3(double a, double b, double c)
        {
            var Q = (a * a - 3 * b) / 9;
            var R = (2 * a * a * a - 9 * a * b + 27 * c) / 54;

            if (R * R < Q * Q * Q)
            {
                var t = Math.Acos(R / Math.Sqrt(Q * Q * Q)) / 3;
                return new List<double> {-2.0 * Math.Sqrt(Q) * Math.Cos(t ) - a / 3,
                                         -2.0 * Math.Sqrt(Q) * Math.Cos(t + (2 * Math.PI / 3)) - a / 3,
                                         -2.0 * Math.Sqrt(Q) * Math.Cos(t - (2 * Math.PI / 3)) - a / 3,
                };
            }
            else
            {
                var F = -R + Math.Sqrt(R * R - Q * Q * Q);
                var A = Math.Sign(F) * Math.Pow(Math.Abs(F), 1.0 / 3.0);
                if (A == 0.0)
                {
                    return new List<double> { -a / 3 };
                }
                else
                {
                    double B = Q / A, x1 = (A + B) - a / 3, x2 = -A - a / 3;
                    return Math.Abs(x2 * (x2 * (x2 + a) + b) + c) < 1e-6 ? new List<double> { x1, x2 } : new List<double> { x1 };
                }
            }
        }


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

            double nu(double i)
            {
                return (Math.Pow(b - A,i - Al) - Math.Pow(a - A, i - Al))/(i - Al);
            }

            double nu0 = nu(1.0);
            double nu1 = nu(2.0) + A * nu0;
            double nu2 = nu(3.0) + 2 * A * nu1 - A * A * nu0;

            double a1 = (nu2 - nu1 * (x2 + x3) + nu0 * x2 * x3) / ((x2 - x1) * (x3 - x1));
            double a2 = -(nu2 - nu1 * (x1 + x3) + nu0 * x1 * x3) / ((x2 - x1) * (x3 - x2));
            double a3 = (nu2 - nu1 * (x2 + x1) + nu0 * x2 * x1) / ((x3 - x2) * (x3 - x1));

            return a1 * Function(x1) + a2 * Function(x2) + a3 * Function(x3);
        }

        public static double GaussIQF(double a = A, double b = B)
        {
            double koeff = 1 / (1 - Al);
            double koeffB = Math.Pow(b - A, 1 - Al);
            double koeffA = Math.Pow(a - A, 1 - Al);
            double[] nu = new double[6];
            double nuGen(int n)
            {
                return koeff * (1.0 / (n / (1.0 - Al) + 1.0)) * (koeffB * Math.Pow(b, n) - koeffA * Math.Pow(a, n) + A * n * nu[n - 1]);
            }

            nu[0] = koeff * (koeffB - koeffA);
            nu[1] = nuGen(1);
            nu[2] = nuGen(2);
            nu[3] = nuGen(3);
            nu[4] = nuGen(4);
            nu[5] = nuGen(5);

            Vector<double> tmpa = Matrix<double>.Build.DenseOfArray(new double[,] { { nu[0], nu[1], nu[2] },
                                                                                    { nu[1], nu[2], nu[3] },
                                                                                    { nu[2], nu[3], nu[4] } })
                               .Solve(Vector<double>.Build.DenseOfArray(new double[] { -nu[3], -nu[4], -nu[5] }));

            var tmpX = SolvePoly3(tmpa[2], tmpa[1], tmpa[0]);

            Vector<double> tmpA = Matrix<double>.Build.DenseOfArray(new double[,] {
                { 1, 1, 1 }, { tmpX[0], tmpX[1],tmpX[2] },{ tmpX[0] * tmpX[0], tmpX[1] * tmpX[1],tmpX[2] * tmpX[2] }})
                               .Solve(Vector<double>.Build.DenseOfArray(new double[] { nu[0], nu[1], nu[2] }));

            return tmpA[0] * Function(tmpX[0]) + tmpA[1] * Function(tmpX[1]) + tmpA[2] * Function(tmpX[2]);
        }


        public static double MethodCQF(int n, bool isNewton = true)
        {
            Func<double,double,double> method = isNewton 
                ? new Func<double, double, double>(NewtonCotesIQF) 
                : new Func<double, double, double>(GaussIQF);

            double h = (B - A) / n;
            double Sum = 0;
            for (int i = 0; i < n; ++i)
            {
                double a = A + i * h;
                double b = A + (i + 1) * h;
                Sum += method(a,b);
            }

            return Sum;
        }

        public static double MethodCQFHalf(int n, int L, double eps, bool isNewton = true)
        {
            double Sum1 = MethodCQF(n,isNewton);
            double Sum2 = MethodCQF(n * L,isNewton);
            double Sum3 = MethodCQF(n * L * L,isNewton);

            double m = -Math.Log((Sum3 - Sum2) / (Sum2 - Sum1)) / Math.Log(L);
            double richardson = (Sum3 - Sum2) / (Math.Pow(L, m) - 1);
            int l = n * L * L;

            while (Math.Abs(richardson) > eps)
            {
                Sum1 = Sum2;
                Sum2 = Sum3;
                l *= L;
                Sum3 = MethodCQF(l,isNewton);

                m = -Math.Log((Sum3 - Sum2) / (Sum2 - Sum1)) / Math.Log(L);
                richardson = (Sum3 - Sum2) / (Math.Pow(L, m) - 1);
            }

            return Sum3 + richardson;
        }

        public static double MethodCQFOptimal(double eps,bool isNewton = true)
        {
            int L = 2;
            double Sum1 = MethodCQF(1,isNewton);
            double Sum2 = MethodCQF(2,isNewton);
            double Sum3 = MethodCQF(4,isNewton);
            double m = isNewton ? 3 : 6;
            double richardson = (Sum3 - Sum2) / (Math.Pow(L, m) - 1);

            while (Math.Abs(richardson) > eps)
            {
                double hopt = Math.Pow((B - A) / L * (eps * (1 - Math.Pow(L, -m))) / (Math.Abs(Sum2 - Sum1)), 1.0 / m);
                int n = (int)Math.Ceiling((B - A) / (hopt));

                Sum1 = Sum2;
                Sum2 = Sum3;
                Sum3 = MethodCQF(n,isNewton); 

                richardson = (Sum3 - Sum2) / (Math.Pow(L, m) - 1);
            }

            return Sum3 + richardson;
        }
    }
}
