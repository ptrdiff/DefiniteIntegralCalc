using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DefiniteIntegralCalc
{
   
    class Program
    {
        static void Main(string[] args)
        {
            //Console.WriteLine(CalcDefInteg.RiemannSum(100));
            Console.WriteLine(CalcDefInteg.RealResult());

            Console.WriteLine();

            Console.WriteLine("NewtonIQF: {0}", CalcDefInteg.NewtonCotesIQF());
            Console.WriteLine("NewtonCQF: {0}", CalcDefInteg.MethodCQF(3));
            Console.WriteLine("NewtonCQFHalf: {0}", CalcDefInteg.MethodCQFHalf(1, 2, 1e-6));
            Console.WriteLine("NewtonCQFOptimal: {0}",CalcDefInteg.MethodCQFOptimal(1e-6));

            Console.WriteLine();

            Console.WriteLine("GaussIQF: {0}", CalcDefInteg.GaussIQF());
            Console.WriteLine("GaussCQF: {0}", CalcDefInteg.MethodCQF(3,false));
            Console.WriteLine("GaussCQFHalf: {0}", CalcDefInteg.MethodCQFHalf(1, 2, 1e-6,false));
            Console.WriteLine("GaussCQFOptimal: {0}", CalcDefInteg.MethodCQFOptimal(1e-6,false));
        }
    }
}
