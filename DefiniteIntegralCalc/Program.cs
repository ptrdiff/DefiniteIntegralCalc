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

            Console.WriteLine(CalcDefInteg.NewtonCotesIQF());
            Console.WriteLine(CalcDefInteg.MethodCQF(3));
            Console.WriteLine(CalcDefInteg.MethodCQFHalf(1, 1e-6));
            Console.WriteLine(CalcDefInteg.MethodCQFOptimal(1e-6));

            Console.WriteLine();

            Console.WriteLine(CalcDefInteg.GaussIQF());
            Console.WriteLine(CalcDefInteg.MethodCQF(3,false));
            Console.WriteLine(CalcDefInteg.MethodCQFHalf(1, 1e-6,false));
            Console.WriteLine(CalcDefInteg.MethodCQFOptimal(1e-6,false));
        }
    }
}
