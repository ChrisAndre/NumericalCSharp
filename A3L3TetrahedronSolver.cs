using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GaussNewton
{
    class A3L3TetrahedronSolver
    {
        public MVNRA approx;
        public static double scalingFactor = 1;
        //strip negatives on i,j,k and a,b,c; may be causing divergence because of borked derivatives
        public A3L3TetrahedronSolver()
        {
            var func = new MVNRA.MultivariableFunction(leastSquares, dLeastSquares);
            approx = new MVNRA(func);
            approx.setMaxIterations(5000);
            approx.setTargetTolerance(1e-10);
            approx.setDeltaAttenuation(1.5);
            //approx.setUseAllIterations(true);
        }
        public int solve(double a, double b, double c, double theta_a, double theta_b, double theta_c, double[] initialGuess)
        {
            double[] constants = new double[6] { a, b, c, theta_a, theta_b, theta_c };
            int iters = 0;
            approx.run(initialGuess, constants, ref iters);
            return iters;
        }
        //Switch to cosines rather than reevaluating
        public static double leastSquares(double[] sideLengths, double[] constants)
        {
            double i = Math.Abs(sideLengths[0]);
            double j = Math.Abs(sideLengths[1]);
            double k = Math.Abs(sideLengths[2]);
            double a = constants[0];
            double b = constants[1];
            double c = constants[2];
            double theta_a = constants[3];
            double theta_b = constants[4];
            double theta_c = constants[5];
            double _return = Math.Pow(lawOfCosines(a, i, j, theta_a), 2) +
                   Math.Pow(lawOfCosines(b, j, k, theta_b), 2) +
                   Math.Pow(lawOfCosines(c, k, i, theta_c), 2);
            return scalingFactor * _return;
            // zero of leastSquares is where all functions have zero values;
        }
        public static double[] dLeastSquares(double[] sideLengths, double[] constants)
        {
            double i = Math.Abs(sideLengths[0]);
            double j = Math.Abs(sideLengths[1]);
            double k = Math.Abs(sideLengths[2]);
            double a = constants[0];
            double b = constants[1];
            double c = constants[2];
            double theta_a = constants[3];
            double theta_b = constants[4];
            double theta_c = constants[5];
            double fa = lawOfCosines(a, i, j, theta_a);
            double fb = lawOfCosines(b, j, k, theta_b);
            double fc = lawOfCosines(c, k, i, theta_c);
            double[] ret = new double[3];
            ret[0] = 4 * fa * (i - j * Math.Cos(theta_a)) + 4 * fc * (i - k * Math.Cos(theta_c));
            ret[1] = 4 * fa * (j - i * Math.Cos(theta_a)) + 4 * fb * (j - k * Math.Cos(theta_b));
            ret[2] = 4 * fb * (k - j * Math.Cos(theta_b)) + 4 * fc * (k - i * Math.Cos(theta_c));
            ret.scale(scalingFactor);
            return ret;
        }
        static double lawOfCosines(double c, double a, double b, double theta_c)
        {
            a = Math.Abs(a);
            b = Math.Abs(b);
            return a * a + b * b - c * c - 2 * a * b * Math.Cos(theta_c);
        }
    }
}
