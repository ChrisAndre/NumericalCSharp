﻿using System;

namespace NumericalCSharp
{
    class A3L3TetrahedronSolver
    {
        private NewtonRaphsonND approx;
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="maxIterations">Iterations the Newton-Raphson solver should use at maximum.</param>
        /// <param name="targetTolerance">Used in least squares; lower generally more accurate. Use caution.</param>
        /// <param name="deltaScaling">Scales Newton-Raphson jumps for faster long-distance convergence. Do not go above 2.</param>
        public A3L3TetrahedronSolver(int maxIterations = 5000, double targetTolerance = 1e-10, double deltaScaling = 1.5)
        {
            var func = new NewtonRaphsonND.MultivariableFunction(leastSquares, dLeastSquares);
            approx = new NewtonRaphsonND(func);
            approx.setMaxIterations(maxIterations);
            approx.setTargetTolerance(targetTolerance);
            approx.setDeltaAttenuation(deltaScaling);
        }
        /// <summary>
        /// Evolves side lengths in initialGuess using newton-raphson and non-linear least squares.
        /// </summary>
        /// <param name="a">Side a of base triangle.</param>
        /// <param name="b">Side b of base triangle.</param>
        /// <param name="c">Side c of base triangle.</param>
        /// <param name="theta_a">Angle between apex, a, and b.</param>
        /// <param name="theta_b">Angle between apex, b, and c.</param>
        /// <param name="theta_c">Angle between apex, c, and a.</param>
        /// <param name="initialGuess">Contains sidelengths to be evolved.</param>
        /// <returns>Iterations used.</returns>
        /// <remarks>Algorithm sometimes fails to converge when initialGuess lengths start out lower than actual lengths; have them longer than anticipated.</remarks>
        public int solve(double a, double b, double c, double theta_a, double theta_b, double theta_c, double[] initialGuess)
        {
            double[] constants = new double[6] { a, b, c, theta_a, theta_b, theta_c };
            int iters = 0;
            approx.run(initialGuess, constants, ref iters);
            return iters;
        }
        private static double leastSquares(double[] sideLengths, double[] constants)
        {
            //TODO: switch to cosines rather than reevaluating
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
            return _return;
            // zero of leastSquares is where all functions have zero values;
        }
        private static double[] dLeastSquares(double[] sideLengths, double[] constants)
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
            return ret;
        }
        private static double lawOfCosines(double c, double a, double b, double theta_c)
        {
            a = Math.Abs(a);
            b = Math.Abs(b);
            return a * a + b * b - c * c - 2 * a * b * Math.Cos(theta_c);
        }
    }
}
