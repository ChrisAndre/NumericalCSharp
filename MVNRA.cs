using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GaussNewton
{
    class MVNRA
    {
        public class MultivariableFunction
        {
            public MultivariableFunction(Function f, Derivative delta)
            {
                this.f = f;
                this.delta = delta;
            }
            public MultivariableFunction(Function f)
            {
                this.f = f;
            }
            public delegate double Function(double[] variables, double[] constants);
            public delegate double[] Derivative(double[] variables, double[] constants);
            Function f;
            Derivative delta = null;
            public double eval(double[] variables, double[] constants)
            {
                return f(variables, constants);
            }
            public double[] gradient(double[] variables, double[] constants, double spread = 0.01)
            {
                if (delta == null)
                {
                    double[] deriv = new double[variables.Length];
                    for (int i = 0; i < variables.Length; i++)
                    {
                        deriv[i] = singleVariableDerivative(variables, constants, i, spread);
                    }
                    return deriv;
                }
                else
                {
                    return delta(variables, constants);
                }
            }
            double singleVariableDerivative(double[] variables, double[] constants, int index, double spread = 0.01)
            {
                variables[index] += spread;
                double fp = f(variables, constants);
                variables[index] -= 2 * spread;
                double fm = f(variables, constants);
                return (fp - fm) / 2 / spread;
            }
            double singleVariableDerivative2(double[] variables, double[] constants, int index, double spread = 0.01)
            {
                double fc = f(variables, constants);
                variables[index] += spread;
                double fp = f(variables, constants);
                variables[index] -= 2 * spread;
                double fm = f(variables, constants);
                double d1 = (fp - fc) / spread;
                double d2 = (fc - fm) / spread;
                return (d1 - d2) / spread;
            }
        }
        public MultivariableFunction func;
        public MVNRA(MultivariableFunction func)
        {
            this.func = func;
        }
        double targetValue = 0;
        double targetTolerance = 1e-8;
        int maxIterations = 1;
        bool cappedIterations = false;
        bool useAllIterations = false;
        double deltaAttenuation = 1;
        double spread = 1e-6;
        public void setUseAllIterations(bool useAllIterations)
        {
            this.useAllIterations = useAllIterations;
        }
        public void setSpread(double spread)
        {
            this.spread = spread;
        }
        public void setMaxIterations(int iters)
        {
            cappedIterations = true;
            maxIterations = iters;
        }
        public void uncapIterations()
        {
            cappedIterations = false;
        }
        public void setTargetTolerance(double targetTolerance)
        {
            this.targetTolerance = targetTolerance;
        }
        public void setTargetValue(double targetValue)
        {
            this.targetValue = targetValue;
        }
        public void setDeltaAttenuation(double deltaAttenuation)
        {
            this.deltaAttenuation = deltaAttenuation;
        }
        public double[] run(double[] initialGuess, double[] constants, ref int iters)
        {
            double[] best = initialGuess;
            double bestError = distance(func.eval(initialGuess, constants));
            while (!isAtTarget(initialGuess, constants) && hasIterationsLeft(iters) || canContinueUsingAllIterations(iters))
            {
                getGuess(initialGuess, constants);
                double error = distance(func.eval(initialGuess, constants));
                if (error < bestError)
                {
                    bestError = error;
                    best = initialGuess.copy();
                }
                iters++;
            }
            //Console.WriteLine("besterror:" + bestError);
            return best;
        }
        void getGuess(double[] initialGuess, double[] constants)

        {
            var grad = func.gradient(initialGuess, constants, spread);
            var dy = grad.dot(grad) / grad.magnitude();
            var f = func.eval(initialGuess, constants) - targetValue;
            grad.normalize();
            grad.scale(-f * deltaAttenuation / dy);
            initialGuess.add(grad);
        }
        bool isAtTarget(double[] guess, double[] constants)
        {
            return Math.Abs(func.eval(guess, constants) - targetValue) <= targetTolerance;
        }
        bool hasIterationsLeft(int iteration)
        {
            return iteration < maxIterations || !cappedIterations;
        }
        bool canContinueUsingAllIterations(int iteration)
        {
            if (!useAllIterations) return false;
            return iteration < maxIterations;
        }
        double distance(double f)
        {
            return Math.Abs(f-targetValue);
        }
    }
}