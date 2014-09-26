using System;

namespace NumericalCSharp
{
    /// <summary>
    /// Base class for single-dimensional iterative root finders
    /// </summary>
    public abstract class RootFinder1D
    {
        public class Solution
        {
            public double Value;
            public bool Satisfies;
            public double RangeDerivative; //TODO
            public int Iterations;
        }
        public delegate double Function(double u);
        public double TargetValue = 0;
        public double TargetTolerance = 1e-8;
        protected bool cappedIterations = false;
        public int MaxIterations
        {
            get
            {
                return _MaxIterations;
            }
            set
            {
                cappedIterations = true;
                _MaxIterations = value;
            }
        }
        int _MaxIterations;
        public bool UseAllIterations = false;
        public void uncapIterations()
        {
            cappedIterations = false;
            UseAllIterations = false;
        }
        protected bool isZero(double x)
        {
            return Math.Abs(x) <= TargetTolerance;
        }
        protected bool isAtTarget(double f)
        {
            return isZero(f - TargetValue);
        }
        protected bool canContinueCappedIterations(int iters)
        {
            return iters < MaxIterations || !cappedIterations;
        }
        protected bool canContinueUsingAllIterations(int iters)
        {
            return iters < MaxIterations && UseAllIterations && cappedIterations;
        }
        public abstract Solution getSolution(double initialGuess);
    }
    /// <summary>
    /// Linear solver.
    /// </summary>
    /// <remarks>
    /// Multidimensional versions exist.
    /// </remarks>
    public class NewtonRaphson1D : RootFinder1D
    {
        public NewtonRaphson1D(Function f, Function d)
            : this(f)
        {
            this.d = d;
        }
        public NewtonRaphson1D(Function f)
        {
            this.f = f;
        }
        Function f, d;
        double derivativeSpread = 1e-6;
        public override Solution getSolution(double initialGuess)
        {
            int iters = 0;
            while (!isAtTarget(f(initialGuess)) && canContinueCappedIterations(iters) || canContinueUsingAllIterations(iters))
            {
                initialGuess += getDeltaX(initialGuess);
                iters++;
            }
            return new Solution()
            {
                Value = initialGuess,
                Iterations = iters,
                Satisfies = isAtTarget(f(initialGuess)),
                RangeDerivative = getDerivative(initialGuess),
            };
        }
        double getDeltaX(double guess)
        {
            double fg = f(guess) - TargetValue;
            return -(fg) / getDerivative(guess);
        }
        double getDerivative(double x)
        {
            if (d != null)
                return d(x);
            else
                return approxDerivative(x);
        }
        double approxDerivative(double x)
        {
            double fx_plus = f(x + derivativeSpread);
            double fx_minus = f(x - derivativeSpread);
            return (fx_plus - fx_minus) / 2 / derivativeSpread;
        }
        double approxSecondDerivative(double x)
        {
            double fx_plus = f(x + derivativeSpread);
            double fx_minus = f(x - derivativeSpread);
            double fx_at = f(x);
            double fpx_plus = (fx_plus - fx_at) / derivativeSpread;
            double fpx_minus = (fx_at - fx_minus) / derivativeSpread;
            return (fpx_plus - fpx_minus) / derivativeSpread;
        }
    }
    /// <summary>
    /// Linear over linear solver of Pade approximation.
    /// </summary>
    /// <remarks>
    /// 1st derivative at zero of f should be non-zero for cubic convergence.
    /// Multidimensional versions exist. (TODO)
    /// </remarks>
    public class Halley1D : RootFinder1D
    {
        public Halley1D(Function f)
        {
            this.f = f;
        }
        public Halley1D(Function f, Function d)
            : this(f)
        {
            this.d = d;
        }
        public Halley1D(Function f, Function d, Function dd)
            : this(f, d)
        {
            this.dd = dd;
        }
        Function f, d, dd;
        double derivativeSpread = 1e-6;
        double secondDerivativeSpread = 1e-3;
        public override Solution getSolution(double initialGuess)
        {
            int iters = 0;
            while (!isAtTarget(f(initialGuess)) && canContinueCappedIterations(iters) || canContinueUsingAllIterations(iters))
            {
                initialGuess += getDeltaX(initialGuess);
                iters++;
            }
            return new Solution()
            {
                Value = initialGuess,
                Satisfies = isAtTarget(f(initialGuess)),
                RangeDerivative = getDerivative(initialGuess),
                Iterations = iters
            };
        }
        double getDeltaX(double guess)
        {
            double fg = f(guess) - TargetValue;
            double fd = getDerivative(guess);
            double fdd = getSecondDerivative(guess);
            return -2 * fg * fd / (2 * fd * fd - fg * fdd);
        }
        double getDerivative(double x)
        {
            if (d != null)
                return d(x);
            return approxDerivative(x);
        }
        double getSecondDerivative(double x)
        {
            if (dd != null)
                return dd(x);
            return approxSecondDerivative(x);
        }
        double approxDerivative(double x)
        {
            double fx_plus = f(x + derivativeSpread);
            double fx_minus = f(x - derivativeSpread);
            return (fx_plus - fx_minus) / 2 / derivativeSpread;
        }
        double approxSecondDerivative(double x)
        {
            double fpx_plus = getDerivative(x + secondDerivativeSpread);
            double fpx_minus = getDerivative(x - secondDerivativeSpread);
            return (fpx_plus - fpx_minus) / secondDerivativeSpread / 2;
        }
    }
}
