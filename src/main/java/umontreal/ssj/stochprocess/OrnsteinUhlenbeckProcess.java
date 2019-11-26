/*
 * Class:        OrnsteinUhlenbeckProcess
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.stochprocess;

import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.*;

/**
 * This class represents an *Ornstein-Uhlenbeck* process @f$\{X(t) : t \geq0 \}@f$, sampled at
 * times @f$0 = t_0 < t_1 < \cdots< t_d@f$. This process obeys the stochastic differential equation
 * 
 * @anchor REF_stochprocess_OrnsteinUhlenbeckProcess_eq_ornstein @f[ dX(t) = \alpha(b - X(t)) dt +
 *         \sigma dB(t) \tag{ornstein} @f] with initial condition @f$X(0)= x_0@f$,
 *         where @f$\alpha@f$, @f$b@f$ and
 * @f$\sigma@f$ are positive constants, and @f$\{B(t), t\ge0\}@f$ is a standard Brownian motion
 *              (with drift 0 and volatility 1). This process is *mean-reverting* in the sense that
 *              it always tends to drift toward its general mean @f$b@f$. The process is generated
 *              using the sequential technique @cite fGLA04a&thinsp; (p. 110)
 * @anchor REF_stochprocess_OrnsteinUhlenbeckProcess_eq_ornstein_seq @f[ X(t_j) = e^{-\alpha(t_j -
 *         t_{j-1})} X(t_{j-1}) + b\left(1 - e^{-\alpha(t_j - t_{j-1})}\right) + \sigma\sqrt{\frac{1
 *         - e^{-2\alpha(t_j - t_{j-1})}}{2\alpha}} Z_j \tag{ornstein-seq} @f] where @f$Z_j \sim
 *         N(0,1)@f$. The time intervals @f$t_j - t_{j-1}@f$ can be arbitrarily large.
 *
 *         <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class OrnsteinUhlenbeckProcess extends StochasticProcess {
	protected NormalGen gen;
	protected double alpha, beta, sigma;
	// Precomputed values
	protected double[] badt, alphadt, sigmasqrdt;

	/**
	 * Constructs a new `OrnsteinUhlenbeckProcess` with parameters
	 * 
	 * @f$\alpha=@f$ `alpha`, @f$b@f$, @f$\sigma=@f$ `sigma` and initial value @f$X(t_0) =@f$ `x0`.
	 *               The normal variates @f$Z_j@f$ will be generated by inversion using the stream
	 *               `stream`.
	 */
	public OrnsteinUhlenbeckProcess(double x0, double alpha, double b, double sigma,
	        RandomStream stream) {
		this(x0, alpha, b, sigma, new NormalGen(stream));
	}

	/**
	 * Here, the normal variate generator is specified directly instead of specifying the stream.
	 * The normal generator `gen` can use another method than inversion.
	 */
	public OrnsteinUhlenbeckProcess(double x0, double alpha, double b, double sigma,
	        NormalGen gen) {
		this.alpha = alpha;
		this.beta = b;
		this.sigma = sigma;
		this.x0 = x0;
		this.gen = gen;
	}

	public double nextObservation() {
		double xOld = path[observationIndex];
		double x = badt[observationIndex] + xOld * alphadt[observationIndex]
		        + sigmasqrdt[observationIndex] * gen.nextDouble();
		observationIndex++;
		path[observationIndex] = x;
		return x;
	}

	/**
	 * Generates and returns the next observation at time @f$t_{j+1} =@f$ `nextTime`, using the
	 * previous observation time @f$t_j@f$ defined earlier (either by this method or by
	 * <tt>setObservationTimes</tt>), as well as the value of the previous observation @f$X(t_j)@f$.
	 * *Warning*: This method will reset the observations time @f$t_{j+1}@f$ for this process to
	 * `nextTime`. The user must make sure that the @f$t_{j+1}@f$ supplied is
	 * 
	 * @f$\geq t_j@f$.
	 */
	public double nextObservation(double nextTime) {
		double previousTime = t[observationIndex];
		double xOld = path[observationIndex];
		observationIndex++;
		t[observationIndex] = nextTime;
		double dt = nextTime - previousTime;
		double tem = Math.exp(-alpha * dt);
		double tem1 = -Math.expm1(-alpha * dt);
		double x = tem * xOld + beta * tem1
		        + sigma * Math.sqrt(tem1 * (1.0 + tem) / (2.0 * alpha)) * gen.nextDouble();
		path[observationIndex] = x;
		return x;
	}

	/**
	 * Generates an observation of the process in `dt` time units, assuming that the process has
	 * value @f$x@f$ at the current time. Uses the process parameters specified in the constructor.
	 * Note that this method does not affect the sample path of the process stored internally (if
	 * any).
	 */
	public double nextObservation(double x, double dt) {
		double tem = Math.exp(-alpha * dt);
		double tem1 = -Math.expm1(-alpha * dt);
		x = tem * x + beta * tem1
		        + sigma * Math.sqrt(tem1 * (1.0 + tem) / (2.0 * alpha)) * gen.nextDouble();
		return x;
	}

	public double[] generatePath() {
		double x;
		double xOld = x0;
		for (int j = 0; j < d; j++) {
			x = badt[j] + xOld * alphadt[j] + sigmasqrdt[j] * gen.nextDouble();
			path[j + 1] = x;
			xOld = x;
		}
		observationIndex = d;
		return path;
	}

	/**
	 * Generates a sample path of the process at all observation times, which are provided in array
	 * `t`. Note that `t[0]` should be the observation time of `x0`, the initial value of the
	 * process, and `t[]` should have at least
	 * 
	 * @f$d+1@f$ elements (see the `setObservationTimes` method).
	 */
	public double[] generatePath(RandomStream stream) {
		gen.setStream(stream);
		return generatePath();
	}

	/**
	 * Resets the parameters @f$X(t_0) =@f$ `x0`, @f$\alpha=@f$ `alpha`,
	 * 
	 * @f$b =@f$ `b` and @f$\sigma=@f$ `sigma` of the process. *Warning*: This method will recompute
	 *      some quantities stored internally, which may be slow if called too frequently.
	 */
	public void setParams(double x0, double alpha, double b, double sigma) {
		this.alpha = alpha;
		this.beta = b;
		this.sigma = sigma;
		this.x0 = x0;
		if (observationTimesSet)
			init(); // Otherwise not needed.
	}

	/**
	 * Resets the random stream of the normal generator to `stream`.
	 */
	public void setStream(RandomStream stream) {
		gen.setStream(stream);
	}

	/**
	 * Returns the random stream of the normal generator.
	 */
	public RandomStream getStream() {
		return gen.getStream();
	}

	/**
	 * Returns the value of @f$\alpha@f$.
	 */
	public double getAlpha() {
		return alpha;
	}

	/**
	 * Returns the value of @f$b@f$.
	 */
	public double getB() {
		return beta;
	}

	/**
	 * Returns the value of @f$\sigma@f$.
	 */
	public double getSigma() {
		return sigma;
	}

	/**
	 * Returns the normal random variate generator used. The `RandomStream` used for that generator
	 * can be changed via `getGen().setStream(stream)`, for example.
	 */
	public NormalGen getGen() {
		return gen;
	}

	protected void initArrays(int d) {
		double dt, tem, tem1;
		for (int j = 0; j < d; j++) {
			dt = t[j + 1] - t[j];
			tem = Math.exp(-alpha * dt);
			tem1 = -Math.expm1(-alpha * dt);
			badt[j] = beta * tem1;
			alphadt[j] = tem;
			sigmasqrdt[j] = sigma * Math.sqrt(tem1 * (1.0 + tem) / (2.0 * alpha));
		}
	}

	// This is called by setObservationTimes to precompute constants
	// in order to speed up the path generation.
	protected void init() {
		super.init();
		badt = new double[d];
		alphadt = new double[d];
		sigmasqrdt = new double[d];
		initArrays(d);
	}

}