package wafomExperiments;

/**
 * The Kahan Summation Algorithm improves the accuracy of floating-point
 * summation by reducing numerical errors introduced during addition. It
 * maintains a running total and a small error term to correct the drift in
 * precision as new elements are added to the sum.
 */

public class KahanSummation {
	private double sum; // The running total
	private double c; // A compensation for lost low-order bits

	public KahanSummation() {
		this.sum = 0.0;
		this.c = 0.0;
	}

	public void add(double number) {
		double y = number - this.c; // So far, so good: c is zero.
		double t = this.sum + y; // Alas, sum is big, y small, so low-order digits of y are lost.
		this.c = (t - this.sum) - y; // (t - sum) recovers the high-order part of y; subtracting y recovers -(low
										// part of y)
		this.sum = t; // Algebraically, c should always be zero. Beware overly-aggressive optimizing
						// compilers!
	}

	public double getSum() {
		return this.sum;
	}
}
