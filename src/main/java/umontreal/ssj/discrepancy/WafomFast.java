package umontreal.ssj.discrepancy;

import java.io.IOException;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.PointSetIterator;
import wafomExperiments.KahanSummation;

/**
 * This class computes the Walsh Figure of Merit (WAFOM) for a Digital Net in
 * Base 2 using a lookup table to accelerate the computation.
 * 
 * For Matsumoto's definition
 * 
 * WAFOM(P) = -1 + \frac{1}{|P|} \left( \sum_{\mathbf{x}_i \in P}
 * \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 + \eta(x_{i,j,l}) 2^{-l} \right)
 * \right)
 * 
 * The same formula applies to Yoshiki's definition; simply replace \( l \) with
 * \( l + 1 \).
 * 
 * WAFOM(P) = -1 + \frac{1}{|P|} \left( \sum_{\mathbf{x}_i \in P}
 * \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 + \eta(x_{i,j,l}) 2^{-(l+1)} \right)
 * \right).
 * 
 * For Goda's definition:
 * 
 * \[ \mathcal{W}(P, \mu) = \sqrt{-1 + \frac{1}{|P|} \left( \sum_{\mathbf{x}_i
 * \in P} \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 + \eta(x_{i,j,l}) 2^{-2l}
 * \right) \right)} \] The same formula applies to Yoshiki's definition; simply
 * replace \( l \) with \( l + 1 \).
 * 
 * The same formula applies to Yoshiki's definition; simply replace \( l \) with
 * \( l + 1 \).
 * 
 * \[ \mathcal{W}(P, \mu + h) = \sqrt{-1 + \frac{1}{|P|} \left(
 * \sum_{\mathbf{x}_i \in P} \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 +
 * \eta(x_{i,j,l}) 2^{-2(l+1)} \right) \right)} \]
 *
 * The method, proposed by Shin Harase in "A search for extensible low-WAFOM
 * point sets," Monte Carlo Methods and Applications 22.4 (2016), pp. 349–357,
 * involves dividing a coordinate of a point with 'w' precision into 'q'
 * subparts, each containing 'l' elements if q is a divisor of 'w'.
 * When 'q' is not a divisor of the number of output digits then there will
 * another be a segment of lehgth t where \(t = w \mod q\) .
 * 
 * A coordinate of a point is expressed as \(X^i = d_1^i, \ldots, d_q^i,
 * d_{q+1}^i\), where for \(1 \leq c \leq q\), the length of each segment is \(l
 * = w \div q\), and the remainder \(d_{q+1}\) is of length \(t = w \mod q\).
 * Each segment \(d_c^i\) is represented as \(x_{i, (c-1)*l+1}, \ldots, x_{i,
 * c*l}\).
 * 
 * Instead of computing \(\prod_{j=1}^{w} [1 + (-1)^{(-1)^{x_{i,j}}} 2^{-j}] -
 * 1\), we use \(\prod_{1 \leq c \leq q} \text{table}_c[d_c^i]\), where
 * \(\text{table}_c[d_c^i]\) are precomputed values: \(\text{table}_c[d_c^i] =
 * \prod_{1 \leq j \leq l} (1 + (-1)^{x_{i, (c-1)*l + j}} \cdot 2^{-((c-1)*l + j
 * + 1)})\). This significantly reduces computation time.
 * 
 * We use Kahan summation algorithm improves numerical accuracy by compensating
 * for floating-point errors, reducing round-off error accumulation.
 **/
public class WafomFast extends WafomBase {

	private int q;
	private LookUpTableWafom tableC;
	int dimension;
	int numPoints;

	/**
	 * Initializes the WafomNUS object
	 * 
	 * @param digitalNet The Digital Net Base 2 instance.
	 * @param outDigits  The precision of the points.
	 * @param q          The positive integer that divides the output digits.
	 * @param tableC     The lookup table for Wafom calculation.
	 */
	public WafomFast(DigitalNetBase2 digitalNet, int outDigits, int q, LookUpTableWafom tableC) {
		super(digitalNet, outDigits);
		this.q = q;
		this.tableC = tableC;
		this.dimension = digitalNet.getDimension();
		this.numPoints = digitalNet.getNumPoints();

	}

	/**
	 * Compute the merit value of the whole digital net
	 */
	@Override
	public double computeWafom() {
		double prod = 1.0, mOne = -1.0;

		PointSetIterator iter = digitalNet.iteratorNoGray();
		KahanSummation kahanSum = new KahanSummation();

		while (iter.hasNextPoint()) {
			prod = 1.0;

			while (iter.hasNextCoordinate()) {

				int u_i_j = iter.nextInt(); // added this in pointSetiterator to return the the coodinate in int format

				for (int c = 1; c <= q + (outDigits % q != 0 ? 1 : 0); ++c) {
					int d_c = compute_d_c(u_i_j, c, outDigits, q); // Compute d_c for each segment
					prod *= tableC.get(c, d_c);
				}

			}
//			KahanSummation kahanSum2 = new KahanSummation();
//			kahanSum2.add(prod);
//			kahanSum2.add(mOne);

//			kahanSum.add(kahanSum2.getSum());
			kahanSum.add(prod + mOne);
			iter.resetToNextPoint();
		}

		return ((kahanSum.getSum() / (double) numPoints));

	}

	/**
	 * Compute the merit value of a subset of points for all dimension of th digital
	 * net
	 */
	public double computeWafomSum(int start, int end) {
		double prod = 1.0, mOne = -1.0;

		// int[] currentPoint;
		KahanSummation kahanSum = new KahanSummation();
		PointSetIterator iter = digitalNet.iteratorNoGray();
		iter.setCurPointIndex(start);
		while (iter.getCurPointIndex() < end) {
			prod = 1.0;
			while (iter.hasNextCoordinate()) {

				int u_i_j = iter.nextInt(); // added this in pointSetiterator to return the the coodinate in int format

				for (int c = 1; c <= q + (outDigits % q != 0 ? 1 : 0); ++c) {
					int d_c = compute_d_c(u_i_j, c, outDigits, q); // Compute d_c for each segment
					prod *= tableC.get(c, d_c);
				}

			}

//
			kahanSum.add(prod + mOne);
			iter.resetToNextPoint();
		}

		return kahanSum.getSum();
	}

	/**
	 * Compute the merit value of a subset of points for a specific dimension j
	 */
	public double[] computeWafom2(int start, int end, int j) {
		double prod = 1.0;

		int numPoints = end - start;
		double[] wafomValues = new double[numPoints];
		PointSetIterator iter = digitalNet.iteratorNoGray();
		int ii = 0;
		for (int i = start; i < end; i++) {
			prod = 1.0;
			iter.setCurPointIndex(i);
			iter.setCurCoordIndex(j);

			int u_i_j = iter.nextInt();

			for (int c = 1; c <= q + (outDigits % q != 0 ? 1 : 0); ++c) {
				int d_c = compute_d_c(u_i_j, c, outDigits, q);
				prod *= tableC.get(c, d_c);
			}

//			wafomValues[i - start] = prod;
			wafomValues[i - start] = prod;

			iter.resetToNextPoint();
		}
		return wafomValues;
	}

	/**
	 * Compute the merit value for the pointset for a specific dimension j
	 */

	public double[] computeWafomForCBC(int j) {

		double prod = 1.0;
		PointSetIterator iter = digitalNet.iteratorNoGray();

		double[] wafomValues = new double[numPoints];

		for (int i = 0; i < numPoints; i++) {
			prod = 1.0;
			iter.setCurPointIndex(i);
			iter.setCurCoordIndex(j);
			int u_i_j = iter.nextInt();
			for (int c = 1; c <= q + (outDigits % q != 0 ? 1 : 0); ++c) {
				int d_c = compute_d_c(u_i_j, c, outDigits, q); // Compute d_c for each segment
				prod *= tableC.get(c, d_c);
			}
			wafomValues[i] = prod;
		}
		return wafomValues;
	}

	/// helper functions

	public double elementWiseProd(double[] wafomValues1, double[] wafomValues) {
		double mOne = -1.0;
		KahanSummation kahanSum = new KahanSummation();
		for (int s = 0; s < wafomValues1.length; s++) {
			wafomValues1[s] *= wafomValues[s];
//			wafomValues1[s]+=mOne;
			kahanSum.add(wafomValues1[s] + mOne);

		}
		return kahanSum.getSum();
	}

	public double elementWiseProd2(double[] wafomValues1, double[] wafomValues) {
		double mOne = -1.0;
		KahanSummation kahanSum = new KahanSummation();
		for (int s = 0; s < wafomValues1.length; s++) {
			wafomValues1[s] = Math.exp(Math.log(wafomValues1[s]) + Math.log(wafomValues[s]));

			kahanSum.add(wafomValues1[s] + mOne);
		}
		return kahanSum.getSum();
	}

//	public double sum2(double[] wafomValues) {
//		double somme = 0;
//
//		for (double val : wafomValues)
//			somme += val;
//		return somme;
//	}
//
	public double sumKahan(double[] array) {

		KahanSummation kahanSum = new KahanSummation();
		for (double val : array) {
			kahanSum.add(val);
		}
		return kahanSum.getSum();
	}

}
