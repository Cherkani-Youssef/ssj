package wafomExperiments;

import umontreal.ssj.discrepancy.LookUpTableWafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;

/**
 * This algorithm, described by Shin Harase in "Monte Carlo Methods and
 * Applications 22.4 (2016), pp. 349–357," efficiently searches for low WAFOM
 * point sets. It extends a point set \(P_d = \{\mathbf{x}_0, \dots,
 * \mathbf{x}_{2^d-1}\}\) for \(0 \leq d \leq k\), where each \(P_d\) is a
 * subset of \(P\) with \(|P|=2^k\). The algorithm evaluates numerous point sets
 * for each \(d\), retaining the best in terms of WAFOM.
 * 
 * This class constructs the columns of the generating matrix \(C\) column by
 * column randomly, ensuring the upper left (\(k \times k\)) submatrices are
 * invertible.
 */

public class ExplicitExtendSize {

	private int w;
	private int k;
	private int nBTries;
	private int dim;
	RandomStream stream;
	private LookUpTableWafom table_c; // for the fast Wafom
	private int q;
	double h;
	double factor;
	String SearchInfos = " ";
	/**
	 * @param k      number of columns
	 * @param dim    dimension of the net
	 * @param w      output digits (precision)
	 * @param q      a positive integer
	 * @param l      length of the segment
	 * @param h      for Matsumoto's definition set to 0, and for Yoshiki's set to 1
	 * @param factor set to 1 for wafom, set to 2 to get Wafom for RMSE
	 */
	public ExplicitExtendSize(int w, int k, int dim, int nbTries, LookUpTableWafom table_c, int q, double h,
			double factor, RandomStream stream) {

		this.w = w;
		this.k = k;
		this.dim = dim;
		this.nBTries = nbTries;
		this.stream = stream;
		this.table_c = table_c;
		this.q = q;
		this.h = h;
		this.factor = factor;

	}

	public int[] randomSearchColumnBYColumn(int verbose) {

		GeneratingMatrix[] bestMatrices = null;
		GeneratingMatrix[] currentMatrices = null;

		double lowestWafom = Double.MAX_VALUE;
		int numPoints = (1 << k);
		int end = 0;

		double lastWafVal = 0.0;
		double wafVal = 0.0;
		for (int colIndex = 1; colIndex <= k; colIndex++) {
			int start = (end);
			end = (1 << colIndex);
//System.out.println(" start " + start + " end " + end);
			 wafVal = (colIndex == 1) ? 0.0 : (lastWafVal);
			lowestWafom = Double.MAX_VALUE;
			numPoints =end;
			for (int attempt = 0; attempt <= nBTries; attempt++) {

				GeneratingMatrix[] tempMatrices = GeneratingMatrix.copyGeneratingMatrixArray(bestMatrices);

				if (colIndex == 1) {
					// Generate a new set of matrices
					tempMatrices = GeneratingMatrix.generateMatrixArray(dim, w, 1, stream);

				} else {
					// add the columns
					GeneratingMatrix.addNewColumnToEachMatrix(tempMatrices, stream);
				}

				// Convert the generated matrices into a format usable by DigitalNetBase2
				int[] C = GeneratingMatrix.generatorMatricesFromStandardFormat(tempMatrices);
				DigitalNetBase2 digitalNet = new DigitalNetBase2(colIndex, w, dim, C);
//				numPoints = digitalNet.getNumPoints();
				WafomFast wafom = new WafomFast(digitalNet, w, q, table_c);

				double somme = wafom.computeWafomSum(start, end) + wafVal;
				double currentWafom =  (somme) / numPoints;
				
				// Update the best matrices if the current Wafom is lower
				if (currentWafom < lowestWafom) {
					lowestWafom = currentWafom;
					currentMatrices = GeneratingMatrix.copyGeneratingMatrixArray(tempMatrices);
					lastWafVal = somme;

				}
			}
			 wafVal =  (lastWafVal);
			if ((verbose > 0)) {
				System.out.println("-----> Column  " + colIndex + " lowest Wafom so far " + lowestWafom + " log10  "
						+ Math.log10(lowestWafom));
			}
			bestMatrices = GeneratingMatrix.copyGeneratingMatrixArray(currentMatrices);

		}

		int[] genMat = GeneratingMatrix.generatorMatricesFromStandardFormat(bestMatrices);
		return genMat; // Return the matrices with the lowest Wafom value
	}

	public static void main(String[] args) {

		int w = 31;
		int q = 3;
		int l = w / q;
		int k = 15;
		int dim = 5;
		int nBTries = 10;

		double h = 1.0;
		double factor = 2.0;
		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);

		int verbose = 1;

		RandomStream stream = new MRG32k3a();
		ExplicitExtendSize searcher = new ExplicitExtendSize(w, k, dim, nBTries, table_c, q, h, factor, stream);

		long startTime = System.currentTimeMillis();
		int[] genMat = searcher.randomSearchColumnBYColumn(verbose);
		long endTime = System.currentTimeMillis();

		DigitalNetBase2 digitalNet = new DigitalNetBase2(k, w, dim, genMat);

		WafomFast wafom = new WafomFast(digitalNet, w, q, table_c);
		double wafomValue = wafom.computeWafom();

		String wafomType = (factor == 1) ? " Matsumoto " : " Goda";
		System.out.println();
		System.out.print(wafomType + " Wafom : " + wafomValue + " | ");
		double factor2 = (factor == 1) ? 2 : 1;

		LookUpTableWafom table_c2 = new LookUpTableWafom(w, q, w / q, h, factor2);
		wafom = new WafomFast(digitalNet, w, q, table_c2);
		wafomValue = wafom.computeWafom();
		wafomType = (factor2 == 1) ? " Matsumoto " : " Goda";
		System.out.println(
				wafomType + " Wafom : " + wafomValue + " Duration " + (endTime - startTime) / 1000.0 + " seconds");

	}

}
