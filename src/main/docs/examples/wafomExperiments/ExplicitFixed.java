package wafomExperiments;

import umontreal.ssj.discrepancy.LookUpTableWafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;

/**
 * In this algorithm we generate the generating matrices of the digital net at
 * random and we retain the matrices with the lowest wafom.
 */

public class ExplicitFixed {

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

	public ExplicitFixed(int w, int k, int dim, int nbTries, LookUpTableWafom table_c, int q, double h, double factor,
			RandomStream stream) {

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

	public int[] randomSearchExplicit(int verbose) {

		double lowestWafom = Double.MAX_VALUE;

		int count = 0;
		int[] bestMatrices = new int[dim * k];

		for (int attempt = 0; attempt < nBTries; attempt++) {

			// Generate a new set of matrices
			GeneratingMatrix[] currentMatrices = GeneratingMatrix.generateMatrixArray(dim, w, k, stream);
//			stream.resetNextSubstream();
			int[] fromStandardFormat = GeneratingMatrix.generatorMatricesFromStandardFormat(currentMatrices);
			DigitalNetBase2 digitalNet = new DigitalNetBase2(k, w, dim, fromStandardFormat);

			WafomFast wafom = new WafomFast(digitalNet, w, q, table_c);
			double currentWafom = wafom.computeWafom(); 

			// If the current Wafom is lower, update best matrices and delete old file
			if (currentWafom < lowestWafom) {
				lowestWafom = currentWafom;
				bestMatrices = digitalNet.getGeneratorMatricesTrans();
			}

			if ((verbose > 0) && ((nBTries > 100 && attempt % 100 == 0) || (attempt % 10 == 0))) {
				System.out.println("Net " + attempt + "/" + nBTries + " currentWafom " + currentWafom
						+ " lowest Wafom  " + Math.log10(Math.abs(lowestWafom)));
			}

		}

		return bestMatrices; // Return the matrices with the lowest Wafom value
	}

	public static void main(String[] args) {

		int w = 31; // when using LMS set to 30
		int q = 3;
		int l = w / q;
		int k = 15;
		int dim = 5;
		int nBTries = 10_000;

		double h = 1.0;
		double factor = 1.0;
		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);

		int verbose = 1;

		RandomStream stream = new MRG32k3a();

		System.out.println("  : w " + w + " k " + k + " dim " + dim + " nBTries " + nBTries + " h " + h + " factor "
				+ factor + " nBTries " + nBTries);

		ExplicitFixed searcher = new ExplicitFixed(w, k, dim, nBTries, table_c, q, h, factor, stream);
		int[] genMat = searcher.randomSearchExplicit(verbose);

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
		System.out.println(wafomType + " Wafom : " + wafomValue);

	}

}
