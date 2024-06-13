package wafomExperiments;

import java.io.IOException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.LookUpTableWafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.LFSR258;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import wafomExperiments.HelperFunctions.SearchStrategy;

/**
 * This algorithm attempts to find extensible scrambled nets, both the LMS and
 * the scrambled generating matrices are constructed gradually. The goal is to
 * find extensible digital nets with low-Wafom for each subset of points
 */
public class LMSExtendSize {
	static String filePath = "/Users/cher/Downloads/new-joe-kuo-6.21201";

	private int[][] bestL;// stores the best LMS of the search
	private int w;
	private int k;
	private int nBTries;
	private int dim;
	private RandomStream stream;
	private LookUpTableWafom table_c;
	private int q;
	double h;
	double factor;
	String SearchInfos = " ";
	int nbTrirage;

	/**
	 * @param k      number of columns
	 * @param dim    dimension of the net
	 * @param w      output digits (precision)
	 * @param q      a positive integer
	 * @param l      length of the segment
	 * @param h      for Matsumoto's definition set to 0, and for Yoshiki's set to 1
	 * @param factor set to 1 for wafom, set to 2 to get Wafom for RMSE
	 */
	public LMSExtendSize(int w, int k, int dim, int nbTries, LookUpTableWafom table_c, int q, double h, double factor,
			RandomStream stream) {

		this.w = w;
		this.k = k;
		this.dim = dim;
		this.nBTries = nbTries;
		this.nbTrirage = nbTries;
		this.stream = stream;
		this.table_c = table_c;
		this.q = q;
		this.h = h;
		this.factor = factor;
		this.bestL = new int[dim][k];

	}

	public int[] ExtendSize(int verbose, int cc) throws IOException {

		final int twoRm1 = (1 << (w - 1));

		int end = 0;
		double lowestWafom = Double.MAX_VALUE;
		double lastWafVal = 0.0;

		for (int d = 0; d < k; d++) {

			int start = (d == 0) ? 0 : (end);
			end = (1 << (d + 1));

			double wafVal = (d == 0) ? 0.0 : (lastWafVal);

			lowestWafom = Double.MAX_VALUE;

			int[][] L = new int[dim][d + 1];
			for (int i = 0; i < dim; i++) {
				// Copy 'd' elements from each row of 'bestL' to 'L'
				System.arraycopy(bestL[i], 0, L[i], 0, d);
			}
			DigitalNetBase2 Sob = new SobolSequence(filePath, d + 1, w, dim);
			int numPoints = Sob.getNumPoints();
			for (int attempt = 0; attempt < nbTrirage; attempt++) {

				Sob.resetGeneratorMatrices();

				for (int j = 0; j < dim; j++) {
					L[j][d] = (twoRm1 + stream.nextInt(0, twoRm1 - 1)) >> (d);

					Sob.leftMultiplyMatBis2(j, L[j]);
				}

				WafomFast wafom = new WafomFast(Sob, w, q, table_c);

				double somme = wafom.computeWafomSum(start, end) + wafVal;// / end
				double currentWafom = (somme) / numPoints;

				if ((currentWafom) < (lowestWafom)) {
					lowestWafom = currentWafom;
					lastWafVal = somme;

					for (int i = 0; i < dim; i++)
						bestL[i][d] = L[i][d];

				}

			}
//			int u = stream.nextInt(100,250);
			if (verbose > 0)
				System.out.println("## Column  " + (d + 1) + " lowest Wafom so far " + lowestWafom + " log10  "
						+ Math.log10(lowestWafom) + " nbTries " + nbTrirage + " u ");

		}

		SearchInfos = String.format("LMS-EXTEND_SIZE, %d, %d, %d, %d, %s, %s ", w, k, dim, nBTries, lowestWafom,
				Math.log10(lowestWafom));
		System.out.println(SearchInfos);

		SobolSequence Sobol = new SobolSequence(filePath, k, w, dim);
		for (int j = 0; j < dim; j++) {
			Sobol.leftMultiplyMatBis2(j, bestL[j]);
		}

		return Sobol.getGeneratorMatricesTrans();
	}

	public static void main(String[] args) throws IOException {

		int w = 31; // when using LMS set to 30
		int q = 3;
		int l = w / q;
		int k = 17;
		int dim = 5;
		int nBTries = 300;

		double h = 1.0;
		double factor = 2.0;
		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);
		int verbose = 0;

		RandomStream stream = new MRG32k3a();
		int cc = (factor == 1.0) ? 108 : 56;

		int[] nbTirage = { 10, 20, 50, 300 };

//		for (int nBTries : nbTirage)
//
//			for (dim = 6; dim < 13; dim += 6)
//
//				for (k = 6; k < 21; k++) 
		{

			SearchStrategy searchStrategy = SearchStrategy.EXTEND_SIZE;

			String wafomType = (factor == 1.0) ? "Matsumoto" : "Goda";
			String searchType = searchStrategy.toString();
			System.out
					.println("          " + "   " + searchType + "  : w " + w + " k " + k + " dim " + dim + " nBTries "
							+ nBTries + " h " + h + " factor " + factor + " wafomType " + wafomType + " cc " + cc);
			System.out.println();

			LMSExtendSize searcher = new LMSExtendSize(w, k, dim, nBTries, table_c, q, h, factor, stream);

			long startTime = System.currentTimeMillis();
			int[] genMat = searcher.ExtendSize(verbose, cc);
			long endTime = System.currentTimeMillis();
			System.out.println("Duration: " + (endTime - startTime) / 1000.0 + " seconds");

			DigitalNetBase2 digitalNet = new DigitalNetBase2(k, w, dim, genMat);

			WafomFast wafom = new WafomFast(digitalNet, w, q, table_c);
			double wafomValue = wafom.computeWafom();

			wafomType = (factor == 1) ? " Matsumoto " : " Goda";
			System.out.println();
			System.out.print(wafomType + "    : " + wafomValue + " log10  " + Math.log10(wafomValue) + " | ");
			double factor2 = (factor == 1) ? 2 : 1;

			LookUpTableWafom table_c2 = new LookUpTableWafom(w, q, w / q, h, factor2);
			wafom = new WafomFast(digitalNet, w, q, table_c2);
			wafomValue = wafom.computeWafom();
			wafomType = (factor2 == 1) ? "Matsumoto" : "Goda";
			System.out.println(wafomType + "    : " + wafomValue);

			searcher.SearchInfos += String.format(" , %f seconds", (endTime - startTime) / 1000.0);

//			for (int i = 0; i < searcher.bestL.length; i++)
//				System.out.println(Arrays.toString(searcher.bestL[i]));
			stream.resetStartStream();
			System.out.println(
					"------------------------------------------------------------------------------------------");

		}

	}

}
