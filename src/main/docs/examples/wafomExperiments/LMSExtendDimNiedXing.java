package wafomExperiments;

import java.io.IOException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.LookUpTableWafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.NiedXingSequenceBase2;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import wafomExperiments.HelperFunctions.SearchStrategy;

/**
 * Niederreiter-Xing
 * 
 * Random Component-By-Component (CBC) approach: This method starts with an
 * empty net, adding coordinates sequentially. For every new coordinate, the
 * algorithm aims to minimize the figure of merit. It does this by generating a
 * left Matrix Scramble and calculating the WAFOM (Walsh Figure of Merit) for
 * the scrambled point set within that dimension. After determining the optimal
 * LMS for the current coordinate, the algorithm progresses to the next
 * coordinate, with no backtracking to previously processed coordinates.
 */
public class LMSExtendDimNiedXing {

	static String filePath = "/Users/cher/Downloads/new-joe-kuo-6.21201";
	private int w;
	private int k;
	private int nBTries;
	private int dim;

	private RandomStream stream;
	private LookUpTableWafom table_c; // for the fast Wafom
	private int q;
	double h;
	double factor;
	String SearchInfos = " ";
	int[][] bestLMS;

	/**
	 * @param k      number of columns
	 * @param dim    dimension of the net
	 * @param w      output digits
	 * @param q      a positive integer
	 * @param l      length of the segment
	 * @param h      for Matsumoto's definition set to 0, and for Yoshiki's set to 1
	 * @param factor set to 1 for wafom, set to 2 to get Wafom for RMSE
	 */

	public LMSExtendDimNiedXing(int w, int k, int dim, int nbTries, LookUpTableWafom table_c, int q, double h,
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
		this.bestLMS = new int[dim][k];
	}

	/**
	 * @param: verbose: verbosity Level
	 */
	public int[] run(int verbose) {

		double lowestWafomForDim = Double.MAX_VALUE;
//		int numPoints = (1 << k);
//		double[] wafomValues = new double[numPoints];

		for (int s = 1; s <= dim; s++) {

			lowestWafomForDim = Double.MAX_VALUE;
//			double[] tempWafomValues = wafomValues.clone();

			int[][] tempLMS = new int[s][dim];
			for (int i = 0; i < s - 1; i++) {
				tempLMS[i] = bestLMS[i].clone();
			}
			DigitalNetBase2 NiedXing = new NiedXingSequenceBase2(k, w, s);

			for (int attempt = 0; attempt < nBTries; attempt++) {
//				double[] wafomValues2 = tempWafomValues.clone();
//				if (s == 1)
//					Arrays.fill(wafomValues2, 1); // Fill the array with ones

				NiedXing.resetGeneratorMatrices();

				// Generate LMS for dimension s
				int[] currentLMS = HelperFunctions.GenerateleftMatrixScrambleBis(stream, k, w);

				tempLMS[s - 1] = currentLMS.clone();
				// Compute L*C up to dimension s
				for (int j = 0; j < s; j++) {
					NiedXing.leftMultiplyMatBisV2(j, tempLMS[j], k);
				}

				WafomFast wafom = new WafomFast(NiedXing, w, q, table_c);

//				double[] wafomValuesForDim = wafom.computeWafomForCBC(s - 1);
//				double sum = wafom.elementWiseProd2(wafomValues2, wafomValuesForDim);
//				double currentWafom2 = sum / numPoints; // wafom.sumKahan(wafomValues2) / numPoints;
				double currentWafom = wafom.computeWafom();
//				System.out.println("currentWafom " + currentWafom + " true " + currentWafom2 + " diff  "
//						+ Math.abs(currentWafom - currentWafom2) );

				// Update the best matrix for this dimension if the current Wafom is lower
				if ((currentWafom) < (lowestWafomForDim)) {

					lowestWafomForDim = currentWafom;
//					wafomValues = wafomValues2.clone();
					bestLMS[s - 1] = currentLMS.clone();

				}
			}
			// Logging based on verbosity level
			if (verbose > 0) {
				System.out.println("Dim " + s + " Lowest Wafom for Dim: " + lowestWafomForDim + " log 10: "
						+ Math.log10((lowestWafomForDim)) + " nbTries " + nBTries);
			}

		}
		System.out.println("lowestWafomForDim " + lowestWafomForDim + " log10  " + Math.log10((lowestWafomForDim)));



		NiedXingSequenceBase2 NiedXing = new NiedXingSequenceBase2(k, w, dim);
		NiedXing.resetGeneratorMatrices();
		for (int j = 0; j < dim; j++) {
			NiedXing.leftMultiplyMatBisV2(j, bestLMS[j], k);
		}

//		WafomFast2 wafom = new WafomFast2(NiedXing, w, q, table_c);
//		double wafomValue = wafom.computeWafom();
//		System.out.println("  wafomValue  : " + wafomValue);

		return NiedXing.getGeneratorMatricesTrans();
	}

	public static void main(String[] args) throws IOException {

		int w = 31;
		int q = 3;
		int l = w / q;
		int k = 15;
		int dim = 5;
		int nBTries = 100;

		double h = 1.0;
		double factor = 1.0;

		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);
		int verbose = 1;

		RandomStream stream = new MRG32k3a();

//		int[] nbTirage = { 10, 20, 50, 300 };
//
//		for (int nBTries : nbTirage)
//
//			for (dim = 6; dim < 13; dim += 6)
//
//				for (k = 6; k < 21; k++) 

		{

			SearchStrategy searchStrategy = SearchStrategy.EXTEND_DIM;

			String searchType = searchStrategy.toString();

			System.out.println("              " + searchType + "  : w " + w + " k " + k + " dim " + dim + " nBTries "
					+ nBTries + " h " + h + " factor " + factor);
			System.out.println();

			String wafomType = (factor == 1) ? "Matsumoto" : "Goda";

//					String filePathResults = "/Users/cher/Desktop/Projet/WafomResults/Results/" + searchType + "/"
//							+ wafomType + "/randomSearch-" + searchType + "-" + dim + "-" + nBTries + ".txt";

			String filePathMatrices = "/Users/cher/Desktop/Projet/WafomResults/GeneratingMatrices/" + searchType + "/"
					+ wafomType + "/randomSearch-" + searchType + "-" + dim + "-" + nBTries + ".txt";

			LMSExtendDimNiedXing searcher = new LMSExtendDimNiedXing(w, k, dim, nBTries, table_c, q, h, factor, stream);

			long startTime = System.currentTimeMillis();
			int[] genMat = searcher.run(verbose);
			long endTime = System.currentTimeMillis();
			System.out.println("Duration: " + (endTime - startTime) / 1000.0 + " seconds");

			DigitalNetBase2 digitalNet = new DigitalNetBase2(k, w, dim, genMat);
			String matricesAsString = digitalNet.getGeneratorMatricesTransAsString(dim);

//			HelperFunctions.writeToFile(filePathMatrices, searchStrategy, dim, k, w, matricesAsString);

//					table_c = new LookUpTableWafom(w, q, l, h, factor);
			WafomFast wafom = new WafomFast(digitalNet, w, q, table_c);
			double wafomValue = wafom.computeWafom();
			System.out.print(wafomType + "    : " + wafomValue + " | ");

			wafomType = (factor == 1) ? " Matsumoto " : " Goda";
			System.out.println();

			double factor2 = (factor == 1) ? 2 : 1;

			LookUpTableWafom table_c2 = new LookUpTableWafom(w, q, w / q, h, factor2);
			wafom = new WafomFast(digitalNet, w, q, table_c2);
			wafomValue = wafom.computeWafom();
			wafomType = (factor2 == 1) ? "Matsumoto" : "Goda";
			System.out.print(wafomType + "    : " + wafomValue);

			searcher.SearchInfos += String.format(" , %f seconds", (endTime - startTime) / 1000.0);
			System.out.println("  ---->  Duration " + ((endTime - startTime) / 1000.0));
			// Save the search results to a CSV file
//					LowWafomPointsetsSearch.saveSearchResults(searcher.SearchInfos, filePathResults);
			stream.resetStartStream();

			System.out.println(
					"------------------------------------------------------------------------------------------");

		}

	}

}
