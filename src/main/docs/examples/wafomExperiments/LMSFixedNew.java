package wafomExperiments;

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import umontreal.ssj.discrepancy.LookUpTableWafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
//import umontreal.ssj.hups.NiedXingSequenceBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.LFSR258;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import wafomExperiments.HelperFunctions.SearchStrategy;

public class LMSFixedNew {
	static String filePath = "/Users/cher/Downloads/new-joe-kuo-6.21201";
	private DigitalNetBase2 digitalNet;
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
	int[][] bestL;
	int[] genMat;

	/**
	 *
	 * @param DigitalNetBase2
	 * @param k               number of columns
	 * @param dim             dimension of the net
	 * @param w               output digits
	 * @param q               a positive integer
	 * @param l               length of the segment
	 * @param h               for Matsumoto's definition set to 0, and for Yoshiki's
	 *                        set to 1
	 * @param factor          set to 1 for wafom, set to 2 to get Wafom for RMSE
	 */

	public LMSFixedNew(DigitalNetBase2 digitalNet, int w, int k, int dim, int nbTries, int q, double h, double factor,
			RandomStream stream) {
		this.digitalNet = digitalNet;
		this.w = w;
		this.k = k;
		this.dim = dim;
		this.nBTries = nbTries;
		this.stream = stream;

		this.q = q;
		this.h = h;
		this.factor = factor;
		this.bestL = new int[dim][k];

		this.table_c = new LookUpTableWafom(w, q, w / q, h, factor);
		this.genMat = digitalNet.getGeneratorMatricesTrans();

	}

	public LMSFixedNew(DigitalNetBase2 digitalNet, int w, int k, int dim, int nbTries, LookUpTableWafom table_c, int q,
			double h, double factor, RandomStream stream) {
		this(digitalNet, w, k, dim, nbTries, q, h, factor, stream);
		this.table_c = table_c;

	}

	public int[][] getLMS() {
		return bestL;
	}

	public int[] SobolLMS(int verbose) {
		ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		double lowestDistance = Double.MAX_VALUE;
//        int[] genMat = new int[k * dim];
//        double[] bestWafomValues = new double[k];
		double[] fixed = (factor == 1.0) ? WafomValues.dataMatsu.get(dim) : WafomValues.dataGoda.get(dim);
	
		List<Future<TrialResult>> futures = new ArrayList<>();
		for (int attempt = 0; attempt <= nBTries; attempt++) {
			int finalAttempt = attempt; // for use in the lambda
			Callable<TrialResult> trialTask = () -> {

				DigitalNetBase2 localNet = new DigitalNetBase2(k, w, dim, genMat);
//                localNet.resetGeneratorMatrices();
				int[][] tempL = localNet.leftMatrixScrambleBis(stream);
				
				
//				double[] wafomValues = computeWafoms(localNet);
//				double distance = (factor == 1.0) ? computeDistance(localNet) : computeDistance(localNet, fixed);
//				double distance = computeDistance(localNet);
				double distance = computeDistance(localNet, fixed);
				// Print intermediate results if verbose mode is enabled
				if ((verbose > 0) && ((nBTries > 100 && finalAttempt % 100 == 0) )) {
					System.out.print(" -> attempt " + finalAttempt + "/" + nBTries);
					stream.resetNextSubstream();
					
				}

				return new TrialResult(tempL, localNet.getGeneratorMatricesTrans(), distance);
			};
			futures.add(executor.submit(trialTask));
		}
		String ANSI_RED = "\u001B[31m";
		String ANSI_RESET = "\u001B[0m";

		// Print best distance

		// Collect results
		for (Future<TrialResult> future : futures) {
			try {
				TrialResult result = future.get();
				if (result.distance < lowestDistance) {
					lowestDistance = result.distance;
					genMat = result.genMat;
					bestL = result.tempL;

					System.out.println(ANSI_RED + " \n Updated lowest distance: " + lowestDistance + " " + ANSI_RESET
							+ "  weights " + weights);
				}
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace(); // handle exceptions appropriately
			}
		}

		executor.shutdown(); // Always remember to shutdown the executor

		if (verbose > 0) {
			System.out.println("Completed all trials. Lowest distance: " + lowestDistance);
		}

		return genMat;
	}

	private class TrialResult {
		int[][] tempL;
		int[] genMat;
//        double[] wafomValues;
		double distance;

		public TrialResult(int[][] tempL, int[] genMat, double distance) {
			this.tempL = tempL;
			this.genMat = genMat;
//            this.wafomValues = wafomValues;
			this.distance = distance;
		}
	}

	public double[] computeWafoms(DigitalNetBase2 localDigitalNet) {
//		int end = 0;

		double[] wafomValues = new double[k];
		int end = 0;
		double wafVal = 0.0;
		for (int d = 1; d <= k; d++) {

			int start = (d == 1) ? 0 : (end);
			end = (1 << (d));
			int numPoints = end;

			WafomFast wafom = new WafomFast(localDigitalNet, w, q, table_c);
			double somme = wafom.computeWafomSum(start, end) + wafVal;// / end
			double currentWafom = (somme) / numPoints;
			wafomValues[d - 1] = Math.log10(Math.abs(currentWafom));

			wafVal = somme;

		}

		return wafomValues;
	}

	static double weights = 1.0;

	public double computeDistance(DigitalNetBase2 localDigitalNet) {
//		int end = 0;

//		double[] wafomValues = new double[k];
		int end = 0;
		double wafVal = 0.0;
		double sum = 0.0;
		WafomFast wafom = new WafomFast(localDigitalNet, w, q, table_c);
		for (int d = 1; d <= k; d++) {

			int start = end;
			end = (1 << (d));
			int numPoints = end;

			double somme = wafom.computeWafomSum(start, end) + wafVal;// / end
			double currentWafom = (somme) / numPoints;

			double bound = -((1. * (d) * (d)) / (dim * 1.));
			if (d < k)
				sum += Math.abs(Math.log10(Math.abs(currentWafom)) - (bound));
			else
				sum += weights * Math.abs(Math.log10(Math.abs(currentWafom)) - (bound));
//			sum += Math.pow(  Math.log10(Math.abs(currentWafom)) - (bound),2);
			wafVal = somme;

		}
//		return Math.sqrt(sum);
		return sum;
	}

	public double computeDistance(DigitalNetBase2 localDigitalNet, double[] fixed) {

		int end = 0;
		double wafVal = 0.0;
		double sum = 0.0;
		WafomFast wafom = new WafomFast(localDigitalNet, w, q, table_c);
		for (int d = 1; d <= k; d++) {

			int start = end;
			end = (1 << (d));
			int numPoints = end;
			double somme = wafom.computeWafomSum(start, end) + wafVal;// / end
			double currentWafom = (somme) / numPoints;

			if (d < 20)
				sum += Math.pow(Math.log10(Math.abs(currentWafom)) - (fixed[d - 1]), 2);
			else
				sum += weights * Math.pow(Math.log10(Math.abs(currentWafom)) - (fixed[d - 1]), 2);

			wafVal = somme;

		}
		return Math.sqrt(sum);
//		return sum;
	}

	public double distance(double[] wafomValues) {
		double sum = 0;

		for (int j = 0; j < k; j++) {
			double bound = -((1. * (j + 1.) * (j + 1.)) / (dim * 1.));
			sum += Math.abs((wafomValues[j]) - (bound));
		}
		return sum;

	}

	public double distance(double[] wafomValues, double[] fixed) {
		double sum = 0;

		for (int j = 0; j < wafomValues.length; j++) {
			sum += Math.pow((wafomValues[j]) - (fixed[j]), 2);
//				
		}
//			
		return Math.sqrt(sum);

	}

	public static void main(String[] args) {
//		[502.8596695762276, 502.6849436213789]      [405.34588863501506, 402.4530766056343]
		// fixed
//		goda abs 392.76476566317285 pow 123.0457918438694
//		matsu  abs 499.7147929764095 pow 149.07395339170435

		int w = 31;
		int q = 3;
		int l = w / q;
		int k = 20;
		int dim = 5;
		int nBTries = 100_000;
		int verbose = 1;

		double h = 1.0;
		double factor = 2.0;
		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);
		double DURATION = 0.0;
		RandomStream stream = new MRG32k3a();

		SearchStrategy searchStrategy = SearchStrategy.LMS_FIX;

		String wafomType = (factor == 1.0) ? "Matsumoto" : "Goda";
		String searchType = searchStrategy.toString();
		System.out.println("          " + "   " + searchType + "  : w " + w + " k " + k + " dim "
				+ dim + " nBTries " + nBTries + " h " + h + " factor " + factor + " wafomType " + wafomType);
		System.out.println();

		DigitalNetBase2 digitalNet = new SobolSequence(filePath, k, w, dim);

		LMSFixedNew searcher = new LMSFixedNew(digitalNet, w, k, dim, nBTries, table_c, q, h, factor, stream);

		long startTime = System.currentTimeMillis();
		int[] genMat = searcher.SobolLMS(verbose);
		long endTime = System.currentTimeMillis();

		System.out.println(", Duration: " + (endTime - startTime) / 1000.0 + " seconds");
		DURATION += (endTime - startTime) / 1000.0;
		DigitalNetBase2 digitalNetFound = new DigitalNetBase2(k, w, dim, genMat);

		double[] wafomValues = searcher.computeWafoms(digitalNetFound);

		double[] fixed = (factor == 1.0) ? WafomValues.dataMatsu.get(dim) : WafomValues.dataGoda.get(dim);
		double distance = searcher.distance(wafomValues, fixed);
		System.out.println("Distance " + distance + "   " + Arrays.toString(wafomValues));
		double[][] data = { fixed, wafomValues };
		String[] labels = new String[] { "LMS-Fix", "New" };
		String chartTitle = searchType + " " + wafomType + " Dimension " + dim + " nBTries " + nBTries;
		String xAxisLabel = "number of Points";
		String yAxisLabel = "log10 Wafom";

		ChartHelper chartHelper = new ChartHelper(chartTitle, xAxisLabel, yAxisLabel, data, labels, 1);
		chartHelper.displayChart();

		WafomFast wafom = new WafomFast(digitalNetFound, w, q, table_c);
		double wafomValue = wafom.computeWafom();

		wafomType = (factor == 1) ? " Matsumoto " : " Goda";
		System.out.println();
		System.out.print(wafomType + "    : " + wafomValue + " | ");

//				wafValues[k - 1] = Math.log10((wafomValue));
//				wafValues2[k - 1] = (wafomValue);
		double factor2 = (factor == 1) ? 2 : 1;

		LookUpTableWafom table_c2 = new LookUpTableWafom(w, q, w / q, h, factor2);
		wafom = new WafomFast(digitalNetFound, w, q, table_c2);
		wafomValue = wafom.computeWafom();
		wafomType = (factor2 == 1) ? "Matsumoto" : "Goda";
		System.out.println(wafomType + "    : " + wafomValue);
//				wafValuesG[k - 1] = Math.log10((wafomValue));
//				wafValues2G[k - 1] = (wafomValue);

		searcher.SearchInfos += String.format(" , %f seconds", (endTime - startTime) / 1000.0);

		// Save the search results to a CSV file
//				LowWafomPointsetsSearch.saveSearchResults(searcher.SearchInfos, filePathResults);
//				CompareWafom.printMatrixInBinary(searcher.bestL, w);

		Search.printJavaArrayDeclaration(searcher.bestL);
		stream.resetStartStream();

		System.out
				.println("------------------------------------------------------------------------------------------");

	}

}
