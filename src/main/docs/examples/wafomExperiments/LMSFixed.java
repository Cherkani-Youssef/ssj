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
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import wafomExperiments.HelperFunctions.SearchStrategy;

/**
 * This Algorithm was cited in @cite{Shin Harase. Quasi-Monte Carlo point sets
 * with small t-values and WAFOM. Applied Mathematics and Computation,
 * 254:318–326, 2015}
 * 
 * This class focuses on randomizing digital nets by applying a left matrix scramble
 * (LMS) on the generating matrices  and retaining the pointset with 
 * the lowest-Wafom out of the N trials.
 * 
 *
 * 
 */
public class LMSFixed {

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

	public LMSFixed(DigitalNetBase2 digitalNet, int w, int k, int dim, int nbTries, int q, double h, double factor,
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
		this.genMat =  digitalNet.getGeneratorMatricesTrans();

	}

	public LMSFixed(DigitalNetBase2 digitalNet, int w, int k, int dim, int nbTries, LookUpTableWafom table_c, int q,
			double h, double factor, RandomStream stream) {
		this(digitalNet, w, k, dim, nbTries, q, h, factor, stream);
		this.table_c = table_c;

	}

	public int[][] getLMS() {
		return bestL;
	}

	public int[] SobolLMS(int verbose) throws InterruptedException, ExecutionException {
		double lowestWafom = Double.MAX_VALUE;
		int[] genMat2 = new int[k * dim];

		ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		List<Future<AttemptResult>> futures = new ArrayList<>();

		for (int attempt = 0; attempt <= nBTries; attempt++) {
//			final int attemptIndex = attempt;
			Callable<AttemptResult> task = () -> {
				DigitalNetBase2 localDigitalNet = new DigitalNetBase2(k,w,dim,genMat);
//				localDigitalNet.resetGeneratorMatrices();

				int[][] tempL = localDigitalNet.leftMatrixScrambleBis(stream);

				WafomFast wafom = new WafomFast(localDigitalNet, w, q, table_c);

				double currentWafom = wafom.computeWafom();

				return new AttemptResult(currentWafom, tempL, localDigitalNet.getGeneratorMatricesTrans());
			};
			futures.add(executor.submit(task));
		}

		for (Future<AttemptResult> future : futures) {
			AttemptResult result = future.get();
			if (result.currentWafom < lowestWafom) {
				lowestWafom = result.currentWafom;
				genMat2 = result.genMat;
				bestL = result.tempL;

				for (int i = 0; i < dim; i++) {
					System.arraycopy(result.tempL[i], 0, bestL[i], 0, result.tempL[i].length);
				}
			}
		}

		executor.shutdown();

		if (verbose > 0 && nBTries > 100) {
			System.out.println("Lowest Wafom: " + Math.log10(lowestWafom));
		}

		return genMat2;
	}

	private static class AttemptResult {
		double currentWafom;
		int[][] tempL;
		int[] genMat;

		AttemptResult(double currentWafom, int[][] tempL, int[] genMat) {
			this.currentWafom = currentWafom;
			this.tempL = tempL;
			this.genMat = genMat;
		}
	}
	
	public static void main(String[] args) throws InterruptedException, ExecutionException {

		int w = 31;
		int q = 3;
		int l = w / q;
		int k = 15;
		int dim = 1;
		int nBTries = 10_000;
		int verbose = 0;

		double h = 1.0;
		double factor = 2.0;
		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);
		double DURATION = 0.0;
		RandomStream stream = new MRG32k3a();
		for (dim = 5; dim < 7; dim++) {

			double[] wafValues = new double[20];
			double[] wafValues2 = new double[20];
			double[] wafValuesG = new double[20];
			double[] wafValues2G = new double[20];

			for (k = 20; k < 21; k++) {
//
				SearchStrategy searchStrategy = SearchStrategy.LMS_FIX;
//
				String wafomType = (factor == 1.0) ? "Matsumoto" : "Goda";
				String searchType = searchStrategy.toString();
				System.out.println("          " + "   " + searchType + "  : w " + w + " k " + k + " dim " + dim
						+ " nBTries " + nBTries + " h " + h + " factor " + factor + " wafomType " + wafomType);
				System.out.println();

//				String filePathResults = "/Users/cher/Desktop/Projet/WafomResults/Results/" + searchType + "/"
//						+ wafomType + "/randomSearch-" + searchType + "-" + nBTries + ".txt";
//				String filePathMatrices = "/Users/cher/Desktop/Projet/WafomResults/GeneratingMatrices/" + searchType
//						+ "/" + wafomType + "/randomSearch-" + searchType + "-" + dim + "-" + nBTries + ".txt";

				DigitalNetBase2 digitalNet = new SobolSequence(filePath, k, w, dim);
//				DigitalNetBase2 digitalNet = new NiedXingSequenceBase2(filePath, k, w, dim);
				LMSFixed searcher = new LMSFixed(digitalNet, w, k, dim, nBTries, table_c, q, h, factor, stream);

//				LMSFixed searcher = new LMSFixed(w, k, dim, nBTries, table_c, q, h, factor, stream);

				long startTime = System.currentTimeMillis();
				int[] genMat = searcher.SobolLMS(verbose);
				long endTime = System.currentTimeMillis();

				System.out.println(", Duration: " + (endTime - startTime) / 1000.0 + " seconds");
				DURATION +=  (endTime - startTime) / 1000.0;
				DigitalNetBase2 digitalNetFound = new DigitalNetBase2(k, w, dim, genMat);
//				String matricesAsString = digitalNet.getGeneratorMatricesTransAsString(dim);
				
				
//				LowWafomPointsetsSearch.writeToFile(filePathMatrices, searchStrategy, dim, k, w, matricesAsString);

				WafomFast wafom = new WafomFast(digitalNetFound, w, q, table_c);
				double wafomValue = wafom.computeWafom();

				wafomType = (factor == 1) ? " Matsumoto " : " Goda";
				System.out.println();
				System.out.print(wafomType + "    : " + wafomValue + " | ");

				wafValues[k - 1] = Math.log10((wafomValue));
				wafValues2[k - 1] = (wafomValue);
				double factor2 = (factor == 1) ? 2 : 1;

				LookUpTableWafom table_c2 = new LookUpTableWafom(w, q, w / q, h, factor2);
				wafom = new WafomFast(digitalNetFound, w, q, table_c2);
				wafomValue = wafom.computeWafom();
				wafomType = (factor2 == 1) ? "Matsumoto" : "Goda";
				System.out.println(wafomType + "    : " + wafomValue);
				wafValuesG[k - 1] = Math.log10((wafomValue));
				wafValues2G[k - 1] = (wafomValue);

				searcher.SearchInfos += String.format(" , %f seconds", (endTime - startTime) / 1000.0);

				// Save the search results to a CSV file
//				LowWafomPointsetsSearch.saveSearchResults(searcher.SearchInfos, filePathResults);
//				CompareWafom.printMatrixInBinary(searcher.bestL, w);

//				Search.printJavaArrayDeclaration(searcher.bestL);
				stream.resetStartStream();

				System.out.println(
						"------------------------------------------------------------------------------------------");

			}
			System.out.println(", Duration: " +DURATION + " seconds");

//			System.out.println();
//			System.out.println(" " + dim + " : " + Arrays.toString(wafValues));
//			System.out.println(" " + dim + " : " + Arrays.toString(wafValues2));
//			System.out.println();
//			System.out.println(" " + dim + " : " + Arrays.toString(wafValuesG));
//			System.out.println(" " + dim + " : " + Arrays.toString(wafValues2G));
		}
		
	}

	
	



}
