package wafomExperiments;

import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

import java.util.stream.IntStream;

import umontreal.ssj.discrepancy.LookUpTableWafom;

import umontreal.ssj.discrepancy.WafomFast;

import umontreal.ssj.hups.SobolSequence;

import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;

import wafomExperiments.HelperFunctions.SearchStrategy;

public class LMSExtendTwoGenetic extends Search {

	public LMSExtendTwoGenetic(int w, int k, int dim, int populationSize, LookUpTableWafom table_c, int q, double h,
			double factor, RandomStream stream) {
		super(w, k, dim, populationSize, table_c, q, h, factor, stream);
	}

	public void run(int verbose, int end, double exponent, double slowExponent) {

		final int twoRm1 = (1 << (w - 1));

//		int[] nPop = nonlinearInterpolationSlow(populationSize, end, 20, exponent);
		int[] nPop = nonlinearInterpolationVariableSmooth(populationSize, end, 20 * dim, exponent, slowExponent, 6*dim);
//		nPop[19] = nPop[18];
//		System.out.println(Arrays.toString(nPop) + " length " + nPop.length);

		Listener listener = new Listener(dim, 1, populationSize);
		populationSize = nPop[0];
//		int h = populationSize;
		int dd = -1;
		int N = 5;
		int index = 0;
		for (int d = 0; d < k; d++) {
			
			System.out.print("  ---->COL " + (d + 1));

			for (int s = 1; s <= dim; s++) {
				SobolSequence Sob = new SobolSequence(filePath, d + 1, 31, s);
//				int[][][] tempL = listener.getPopulation();
				Listener listener2 = new Listener(dim, d + 1, populationSize);
				listener.currentIndex = 0;

				if (d >= dd)
					listener2.resetMeritValuesArray();

				for (int i = 0; i < populationSize; i++) {
					int[][] tempL = listener.getLeftMatrixScramble(i);
					if (d < dd) {

						tempL[s - 1][d] = (twoRm1 + stream.nextInt(0, twoRm1 - 1)) >> (d);
						listener2.add(tempL);
					} else {

						for (int attempt = 0; attempt < N; attempt++) {
							Sob.resetGeneratorMatrices();
							tempL[s - 1][d] = (twoRm1 + stream.nextInt(0, twoRm1 - 1)) >> (d);

							

							for (int j = 0; j < s; j++) {
								Sob.leftMultiplyMatBisV2(j, tempL[j], d + 1);
							}

							WafomFast wafom = new WafomFast(Sob, w, q, table_c);
							double currentWafom = wafom.computeWafom();
							listener2.update(tempL, currentWafom, i);
						}

					}
				}
				
//				System.out.print(" =>Dim " + (s)  + "  populationSize " + populationSize);

//				int[][][] tempL = listener2.getPopulation();
//				int[] bestIndices = selectBest(populationSize, nPop[index++], tempL, s, d + 1, false);
				Integer[] indices = IntStream.range(0, populationSize).boxed().toArray(Integer[]::new);
				double[] wafomValues = listener2.getMeritValues();
				Arrays.sort(indices, Comparator.comparingDouble(i -> wafomValues[i]));
				int[] bestIndices = Arrays.stream(indices).limit(nPop[index++]).mapToInt(Integer::intValue).toArray();

				populationSize = bestIndices.length;

				for (int i = 0; i < populationSize; i++) {

					listener.leftMatrixScrambles[i] = listener2
							.copyMatrix(listener2.leftMatrixScrambles[bestIndices[i]]);
				}

			}
			if (d < k - 1) {
				int[][][] tempL = listener.getPopulation();
				listener = new Listener(tempL, dim, populationSize);
			}

		}
		System.out.println();

		bestLs = listener.getPopulation();

	}

	public static void main(String[] args) throws IOException {
		double factor = 1.0;

		int w = 31; // when using LMS set to 30
		int q = 3;
		int l = w / q;
		int k = 20;
		int dim = 5;
//		int nBTries = 1;// 10
//		int r = 1;

		double slowExponent = 1;
		int populationSize = 100;

		double exponent = 0.05;
//		int n = (factor == 2.0) ? 2 : 3;
		int n = 3;
		double h = 1.0;

		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);
		int verbose = 1;
		RandomStream stream = new MRG32k3a();
//
		int[] pop = new int[] { 1000, 2000, 5000, 10_000, 20_000, 30_000 };// , 40_000, 50_000, 60_000
		double[] expo = new double[] { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };
//		double[] expo2 = new double[] { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };
////////
		for (int p : pop)
			for (double x : expo)
//				for (double x2 : expo2)
			{
				System.out.println(
						"                                   ************************************************************************");
				populationSize = p;
//// 				exponent = 0.8;
				exponent = x;
////			slowExponent = x2;
				stream.resetStartStream();

				SearchStrategy searchStrategy = SearchStrategy.EXTEND_TWO;
				String searchType = searchStrategy.toString();

//				String wafomType = (factor == 1) ? "Matsumoto" : "Goda";
				System.out.println("              " + searchType + "  : w " + w + " k " + k + " dim " + dim + " h " + h
						+ " factor " + factor + " populationSize " + populationSize + " n " + n + " exponent "
						+ exponent + " slowExponent " + slowExponent);
//					System.out.println();

				LMSExtendTwoGenetic searcher = new LMSExtendTwoGenetic(w, k, dim, populationSize, table_c, q, h, factor,
						stream);
				long startTime = System.currentTimeMillis();
				searcher.run(verbose, n, exponent, slowExponent);
				long endTime = System.currentTimeMillis();

//				double[] fixed = (factor == 1.0) ? lmsFixedMatsu : lmsFixedGoda;
				int[] indices = searcher.selectBest(searcher.populationSize, n, searcher.bestLs, dim, k, true);
//					System.out.println(
//							Arrays.toString(indices) + " Duration " + (endTime - startTime) / 1000.0 + " seconds");
//				double[][] wafomValues = searcher.calcWafomPopulation(fixed, indices);
//
//			for (int i = 0; i < wafomValues.length; i++) {
//				System.out.println(Arrays.toString(wafomValues[i]));
//				System.out.println();
//			}

//				String chartTitle = searchType + " " + wafomType + " Dimension " + dim + " Population " + populationSize
//						+ " n " + n + " exponent " + exponent + " slowExponent " + slowExponent;
//				String xAxisLabel = "number of Points";
//				String yAxisLabel = "log10 Wafom";
//
//				String[] labels = new String[indices.length + 1]; // Create a String array of the same length
//
//				for (int i = 0; i < indices.length; i++) {
//					labels[i] = String.valueOf(indices[i]); // Convert each int to String
//				}
//				labels[indices.length] = "FIX";
//
//				ChartHelper chartHelper = new ChartHelper(chartTitle, xAxisLabel, yAxisLabel, wafomValues, labels, 1);
//				chartHelper.displayChart();

//			System.out.println("              " + searchType + "  : w " + w + " k " + k + " dim " + dim + " h " + h
//					+ " factor " + factor + " populationSize " + populationSize + " n " + n + " exponent " + exponent
//					+ " slowExponent " + slowExponent

//				System.out.println("************************************************************************");
//				for (int index : indices)
//					printJavaArrayDeclaration(searcher.bestLs[index]);
//				System.out.println("************************************************************************");

			}

	}

}
