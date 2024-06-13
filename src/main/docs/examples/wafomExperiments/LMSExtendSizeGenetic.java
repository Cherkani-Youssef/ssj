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





public class LMSExtendSizeGenetic extends Search {

	public LMSExtendSizeGenetic(int w, int k, int dim, int populationSize, LookUpTableWafom table_c, int q, double h,
			double factor, RandomStream stream) {
		super(w, k, dim, populationSize, table_c, q, h, factor, stream);
	}

	public void run(int verbose, int end, double exponent, double slowExponent) {

		final int twoRm1 = (1 << (w - 1));

//		int[] nPop = nonlinearInterpolationSlow(populationSize, end, 20, exponent);
		int[] nPop = nonlinearInterpolationVariableSmooth(populationSize, end, 20, exponent, slowExponent, 6);

		System.out.println(Arrays.toString(nPop) + " length " + nPop.length);
		SobolSequence Sob = null;
		Listener listener = new Listener(dim, 1, populationSize);
		populationSize = nPop[0];
//		int h = populationSize;
		int dd = 21;
		int N = 2;

		for (int d = 0; d < k; d++) {

			for (int s = 1; s <= dim; s++) {
				listener.currentIndex = 0;

				if (d >= dd) {
					listener.resetMeritValuesArray();
					Sob = new SobolSequence(filePath, d + 1, 31, s);
				}
				for (int i = 0; i < populationSize; i++) {
					int[][] tempL = listener.getLeftMatrixScramble(i);
					if (d < dd) {

						tempL[s - 1][d] = (twoRm1 + stream.nextInt(0, twoRm1 - 1)) >> (d);
						listener.add(tempL);
					} else {

						for (int attempt = 0; attempt < N; attempt++) {

							tempL[s - 1][d] = (twoRm1 + stream.nextInt(0, twoRm1 - 1)) >> (d);
							Sob.resetGeneratorMatrices();

							for (int j = 0; j < s; j++) {
								Sob.leftMultiplyMatBisV2(j, tempL[j], d + 1);
							}

							WafomFast wafom = new WafomFast(Sob, w, q, table_c);
							double currentWafom = wafom.computeWafom();
							listener.update(tempL, currentWafom, i);
						}

					}
				}

			}
			
			

			if (d < k - 1) {
				if (d < dd) {

					System.out.print(" =>COL " + (d + 1)+ " populationSize " + populationSize + "  "); // + " populationSize " + populationSize

					int[][][] tempL = listener.getPopulation();
					int[] bestIndices = selectBest(populationSize, nPop[d + 1], tempL, dim, d + 1, true);
					populationSize = bestIndices.length;

					listener = new Listener(tempL, dim, populationSize, bestIndices);
				} else {
					System.out.print(" -> COL " + (d + 1)+ " populationSize " + populationSize + " "); // + " populationSize " + populationSize

					int[][][] tempL = listener.getPopulation();
//					Integer[] indices = IntStream.range(0, populationSize).boxed().toArray(Integer[]::new);
//					double[] wafomValues = listener.getMeritValues();
//					Arrays.sort(indices, Comparator.comparingDouble(i -> wafomValues[i]));
//					int[] bestIndices = Arrays.stream(indices).limit(nPop[d + 1]).mapToInt(Integer::intValue).toArray();
					int[] bestIndices = selectBest(populationSize, nPop[d + 1], tempL, dim, d + 1, true);
					populationSize = bestIndices.length;
					listener = new Listener(tempL, dim, populationSize, bestIndices);

				}
			}

		}
		System.out.println();

		bestLs = listener.getPopulation();

	}

	public static void main(String[] args) throws IOException {
		double factor = 2.0;

		int w = 31; // when using LMS set to 30
		int q = 3;
		int l = w / q;
		int k = 20;
		int dim = 6;
//		int nBTries = 1;// 10
//		int r = 1;

		double slowExponent = 1;
		int populationSize = 10_000;

		double exponent = 0.2;
//		int n = (factor == 2.0) ? 2 : 3;
		int n = 4;
		double h = 1.0;

		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);
		int verbose = 1;
		RandomStream stream = new MRG32k3a();

		
		{
			System.out.println(
					"                                   ************************************************************************");


			stream.resetStartStream();

			SearchStrategy searchStrategy = SearchStrategy.EXTEND_SIZE;
			String searchType = searchStrategy.toString();

			String wafomType = (factor == 1) ? "Matsumoto" : "Goda";
			System.out.println("              " + searchType + "  : w " + w + " k " + k + " dim " + dim + " h " + h
					+ " factor " + factor + " populationSize " + populationSize + " n " + n + " exponent " + exponent
					+ " slowExponent " + slowExponent);
//					System.out.println();

			LMSExtendSizeGenetic searcher = new LMSExtendSizeGenetic(w, k, dim, populationSize, table_c, q, h, factor,
					stream);
			long startTime = System.currentTimeMillis();
			searcher.run(verbose, n, exponent, slowExponent);
			long endTime = System.currentTimeMillis();

			double[] fixed = (factor == 1.0) ? lmsFixedMatsu : lmsFixedGoda;
			int[] indices = searcher.selectBest(searcher.populationSize, n, searcher.bestLs, dim, k, true);
			System.out.println(Arrays.toString(indices) + " Duration " + (endTime - startTime) / 1000.0 + " seconds");
			double[][] wafomValues = searcher.calcWafomPopulation(fixed, indices);

			System.out.println(Arrays.toString(wafomValues[indices[0]]));
//			for (int i = 0; i < wafomValues.length; i++) {
//				System.out.println(Arrays.toString(wafomValues[i]));
//				System.out.println();
//			}

				String chartTitle = searchType + " " + wafomType + " Dimension " + dim + " Population " + populationSize
						+ " n " + n + " exponent " + exponent + " slowExponent " + slowExponent;
				String xAxisLabel = "number of Points";
				String yAxisLabel = "log10 Wafom";

				String[] labels = new String[indices.length + 1]; // Create a String array of the same length

				for (int i = 0; i < indices.length; i++) {
					labels[i] = String.valueOf(indices[i]); // Convert each int to String
				}
				labels[indices.length] = "FIX";

				ChartHelper chartHelper = new ChartHelper(chartTitle, xAxisLabel, yAxisLabel, wafomValues, labels, 1);
				chartHelper.displayChart();

			System.out.println("************************************************************************");
			for (int index : indices)
				printJavaArrayDeclaration(searcher.bestLs[index]);
			System.out.println("********************************************************************************************************");
			
			System.out.println("              " + searchType + "  : w " + w + " k " + k + " dim " + dim + " h " + h
					+ " factor " + factor + " populationSize " + populationSize + " n " + n + " exponent " + exponent
					+ " slowExponent " + slowExponent);

		}

	}

}
