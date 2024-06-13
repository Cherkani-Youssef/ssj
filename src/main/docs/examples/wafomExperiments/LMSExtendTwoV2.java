package wafomExperiments;

import java.io.IOException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.LookUpTableWafom;
import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.mcqmctools.RQMCExperiment;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import wafomExperiments.HelperFunctions.SearchStrategy;

/**
 * We attempt to find extensible low-Wafom digital nets in both dimension and
 * number of points
 * 
 * We populate we sample N random choices for each column of the first dimension
 * and select the column that gives the best Wafom update, then we move to the
 * next dimensions without ever looking back.
 */

public class LMSExtendTwoV2 {

	private int w;
	private int k;
	int nBTries;
	private int dim;
	private RandomStream stream;
	private LookUpTableWafom table_c; // for the fast Wafom
	private int q;
	double h;
	double factor;
	String SearchInfos = " ";

	static String filePath = "/Users/cher/Downloads/new-joe-kuo-6.21201";

	int[][] bestL;
	double[] wafomValuesStorage;// = new double[1 << k];

	public LMSExtendTwoV2(int w, int k, int dim, int nbTries, LookUpTableWafom table_c, int q, double h, double factor,
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
		this.bestL = new int[dim][k];
		this.wafomValuesStorage = new double[1 << k];
		Arrays.fill(wafomValuesStorage, 1);

	}

	public void copyBestL(int[][] tempL, int s, int d) {

		for (int i = 0; i < s; i++) {

			System.arraycopy(bestL[i], 0, tempL[i], 0, d);
		}

	}

	public int[] ExtendTwo(int verbose) throws IOException {
		final int twoRm1 = (1 << (w - 1));
		RandomStream stream2 = new MRG32k3a();
		int nbTirage = nBTries;
		double lowestWafomForDim = Double.MAX_VALUE;
		final int twoR = (1 << (w));
//		SobolSequence Sobol = new SobolSequence(filePath, k, w, dim);
//		int[] genMat = Sobol.getGeneratorMatricesTrans();
		int end = 0;
//		int[][] tempL = new int[dim][k];
		KahanSummation kahanSum = new KahanSummation();

		for (int s = 1; s <= dim; s++) {
			double lastSum = 0.0;
			double lastSumTemp = 0.0;
			kahanSum = new KahanSummation(); // Reset for each dimension
			nbTirage = nBTries;

			for (int d = 1; d <= k; d++) {

				int start = (d == 1) ? 0 : (end);
				end = (1 << d);

				lowestWafomForDim = Double.MAX_VALUE;
				int longueur = end - start;
				int[][] tempL = new int[s][d];
				copyBestL(tempL, s, d);

				double[] tempWafomValues = new double[longueur];
				System.arraycopy(wafomValuesStorage, start, tempWafomValues, 0, longueur);

				int numPoints = end;
				SobolSequence Sob = new SobolSequence(filePath, d, w, s);
				for (int attemptCol = 0; attemptCol < nbTirage; attemptCol++) {
					double[] wafomValues2 = tempWafomValues.clone();
//					if (s == 1)
//						Arrays.fill(wafomValues2, 1);

					Sob.resetGeneratorMatrices();
					tempL[s - 1][d - 1] = ((stream.nextInt(twoRm1, twoR - 1)) >> (d - 1));

					for (int j = 0; j < s; j++) {
						Sob.leftMultiplyMatBisV2(j, tempL[j], d);
					}
					WafomFast wafom = new WafomFast(Sob, w, q, table_c);

					double[] currentWafomValues = wafom.computeWafom2(start, end, s - 1);

//					wafom.elementWiseProd(wafomValues2, currentWafomValues);
					double somme = wafom.elementWiseProd2(wafomValues2, currentWafomValues); // wafom.sumKahan(wafomValues2);

//					double currentWafom =  -1+(lastSum + somme ) / numPoints;

					KahanSummation kahanSum2 = new KahanSummation();
					kahanSum2.add(kahanSum.getSum());
					kahanSum2.add(somme);
//					double currentWafom = -1 + kahanSum2.getSum() / numPoints;
					double currentWafom = kahanSum2.getSum() / numPoints;
//					double currWafom = wafom.computeWafom();
//					System.out.println(" currentWafom Cumulative " + currentWafom + " currWafom Original " + currWafom
//							+ " error " + Math.abs(currWafom - currentWafom));

					if ((currentWafom) < (lowestWafomForDim)) {
						lowestWafomForDim = currentWafom;
						bestL[s - 1][d - 1] = tempL[s - 1][d - 1];
						System.arraycopy(wafomValues2, 0, wafomValuesStorage, start, longueur);
						lastSumTemp = somme;
					}

				}
				lastSum += lastSumTemp;
				kahanSum.add(lastSumTemp); // U

			}

			if (verbose > 0)
				System.out.println("## Dimension  " + (s) + " lowest Wafom so far " + lowestWafomForDim + " log10  "
						+ Math.log10(lowestWafomForDim) + " nbTirage " + nbTirage);
		}

		System.out.println();
		SobolSequence Sobol = new SobolSequence(filePath, k, w, dim);
		for (int j = 0; j < dim; j++) {
			Sobol.leftMultiplyMatBisV2(j, bestL[j], k);
		}
		LookUpTableWafom table_c = new LookUpTableWafom(w, q, w / q, h, factor);
		WafomFast wafom = new WafomFast(Sobol, w, q, table_c);
		double wafomValue = wafom.computeWafom();
		System.out.println(" ---> Wafom : " + wafomValue + " log10 " + Math.log10(wafomValue));
		System.out.println();

		return Sobol.getGeneratorMatricesTrans();
	}

	public static void main(String[] args) throws IOException { // 77 126 146
		int w = 31;
		int k = 15;
		int dim = 5;
		int q = 3;
		int l = w / q;
		double factor = 2.0;
		double h = 1.0;
		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);

		int nBTries = 10;
		int verbose = 1;
		RandomStream stream = new MRG32k3a();// MRG32k3a LFSR258

		int cc = 36;
		{
			System.out.println(
					"**********************************************************************************************************************************");

			{
				SearchStrategy searchStrategy = SearchStrategy.EXTEND_TWO;
				String searchType = searchStrategy.toString();
				String wafomType = (factor == 1) ? "Matsumoto" : "Goda";

				System.out.println("          " + "   " + searchType + "  : w " + w + " k " + k + " dim " + dim
						+ " nbTries " + nBTries + " h " + h + " factor " + factor + " ");
				System.out.println();
				stream.resetStartStream();

				System.out.println();

				LMSExtendTwoV2 searcher = new LMSExtendTwoV2(w, k, dim, nBTries, table_c, q, h, factor, stream);
				long startTime = System.currentTimeMillis();
				int[] genMat = searcher.ExtendTwo(verbose);
				long endTime = System.currentTimeMillis();

				DigitalNetBase2 digitalNet = new DigitalNetBase2(k, w, dim, genMat);

//			LowWafomPointsetsSearch.writeToFile(filePathMatrices, searchStrategy, dim, k, w, matricesAsString);

				System.out.println("Duration: " + (endTime - startTime) / 1000.0 + " seconds");
				h = 1.0;
				table_c = new LookUpTableWafom(w, q, l, h, factor);
				WafomFast wafom = new WafomFast(digitalNet, w, q, table_c);
				double wafomValue = wafom.computeWafom();

				wafomType = (factor == 1) ? " Matsumoto " : " Goda";
				System.out.println();
				System.out.print(wafomType + " Wafom : " + wafomValue + " | ");
				double factor2 = (factor == 1) ? 2 : 1;

				LookUpTableWafom table_c2 = new LookUpTableWafom(w, q, w / q, h, factor2);
				wafom = new WafomFast(digitalNet, w, q, table_c2);
				wafomValue = wafom.computeWafom();
				wafomType = (factor2 == 1) ? " Matsumoto " : " Goda";
				System.out.println(wafomType + " Wafom : " + wafomValue);

				searcher.SearchInfos += String.format(" , %f seconds", (endTime - startTime) / 1000.0);

				// Save the search results to a CSV file
//			LowWafomPointsetsSearch.saveSearchResults(searcher.SearchInfos, filePathResults);
//				CompareWafom.printMatrixInBinary(searcher.bestL, w)
				Search.printJavaArrayDeclaration(searcher.bestL);
				;

				System.out.println(
						"------------------------------------------------------------------------------------------");

			}

//			
		}

	}

}
