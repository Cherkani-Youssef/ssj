package wafomExperiments;

import java.io.IOException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.LookUpTableWafom;
import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import wafomExperiments.HelperFunctions.SearchStrategy;

/**
 * We attempt to find extensible low-Wafom digital nets in both dimension and
 * number of points
 * 
 * We populate we sample N random choices for the first column of the first
 * dimension then we move to the next dimensions without ever looking back.
 * then we move to the next column
 */ 
public class LMSExtendTwo {

	private int w;
	private int k;
	int nbTries;
	private int dim;

	private RandomStream stream;
	private LookUpTableWafom table_c; // for the fast Wafom
	private int q;
	double h;
	double factor;
	String SearchInfos = " ";
	static String filePath = "/Users/cher/Downloads/new-joe-kuo-6.21201";
	int[][] bestL;

	public LMSExtendTwo(int w, int k, int dim, int nbTries, LookUpTableWafom table_c, int q, double h, double factor,
			RandomStream stream) {

		this.w = w;
		this.k = k;
		this.dim = dim;
		this.nbTries = nbTries;

		this.stream = stream;
		this.table_c = table_c;
		this.q = q;
		this.h = h;
		this.factor = factor;
		this.bestL = new int[dim][k];

	}

	public void copyBestL(int[][] tempL, int s, int d) {

		for (int i = 0; i < s; i++) {

			System.arraycopy(bestL[i], 0, tempL[i], 0, d);
		}

	}

	public int[] ExtendTwo(int verbose) throws IOException {
		final int twoRm1 = (1 << (w - 1));
		int end = 0;
//		double[] wafomValuesForDim = new double[dim];
		KahanSummation[] kahanSums = new KahanSummation[dim];
		// Initialize each element of the array
		for (int i = 0; i < kahanSums.length; i++) {
			kahanSums[i] = new KahanSummation();
		}
		double lowestWafomForDim = Double.MAX_VALUE;
//		double[] xx = new double[k];
//		int nBtirage = 0;
		for (int d = 0; d < k; d++) {
//			System.out
//					.println(" -------------------------Column  " + (d + 1) + " -------------------------------------");

			int start = (d == 0) ? 0 : (end);
			end = (1 << (d + 1));

			int[][] tempL = new int[dim][d + 1];
			copyBestL(tempL, dim, d + 1);

			int numPoints = end;
			int longueur = end - start;
			double[] wafomValues = new double[longueur];
			double lastSumTemp = 0.0;

			for (int s = 1; s <= dim; s++) {
//				System.out
//						.println(" -------------------------Dimension " + s + " -------------------------------------");

				double[] tempWafomValues = wafomValues.clone();

				lowestWafomForDim = Double.MAX_VALUE;

				for (int attemptDim = 0; attemptDim < nbTries; attemptDim++) {
					double[] wafomValues2 = tempWafomValues.clone();
					if (s == 1)
						Arrays.fill(wafomValues2, 1.0);

					DigitalNetBase2 Sob = new SobolSequence(filePath, d + 1, w, s);

					tempL[s - 1][d] = (twoRm1 + stream.nextInt(0, twoRm1 - 1)) >> (d);

					for (int j = 0; j < s; j++) {
						Sob.leftMultiplyMatBisV2(j, tempL[j], d + 1);
					}
					WafomFast wafom = new WafomFast(Sob, w, q, table_c);

					double[] currentWafomValues = wafom.computeWafom2(start, end, s - 1);

					double somme = wafom.elementWiseProd(wafomValues2, currentWafomValues);
//					double somme2 =  wafom.sumKahan(wafomValues2);
//					wafom.elementWiseProd2(wafomValues2, currentWafomValues);
//					double somme = wafom.sumKahan(wafomValues2);

					double currentWafom = (somme + kahanSums[s - 1].getSum()) / (double) numPoints;
//					double currentWafom = -1.0 + (somme + kahanSums[s - 1].getSum()) / (double) numPoints;
//					double currWafom = wafom.computeWafom();
//					System.out.println(" currWafom Original " + currWafom + " currentWafom Cumulative " + currentWafom
//							+ " error " + Math.abs(currWafom - currentWafom));
					if ((currentWafom) < (lowestWafomForDim)) {
						lowestWafomForDim = currentWafom;

						bestL[s - 1][d] = tempL[s - 1][d];
						wafomValues = wafomValues2.clone();
						lastSumTemp = somme;
//						if(s==3)
//						xx[d] = Math.log10(lowestWafomForDim);
					}
				}
				kahanSums[s - 1].add(lastSumTemp);
//				wafomValuesForDim[s - 1] += lastSumTemp;
			}

			if (verbose > 0)
				System.out.println("## Column  " + (d + 1) + " lowest Wafom so far " + lowestWafomForDim + " log10 "
						+ Math.log10(lowestWafomForDim));
//			xx[d] = Math.log10(lowestWafomForDim);
//			tries += 10;
		}

		System.out.println();
		SobolSequence Sobol = new SobolSequence(filePath, k, w, dim);
		for (int j = 0; j < dim; j++) {
			Sobol.leftMultiplyMatBisV2(j, bestL[j], k);
//			System.out.println(Arrays.toString(bestL[j]));
		}
		LookUpTableWafom table_c = new LookUpTableWafom(w, q, w / q, h, factor);
		WafomFast wafom = new WafomFast(Sobol, w, q, table_c);
		double wafomValue = wafom.computeWafom();
		System.out.println(" Wafom : " + wafomValue + " log10 " + Math.log10(wafomValue));
		System.out.println();
//		System.out.println();
//		System.out.println(Arrays.toString(xx));
//		System.out.println();
		return Sobol.getGeneratorMatricesTrans();

	}

	public static void main(String[] args) throws IOException {

		int w = 31; // when using LMS set to 30 1.8048372250945022E-5 000 000 000 000 00
		int q = 3;
		int l = w / q;
		int k = 20;
		int dim = 5;
		int nbTries = 3;

		int verbose = 1;

		double h = 1.0;
		double factor = 1.0;
		RandomStream stream = new MRG32k3a();
		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);

		{

			SearchStrategy searchStrategy = SearchStrategy.EXTEND_TWO;
			String searchType = searchStrategy.toString();
			String wafomType = (factor == 1) ? "Matsumoto" : "Goda";

			System.out.println("          " + "   " + searchType + "  : w " + w + " k " + k + " dim " + dim
					+ " nbTries " + nbTries + " h " + h + " factor " + factor + " ");
			System.out.println();

			
			LMSExtendTwo searcher = new LMSExtendTwo(w, k, dim, nbTries, table_c, q, h, factor, stream);

			long startTime = System.currentTimeMillis();
			int[] genMat = searcher.ExtendTwo(verbose);
			long endTime = System.currentTimeMillis();

			DigitalNetBase2 digitalNet = new DigitalNetBase2(k, w, dim, genMat);
			

//		LowWafomPointsetsSearch.writeToFile(filePathMatrices, searchStrategy, dim, k, w, matricesAsString);

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
//		LowWafomPointsetsSearch.saveSearchResults(searcher.SearchInfos, filePathResults);
//			CompareWafom.printMatrixInBinary(searcher.bestL, w)
			Search.printJavaArrayDeclaration(searcher.bestL);
			;

			System.out.println(
					"------------------------------------------------------------------------------------------");

		}
	}

}
