package wafomExperiments;

import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.IntStream;

import umontreal.ssj.discrepancy.LookUpTableWafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.RandomStream;

public abstract class Search {
	protected String filePath = "/Users/cher/Downloads/new-joe-kuo-6.21201";
	protected int w;
	protected int k;
	protected int dim;
	protected RandomStream stream;
	protected LookUpTableWafom table_c; // for the fast Wafom
	protected int q;
	protected double h;
	protected double factor;
	protected String SearchInfos = " ";
	protected int[][] bestL;
	protected int indexbestL;
	protected double[] lowestWafomValues;
	protected int[][][] bestLs;
	protected int populationSize;
	protected static double[] lmsFixedMatsu;

	public Search(int w, int k, int dim, int populationSize, LookUpTableWafom table_c, int q, double h, double factor,
			RandomStream stream) {
		this.w = w;
		this.k = k;
		this.dim = dim;

		this.stream = stream;
		this.table_c = table_c;
		this.q = q;
		this.h = h;
		this.factor = factor;
		this.bestL = new int[dim][k];
		this.populationSize = populationSize;
		this.lowestWafomValues = new double[populationSize];
		this.bestLs = new int[populationSize][dim][k];
		lmsFixedMatsu = WafomValues.dataMatsu.get(dim);
		lmsFixedGoda =  WafomValues.dataGoda.get(dim);
	}

//	public Search(int w, int k, int dim, int populationSize, LookUpTableWafom table_c, int q,
//			double h, double factor, RandomStream stream) {
//		this(w, k, dim, populationSize, table_c, q, h, factor, stream);
//		
//	}

	public abstract void run(int verbose, int end, double exponent, double slowExponent);

	
	public int[] selectBest2(int populationSize, int n, int[][][] bestLs, int dim, int k, boolean flag) {
		double[] fixed = (factor == 1) ? lmsFixedMatsu : lmsFixedGoda;
		System.out.println(Arrays.toString(fixed));
//		double[] fixed = (factor == 1) ? data.get(dim) : dataGoda.get(dim);
		double[][] wafomValues = new double[populationSize][k];
		double[] distances = new double[populationSize];
		double[] distances2 = new double[populationSize];
		double[] distances3 = new double[populationSize];
		SobolSequence[] digitalNet = new SobolSequence[k];
		for (int d = 0; d < k; d++) {
			digitalNet[d] = new SobolSequence(filePath, d + 1, w, dim);
		}
		// Compute WAFOM values
		for (int i = 0; i < populationSize; i++) {
			double[] wafomValues2 = new double[k];
			for (int d = 0; d < k; d++) {
				digitalNet[d].resetGeneratorMatrices();
//				SobolSequence Sob = new SobolSequence(filePath, d + 1, w, dim);
				for (int j = 0; j < dim; j++) {
					digitalNet[d].leftMultiplyMatBisV2(j, bestLs[i][j], d + 1);
				}
				WafomFast wafom = new WafomFast(digitalNet[d], w, q, table_c);
				double currentWafom = wafom.computeWafom();
//				wafomValues2[d] = Math.log10(currentWafom);
				wafomValues2[d] = Math.log10(Math.abs(currentWafom));
			}
			wafomValues[i] = wafomValues2;
		}

		// Calculate Euclidean distances between each sample and the fixed array
		for (int i = 0; i < populationSize; i++) {
			double sum = 0;
			double sum2 = 0;
			double sum3 = 0;
			for (int j = 0; j < k; j++) {
				double bound = -((1. * (j + 1) * (j + 1)) / (dim * 1.));
//				sum += Math.abs((wafomValues[i][j]) - (bound));
				sum += Math.pow((wafomValues[i][j]) - (bound), 2.0);
				sum2 += Math.pow((wafomValues[i][j] - fixed[j]), 2.0);
				sum3 += Math.abs((wafomValues[i][j]) - (fixed[j]));
//				sum3 += Math.pow((wafomValues[i][j]) - (fixed[j]), 2.0);
			}
//			distances[i] = (sum);
			distances[i] = Math.sqrt(sum);
			distances2[i] = Math.sqrt(sum2);
			distances3[i] = (sum3);
//			distances3[i] = Math.sqrt(sum3);
		}

		Integer[] indices = IntStream.range(0, populationSize).boxed().toArray(Integer[]::new);
//		Arrays.sort(indices, Comparator.comparingDouble(i -> distances[i]));// matsu
		if (factor == 1.0)
			Arrays.sort(indices, Comparator.comparingDouble(i -> distances[i]));
		else
			Arrays.sort(indices, Comparator.comparingDouble(i -> distances3[i]));

//
//		// Select the top n indices
		int[] bestIndices = Arrays.stream(indices).limit(n).mapToInt(Integer::intValue).toArray();

		if (flag) {
//			System.out.print("Best distance " + distances2[bestIndices[0]] + " index " + bestIndices[0]+ "  |");
////			System.out.println();// Arrays.toString(wafomValues[bestIndices[0]]) +
//			System.out.println(" Worst  " + distances2[bestIndices[bestIndices.length - 1]] + " index "
//					+ (bestIndices.length - 1));

			String ANSI_RED = "\u001B[31m";
			String ANSI_RESET = "\u001B[0m";
			double goal = 3.3;
			// Print best distance
			if (distances2[bestIndices[0]] < goal) {
				System.out.print(ANSI_RED + "Best distance " + distances2[bestIndices[0]] + " index " + bestIndices[0]
						+ ANSI_RESET + "  |");
			} else {
				System.out.print("Best distance " + distances2[bestIndices[0]] + " index " + bestIndices[0] + "  |");
			}

			// Print worst distance
			if (distances2[bestIndices[bestIndices.length - 1]] < goal) {
				System.out.println(ANSI_RED + " Worst  " + distances2[bestIndices[bestIndices.length - 1]] + " index "
						+ bestIndices.length + ANSI_RESET);
			} else {
				System.out.println(
						" Worst  " + distances2[bestIndices[bestIndices.length - 1]] + " index " + bestIndices.length);
			}
//			System.out.println();
		}

		// If indexbestL is already in bestIndices, return bestIndices directly
		return bestIndices;

	}

	public int[] selectBest(int populationSize, int n, int[][][] bestLs, int dim, int k, boolean flag) {
		double[] fixed = (factor == 1) ? lmsFixedMatsu : lmsFixedGoda;
		double[][] wafomValues = new double[populationSize][k];
		double[] distances = new double[populationSize];
		double[] distances2 = new double[populationSize];
		double[] distances3 = new double[populationSize];
	    SobolSequence[] digitalNet = new SobolSequence[k];
	    for (int d = 0; d < k; d++) {
	        digitalNet[d] = new SobolSequence(filePath, d + 1, w, dim);
	    }
	   
		// Compute WAFOM values in parallel
		IntStream.range(0, populationSize).parallel().forEach(i -> {
			double[] wafomValues2 = new double[k];
			 int end=0;
			 double wafVal=0;
			for (int d = 0; d < k; d++) {
				int start = (d == 0) ? 0 : (end);
				end = (1 << (d + 1));
				int numPoints =end;
//				SobolSequence digitalNet = new SobolSequence(filePath, d + 1, w, dim);
				DigitalNetBase2 dn = new DigitalNetBase2( d + 1, w, dim, digitalNet[d].getGeneratorMatricesTrans());
//	            digitalNet[d].resetGeneratorMatrices();
				for (int j = 0; j < dim; j++) {
					dn.leftMultiplyMatBisV2(j, bestLs[i][j], d + 1);
				}
				
				WafomFast wafom = new WafomFast(dn, w, q, table_c);
//				double currentWafom2 = wafom.computeWafom();
//				wafomValues2[d] = Math.log10((currentWafom));
				
				double somme = wafom.computeWafomSum(start, end) + wafVal;// / end
				double currentWafom = (somme) / numPoints;
				wafomValues2[d] = Math.log10((currentWafom));
				
				wafVal =somme;
//				System.out.println("currentWafom2 " + currentWafom + " currentWafom2 " + currentWafom2 + "  abs " + Math.abs(currentWafom2- currentWafom));
			}
			wafomValues[i] = wafomValues2;
		});

		// Calculate Euclidean distances between each sample and the fixed array in
		// parallel
		IntStream.range(0, populationSize).parallel().forEach(i -> {

			double sum = 0;
			double sum2 = 0;
			double sum3 = 0;
			for (int j = 0; j < k; j++) {
				double bound = -((1. * (j + 1.) * (j + 1.)) / (dim * 1.));
				sum += Math.abs((wafomValues[i][j]) - (bound));
//				sum += Math.pow((wafomValues[i][j]) - (bound), 2.0);
				sum2 += Math.pow((wafomValues[i][j] - fixed[j]), 2.0);
				sum3 += Math.abs((wafomValues[i][j]) - (fixed[j]));

			}
			distances[i] = (sum);
//			distances[i] = Math.sqrt(sum);
			distances2[i] = Math.sqrt(sum2);
			distances3[i] = (sum3);
		});

		Integer[] indices = IntStream.range(0, populationSize).boxed().toArray(Integer[]::new);
		if (factor == 1.0)
			Arrays.sort(indices, Comparator.comparingDouble(i -> distances[i]));
		else
			Arrays.sort(indices, Comparator.comparingDouble(i -> distances[i]));

		int[] bestIndices = Arrays.stream(indices).limit(n).mapToInt(Integer::intValue).toArray();

		if (flag) {


			String ANSI_RED = "\u001B[31m";
			String ANSI_RESET = "\u001B[0m";
			double goal = (factor == 1) ? 1 : 3;
			// Print best distance
			if (distances2[bestIndices[0]] < goal) {
				System.out.print(ANSI_RED + "Best distance " + distances2[bestIndices[0]] + " index " + bestIndices[0]
						+ ANSI_RESET + "  |");
			} else {
				System.out.print("Best distance " + distances2[bestIndices[0]] + " index " + bestIndices[0] + "  |");
			}

			// Print worst distance
			if (distances2[bestIndices[bestIndices.length - 1]] < goal) {
				System.out.println(ANSI_RED + " Worst  " + distances2[bestIndices[bestIndices.length - 1]] + " index "
						+ bestIndices.length + ANSI_RESET);
			} else {
				System.out.println(
						" Worst  " + distances2[bestIndices[bestIndices.length - 1]] + " index " + bestIndices.length);
			}
//			System.out.println();
		}

		return bestIndices;
	}

	public double[] computesWafoms(int i) {
//		int end = 0;

		double[] wafomValues = new double[k];
		int end =0;
		double wafVal =0.0;
		for (int d = 1; d <= k; d++) {
			
			int start = (d == 1) ? 0 : (end);
			end = (1 << (d ));
			int numPoints =end;
			DigitalNetBase2 Sob = new SobolSequence(filePath, d, w, dim);
			for (int j = 0; j < dim; j++) {
				Sob.leftMultiplyMatBisV2(j, bestLs[i][j], d);
			}
			WafomFast wafom = new WafomFast(Sob, w, q, table_c);
//			double currentWafom = wafom.computeWafom();
//			wafomValues[d - 1] = Math.log10(currentWafom);
			double somme = wafom.computeWafomSum(start, end) + wafVal;// / end
			double currentWafom = (somme) / numPoints;
			wafomValues[d-1] = Math.log10((currentWafom));
			
			wafVal = somme;

		}

		return wafomValues;
	}

	public double[][] calcWafomPopulation(double[] fixed, int[] indices) {

		double[][] wafomValues = new double[indices.length + 1][k];
		
//			for (int i = 0; i < bestLs.length; i++) 
		int j = 0;
		for (int i : indices) {
			wafomValues[j] = computesWafoms(i);
			j++;
		}
		wafomValues[indices.length] = fixed;
		return wafomValues;
	}

	public static int[] nonlinearInterpolationSlow(double xStart, double xEnd, int numPoints, double exponent) {
		int[] xx = new int[numPoints];
		double h = (xStart - xEnd);
		for (int i = 0; i < numPoints; i++) {
			double t = (double) i / (numPoints - 1);
			xx[i] = (int) (xStart - h * Math.pow(t, exponent));
		}
		return xx;
	}

	public static int[] nonlinearInterpolation(double xStart, double xEnd, int numPoints) {
		int[] xx = new int[numPoints];
		double ratio = Math.pow(xEnd / xStart, 1.0 / (numPoints - 1));
		for (int i = 0; i < numPoints; i++) {
			xx[i] = (int) (xStart * Math.pow(ratio, i));
		}
		return xx;
	}

	public static int[] nonlinearInterpolationVariableSmooth(double xStart, double xEnd, int numPoints,
			double fastExponent, double slowExponent, int switchPoint) {
		int[] xx = new int[numPoints];
		double transitionValue = xStart
				- (xStart - xEnd) * Math.pow((switchPoint - 1) / (double) (numPoints - 1), fastExponent);
		for (int i = 0; i < numPoints; i++) {
			double t = i / (double) (numPoints - 1);
			if (i < switchPoint) {
				xx[i] = (int) (xStart - (xStart - xEnd) * Math.pow(t, slowExponent));
			} else {
				double tAdjusted = (i - switchPoint + 1) / (double) (numPoints - switchPoint);
				xx[i] = (int) (transitionValue - (transitionValue - xEnd) * Math.pow(tAdjusted, fastExponent));
			}
		}
		return xx;
	}

	public static void printJavaArrayDeclaration(int[][] array) {
		StringBuilder javaArrayStr = new StringBuilder("static int[][] data = {\n");
		for (int[] row : array) {
			javaArrayStr.append("    {");
			for (int i = 0; i < row.length; i++) {
				javaArrayStr.append(row[i]);
				if (i < row.length - 1) {
					javaArrayStr.append(", ");
				}
			}
			javaArrayStr.append("},\n");
		}
		// Remove the last comma and newline
		javaArrayStr.setLength(javaArrayStr.length() - 2);
		javaArrayStr.append("\n};");
		System.out.println(javaArrayStr.toString());
	}

	
	static Map<Integer, double[]> dataGoda = new HashMap<>();
	static {
		dataGoda.put(6, new double[] { -0.93765655, -1.43556268, -1.71077301, -2.27080655, -3.58914217, -4.15003855,
				-4.98318009, -5.70503002, -6.5526066, -7.37683742, -8.13608116, -9.08483064, -9.99496139, -10.95205464,
				-11.91803976, -12.98240962, -13.99330764, -15.0234905, -15.9454432, -17.04319817 });

		dataGoda.put(5,
				new double[] { -1.1111838329815749, -1.6120324532949362, -1.9827540689756675, -3.0254303782844363,
						-4.003218981139232, -4.655708520091044, -5.525882242442834, -6.388975564294801,
						-7.3280367577973955, -8.370933613641608, -9.224427120434658, -10.246049260793168,
						-11.32492854974176, -12.301234358825228, -13.44072249732897, -14.64328223571436,
						-15.744254624458907, -16.886146954950906, -18.175022693269263, -18.54382591931164 });
//				new double[] { -1.1112776855699573, -1.6326283939910866, -1.991473907379913, -3.0347766120312207,
//						-4.036271044595612, -4.722051563515058, -5.7169030884427405, -6.604192919749763,
//						-7.515154629949739, -8.455937478773434, -9.369540996447157, -10.417459100484201,
//						-11.358104475653453, -12.451835684083779, -13.742416476675452, -14.699743799116213,
//						-15.781795683397236, -17.208866084120004, -18.331145190655878, -18.54382591931164 });

		dataGoda.put(4,
				new double[] { -1.3199777, -1.80205477, -2.30077231, -3.5897983, -4.57259739, -5.50064598, -6.44583793,
						-7.35764631, -8.58601552, -9.57518864, -10.61012525, -11.80910192, -12.93224841, -14.2690328,
						-15.6597198, -16.87470886, -18.58354898, -18.18476377, -18.0752598, -17.73576507 });

		dataGoda.put(3,
				new double[] { -1.59725034, -2.21772021, -3.42439805, -4.4861429, -5.49465778, -6.66634812, -7.80705545,
						-8.92481843, -10.46572406, -11.63359286, -12.88642675, -14.37922036, -16.05167211, -17.36829263,
						-17.61090141, -17.13649168, -17.74555903, -17.47970085, -17.65049577, -17.8249246 });

		dataGoda.put(2,
				new double[] { -2.01205825, -3.10931631, -4.42572, -5.90046345, -7.33801485, -8.67932797, -10.40083486,
						-11.99409168, -13.56307921, -15.23185304, -17.05107587, -16.59105858, -16.77037538, -16.7057141,
						-16.74469989, -17.08982846, -17.5157972, -17.13987968, -17.65471105, -17.37593996 });

		dataGoda.put(1, new double[] { -2.85733052, -4.9746888, -7.22538946, -9.53800845, -11.94727416, -15.05576178,
				-16.85767976, -16.54991638, -16.52084493, -16.52240996, -16.53433096, -16.55453434, -16.69479345,
				-16.70229034, -16.72146409, -17.06786433, -17.52395606, -17.72613815, -18.19501644, -17.38741139 });
	}

	// Print the map (optional)

	

	protected static double[] lmsFixedGoda;
	protected static Map<Integer, double[]> data = new HashMap<>();
	static {
		data.put(1,
				new double[] { -1.0771947333189962, -2.0376234312757305, -3.1374219789278285, -4.369954202502529,
						-5.541681512217073, -7.034071558988664, -8.417357206706108, -9.920208615155643,
						-11.574360869990832, -13.040480765437148, -14.04320817331316, -16.46713783372254,
						-17.58778157314261, -16.846742675486773, -16.99267686779873, -15.541768232264218,
						-17.09651007391606, -17.16205996441246, -15.555365609174643, -18.388014576312174 });
		data.put(2,
				new double[] { -0.3657805213190567, -0.9095280242258541, -1.5749064392283134, -2.2693310160455162,
						-3.044663150281978, -3.7546831062898143, -4.5402147875973045, -5.39278758540023,
						-6.287238404846141, -7.13219632029019, -8.048377561101754, -9.162019967283102,
						-10.298239857311296, -11.336817745559815, -12.488459237243768, -13.558600775378727,
						-14.764140942954903, -16.088101509498724, -17.336760333334748, -18.72185167814452 });

		data.put(3,
				new double[] { 0.04529018646613196, -0.39417490942650985, -0.8925371461806126, -1.4239279599925485,
						-1.962126081113105, -2.533593876518451, -3.175983080778369, -3.7142281536751853,
						-4.359238619406787, -5.132962616711819, -5.780372272288737, -6.524572234572103,
						-7.378797469916112, -8.19249642895425, -8.83727327852428, -9.665866360659445,
						-10.524828251605168, -11.781288479296101, -12.344965758297958, -13.267235556448446 });
		data.put(4,
				new double[] { 0.3530389995872257, -0.03257939688643641, -0.4211039626712915, -0.8915495022594677,
						-1.3602599117551963, -1.847870082538669, -2.339243213353186, -2.8495981383881412,
						-3.4444757165486215, -3.937631699395851, -4.577038061519791, -5.126896544600654,
						-5.752344942637669, -6.3827790186698286, -7.093529131126328, -7.771607519506904,
						-8.389623633209816, -9.160678299637256, -9.960464953559208, -10.587276328048178 });

		data.put(5,
				new double[] { 0.613828734064777, 0.2560450908235985, -0.09010192226414528, -0.5064484643293925,
						-0.9299864816872362, -1.3434424519854284, -1.7982806648204737, -2.236462994276002,
						-2.730066204798715, -3.256202558402463, -3.7465453426349766, -4.225576353037069,
						-4.75467995065147, -5.306310654154074, -5.878611987018624, -6.449929319780994,
						-7.059955651365665, -7.7177367577377245, -8.332833533822935, -8.90124313407699 });
//				new double[] { 0.6137993501692909, 0.25521497838575036, -0.0984302950420468, -0.5064484643293925,
//						-0.9425345442726503, -1.356940887693153, -1.8189633921014587, -2.267077031611006,
//						-2.758023630462673, -3.256202558402463, -3.7465453426349766, -4.2638237271395925,
//						-4.782013292628121, -5.352782585305871, -5.959097895944288, -6.489135591573625,
//						-7.114080441035466, -7.759763709285999, -8.341333787043553, -8.975662107905057 });
		data.put(6,
				new double[] { 0.8505095178058377, 0.5117559165621622, 0.18443646063127495, -0.17792101623009046,
						-0.5863103442479798, -0.962384564392268, -1.3569183979996995, -1.777556594735675,
						-2.199000238289333, -2.687432566451734, -3.0876299349842884, -3.53590495122952,
						-4.020930942433731, -4.5326627490637605, -4.991108310125201, -5.50797990813703,
						-6.063708846101859, -6.6031601236331765, -7.163815052660158, -7.682028345109087 });

		data.put(7,
				new double[] { 1.0731319369908072, 0.7474415108035919, 0.42373682079901537, 0.09119804731103145,
						-0.28912936747446116, -0.6403936171879462, -1.0265624670413083, -1.4069129132837888,
						-1.7881736816150389, -2.2029596866137218, -2.602440397335511, -3.027309248661043,
						-3.4995158899376415, -3.9087203907393797, -4.37422293727433, -4.834644549163492,
						-5.302169484491636, -5.8173218067706305, -6.3354407527795, -6.777799816648749 });

		data.put(8,
				new double[] { 1.2874470762860533, 0.970241189174632, 0.6507808193252246, 0.3316616536597601,
						-0.012122579640062383, -0.3559953090382231, -0.7187574606199624, -1.0782020722486927,
						-1.4475958435523526, -1.8274224729200077, -2.2000481850921028, -2.599144004066138,
						-3.0254803368835965, -3.40691602084168, -3.8466563775135243, -4.265787214878424,
						-4.678701685521014, -5.158148545476161, -5.655435438836999, -6.0541174435668275 });

		data.put(9,
				new double[] { 1.4967754950449323, 1.1851020924078188, 0.8719867756274303, 0.561376227428838,
						0.22936574287065625, -0.10621421508405146, -0.4440286838170533, -0.7861846273184335,
						-1.1363692028489003, -1.5031832078937863, -1.8595966713048127, -2.2377484922840454,
						-2.612913136154551, -2.9791667229195102, -3.41934225813117, -3.7882880366424883,
						-4.21563776850349, -4.616863189681715, -5.061335182785344, -5.4902076428222735 });

		data.put(10,
				new double[] { 1.7030639420166074, 1.3952171266759978, 1.0856451797418736, 0.7757636338301119,
						0.45402149179957785, 0.13105596451922547, -0.19825053805141762, -0.5233520139050871,
						-0.8584571514120322, -1.214914318747674, -1.557315010825674, -1.9200845087112144,
						-2.2778773335488576, -2.635438126049935, -3.005223209129495, -3.393229715959083,
						-3.789277609521868, -4.176964645555049, -4.600208270145954, -4.997140924166988 });

		data.put(11,
				new double[] { 1.9074782901611853, 1.6020084777578631, 1.2952326650661348, 0.9891927672502602,
						0.6742582115253081, 0.35491628208156095, 0.034535101230651455, -0.2827235280505435,
						-0.6062381932146864, -0.9518451332769117, -1.2778490796803326, -1.625477640569048,
						-1.969983120102105, -2.321952604071321, -2.681655138151792, -3.0450654075502217,
						-3.4130725526010255, -3.7953902462897497, -4.177096382239556, -4.550697055416168 });

		data.put(12,
				new double[] { 2.1107120264533012, 1.8068720776496967, 1.5021742372682874, 1.1974912227685868,
						0.8881060282544806, 0.5725056304631957, 0.25739871030037337, -0.05574746468742,
						-0.3716030895744198, -0.6990873624948508, -1.0154882671060248, -1.3615127823505513,
						-1.69768908383784, -2.0419935831013487, -2.38136904824972, -2.7367762284692336,
						-3.087535347625153, -3.4346332262459853, -3.8147147611241348, -4.200892473477781 });

		data.put(13,
				new double[] { 2.313208979490564, 2.010366241982424, 1.706899384743287, 1.4030746895147428,
						1.0958682081945592, 0.7847760756944349, 0.47546004918286744, 0.1652149755422599,
						-0.14724972310272036, -0.46256618549318346, -0.7829438529922964, -1.1151436329900064,
						-1.4465883373989168, -1.775735704518679, -2.105431344827432, -2.457397779284936,
						-2.8039732518129594, -3.1378709566438023, -3.483236587477671, -3.8583646036210673 });

		data.put(14,
				new double[] { 2.515246674059675, 2.213040213122262, 1.9104158217907756, 1.6072457064499273,
						1.3023235333117273, 0.9941910158669789, 0.6868002464311127, 0.38143973575688306,
						0.07168512169223772, -0.24177046779014094, -0.5562238494386257, -0.8822780674223969,
						-1.2008059718863957, -1.5255949707738619, -1.8526948482421401, -2.196129610080445,
						-2.51918111476966, -2.8512889691789867, -3.1809865185793638, -3.5376606983629575 });

		data.put(15,
				new double[] { 2.7169942997883805, 2.4152137335266795, 2.1131285708016754, 1.810640633370065,
						1.5069996278992077, 1.2013565369730337, 0.8961344954632372, 0.5908232855419872,
						0.28560394892133106, -0.024644909772507717, -0.3362852070634835, -0.6523004297268812,
						-0.9675291917797572, -1.2888083481626125, -1.6054545695925038, -1.9363258850255034,
						-2.263640011003441, -2.584695816867404, -2.915098433319791, -3.2570559905547967 });

		data.put(16,
				new double[] { 2.9185599916017995, 2.6170469846940576, 2.3153237736305754, 2.0133993610691077,
						1.7107152080351853, 1.4064613472494925, 1.1024420566576032, 0.7992012149067684,
						0.4946111840014805, 0.1865650948486217, -0.12000384541939174, -0.43116828450757316,
						-0.7428389566804021, -1.057018998373035, -1.3683221376958794, -1.694330483752613,
						-2.0148043123605133, -2.333680740275048, -2.657812004462853, -2.987074685484651 });

	}

}
