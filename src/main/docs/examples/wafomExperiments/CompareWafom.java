package wafomExperiments;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import umontreal.ssj.discrepancy.LookUpTableWafom;

import umontreal.ssj.discrepancy.WafomFast;

import umontreal.ssj.hups.SobolSequence;





// Helper class to compare the wafom values of the new methods with LMS-Fixed
public class CompareWafom {
	static String filePath = "/Users/cher/Downloads/new-joe-kuo-6.21201";

	

	

	public static void main(String[] args) throws IOException {
//		int r = 10;
		int dim = 5;
		int[][][] LMS = {  LMSExtendStorage.LMSMatsu.get(dim),LMSExtendStorage.LMSGoda.get(dim)};
		double[][] wafomValues = new double[LMS.length + 1][];

		double[] distances = new double[LMS.length];
		double[] distances2 = new double[LMS.length]; 
		int k = 20;
//		int dim = 15;
		int w = 31;
		int q = 3;
		int l = w / q;

		double factor = 2.0;
		double[] fixed = (factor == 1) ? WafomValues.dataMatsu.get(dim) : WafomValues.dataGoda.get(dim);
		wafomValues[LMS.length] =  fixed;
		
//		wafomValues[LMS.length + 1] = lmsFixedMatsu;
						
		double h = 1.0;
		LookUpTableWafom table_c = new LookUpTableWafom(w, q, l, h, factor);
		for (int i = 0; i < LMS.length; i++) {

			double[] wafomValues2 = new double[k];
			for (int d = 0; d < k; d++) {
				SobolSequence Sob = new SobolSequence(filePath, d + 1, 31, dim);
				for (int j = 0; j < dim; j++) {
					Sob.leftMultiplyMatBisV2(j, LMS[i][j], d + 1);
				}
				WafomFast wafom = new WafomFast(Sob, w, q, table_c);
				double currentWafom = wafom.computeWafom();
				wafomValues2[d] = Math.log10(Math.abs(currentWafom));
			}
			wafomValues[i] = wafomValues2;

		}
		printJavaArrayDeclaration(wafomValues);

		String chartTitle = " Dimension " + dim;
		String xAxisLabel = " number of Points ";
		String yAxisLabel = " log10 Wafom ";
//	double[][] wafomValues = { wafomValuesFix, wafomValuesTwo };

		String[] labels = new String[wafomValues.length + 1]; // Create a String array of the same length

		labels[0] = "Matsu";
		labels[1] = "Matsu";
		labels[wafomValues.length] = "FIX";

//		ChartHelper chartHelper = new ChartHelper(chartTitle, xAxisLabel, yAxisLabel, wafomValues, labels, 1);
//		chartHelper.displayChart(); [502.8596695762276, 502.6849436213789][405.34588863501506, 402.4530766056343]
		double sum3 = 0;
		double sum4 = 0;
		for (int i = 0; i < LMS.length; i++) {
			double sum = 0;
			double sum2=0;
			 sum3=0;
			sum4=0;
//		
			for (int j = 0; j < k; j++) {

				sum += Math.pow((wafomValues[i][j] - fixed[j]), 2);
				double bound = -((1. * (j + 1.) * (j + 1.)) / (dim * 1.));
				sum2 += Math.abs((wafomValues[i][j] - bound));
				sum3 += Math.abs((fixed[j] - bound));
				sum4 += Math.pow((fixed[j] - bound),2);
			}
			distances2[i] = (sum2);
			distances[i] = Math.sqrt(sum);
		}
		System.out.println();
		System.out.println(Arrays.toString(distances));
		System.out.println(Arrays.toString(distances2));
		System.out.println(" abs " + sum3 + " pow " + Math.sqrt(sum4));

	}

	public static void printJavaArrayDeclaration(double[][] array) {
		StringBuilder javaArrayStr = new StringBuilder("int[][] data = {\n");
		for (double[] row : array) {
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

}
