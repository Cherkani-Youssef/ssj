package wafomExperiments;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.LookUpTableWafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.RandomStream;


public class HelperFunctions {
	
	public enum SearchStrategy {
//		NAIVE_NIED_XING,  COLUMN_BY_COLUMN, NAIVE_EXPLICIT, CBC_NIED_XING, 
		LMS_FIX, EXTEND_DIM, EXTEND_SIZE, EXTEND_TWO
		// Add more strategies as needed
	}

	/**
	 * Writes the specified parameters and generator matrices to a file in the same
	 * format as in latnetbuilder.
	 * 
	 * @param filePath The path of the file to write to.
	 * @param strategy The search strategy used.
	 * @param s        The dimensionality (s).
	 * @param k        The value of k, where n = 2^k.
	 * @param r        The number of digits (r).
	 * @param matrices The generator matrices, where each row represents a matrix.
	 * @throws IOException If an error occurs during file writing.
	 */
	public static void writeToFile(String filePath, SearchStrategy strategy, int s, int k, int r,
			String matricesAsString) throws IOException {
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath, true))) {
			// Write the header information
			writer.write("# dnet\n");
			writer.write("# " + strategy.toString() + "\n");
			writer.write("# basis b = 2\n");
			writer.write(s + " # s = " + s + " dimensions\n");
			writer.write(k + " # k = " + k + ", so n = 2^" + k + " = " + (int) Math.pow(2, k) + " points\n");
			writer.write(r + " # r = " + r + " digits\n");
			writer.write("# The columns of gen. matrices C_1, ..., C_s, one matrix per line:\n");

			// Write the generator matrices
			writer.write(matricesAsString);
			writer.write("___________________________\n");
		}
	}
	
	/**
	 * Generates one LMS where each cell is a column
	 * @param stream
	 * @param numCols
	 * @param outDigits
	 * @return
	 */

	public static int[] GenerateleftMatrixScrambleBis(RandomStream stream, int numCols, int outDigits) {
		int c;
		int[] L = new int[numCols];
		int twomr1 = 1 << (outDigits - 1);

		for (c = 0; c < numCols; c++) {

			L[c] = (twomr1 + stream.nextInt(0, twomr1 - 1)) >> c;

		}
		return L;
	}

	
}
