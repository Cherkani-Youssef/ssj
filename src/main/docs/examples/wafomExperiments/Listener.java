
package wafomExperiments;

import java.util.AbstractMap.SimpleEntry;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.RandomStream;

public class Listener {

	int[][][] leftMatrixScrambles;
//	private int[][][] leftMatrixScramblesTemp;
	private int populationSize;
	int currentIndex;
	double[] meritValues;
	int s, k;
	int[] indices;

	// Constructor
	public Listener(int s, int k, int populationSize) {
		this.leftMatrixScrambles = new int[populationSize][s][k];
		this.meritValues = new double[populationSize];
		this.indices = new int[populationSize];
//		Arrays.fill(meritValues, Double.MAX_VALUE);
		this.currentIndex = 0;
		this.populationSize = populationSize;
		this.s = s;
		this.k = k;

	}

	// Constructor
	public Listener(int[][][] SeqL, int s, int populationSize) {
		this.leftMatrixScrambles = new int[populationSize][SeqL[0].length][SeqL[0][0].length + 1];
		this.meritValues = new double[populationSize];
//		Arrays.fill(meritValues, Double.MAX_VALUE);
		this.currentIndex = 0;
		this.populationSize = populationSize;

		for (int i = 0; i < populationSize; i++) {
			for (int j = 0; j < s; j++) {

				for (int k = 0; k < SeqL[i][j].length; k++) {
					// Copy each value from SeqL to the corresponding place in leftMatrixScrambles
					leftMatrixScrambles[i][j][k] = SeqL[i][j][k];
				}

			}
		}
	}


	public Listener(int[][][] SeqL, int s, int populationSize, int[] indices) {
		this.leftMatrixScrambles = new int[populationSize][s][SeqL[0][0].length + 1];
		this.meritValues = new double[populationSize];
//		Arrays.fill(meritValues, Double.MAX_VALUE);
		this.currentIndex = 0;
		this.populationSize = populationSize;

		int i = 0;
//		for (int i = 0; i < populationSize; i++) 
		for (int ii : indices) {
			for (int j = 0; j < s; j++) {

				for (int k = 0; k < SeqL[i][j].length; k++) {
					// Copy each value from SeqL to the corresponding place in leftMatrixScrambles
					leftMatrixScrambles[i][j][k] = SeqL[ii][j][k];
				}

			}
			i++;
		}
	}
	
	public void add(int[][] newMatrixScramble) {
		leftMatrixScrambles[currentIndex++] = copyMatrix(newMatrixScramble);

	}


	public int[][] getLeftMatrixScramble(int i) {
		if (i >= populationSize)
			throw new IllegalArgumentException("Index out of bounds ");

		return copyMatrix(leftMatrixScrambles[i]);

	}

//	// Getter for the merit value
	public double[] getMeritValues() {
		return meritValues.clone();
	}

	public int[][][] getPopulation() {
		if (leftMatrixScrambles == null) {
			return null; // Or you could return an empty array depending on your use case
		}

		// Create a new 3D array with the same dimensions
		int[][][] copy = new int[leftMatrixScrambles.length][][];

		for (int i = 0; i < leftMatrixScrambles.length; i++) {
			if (leftMatrixScrambles[i] == null) {
				copy[i] = null;
				continue;
			}
			copy[i] = new int[leftMatrixScrambles[i].length][];
			for (int j = 0; j < leftMatrixScrambles[i].length; j++) {
				if (leftMatrixScrambles[i][j] == null) {
					copy[i][j] = null;
					continue;
				}
				copy[i][j] = new int[leftMatrixScrambles[i][j].length];
				System.arraycopy(leftMatrixScrambles[i][j], 0, copy[i][j], 0, leftMatrixScrambles[i][j].length);
			}
		}

		return copy;
	}

	public void resetMeritValuesArray() {

		Arrays.fill(meritValues, Double.MAX_VALUE);

	}


	public void update(int[][] newMatrixScramble, double newMeritValue, int currentIndex) {

		if ((newMeritValue < meritValues[currentIndex])) {
			meritValues[currentIndex] = newMeritValue;
			leftMatrixScrambles[currentIndex] = copyMatrix(newMatrixScramble);

		}

	}

	public void add(int[][] newMatrixScramble, double newMeritValue, int attempt) {

		if (attempt < populationSize) {

			meritValues[currentIndex] = newMeritValue;
			leftMatrixScrambles[currentIndex] = copyMatrix(newMatrixScramble);
			currentIndex++;
			return;
		}
		int index = findMax(meritValues);
		if (index != -1) {

			if ((meritValues[index] > newMeritValue)) {
				meritValues[index] = newMeritValue;
				leftMatrixScrambles[index] = copyMatrix(newMatrixScramble);
			}
		}

	}

	public void update(int[] newMatrixScramble, int dimIndex, double newMeritValue) {

		int index = findMax(meritValues);

		if (index != -1) {
			leftMatrixScrambles[index][dimIndex] = newMatrixScramble;
			currentIndex++;
		}
	}

	public int[][] copyMatrix(int[][] matrix) {
		if (matrix == null) {
			return null;
		}
		int[][] newMatrix = new int[matrix.length][];
		for (int i = 0; i < matrix.length; i++) {
			newMatrix[i] = new int[matrix[i].length];
			System.arraycopy(matrix[i], 0, newMatrix[i], 0, matrix[i].length);
		}
		return newMatrix;
	}

	public static double[] findMin(double[] array) {
		if (array == null || array.length == 0) {
			throw new IllegalArgumentException("Array must not be empty");
		}
		int index = 0;
		double min = array[0];
		for (int i = 1; i < array.length; i++) {
			if (array[i] < min) {
				min = array[i];
				index = i;
			}
		}
		return new double[] { min, index };
	}

	public static int findMax(double[] array) {
		if (array == null || array.length == 0) {
			throw new IllegalArgumentException("Array must not be empty");
		}
		int index = 0; // Start with the first element
		double max = array[0];
		for (int i = 1; i < array.length; i++) { // Start loop from the second element
			if (array[i] > max) {
				max = array[i];
				index = i;
			}
		}
		return index;
	}

}
