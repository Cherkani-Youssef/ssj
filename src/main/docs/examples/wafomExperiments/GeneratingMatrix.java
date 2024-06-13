package wafomExperiments;

import java.util.Arrays;

import umontreal.ssj.rng.LFSR258;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;

/**
 * This class is designed to create and manipulate generating matrices with a
 * specific focus on ensuring that the upper-left submatrix is invertible. It
 * supports generating matrices column by column.
 */
public class GeneratingMatrix {
	private int nRows;
	private int nCols;
	private int[][] matrix;
	private int rank;
	private int[][] rrefMatrix;
	private RandomStream stream;
	private int[][] subMatrix;

	/**
	 * Constructor for creating a generating matrix with a specified number of rows
	 * and columns. If nCols is not specified or set to null, the matrix is
	 * initialized without columns, to be generated later.
	 *
	 * @param nRows  The number of rows in the matrix.
	 * @param nCols  The number of columns in the matrix. Can be null if columns are
	 *               to be generated later.
	 * @param stream The random stream used for generating random elements in the
	 *               matrix.
	 */

	public GeneratingMatrix(int nRows, Integer nCols, RandomStream stream) {
		if (nCols != null && nCols > nRows) {
			throw new IllegalArgumentException("Number of columns must be less than number of nRows.");
		}
		this.nRows = nRows;
		this.nCols = 0; // Initial number of columns
		this.stream = stream;
		this.matrix = new int[nRows][0]; // Initially empty matrix
		this.rank = 0;
		this.rrefMatrix = new int[nRows][0]; // Initially empty matrix
		this.subMatrix = new int[0][]; // Initialize the subMatrix as empty

		if (nCols != null) {
			generateGeneratingMatrix(nCols,stream);
		}
	}

	/**
	 * Constructor for creating a generating matrix with a specified number of rows.
	 * Columns columns added when calling generate generateOneColumn()
	 *
	 * @param nRows  The number of rows in the matrix.
	 * @param stream The random stream used for generating random elements in the
	 *               matrix.
	 */
	public GeneratingMatrix(int nRows, RandomStream stream) {

		this.nRows = nRows;
		this.nCols = 0; // Initial number of columns
		this.stream = stream;
		this.matrix = new int[nRows][0]; // Initially empty matrix
		this.rank = 0;
		this.rrefMatrix = new int[nRows][0]; // Initially empty matrix
		this.subMatrix = new int[0][]; // Initialize the subMatrix as empty

	}

	/**
	 * Copy constructor for creating a deep copy of another GeneratingMatrix.
	 *
	 * @param original The original GeneratingMatrix to copy.
	 */
	public GeneratingMatrix(GeneratingMatrix original) {
		this.nRows = original.nRows;
		this.nCols = original.nCols;
		this.matrix = CopyMatrix(original.matrix);
		this.rank = original.rank;
		this.rrefMatrix = CopyMatrix(original.rrefMatrix);
		this.stream = original.stream; // Assuming the stream should be shared; if not, you'll need a clone or
										// equivalent for the stream as well.
		this.subMatrix = CopyMatrix(original.subMatrix);
	}

	// Utility method for deep copying a matrix
	private static int[][] CopyMatrix(int[][] original) {
		if (original == null)
			return null;
		int[][] copy = new int[original.length][];
		for (int i = 0; i < original.length; i++) {
			copy[i] = Arrays.copyOf(original[i], original[i].length);
		}
		return copy;
	}

	/**
	 * Generates an array of GeneratingMatrix objects.
	 *
	 * @param s      The size of the array to generate in other words dimension.
	 * @param nRows  The number of rows in each matrix.
	 * @param nCols  The number of columns in each matrix.
	 * @param stream The random stream used for generating random elements in the
	 *               matrices.
	 * @return An array of GeneratingMatrix objects.
	 */
	public static GeneratingMatrix[] generateMatrixArray(int s, int nRows, int nCols, RandomStream stream) {
		GeneratingMatrix[] matrices = new GeneratingMatrix[s];
		// RandomStream stream = new MRG32k3a();
		for (int i = 0; i < s; i++) {

			matrices[i] = new GeneratingMatrix(nRows, nCols, stream);
		}
		return matrices;
	}

	/**
	 * Adds a new column to each matrix in an array of GeneratingMatrix objects.
	 *
	 * @param matrices The array of GeneratingMatrix objects to modify.
	 */
	public static void addNewColumnToEachMatrix(GeneratingMatrix[] matrices, RandomStream stream) {
		for (GeneratingMatrix matrix : matrices) {
			matrix.generateOneColumn(stream);
		}
	}
	/**
     * Creates a deep copy of an array of GeneratingMatrix objects.
     * 
     * @param originalArray The original array to copy.
     * @return A new array containing deep copies of the original GeneratingMatrix objects.
     */
	public static GeneratingMatrix[] copyGeneratingMatrixArray(GeneratingMatrix[] originalArray) {
		if (originalArray == null)
			return null;
		GeneratingMatrix[] copyArray = new GeneratingMatrix[originalArray.length];
		for (int i = 0; i < originalArray.length; i++) {
			if (originalArray[i] != null) {
				copyArray[i] = new GeneratingMatrix(originalArray[i]);
			} else {
				copyArray[i] = null; // Or however you wish to handle null entries
			}
		}
		return copyArray;
	}

	public static int[] generatorMatricesFromStandardFormat(GeneratingMatrix[] generatingMatrices) {
		int dim = generatingMatrices.length;
		if (dim == 0)
			return new int[0]; // Return an empty array if there are no matrices

		int numRows = generatingMatrices[0].getnRows();
		int numCols = generatingMatrices[0].getnCols();

		int[] genMat = new int[dim * numCols];
		int r, c, j; // Row r, column c, dimension j.
		for (j = 0; j < dim; ++j) {
			int[][] matrix = generatingMatrices[j].getMatrix();
			for (r = 0; r < numRows; ++r) {
				for (c = 0; c < numCols; ++c) {
					if (matrix[r][c] > 0)
						genMat[j * numCols + c] += (1 << (numRows - 1 - r));
				}
			}
		}
		return genMat;
	}

	// Getter for the number of nRows
	public int getnRows() {
		return this.nRows;
	}

	// Getter for the number of columns
	public int getnCols() {
		return this.nCols;
	}

	// Getter for the matrix
	public int[][] getMatrix() {
		return this.matrix;
	}

	// Getter for the rank
	public int getRank() {
		return this.rank;
	}

	// Getter for the RREF matrix
	public int[][] getRrefMatrix() {
		return this.rrefMatrix;
	}

	// Getter for the subMatrix
	public int[][] getSubMatrix() {
		return this.subMatrix;
	}
	/**
     * Converts a given matrix to its Reduced Row Echelon Form (RREF).
     * 
     * @param matrix The matrix to be converted.
     * @return The RREF of the given matrix.
     */

	public static int[][] rref(int[][] matrix) {
		int m = matrix.length;
		if (m == 0)
			return new int[0][0];
		int n = matrix[0].length;

		int[][] M = new int[m][n];
		for (int i = 0; i < m; i++) {
			M[i] = Arrays.copyOf(matrix[i], n);
		}

		int i = 0, j = 0;
		while (i < m && j < n) {
			int maxIndex = i;
			for (int k = i + 1; k < m; k++) {
				if (Math.abs(M[k][j]) > Math.abs(M[maxIndex][j])) {
					maxIndex = k;
				}
			}

			// Swap nRows
			int[] temp = M[i];
			M[i] = M[maxIndex];
			M[maxIndex] = temp;

			for (int row = 0; row < m; row++) {
				if (row != i) {
					int scalar = M[row][j];
					for (int col = j; col < n; col++) {
						M[row][col] ^= M[i][col] * scalar;
					}
				}
			}
			i++;
			j++;
		}

		return M;
	}

	private boolean isLinearlyIndependent(int[] newColumn) {
		if (this.matrix[0].length == 0) {
			newColumn[0] = 1;
			for (int i : newColumn) {
				if (i != 0)
					return true; // Any non-zero column is independent if the matrix is empty
			}
			return false;
		}

		// Combine and compute RREF to check for linear independence.

		int[][] combinedMatrix = (addColumn(this.matrix, newColumn));
		
		int k = this.nCols + 1;
		
		int[][] upperLeftSubMatrix = upperLeftSubMatrix(combinedMatrix, k, k);



		int[][] upperLeftSubMatriRRef = rref(upperLeftSubMatrix);
		int rankAfter = computeRank(upperLeftSubMatriRRef);

		return rankAfter == this.nCols + 1;
	}

	private void updateMatrix(int[] newColumn) {
		this.matrix = addColumn(this.matrix, newColumn);
		this.nCols += 1;
		this.subMatrix = upperLeftSubMatrix(this.matrix, nCols, nCols);
		this.rrefMatrix = rref(this.subMatrix);
		this.rank = computeRank(this.rrefMatrix);
	}

	private void generateGeneratingMatrix(int m, RandomStream stream) {
		if (m > this.nRows) {
			throw new IllegalArgumentException("Number of columns must be less than number of nRows.");
		}
		for (int i = 0; i < m; i++) {
			generateOneColumn(stream);
		}
	}
	/**
     * Generates a new column and adds it to the matrix if it is linearly independent.
     */
	public void generateOneColumn(RandomStream stream) {

		if (this.nCols + 1 > this.nRows) {
			throw new IllegalArgumentException("Number of columns must be less than number of nRows.");
		}
		int[] newColumn = new int[this.nRows];
		for (int i = 0; i < this.nRows; i++) {
			newColumn[i] = stream.nextInt(0, 1); // (stream.nextDouble() < 0.51) ? 1 : 0; //
		}
		// Check if Invertibel
		if (isLinearlyIndependent(newColumn)) {
			updateMatrix(newColumn);

		} else {// if not flip the k-th bit of the column generated
			newColumn[this.nCols] ^= 1;
			updateMatrix(newColumn);
		}
	}

	/**
     * Adds a column to the matrix.
     * 
     * @param matrix The matrix to which the column is to be added.
     * @param newColumn The new column to be added.
     * @return The new matrix with the added column.
     */

	private int[][] addColumn(int[][] matrix, int[] newColumn) {

		int[][] newMatrix = new int[matrix.length][matrix[0].length + 1];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				newMatrix[i][j] = matrix[i][j];
			}
			newMatrix[i][matrix[0].length] = newColumn[i];
		}
		return newMatrix;
	}

	/**
     * Adds a new row to the subMatrix of this GeneratingMatrix.
     * 
     * @param newRow The new row to be added.
     */
	public void addRow(int[] newRow) {
		// Create a new matrix with one more row than subMatrix
		int[][] newSubMatrix = new int[this.subMatrix.length + 1][];

		// Copy existing nRows from subMatrix to newSubMatrix
		for (int i = 0; i < this.subMatrix.length; i++) {
			newSubMatrix[i] = this.subMatrix[i];
		}

		// Add the new row at the end
		newSubMatrix[this.subMatrix.length] = newRow;

		// Replace the old subMatrix with newSubMatrix
		this.subMatrix = newSubMatrix;
	}

	// returns the upper Left SubMatrix
	public static int[][] upperLeftSubMatrix(int[][] matrix, int nnRows, int nnCols) {
		// Ensure the input dimensions are valid for the given matrix.
		if (matrix == null || matrix.length < nnRows || matrix[0].length < nnCols) {
			throw new IllegalArgumentException("Invalid dimensions for the given matrix.");
		}

		// Initialize the submatrix with the specified dimensions.
		int[][] subMatrix = new int[nnRows][nnCols];

		// Copy the relevant portion of the original matrix to the submatrix.
		for (int i = 0; i < nnRows; i++) {
			for (int j = 0; j < nnCols; j++) {
				subMatrix[i][j] = matrix[i][j];
			}
		}

		return subMatrix;
	}

	// Utility method to print a matrix
	static void printMatrix(int[][] matrix) {
		for (int[] row : matrix) {
			for (int element : row) {
				System.out.print(element + " ");
			}
			System.out.println();
		}
	}

	// Helper method to count non-zero nRows in a matrix.
	private static int computeRank(int[][] matrix) {
		int count = 0; // Initialize counter for non-zero nRows
		for (int[] row : matrix) { // Iterate through each row of the matrix
			boolean isNonZeroRow = false; // Flag to check if the current row is non-zero
			for (int element : row) { // Check each element in the row
				if (element != 0) { // If any element is non-zero, mark the row as non-zero
					isNonZeroRow = true;
					break; // Break the inner loop as we found a non-zero element
				}
			}
			if (isNonZeroRow) { // If the row is non-zero, increment the count
				count++;
			}
		}
		return count; // Return the total count of non-zero nRows
	}
	/**
     * Modifies the elements of the lower rows of the matrix based on a given probability.
     * 
     * @param probability The probability with which each element is modified.
     * @param fromcol The starting column index from which modifications begin.
     * @param stream The random stream used for generating the probability of modification.
     */
	public void modifyLowerRowsWithProbability(double probability, int fromcol, RandomStream stream) {
		// Check if the matrix is square (m = n)
		if (this.nRows == this.nCols) {
			throw new IllegalArgumentException("No changes allowed! The matrix must be invertible.");
		}

		// Iterate over the last (m - n) rows
		for (int i = this.nCols; i < this.nRows; i++) {
			for (int j = fromcol; j < this.nCols; j++) {
				// Flip the element based on the probability
				if (stream.nextDouble() < probability) {
					this.matrix[i][j] ^= 1; // Using XOR to flip the bit (assuming binary elements)
				}
			}
		}
	}
	/**
     * Applies the modification of lower rows with a given probability to an array of GeneratingMatrix objects.
     * 
     * @param matrices The array of GeneratingMatrix objects to be modified.
     * @param probability The probability with which each element is modified.
     * @param fromcol The starting column index from which modifications begin.
     * @param stream The random stream used for generating the probability of modification.
     */
	public static void modifyLowerRowsWithProbabilityForArray(GeneratingMatrix[] matrices, double probability,
			int fromcol, RandomStream stream) {
		for (GeneratingMatrix matrix : matrices) {
			// Apply the modification on each matrix in the array
			matrix.modifyLowerRowsWithProbability(probability, fromcol, stream);
		}
	}

}