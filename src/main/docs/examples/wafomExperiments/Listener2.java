package wafomExperiments;


public class Listener2 {
	// Attributes to store the matrix and the merit value
	private int[][] leftMatrixScramble;
	private double meritValue;
	private double somme;

	// Constructor
	public Listener2() {
		// Initialize meritValue with the maximum possible value so that any new value
		// will be lower
		meritValue = Double.MAX_VALUE;
		setSomme(Double.MAX_VALUE);
		leftMatrixScramble = null;
	}

	public void update(int[][] newMatrixScramble, double newMeritValue) {
	    // Check if the new merit value is lower than the current merit value
	    if (newMeritValue < meritValue) {
	        meritValue = newMeritValue;
	        // Perform a deep copy of the matrix
	        leftMatrixScramble = copyMatrix(newMatrixScramble);
	    }
	}
	
	public void update(int[][] newMatrixScramble, double newMeritValue, double sum) {
	    // Check if the new merit value is lower than the current merit value
	    if (newMeritValue < meritValue) {
	        meritValue = newMeritValue;
	        setSomme(sum);
	        somme  = sum;
	        // Perform a deep copy of the matrix
	        leftMatrixScramble = copyMatrix(newMatrixScramble);
	    }
	}

	// Helper method to perform a deep copy of a 2D array
	private int[][] copyMatrix(int[][] matrix) {
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


	// Getter for the merit value
	public double getMeritValue() {
		return meritValue;
	}

	// Getter for the matrix scramble
	public int[][] getLeftMatrixScramble() {
		return leftMatrixScramble;
	}

	public double getSomme() {
		return somme;
	}

	public void setSomme(double somme) {
		this.somme = somme;
	}
}