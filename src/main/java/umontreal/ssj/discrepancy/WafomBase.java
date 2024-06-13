package umontreal.ssj.discrepancy;

import java.io.IOException;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.rng.RandomStream;
import wafomExperiments.QmcPointsFormatFileWriter;

public abstract class WafomBase {
	protected DigitalNetBase2 digitalNet;
	protected int outDigits;
	protected RandomStream stream;
	protected boolean writeFile;
	protected String filePath;
	//protected int[][] points;

	public WafomBase(DigitalNetBase2 digitalNet, int outDigits) {
		this.digitalNet = digitalNet;
		this.outDigits = outDigits;

	}

	public WafomBase(DigitalNetBase2 digitalNet, int outDigits, RandomStream stream, boolean writeFile,
			String filePath) {
		this.digitalNet = digitalNet;
		this.outDigits = outDigits;
		this.stream = stream;
		this.writeFile = writeFile;
		this.filePath = filePath;

	}



	// Abstract method for computing Wafom. Subclasses will provide the
	// implementation.
	public abstract double computeWafom() ;


	protected int compute_d_c(int x_i, int c, int w, int q) {
		int l = w / q; // Length of each full segment
		int remainder = w % q; // Length of the remainder segment, if any
		int segments = q + (remainder > 0 ? 1 : 0); // Total number of segments including remainder if exists

		// Boundary check for c
		if (c < 1 || c > segments) {
			throw new IllegalArgumentException("Segment index c is out of bounds.");
		}

		if (c <= q) {
			int startBitFromMsb = (c * l) - 1;
			int shiftAmount = w - startBitFromMsb - 1;
			return (x_i >> shiftAmount) & ((1 << l) - 1); //(u_i_j >> ((q  -1- c) * q)) & (1 << q);// 
		} else {
			// Handle the remainder segment correctly
			return x_i & ((1 << remainder) - 1); // Extract the remainder segment value
		}
	}

}
