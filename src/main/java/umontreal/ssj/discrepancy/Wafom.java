package umontreal.ssj.discrepancy;

import java.io.IOException;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.PointSetIterator;
import wafomExperiments.KahanSummation;

public class Wafom extends WafomBase {
	private double h;
	private double factor;
	private double c = 0.0;

	/**
	 * 
	 * 
	 * @param digitalNet The Digital Net Base 2 instance.
	 * @param outDigits  The precision of the points.
	 * @param h          The parameter defining Matsumoto's or Yoshiki's definition.
	 * @param factor     set to to 1 for wafom set to 2 to get Wafom for RMSE
	 */
	public Wafom(DigitalNetBase2 digitalNet, int outDigits, double h, double factor) {
		super(digitalNet, outDigits);
		this.h = h;
		this.factor = factor;
	}

	@Override
	public double computeWafom() throws IOException {

		double prod = 1.0, somme = 0.0;
		long numPoints = digitalNet.getNumPoints();

		PointSetIterator iter = digitalNet.iterator();
		KahanSummation kahanSum = new KahanSummation();
		while (iter.hasNextPoint()) {
			prod = 1.0;

			while (iter.hasNextCoordinate()) {

				int u_i_j = iter.nextInt(); // added this in pointSetiterator to return the the coodinate in int format
				for (int l = 1; l <= outDigits; l++) {
					int u_i_j_l = ((u_i_j >> (outDigits - l)) & 1);

//					int exponent = (int)(factor * (l  + h));
					double two_exponent = Math.pow(2.0, -(factor * (l + h)));
					prod *= (1 + (1 - 2 * u_i_j_l) * two_exponent); // Math.pow(2.0, -exponent)); ( 1.0 / (double) (1L
																	// << exponent))
				}
			}
			KahanSummation kahanSum2 = new KahanSummation();
			kahanSum2.add(prod);
			kahanSum2.add(-1.0);

			kahanSum.add(kahanSum2.getSum());
			iter.resetToNextPoint();
		}
		return((kahanSum.getSum() / (double) numPoints));
	}

	public double computeWafom2() throws IOException {
		double logProd = 0.0, somme = 0.0;
		long numPoints = digitalNet.getNumPoints();
		double mOne = -1.0; 

		PointSetIterator iter = digitalNet.iterator();
		KahanSummation kahanSum = new KahanSummation();
		while (iter.hasNextPoint()) {
			logProd = 0.0; // Reset logProd for each new point
			KahanSummation kahanSumProd = new KahanSummation();
			while (iter.hasNextCoordinate()) {
				int u_i_j = iter.nextInt(); // Extract the coordinate in int format
				for (int l = 1; l <= outDigits; l++) {
					int u_i_j_l = ((u_i_j >> (outDigits - l)) & 1);
					double two_exponent = Math.pow(2.0, -(factor * (l + h)));
					double term = 1 + (1 - 2 * u_i_j_l) * two_exponent;

					logProd += Math.log(term); // Summing logs instead of multiplying
//					kahanSumProd.add( Math.log(term));

				}
			}
			// Convert logProd back to linear scale using exponential
			double prod = Math.exp(logProd);

			KahanSummation kahanSum2 = new KahanSummation();
			kahanSum2.add(prod);
			kahanSum2.add(mOne);

			kahanSum.add(kahanSum2.getSum());
			iter.resetToNextPoint();
		}
		return((kahanSum.getSum() / (double) numPoints));
	}

	

}
