//package umontreal.ssj.discrepancy;
//
//public class LookUpTableWafom {
//
//	private final int q;
//	private final int l;
//	private final double[] lookupTable; // Use a primitive array
//	private int h; // set to 1 to get Yoshiki's bouds
//	private int factor;// set to 2 to get Walsh figure of merit for root mean square error
//
//	/**
//	 * @param q a positive integer that divides the output digits
//	 * @param l length of the segment
//	 * @param h for Matsumoto's definition set to 0, and for Yoshiki's set to 1
//	 * @param factor set to to 1 for wafom set to 2 to get Wafom for RMSE
//	 */
// 
//	public LookUpTableWafom(int q, int l, int h, int factor) {
//		if (h != 0 && h != 1)  {
//			throw new IllegalArgumentException("h must be either 0 or 1 : for Matsumoto's definition set to 0, and for Yoshiki's set to 1");
//		}
//    	if (factor != 1 && factor != 2) {
//			throw new IllegalArgumentException("factor must be set to 1 for wafom, and set to 2 for Wafom for RMSE");
//		}
//		
//		this.q = q;
//		this.l = l;
//		this.lookupTable = new double[(int) (q * (1L << l))]; // Allocate array
//		this.h = h;
//		this.factor = factor;
//
//		generate();
//	}
//
//	public double get(int c, int e) {
//
//		if (c == 0 || c > q)
//			throw new IndexOutOfBoundsException("c is out of range");
//		if (e > (1L << l))
//			throw new IndexOutOfBoundsException("e is out of range");
//		return lookupTable[(int) index(c, e)];
//
//	}
//
//	private int index(int c, int e) {
//		return (c - 1) * (1L << l) + e;
//	}
//
//	private void generate() {
//		for (int c = 1; c <= q; ++c) {
//			for (int e = 0; e < (1L << l); ++e) {
//				double product = 1.0;
//				for (int j = 1; j <= l; ++j) {
//					int e_j = (int) ((e >> (l - j)) & 1); // Extract the j-th bit of e
//					int exponenent = (int) (factor * ( (c - 1) * l +  (j + h)  )) ;
//					product *= (1 + (1 - 2 * e_j) * (1.0 / (1L << (exponenent))) );
//				}
//				lookupTable[(int) index((int) c, e)] = product; // Set value directly in the array
//			}
//		}
//	}
//}
package umontreal.ssj.discrepancy;

public class LookUpTableWafom {

    private final int q;
    private final int l;
    private final double[] lookupTable; // Use a primitive array for efficiency
    private double h; // set to 1 to get Yoshiki's bounds
    private double factor; // set to 2 to get Walsh figure of merit for root mean square error
    private int n; // The number of bits

    /**
     * Constructor for LookUpTableWafom.
     *
     * @param n      the number of bits
     * @param q      a positive integer that divides n into l segements
     * @param l      length of the segment
     * @param h      for Matsumoto's definition set to 0, and for Yoshiki's set to 1
     * @param factor set to 1 for wafom, set to 2 to get Wafom for RMSE
     */
    public LookUpTableWafom(int n, int q, int l, double h, double factor) {
        if (h != 0 && h != 1) {
            throw new IllegalArgumentException("h must be either 0 or 1: for Matsumoto's definition set to 0, and for Yoshiki's set to 1");
        }
        if (factor != 1 && factor != 2) {
            throw new IllegalArgumentException("factor must be set to 1 for wafom, and set to 2 for Wafom for RMSE");
        }
        if (n < q) {
            throw new IllegalArgumentException("n must be larger than or equal to q");
        }
 
        this.n = n;
        this.q = q;
        this.l = l;
        // Calculate size needed for the lookup table, including the remainder if it exists.
        long tableSize = q * (1L << l) + (n % q != 0 ? (1L << (n % q)) : 0);
        this.lookupTable = new double[(int) tableSize];
        this.h = h;
        this.factor = factor;

        generate();
    }

    public double get(int c, int e) {
        if (c == 0 || c > q + (n % q != 0 ? 1 : 0))
            throw new IndexOutOfBoundsException("c is out of range");
        if (e >= (1L << l) || (c == q + 1 && e >= (1L << (n % q))))
            throw new IndexOutOfBoundsException("e is out of range");
        return lookupTable[(int) index(c, e)];
    }

    private int index(int m, int e) {
        if (m <= q) {
            return (int) ((m - 1) * (1L << l) + e);
        } else {
            // Index for the remainder segment
            return (int) (q * (1L << l) + e);
        }
    }


    private void generate() {
    	int length= 1<<l;
        // Compute products for the full segments
        for (int c = 1; c <= q; ++c) {
            for (int e = 0; e < length; ++e) {
                lookupTable[ index( c, e)] = computeProduct(c, e, l);
            }
        }
        // Compute product for the remainder segment, if there is a remainder
        if (n % q != 0) {
            int remainderLength = n % q;
            for (int e = 0; e < (1L << remainderLength); ++e) {
                lookupTable[ index(q + 1, e)] = computeProduct(q + 1, e, remainderLength);
            }
        }
    }

    private double computeProduct(int c, int e, int length) {
        double product = 1.0;
        for (int j = 1; j <= length; ++j) {
            int e_j = (int) ((e >> (length - j)) & 1);
            
            double two_exponent =Math.pow(2.0, -factor * ((c - 1) * l  + (j +  h)) ); //(int) ( factor *  ((c - 1) * l + (j + h)) );
            product *= (1 + (1 - 2 * e_j) *two_exponent  ); // * Math.pow(2, -two_exponent ));// ( 1.0 / (double) (1L << two_exponent)  )
        }
        return product;
    }
}
