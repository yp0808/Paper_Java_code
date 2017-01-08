package SALPMPE;

/**
 * Defines the separable approximating architecture.
 * @version 2.0
 */
public class sol_sep{
	/** weight of separable functions in approximation  */ double[][][] r_sep;
	/** intercept value in approximation  */double r_const;
	
	/**
	 * Defines a new instance of the coefficients for approximating the value function.
	 */
	public sol_sep(){
		r_sep   = new double[Parameters.K+1][Parameters.K+1][Parameters.N];
		r_const = 0;
	}
	
	/**
	 * Set the values of the coefficients associated to basis functions.
	 * @param q values for approximating coefficients.
	 */
	public void setr(sol_sep q){
		r_sep   = (double[][][]) q.r_sep.clone();
		r_const = q.r_const;
	}

	
	/**
	 * Smooth update for the coefficients in the approximating architecture. 
	 * @param v     reference level for the approximating coefficients.
	 * @param niter Iteration number
	 */
	public void updater(sol_sep v, int niter){
		
		// Apply smooth update only after a few periods
		double w = 1/Math.pow(Math.max(niter-Parameters.smoothlag,1), Parameters.smooth);			
		r_const  = w*r_const + (1-w)*v.r_const;
		for(int x=1;x<r_sep.length;x++){
			for(int k=1;k<r_sep[0].length;k++){
				for(int n=0;n<r_sep[0][0].length;n++){
					r_sep[x][k][n] = w*r_sep[x][k][n] + (1-w)*v.r_sep[x][k][n];
				}
			}
		}
	}			

}