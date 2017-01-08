package SALPMPE;

/**
 * Defines the single period profit function for the quality ladder competition model. 
 * @version 2.0
 */
public class Profit99 {
	
	/** Array of quality levels */static double[] gm;
	static{
	
	// Quality dependent part of the adjusted mean utility for the profit demand model
		gm	= new double[Parameters.K+1];
		for(int k=1;k<=Parameters.K;k++){
			gm[k] = Math.pow((k-1)/Parameters.psi +1 , Parameters.theta1); 
		}
	}
	
	
	/**
	 * Defines a new instance of the profit function.
	 */
	public Profit99(){		
	}
	
	
	/**
	 * Computes industry equilibrium profits for the single-period price competition game. 
	 * @param s Industry state.
	 * @return  Profit vector.
	 */
	public double[] profit(int[] s){	 
		double[] price  = new double[Parameters.K+1];
		double[] demand = new double[Parameters.K+1];
		double[] pi     = new double[Parameters.K+1];
		double Nsum	    = 1;
		
		// Solving for price equilibrium
		if (s[0]>0){
			price = fsolve(s,1e-6);
			for(int k=1 ; k<=Parameters.K;k++){
				demand[k]=0;
				if (s[k]!=0){
					demand[k] = gm[k]*Math.pow(Parameters.Y-price[k],Parameters.theta2);
					Nsum     += s[k]*demand[k];
				}	
			}
		
			// Profit computation
			for(int k=1;k<=Parameters.K;k++){	
				pi[k] = Parameters.m*demand[k]/Nsum*(price[k]-Parameters.c);
			}
		}else {
			for(int k=1;k<=Parameters.K;k++){
				pi[k] =0;
			}
		}
		return pi;		
	}
		
	
	/**
	 * Computes industry equilibrium market shares for the single-period price competition game. 
	 * @param s Industry state.
	 * @return  Market share vector.
	 */
	public double[] mshare(int[] s){	 
		double[] price  = new double[Parameters.K+1];
		double[] demand = new double[Parameters.K+1];
		double[] share  = new double[Parameters.K+1];
		double Nsum	    = 0;			// Set =1 for report considering outside good
		
		if (s[0]>0){
			// Solving for price equilibrium
			price = fsolve(s,1e-6);
			for(int k=1 ; k<=Parameters.K;k++){
				demand[k]=0;
				if (s[k]!=0){
					demand[k] = gm[k]*Math.pow(Parameters.Y-price[k],Parameters.theta2);
					Nsum     += s[k]*demand[k];
				}	
			}
			
			// Share computation
			for(int k=1;k<=Parameters.K;k++){
				share[k] = demand[k]/Nsum;
			}
		}else{
			for(int k=1;k<=Parameters.K;k++){
				share[k]=0;
			}
		}
		return share;		
	}
	
	
	/**
	 * Computes industry equilibrium consumer surplus for the single-period price competition game. 
	 * @param s Industry state.
	 * @return  Consumer surplus vector.
	 */
	public double[] surplus(int[] s){	 
		double[] price  = new double[Parameters.K+1];
		double[] demand = new double[Parameters.K+1];
		double[] surplus = new double[2];
		double Nsum	    = 1;
		
		if (s[0]>0){
			// Solving for price equilibrium
			price = fsolve(s,1e-6);
			for(int k=1 ; k<=Parameters.K;k++){
				demand[k]=0;
				if (s[k]!=0){
					demand[k] = gm[k]*Math.pow(Parameters.Y-price[k],Parameters.theta2);
					Nsum     += s[k]*demand[k];
				}
			}
			surplus[0] = Parameters.m*Math.log(Nsum); 		// Consumer surplus
			surplus[1] = 0;									// Produce surplus
			for(int k=1 ; k<=Parameters.K;k++){
				if (s[k]!=0){ 
					surplus[1]+= s[k]*Parameters.m*demand[k]/Nsum*(price[k]-Parameters.c);
				}
			}
		}else{
			surplus[0] = 0;
			surplus[1] = 0;
		}
		return surplus;
	}
	

	/**
	 * Newton method implementation to solve for equilibrium prices for quality ladder competition model.	
	 * @param s    Industry state.
	 * @param Xtol Price tolerance. 
	 * @return     Vector of equilibrium prices.
	 */	
	public double[] fsolve(int[] s, double Xtol){
		// Newton method implementation for Profit function 99							
		int n          = 0;						// Problem dimension
		int k 		   = 0;						// Iteration counter
		int aux        = 0;						// auxiliary variable		
		double[] price = new double[s.length];	// Equilibrium prices
		double dX 	   = Xtol + 1.0;			// Step metric		
		double[] X;								// Unknown		
		double[][] Jac;							// Jacobian matrix
		int[] indice;							// auxiliary array
		int[] indeX;							// auxiliary array
		double[] fi;							// auxiliary variables for computing Jacobian
		double den;								// auxiliary variable  for computing Jacobian
		double[] rhs;							// right hand side on Newton step equation
		
		// Getting problem dimension (n)
		for(int x=1;x<s.length;x++){if (s[x]>0){ n++;}}		
		X      = new double[n];
		fi     = new double[n];
		rhs    = new double[n];
		indice = new int[n];
		indeX  = new int[n];
		Jac    = new double[n][n];
		
		// Setting initial prices
		for(int x=1;x<s.length;x++){
			if (s[x]>0){
				indice[aux] = x;
				X[aux]      = Parameters.Y/2; 								
				price[x]    =0;
				aux++;				
			} else{ price[x] = 0;
			}
		}
		
		// Iterate until meet convergence criteria or reaching max iterations
		while ((dX > Xtol) && (k<Parameters.NiterNew)){   								
			den = 1.0;
			for(int i=0;i<n;i++){
				fi[i] = gm[indice[i]]*Math.pow(Parameters.Y-X[i],Parameters.theta2);
				den  += s[indice[i]]*fi[i];
			}
			
			// First order conditions	
			for(int i=0;i<n;i++){
				fi[i]  = fi[i]/den;
				rhs[i] = X[i]-Parameters.Y + Parameters.theta2*(X[i]-Parameters.c)*(1-fi[i]);			
			}
			
			// Computing Jacobian matrix
			for(int i=0;i<n;i++){
				for(int j=0;j<n;j++){
					if (j==i){
						if (X[i]!=Parameters.Y){
							Jac[i][j] = -1 + Parameters.theta2*(fi[i]-1)*(1+((X[i]-Parameters.c)/(Parameters.Y-X[i]))*Parameters.theta2*fi[i]);
						} else{
							Jac[i][j] = -1 - Parameters.theta2;
						}
					} else{
						if (X[j]!=Parameters.Y){
							Jac[i][j] = ((X[i]-Parameters.c)/(Parameters.Y-X[j]))*Parameters.theta2*Parameters.theta2*fi[i]*fi[j];
						} else{
							Jac[i][j] = 0;
						}
					}					
				}
			}
			
			// Computing new prices
			double d[] = solve(Jac, rhs, indeX);
			dX = 0;
		      for (int i=0; i<n; i++) {
		        dX   += d[i]*d[i];
		        X[i] += d[i];
		      }
		      dX = Math.sqrt(dX/n);
		      k++;						
		}
		
		// Recovering equilibrium prices from Newton routine
		for (int i=0; i<n; i++) {
			price[indice[i]] = X[i];		
		}
		if (k == Parameters.NiterNew){ System.out.println("Equilibrium price not found for state S");	
		}
		return price;		
	}
	
	/**
	 *  Method to solve the equation a[][] x[] = b[] with
	 *  the partial-pivoting Gaussian elimination scheme.
	 *  The code is part of the book, "An Introduction to        
	 *  Computational Physics, 2nd Edition," written by Tao Pang and      
	 *  published by Cambridge University Press on January 19, 2006.
	 *  @param a Left hand side matrix.
	 *  @param b Right hand side vector.
	 *  @return x Solution to linear system.
	 */
	  public static double[] solve(double a[][], double b[],
	    int index[]) {
	    int n = b.length;
	    double x[] = new double[n];

	 // Invoke the partial-pivoting Gaussian elimination
	    gaussian(a, index);

	 // Rescale array b[i] with the ratios stored
	    for(int i=0; i<n-1; ++i) {
	      for(int j =i+1; j<n; ++j) {
	        b[index[j]] -= a[index[j]][i]*b[index[i]];
	      }
	    }

	 // Perform the backward substitutions
	    x[n-1] = b[index[n-1]]/a[index[n-1]][n-1];
	    for (int i=n-2; i>=0; --i) {
	      x[i] = b[index[i]];
	      for (int j=i+1; j<n; ++j) {
	        x[i] -= a[index[i]][j]*x[j];
	      }
	      x[i] /= a[index[i]][i];
	    }
	    return x;
	  }

	  
	/**
	 * Method to perform the partial-pivoting Gaussian
	 * elimination.  
	 * @param index Records the pivoting order.
	 * @param a Left hand side matrix.
	 */
	  public static void gaussian(double a[][], int index[]) {
	    int n = index.length;
	    double c[] = new double[n];

	 // Initialize the index
	    for (int i=0; i<n; ++i) index[i] = i;

	 // Find the rescaling factors, one from each row
	    for (int i=0; i<n; ++i) {
	      double c1 = 0;
	      for (int j=0; j<n; ++j) {
	        double c0 = Math.abs(a[i][j]);
	        if (c0 > c1) c1 = c0;
	      }
	      c[i] = c1;
	    }

	 // Search the pivoting element from each column
	    int k = 0;
	    for (int j=0; j<n-1; ++j) {
	      double pi1 = 0;
	      for (int i=j; i<n; ++i) {
	        double pi0 = Math.abs(a[index[i]][j]);
	        pi0 /= c[index[i]];
	        if (pi0 > pi1) {
	          pi1 = pi0;
	          k = i;
	        }
	      }

	   // Interchange rows according to the pivoting order
	      int itmp = index[j];
	      index[j] = index[k];
	      index[k] = itmp;
	      for (int i=j+1; i<n; ++i) {
	        double pj = a[index[i]][j]/a[index[j]][j];

	     // Record pivoting ratios below the diagonal
	        a[index[i]][j] = pj;

	     // Modify other elements accordingly
	        for (int l=j+1; l<n; ++l)
	          a[index[i]][l] -= pj*a[index[j]][l];
	      }
	    }
	  }	

}
