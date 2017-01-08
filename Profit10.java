package SALPMPE;

/**
 * Defines the single period profit function for the capacity competition model. 
 * @version 2.0
 */
public class Profit10 {
	/** Array of capacities associated to quality levels */static double[] qbar;
	static{
	
	// Capacity of firm in quality state k
		qbar	= new double[Parameters.K+1];
		for(int k=1;k<=Parameters.K;k++){
			qbar[k] = Parameters.qmin + (Parameters.qmax-Parameters.qmin)*(k-1)/(Parameters.K-1); 
		}
	}
	
	
	/**
	 * Defines a new instance of the profit function.
	 */
	public Profit10(){		
	}
	
	
	/**
	 * Computes industry equilibrium profits for the single-period quantity competition game. 
	 * @param s Industry state.
	 * @return  Profit vector.
	 */
	public double[] profit(int[] s){	 
		double[] demand = new double[Parameters.K+1];
		double[] pi     = new double[Parameters.K+1];
		double Tdemand	= 0;
		double Price;
		
		// Solving for quantity equilibrium
		if (s[0]>0){
			demand = fsolve(s);
			for(int k=1 ; k<=Parameters.K;k++){
				if (s[k]!=0){
					Tdemand += s[k]*demand[k];
				}
			}
			// Profit computation
			Price 			 = (Parameters.a - Tdemand)/Parameters.b;
			for(int k=1;k<=Parameters.K;k++){
				if (s[k]!=0){
					pi[k] 	 = Price*demand[k];
				}
			}
		}else {
			for(int k=1;k<=Parameters.K;k++){
				pi[k] =0;
			}
		}
		return pi;		
	}
	
	
	/**
	 * Computes industry equilibrium market shares for the single-period quantity competition game. 
	 * @param s Industry state.
	 * @return  Market share vector.
	 */
	public double[] mshare(int[] s){	 
		double[] demand = new double[Parameters.K+1];
		double[] share  = new double[Parameters.K+1];
		double Tdemand	= 0;
		
		if (s[0]>0){
			// Solving for price equilibrium
			demand = fsolve(s);
			for(int k=1 ; k<=Parameters.K;k++){
				if (s[k]!=0){
					Tdemand += s[k]*demand[k];
				}	
			}
			
			// Share computation
			for(int k=1;k<=Parameters.K;k++){
				if(s[k]!=0){
					share[k] = demand[k]/Tdemand;
				}
			}
		}else{
			for(int k=1;k<=Parameters.K;k++){
				share[k]=0;
			}
		}
		return share;		
	}
	
	
	/**
	 * Computes industry equilibrium consumer surplus for the single-period quantity competition game. 
	 * @param s Industry state.
	 * @return  Consumer surplus vector.
	 */
	public double[] surplus(int[] s){	 
		double[] demand = new double[Parameters.K+1];
		double[] surplus = new double[2];
		double Tdemand = 0;
		double Price;
		
		if (s[0]>0){
			// Solving for price equilibrium
			demand = fsolve(s);
			for(int k=1 ; k<=Parameters.K;k++){
				if (s[k]!=0){
					Tdemand += s[k]*demand[k];
				}	
			}		
			
			// Consumer surplus computation
			Price 			 = (Parameters.a - Tdemand)/Parameters.b;
			surplus[0] = (Parameters.a/Parameters.b- Price)*Tdemand/2;  	// Consumer surplus
			surplus[1] = Price*Tdemand;										// Produce surplus
		}else{
			surplus[0] = 0;
			surplus[1] = 0;
		}
		return surplus;
	}
	
	
	/**
	 * Newton method implementation to solve for equilibrium quantities for capacity competition model.	
	 * @param s    Industry state.
	 * @param Xtol Quantity tolerance. 
	 * @return     Vector of equilibrium quantities.
	 */
	public double[] fsolve(int[] s){
		// Newton method implementation for Profit function 99		
		boolean unfeas   = true;					// while condition
		double res;						  			// residual capacity
		int rem;									// remaining firms
		double[]  demand = new double[s.length];	// Equilibrium prices				
		boolean[] active = new boolean[s.length];	// auxiliary array
		double quantity;							// auxiliary variable
		
		// Getting problem dimension (n)
		for(int x=1;x<s.length;x++){
			if (s[x]>0){
				active[x] = true;
			} else{
				active[x] = false;
				demand[x]  = 0;
			}
		}
		// Setting initial quantities
		while (unfeas){
			res = Parameters.a;
			rem = 0;
			for(int x=1;x<s.length;x++){
				if (active[x]){
					rem += s[x];
				} else{
					res -= s[x]*demand[x];					
				}
			}
			quantity = res/(rem+1);
			unfeas = false;
			for(int x=1;x<s.length;x++){
				if (active[x]){
					demand[x] = quantity;
					if (quantity>qbar[x]){
						demand[x] = qbar[x];
						active[x] = false;
						unfeas = true;
					}
				}
			}
		}
		return demand;		
	}	
}
