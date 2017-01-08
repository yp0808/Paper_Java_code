package SALPMPE;


/**
 * Defines the single period profit function. 
 * @version 2.0
 */
public class Profit {
	/** Profit function associated to quality ladder competition model */Profit99 prof99;
	/** Profit function associated to capacity competition model */Profit10 prof10;


	/**
	 * Defines a new instance of the profit function.
	 */
	public Profit(){
		prof99 = new Profit99();
		prof10 = new Profit10();
	}

	
	/**
	 * Computes industry equilibrium profits for the single-period game. 
	 * @param s Industry state.
	 * @return  Profit vector.
	 */
	public double[] profit(int[] s){
	     if(Parameters.useprofit99){return prof99.profit(s);}
	     else if(Parameters.useprofit10){return prof10.profit(s);}
	     else{return null;}
	}

	/**
	 * Computes industry equilibrium market shares for the single-period game. 
	 * @param s Industry state.
	 * @return  Market share vector.
	 */
	public double[] mshare(int[] s){
	     if(Parameters.useprofit99){return prof99.mshare(s);}
	     else if(Parameters.useprofit10){return prof10.mshare(s);}
	     else{return null;}
	}

	
	/**
	 * Computes industry equilibrium consumer surplus for the single-period game. 
	 * @param s Industry state.
	 * @return  Consumer surplus vector.
	 */
	public double[] surplus(int[] s){
	     if(Parameters.useprofit99){return prof99.surplus(s);}
	     else if(Parameters.useprofit10){return prof10.surplus(s);}
	     else{return null;}
	}	
}
