package SALPMPE;

/**
 * Defines the representation of a state's strategy.
 * @version 2.0
 */
public class strategy {
	
	/** Investment policy */int[] policy;
	/** Exit threshold */   double[] exit;
	/** Entry threshold */  double entry;
	
	/**
	 * Creates new strategy object   
	 */
	public strategy(){
		policy = new int[Parameters.K+1];
		exit   = new double[Parameters.K+1];
		entry  = 0;
	}
	
	
	/**
	 * Sets investment policy, and entry and exit thresholds.  
	 * @param inputStrategy Strategy thresholds.
	 */
	public void setStrategy(strategy inputStrategy){
		for(int k=1;k<=Parameters.K;k++){
			policy[k] = inputStrategy.policy[k];
			exit[k]   = inputStrategy.exit[k];
		}
		entry = inputStrategy.entry;
	}
	
	
	/**
	 * Sets investment policy.  
	 * @param inputPolicy Policy per state.
	 */
	public void setPolicy(int[] inputPolicy){
		for(int k=1;k<=Parameters.K;k++){
			policy[k] = inputPolicy[k];
		}
	}
	
	
	/**
	 * Sets exit threshold.  
	 * @param inputExit Exit threshold.
	 */
	public void setExit(double[] inputExit){
		for(int k=1;k<=Parameters.K;k++){
			exit[k] = inputExit[k];
		}
	}
	
	/**
	 * Sets entry threshold.  
	 * @param inputEntry Entry threshold.
	 */
	public void setEntry(double inputEntry){
		entry = inputEntry;
	}

}
