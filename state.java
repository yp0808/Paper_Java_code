package SALPMPE;


/**
 * Defines the representation of the industry state and associated operations.
 * @version 2.0
 */
public class state {

	/**State frequency relative to the sampled set*/ double frequency; 
	/**Profit vector associated to the state */double[] profit;					
	/**Strategy for each firm outside industry */strategy Strat;
	/**Industry state */int[] index;
	

	/**
	 * Defines a new instance of an industry state.
	 * @param inputindex Industry state.
	 */
	public state(int[] inputindex){
		
		// Creates state using input industry state
		index     = (int[]) inputindex.clone();
		Strat     = new strategy() ;
		profit    = new double[Parameters.K+1];
		frequency = 1;
	}	
	
	
	/**
	 * Assigns profit vector to an industry state.
	 * @param inputprofit Profit vector.
	 */
	public void setProfit(double[] inputprofit){
		
		// Sets profits using input profit vector
		for(int k=1;k<=Parameters.K;k++){
			profit[k] = inputprofit[k];
		}
	}
	
	
	/**
	 * Increases the frequency of the state by one. 
	 */
	public void increaseFrequency(){
		
		// Increase state frequency when sampled
		frequency = frequency + 1;
	}
	
	
	/**
	 * Increases the frequency of the state.
	 * @param newfreq Frequency increment.
	 */
	public void setFrequency(double newfreq){
		
		// Increase state frequency when sampled
		frequency = newfreq;
	}
	
	
	/**
	 * Assigns strategy to the industry state.
	 * @param inputpolicy Strategy vector.
	 */
	public void setStrategy(strategy inputStrategy){

		//Set strategies using input strategy vector
		Strat.setStrategy(inputStrategy);
	}
	
	public void setStrategy(int[] investment, double[] exit, double entry){
		Strat.setPolicy(investment);
		Strat.setExit(exit);
		Strat.setEntry(entry);
	}
	
	
	/**
	 * Checks if the industry state is equal to a given industry state.
	 * @param s benchmark industry state.
	 * @return "true" states are equal, "false" states are not equal.
	 */
	public boolean isEqual(int[] s){
		
		// Compare number of firms on each quality level
		for(int k=1;k<=Parameters.K;k++){
			if(s[k]!=index[k]){
				return false;
			}
		}
		return true;
	}

}