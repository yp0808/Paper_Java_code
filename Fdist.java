package SALPMPE;

/**
 * Defines distribution function of scrap value and entry cost random variables.
 * @version 2.0
 */
public class Fdist {

	/**
	 * Defines a new instance of distribution function.
	*/
	public Fdist(){
		
	}
	
	/**
	 * probability of entry for a given threshold entry value.
	 * @param x entry threshold.
	 * @return entry probability.
	 */
	public static double FdistEntry(double x){
		// P[Setup <= x] = P[entry]
		double aux;
		if (Parameters.entryexit){
			aux = 1.0 - Math.exp(-x/Parameters.epsne);				// exponential distribution of the setup cost
			return aux;	
		} else{
			return 0;
		}
	}
	
	
	/**
	 * probability of staying on industry for a given threshold exit value.
	 * @param x exit threshold.
	 * @return probability of not exiting.
	 */
	public static double FdistExit(double x){
		// P[Scrap value <= x] = P[do not exit]		
		double aux; 
		if (Parameters.entryexit){
			aux = 1.0 - Math.exp(-x/(Parameters.epsn));// + Parameters.phii*(i-1)));				// exponential distribution of the scrap value
			return aux;
		} else{
			return 1;
		}
	}
	
	/**
	 * expected value of scrap value conditional of leaving industry for a given exit threshold.
	 * @param x exit threshold.
	 * @return conditional expected scrap value.
	 */
	public static double EdistExit(double x){
		double aux=0;
		if (Parameters.entryexit){
			aux = (Parameters.epsn + x)*Math.exp(-x/Parameters.epsn);			// exponential distribution of the scrap value
			return aux;
		} else{
			return 0;
		}
	}
	
}