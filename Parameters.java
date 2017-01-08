package SALPMPE;

/**
 * Defines the running parameters for the SALP algorithm. 
 * @version 2.0
 */
public class Parameters {
	
	// General algorithm's parameters
	/** Maximum best response iterations */static int     	maxIter        = 10;  
	/** Minimum best response iterations */static int     	minIter        = 10;  	
	/** Maximum iterations best response computation */static int      InnermaxIt     = 5; 
	/** Initial lag best response computation */static int 		InnerLag       = 1; 	
	/** Error tolerance for convergence criterion */static double  	errTol         = 0.005; 
	/** Small double for building objective function */static double   eps            = 0.005; 
	/** Small instance? (only report full solutions for small instances) */static boolean Small           = false;
	
	// Approximating
	/** Reducing approximating architecture through linear interpolation */static boolean reducearch      = false;  		
	/** List all states in constraints (instead of sampling) */static boolean allstates       = true; 
	/** When true sets c vector (LP objective function) to one for every state */static boolean unifCvector     = false; 
			
	// Linear program parameters
	/** Maximum number of constraints on the resulting approximate linear program */static int 	 	maxStates      = 1500; 
	/** Maximum number of constraints for convergence criteria computation */static int 		maxCriteria    = 100; 
	/** Bound for the value of the LP variables */static double  	varBound       = 1000; 
	
	// Industry parameters
	/** Maximum number of firms in industry */static int 	 	N				= 3; 
	/** Number of quality levels (starts at 1) */static int 	 	K				= 10; 
	/** Allow entry and exit? */static boolean entryexit  		= true; 
	/** Scrap mean  value*/static final double epsn    	= 30.0; 	
	/** Entry cost mean value*/static final double epsne   	= 300.0;
	/** Entry state */static int   x_e 				= 2; 
	/** Initial industry state for simulations (currently everyone starts at entry state)*/static int[] s_0; 
		
	// Economic parameters
	/** One period discount factor*/static double  	beta    		= 0.925; 
	/** Marginal cost of investment*/static double  	invCost 		= 1.00; 
	/** Lagrangian penalty for slack variables*/static double 	lambda 			= 2/(1-beta);  
	/** Budget for constraint violation in SALP implementation*/static double   budget 			= 0.0; 
	/** True if relaxing budget constraint*/static boolean relaxbudget      = false;
	 
	//Profit function to use (exactly one must be true)

	/** True if running Logit price competition model */static boolean useprofit99      = false;
	/** True if running quantity competition*/static boolean useprofit10		= true;

	// Profit function 99 parameters
	/** Market size */static double  	m       		= 100.0;
	/** Marginal cost of production*/static double  	c 	    		= 0.5; 
	/** Normalization factor*/static double  	psi     		= 1.0; 
	/** Consumer income*/static double  	Y       		= 1.0;  
	/** Quality elasticity*/static double  	theta1  		= 0.5; 
	/** Price elasticity*/static double  	theta2  		= 0.5; 
	/** Maximum number of iterations for Newton method (price computation)*/static final int NiterNew  		= 100; 
	
	// Profit function 10 parameters
	/** Market size */static double   a				= 40; 
	/** Price elasticity*/static double   b 				= 10; 
	/** Minimum capacity*/static final double qmin		= 5; 
	/** Maximum capacity*/static final double qmax		= 40;
	
	// Transition dynamics parameters
	/** Depreciation factor */static double 	delta 			= 0.7; 
	/** Return to investment */static double 	alpha 			= 3; 
	/** Appreciation factor in model transition dynamics */static double 	gamma 			= 0.0; 
	/** Appreciation factor for sampling routines */static double 	simulGamma		= 0.0; 
	/** Probability of firms using initial policy (additional noise in transitions) */static double   obli	 		= 0.0; 		
	
	// Smooth update parameter
	/** Smooth update parameter */static double  	smooth 			= 0.666666;
	/** Smooth update lag parameter */static int     	smoothlag 		= 2; 
	
	// Sampling procedure tuning parameters
	/** Duration of transient for sampling states under greedy policy */static int   transcient			= 1000; 	
	/** Spacing for sampling states while simulating industry evolution */static int 	 spacing			= 1; 
	/** Sample size for simulating industry evolution*/static int 	 samplesize 		= 10000;
	/** Replications for simulation for objective formulation and sample constraints*/static int 	 nreplic			= 1; 
	/** Sample size for computing response indicator */static int 	 Vreplic    		= 30;
	/** Number of periods for computing value function for candidate best response */static int 	 Vperiods 			= 30; 

	// Possible investment decisions parameters
	/** Step size for investment grid */static double  	di				= 0.04;
	/** Size of investment grid */static int 	   	invLength 		= 50; 
	/** Size of dynamic investment grid*/static int 		localLength     = 20;
	/** Reduction factor for local (dynamic) investment grid*/static double   reduction       = 3.0;
	/** Maximum investment*/static double maxInv		    = (double) (invLength-1)*di;
	/** Investment array*/static double[]	invArray;
	/** Name of output file */static String profmodel			= "-N" + N + "-";
	/** Name of output log file */static String excelname         = "-indlogN" + N + "-";
	
	static{
		s_0 = new int[K+1];
		for(int k=1;k<=K;k++){
			s_0[k]=0;
		}
		s_0[x_e]=N; s_0[0] = N; 						// Initial state for simulations
		
		// output files' names
		if (useprofit99){ profmodel += "99"; excelname +="Q";}
		if (useprofit10){ profmodel += "10"; excelname +="C";}
		if (allstates){profmodel += "All";}
	}
}