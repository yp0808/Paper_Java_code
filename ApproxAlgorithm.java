package SALPMPE;

import ilog.concert.IloException;
import ilog.cplex.IloCplex.CplexStatus;

import java.io.*;
import java.text.DecimalFormat;

/**
 * Defines the approximate best response iterative algorithm. Computes convergence criteria, report solutions and long run indicators.
 * @version 2.0
 */
public class ApproxAlgorithm {
	/** combinatorial numbers, use to index states in small intances */
	static double[][] Combinatorial;
		
	/**
	 * Defines a new instance of the approximate best response iterative algorithm.
	 */
	public ApproxAlgorithm(){		
	}

	
	/**
	 * Class main routine. Declares a new algorithm instance and solves for the equilibrium.  
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException{
		ApproxAlgorithm approx = new ApproxAlgorithm();
		getCombinatorial();
		approx.solveEquilibrium();
		System.gc();	
	}
	
	
	/**
	 * Solves for the approximate equilibrium.
	 * @throws IOException
	 */
	public void solveEquilibrium() throws IOException{
		try{
			long timeini = System.currentTimeMillis();						// Read initial time
			FileWriter fstream_ind = new FileWriter("indicator" + Parameters.profmodel + ".txt");	// Output file for long run indicators
			BufferedWriter out_ind = new BufferedWriter(fstream_ind);			
			FileWriter fstream_sum = new FileWriter("log" + Parameters.profmodel + ".txt");	// Output file for algorithm evolution
			BufferedWriter out_sum = new BufferedWriter(fstream_sum);
			ApproxLP ALP;
			CplexStatus solverStatus;
			int convaux 		   = 2;										// Help in checking convergence in numerical exp
			int niter 			   = 0;										// Iteration counter			
			sol_sep[] q_sep 	   = new sol_sep[1];						// History of solution to approximate LPs
			sol_sep q_firm         = new sol_sep();							// auxiliary solution for update
			sol_sep q_aux          = new sol_sep();							// auxiliary solution for update
			q_sep[0] 		       = new sol_sep();							// solution to approximate LP for auxiliary iteration -1
			state[] sample;
			state[] altsample;
			boolean[] constraint;
			boolean[] exitCriteria;
			double[][] exit;
			double[][] altexit;
			double JBase; 													// Simulated value function using incumbent strategy 
			double JRresponse;												// Simulated value function under unilateral deviation from incumbent strategy						
						
			while (niter<=Parameters.maxIter){
				ALP     	= new ApproxLP();
				altsample 	= new state[1];
				constraint  = new boolean[1];
				if(Parameters.allstates){
					sample      = allStates(ALP, q_sep);
				} else{
					sample  	= ALP.sampleStates(q_sep,Parameters.simulGamma,Parameters.obli);									
					if (Parameters.reducearch){altsample   = ALP.sampleRest(q_sep,sample,constraint);}
				}
				constraint  = ApproxLP.getConstraints(sample,Parameters.maxStates);
				TotalNumberStates(sample.length-1,altsample.length-1);
				q_firm.setr(q_sep[q_sep.length-1]);				
				for(int i=1;i<=maxInner(niter);i++){	
					if(i>1){ALP = new ApproxLP();}					
					exit      = ALP.getExit(sample,q_firm,i);
					altexit   = ALP.getExit(altsample,q_firm,i);
					ALP.buildModel(sample,altsample,exit,altexit,constraint);
					ALP.reportMemory(i,niter);
					solverStatus = ALP.model.getCplexStatus();
					System.out.println(solverStatus.toString());					
					q_aux.setr(q_firm);
					q_firm.setr(ALP.getSolution());
					if(ALP.CompareExit(sample,q_firm,q_aux,i,out_sum)){
						break;
					} else{
						q_firm.updater(q_aux, niter);
						if(i<maxInner(niter)){
							try{
								killCplex(ALP);
								exit    = null;
								altexit = null;
								System.gc();
							}catch(Exception e){
								System.out.println(e.toString());
							}
						}
					}
				}
				
				q_sep            	  = (sol_sep[])resizeArray(q_sep,q_sep.length+1);
				q_sep[q_sep.length-1] = new sol_sep();
				q_sep[q_sep.length-1].setr(q_firm);
				q_sep[q_sep.length-1].updater(q_sep[q_sep.length-2], niter);
				if (niter>Parameters.minIter){
					exitCriteria = ApproxLP.getConstraints(sample, Parameters.maxCriteria);
					JBase      = ALP.ComputeV(sample, q_sep, "Base",exitCriteria);
					JRresponse = ALP.ComputeV(sample, q_sep, "R",exitCriteria);
					exitCriteria = null;
				} else {
					JBase      = 0;
					JRresponse = 1;
				}
				if (JBase > JRresponse){convaux++;}
				long runtime = System.currentTimeMillis()-timeini;
				solverStatus = ALP.model.getCplexStatus();
				System.out.println(solverStatus.toString());
				excelIndicator(sample, niter, runtime,solverStatus.toString());
				if (convergenceCriteria(JBase,JRresponse,Parameters.errTol*(convaux/2),niter)){
					System.out.println(Long.toString(runtime/1000) + " secs");
					out_sum.write("Runtime = " + runtime/1000 + "secs \r\n");
					System.out.println(Long.toString(runtime/1000) + " secs");
					System.out.println("Equilibrium Found at iteration "+ Integer.toString(niter));
					double budslack = ALP.getSlack(sample,altsample, constraint); 
					sample = ALP.sampleStates(q_sep, Parameters.gamma,0.0);
					computeIndicator(sample, out_ind, niter, budslack,runtime);
					reportCoeff(ALP);
					if(Parameters.Small){
						reportValue(ALP, q_sep);
						reportPolicy(ALP, q_sep);	
						System.out.println("Computing best response");
						ALP     	= new ApproxLP();			
						sample  	= ALP.sampleStates(q_sep,Parameters.simulGamma,0.0);												
						altsample 	= new state[1];
						constraint  = new boolean[1];
						if(Parameters.allstates){
							sample      = allStates(ALP, q_sep);
						} else{
							sample  	= ALP.sampleStates(q_sep,Parameters.simulGamma,Parameters.obli);									
							if (Parameters.reducearch){altsample   = ALP.sampleRest(q_sep,sample,constraint);}
						}
						constraint  = ApproxLP.getConstraints(sample,Parameters.maxStates);						
						TotalNumberStates(sample.length-1,altsample.length-1);
						q_firm.setr(q_sep[q_sep.length-1]);				
						for(int i=1;i<=maxInner(niter);i++){	
							if(i>1){ALP = new ApproxLP();}						
							exit      = ALP.getExit(sample,q_firm,i);
							altexit   = ALP.getExit(altsample,q_firm,i);
							ALP.buildModel(sample,altsample,exit,altexit,constraint);			
							System.out.println(ALP.model.getCplexStatus());
							q_aux.setr(q_firm);
							q_firm.setr(ALP.getSolution());
							if(ALP.CompareExit(sample,q_firm,q_aux,i,out_sum)){
								break;
							} else{
								q_firm.updater(q_aux, niter);
								if(i<maxInner(niter)){
									try{
										killCplex(ALP);
										System.gc();
									}catch(Exception e){
										System.out.println(e.toString());
									}
								}
							}
						}					
						q_sep            	  = (sol_sep[])resizeArray(q_sep,q_sep.length+1);
						q_sep[q_sep.length-1] = new sol_sep();
						q_sep[q_sep.length-1].setr(q_firm);
						q_sep[q_sep.length-1].updater(q_sep[q_sep.length-2], niter);					
					}
					niter = Parameters.maxIter+1;
				} else {			
					System.out.println("Iteration "+ Integer.toString(niter) + " finished.\r\n");
				}
				
				ReportProgress(JBase, JRresponse,Parameters.errTol*(convaux/2), niter, out_sum);
				niter++;
				try{
					killCplex(ALP);
					sample = null;
					constraint = null;
					altsample = null;	
					exit    = null;
					altexit = null;
				} catch(Exception e){
					System.out.println(e.toString());
				}
				System.gc();
			}
			out_ind.close();
			out_sum.close();
			if (niter == Parameters.maxIter){System.out.println("Maximum number of iterations reached.");}
		}catch(IloException e){
			System.out.println(e.toString());	
		}
	}
	
	
	/**
	 * Kills CPLEX instances and tries to release memory.
	 * @param ALP Instance of CPLEX model.
	 */
	public void killCplex(ApproxLP ALP){
		ALP.model.end();
		ALP = null;		
	}
	
	
	/**
	 * Checks if convergence criteria is met.
	 * @param JBase	Simulated value function for a firm using the approximate best response to the incumbent strategy when competitors use the incumbent strategy.
	 * @param JRresponse Simulated value function for a firm when all industry uses the incumbent strategy.
	 * @param Tolerr Tolerance error for the convergence criteria.
	 * @param niter Iteration number.
	 * @return "true"- equilibrium found, "false" iterate.
	 */
	public boolean convergenceCriteria(double JBase,double JRresponse, double Tolerr, int niter){
		
		// Compute simulated absolute relative increase in profit by unilateraly deviating from incumbent strategy
		double relErr = Math.min(Math.abs((JRresponse-JBase)/JBase),Math.abs(JRresponse-JBase));
		if (relErr> Tolerr){ 
			if (niter==Parameters.maxIter){
				System.out.println("Max iterations completed. Finishing.");
				return true;
			} else {return false;}
		}
		else { 
			if (niter<Parameters.minIter){				
				return false;
			}else{
				return true;
			}
		}
	}	
	
	
	/**
	 * Reports solution coefficients.
	 * @param ALP Instance of CPLEX model. 
	 */
	public void reportCoeff(ApproxLP ALP){
		try{
			FileWriter fstream = new FileWriter("coefficients" + Parameters.profmodel + ".txt");	// Output file for reporting equilibrium coefficients
			BufferedWriter out = new BufferedWriter(fstream);
			String ln;
			
			// Get last approximate LP solution
			sol_sep u = ALP.getSolution();						
			ln        = " r_sep(x,k,n) coefficients \r\n";
			out.write(ln);
			for(int x=1;x<=Parameters.K;x++){
				for(int k=1;k<=Parameters.K;k++){
					for(int n=0;n<Parameters.N;n++){
						
						// Write equilibrium coefficients
						ln  = Integer.toString(x) + " ; " +	Integer.toString(k) + " ; ";
						ln += Integer.toString(n) + " ;" + Double.toString(u.r_sep[x][k][n]) + " \r\n";
						out.write(ln);
					}					
				}	
			}						
			ln = " r_const coefficient \r\n";				
			ln = Double.toString(u.r_const)+ " \r\n";
			out.write(ln);
			out.close();
		} catch(Exception e){
			System.out.println(e.toString());
		}
	}
	
	
	/**
	 * Reports entry and exit equilibrium policies.
	 * @param ALP Instance of approximate algorithm.
	 * @param q_sep History of approximate linear program solutions.
	 */
	public void reportPolicy(ApproxLP ALP, sol_sep[] q_sep){ 	
		try{
			int[] s;
			//FileWriter fstreamP  = new FileWriter("AppPolicy" + Parameters.profmodel + ".txt");
			FileWriter fstreamEx = new FileWriter("exit"   + Parameters.profmodel + ".txt");
			FileWriter fstreamEn = new FileWriter("entry"  + Parameters.profmodel + ".txt");			
			BufferedWriter outEx = new BufferedWriter(fstreamEx);
			BufferedWriter outEn = new BufferedWriter(fstreamEn); 
			String lnEx;
			String lnEn;
			strategy Strat		= new strategy();
			// Go over all possible values for pair firm state-competitors state
			for(int j = 0;j<(int) Combinatorial[Parameters.N-1][Parameters.K];j++){
				s      = ALP.inverseIndex(j,Combinatorial);
				for(int x=1;x<=Parameters.K;x++){					
				// Recovering industry state										
					s[x]++;s[0]++;					
					// Compute profit and strategy vectors for incumbent industry state												
					Strat  = ALP.IterativeStrategy(s, q_sep, q_sep.length-1);	
					s[x]--;s[0]--;	
					
					// Write equilibrium strategy
					
					lnEx = Double.toString(Strat.exit[x]) + "\r\n";					
					outEx.write(lnEx);					
				}
				Strat= ALP.IterativeStrategy(s, q_sep, q_sep.length-1);
				lnEn = Double.toString(Strat.entry) + "\r\n";
				outEn.write(lnEn);
			}
			outEx.close(); outEn.close();
		} catch(Exception e){
			System.out.println(e.toString());
		}
	}
	
	
	/**
	 * Reports value function approximation and equilibrium investment policy.
	 * @param ALP Instance of approximate algorithm.
	 * @param q_sep History of approximate linear program solutions.
	 */
	public void reportValue(ApproxLP ALP, sol_sep[] q_sep){ 	
		try{
			double value;
			int[] s;
			FileWriter fstreamP  = new FileWriter("investment" + Parameters.profmodel + ".txt");
			FileWriter fstreamV  = new FileWriter("values" + Parameters.profmodel + ".txt");
			BufferedWriter outP  = new BufferedWriter(fstreamP);
			BufferedWriter outV  = new BufferedWriter(fstreamV);
			String lnP;
			String lnV; 
			strategy Strat		 = new strategy();
			sol_sep u            = ALP.getSolution();
			// Go over all possible values for pair firm state-competitors state
			for(int x=1;x<=Parameters.K;x++){				
				for(int j = 0;j<(int) Combinatorial[Parameters.N-1][Parameters.K];j++){					
					value  = u.r_const;
					s      = ALP.inverseIndex(j,Combinatorial);
					// Recovering industry state										
					s[x]++;s[0]++;					
					// Compute profit and strategy vectors for incumbent industry state												
					Strat  = ALP.IterativeStrategy(s, q_sep, q_sep.length-1);	
					s[x]--;s[0]--;
					// Recovering industry state
					// Write state	
					for(int k=1;k<=Parameters.K;k++){
						value += u.r_sep[x][k][s[k]];		// forming value function approximation
					}
					// Write value function approximation
					lnP  = Double.toString(Parameters.invArray[Strat.policy[x]]) + "\r\n";
					lnV = Double.toString(value) +" \r\n";
					outV.write(lnV);
					outP.write(lnP);					
					// Write equilibrium strategy
				}
			}
			outV.close();
			outP.close(); 
		} catch(Exception e){
			System.out.println(e.toString());
		}
	}
	
	
	/**
	 * Computes simulated long run indicators of interest.
	 * @param sample Sampled set of states.
	 * @param out_ind Buffered writer for reporting indicators.
	 * @param niter Iteration number.
	 * @param budslack Value of slack variable in SALP implementation.
	 * @param runningtime Algorithm's total running time.
	 */
	public void computeIndicator(state[] sample , BufferedWriter out_ind, int niter, double budslack, double runningtime){ 	
		try{
			String name = Parameters.excelname;
			if (Parameters.relaxbudget){ name +="relax";}
			FileWriter fstream_excel = new FileWriter(name + ".txt", true);	// Output file for long run indicators
			BufferedWriter out_excel = new BufferedWriter(fstream_excel);
			double[] statdist 	   = new double[Parameters.K+1];
			double[] sigmastatdist = new double[Parameters.K+1];
			double[] statpro 	   = new double[Parameters.K+1];
			double[] sigmastatpro  = new double[Parameters.K+1];
			int[] s;								// Incumbent sampled state
			double TotalInvestment = 0;				// Mean total investment
			double sigmaTotalInvestment = 0;		// Mean total investment
			double ProducerSurplus = 0;				// Mean producer surplus
			double sigmaProducerSurplus = 0;		// Mean producer surplus
			double ConsumerSurplus = 0;				// Mean consumer surplus
			double sigmaConsumerSurplus = 0;		// Mean consumer surplus
			double EntryRate	   = 0;				// Mean entry rate
			double sigmaEntryRate  = 0;				// Mean entry rate
			double total_freq	   = 0;				// Normalizing factor
			double[] Conce		   = new double[Parameters.N+1]; // compute industry concentration			
			double[] FreqConce	   = new double[Parameters.N+1]; // compute industry concentration
			double Conceaux 	   = 0;				// compute industry concentration
			int Indaux = 0; int ncount = 0;			// compute industry concentration
			double[] aux_freq_x	   = new double[Parameters.K+1];
			double aux_freq;  
			double[] surplus 	   = new double[2]; // Consumer and Producer surplus
			double[] mshare;						// Market shares vector for incumbent sampled state
			Profit profit = new Profit();
			String ln ="";
			//Write parameters
			
			if (Parameters.useprofit99){
				out_ind.write("Profit : Logit \r\n");
				out_ind.write("Market size : " + Double.toString(Parameters.m)+ "\r\n");
				out_ind.write("Consumer Income : " + Double.toString(Parameters.Y)+ "\r\n");
				out_ind.write("theta1 : " + Double.toString(Parameters.theta1) + "\r\n");
				out_ind.write("theta2 : " + Double.toString(Parameters.theta2)+ "\r\n");
				out_ind.write("Psi : " + Double.toString(Parameters.psi)+ "\r\n");
				out_ind.write("Marginal Cost : " + Double.toString(Parameters.c)+ "\r\n");
			}
			if (Parameters.useprofit10){
				out_ind.write("Profit : Quantity competition \r\n");
				out_ind.write("Max demand: " + Double.toString(Parameters.a)+ "\r\n");
				out_ind.write("Price elasticity: " + Double.toString(Parameters.b)+ "\r\n");
				out_ind.write("q_min: " + Double.toString(Parameters.qmin)+ "\r\n");
				out_ind.write("q_max: " + Double.toString(Parameters.qmax)+ "\r\n");
			}
			if (Parameters.entryexit){
				out_ind.write(" Entry and exit allowed \r\n");
			}else{
				out_ind.write(" No entry/exit \r\n");
			}
			out_ind.write("Grid space : " + Double.toString(Parameters.di)+ "\r\n");
			out_ind.write("Grid size : " + Double.toString(Parameters.invLength)+ "\r\n");
			out_ind.write("alpha : " + Double.toString(Parameters.alpha)+ "\r\n");
			out_ind.write("delta : " + Double.toString(Parameters.delta)+ "\r\n");
			out_ind.write("up : " + Double.toString(Parameters.gamma)+ "\r\n");
			out_ind.write("Discount factor : " + Double.toString(Parameters.beta)+ "\r\n");
			out_ind.write("Investment Cost : " + Double.toString(Parameters.invCost)+ "\r\n");
			out_ind.write("Maximum Number of Firms : " + Double.toString(Parameters.N)+ "\r\n");
			out_ind.write("xmax : " + Double.toString(Parameters.K) + "\r\n");			
			out_ind.write("espn : " + Double.toString(Parameters.epsn) + "\r\n");
			out_ind.write("espne : " + Double.toString(Parameters.epsne) + "\r\n");
			out_ind.write("inner it max : " + Integer.toString(Parameters.InnermaxIt) + "\r\n");
			out_ind.write("inner it lag : " + Integer.toString(Parameters.InnerLag) + "\r\n");
			
			// Go over all sampled states
			for(int i =1 ;i<sample.length;i++){
				aux_freq = (double)sample[i].frequency;
				total_freq += aux_freq;
				s = (int[]) sample[i].index.clone();
				
				// Compute consumer surplus and market shares for incumbent state
				if(s[0]<Parameters.N){
					EntryRate +=  (Parameters.N-s[0])*aux_freq*Fdist.FdistEntry(sample[i].Strat.entry);
				}
				if(s[0]>0){
					statdist[0] +=s[0]*aux_freq;
					surplus =profit.surplus(s);
					mshare =profit.mshare(s);
					for(int k=1;k<s.length;k++){
						statdist[k] += s[k]*aux_freq;
						if(s[k]>0){
							aux_freq_x[k]+=aux_freq;
							statpro[k] +=sample[i].profit[k]*aux_freq;
							// Computing mean simulated indicators
							TotalInvestment += aux_freq*s[k]*Parameters.invArray[sample[i].Strat.policy[k]];
						}
					}
					ConsumerSurplus += aux_freq*surplus[0];
					ProducerSurplus += aux_freq*surplus[1];
					
					// Compute industry concentration
					Conceaux	= 0;
					Indaux	 	= Parameters.K;
					ncount      = 1;
					while (s[0]>0){
						if (s[Indaux]>0){
							Conceaux 		  += mshare[Indaux];
							Conce[ncount] 	  += Conceaux*aux_freq;
							FreqConce[ncount] += aux_freq;
							ncount ++;
							s[0]--; s[Indaux]--;
						} else{
							Indaux--;
						}
					}
				}
			}
			for(int n=1;n<=Parameters.N;n++){
				if(FreqConce[n]>0){
					Conce[n] = Conce[n]/FreqConce[n];
				} else {
					Conce[n] = 0;
				}
			}
			TotalInvestment		= TotalInvestment/total_freq;
			EntryRate 			= EntryRate/total_freq;
			ProducerSurplus 	= ProducerSurplus/total_freq;
			ConsumerSurplus 	= ConsumerSurplus/total_freq;
			for(int k=0;k<=Parameters.K;k++){
				statdist[k] 	= statdist[k]/total_freq; 
			}
			for(int k=1;k<=Parameters.K;k++){
				if (aux_freq_x[k]>0){
					statpro[k] 	= statpro[k]/aux_freq_x[k];
				}
			}
			// compute variance of indicators
			int auxp;
			for(int i =1 ;i<sample.length;i++){
				aux_freq = (double)sample[i].frequency;				
				s = (int[]) sample[i].index.clone();				
				// Compute consumer surplus and market shares for incumbent state
				if(s[0]<Parameters.N){
					sigmaEntryRate +=  aux_freq*Math.pow((Parameters.N-s[0])*Fdist.FdistEntry(sample[i].Strat.entry)- EntryRate,2);
				}
				if(s[0]>0){
					sigmastatdist[0] +=Math.pow(s[0]-statdist[0],2)*aux_freq;
					surplus =profit.surplus(s);
					auxp=0;
					for(int k=1;k<s.length;k++){
						sigmastatdist[k] +=Math.pow(s[k]-statdist[k],2)*aux_freq;
						if(s[k]>0){
							sigmastatpro[k] +=Math.pow(sample[i].profit[k]-statpro[k],2)*aux_freq;
							// Computing mean simulated indicators
							auxp += s[k]*Parameters.invArray[sample[i].Strat.policy[k]];							
						}
					}
					sigmaTotalInvestment += aux_freq*Math.pow(auxp-TotalInvestment,2);
					sigmaConsumerSurplus += aux_freq*Math.pow(surplus[0]-ConsumerSurplus,2);
					sigmaProducerSurplus += aux_freq*Math.pow(surplus[1]-ProducerSurplus,2);					
				}
			}
			sigmaTotalInvestment	= sigmaTotalInvestment/total_freq;
			sigmaEntryRate 			= sigmaEntryRate/total_freq;
			sigmaProducerSurplus 	= sigmaProducerSurplus/total_freq;
			sigmaConsumerSurplus 	= sigmaConsumerSurplus/total_freq;
			for(int k=0;k<=Parameters.K;k++){
				sigmastatdist[k] 	= sigmastatdist[k]/total_freq; 
			}
			for(int k=1;k<=Parameters.K;k++){
				if (aux_freq_x[k]>0){
					sigmastatpro[k] 	= sigmastatpro[k]/aux_freq_x[k];
				}
			}
			for(int k=1;k<=Parameters.K;k++){
				ln += ":" +Double.toString(statdist[k]) +": stdv : " +Double.toString(Math.sqrt(sigmastatdist[k])) + "\r\n"; 
			}
			ln += ":" + Double.toString(statdist[0]) + ": stdv :" + Double.toString(Math.sqrt(sigmastatdist[0])) + "\r\n";
			// Write mean simulated indicators normalizing by total frequency
			out_ind.write("Iterations : " + Integer.toString(niter)+ "\r\n");
			out_ind.write("Total Investment : " + Double.toString(TotalInvestment) + ": stdv:" + Double.toString(Math.sqrt(sigmaTotalInvestment)) +  "\r\n");
			out_ind.write("Entry Rate : " + Double.toString(EntryRate) + ": stdv:" + Double.toString(Math.sqrt(sigmaEntryRate)) + "\r\n");
			out_ind.write("Producer Surplus : " + Double.toString(ProducerSurplus)+ ": stdv:" + Double.toString(Math.sqrt(sigmaProducerSurplus)) + "\r\n");
			out_ind.write("Consumer Surplus : " + Double.toString(ConsumerSurplus)+ ": stdv:" + Double.toString(Math.sqrt(sigmaConsumerSurplus)) + "\r\n");
			out_ind.write("Industry concentration: \r\n");
			for (int n=1;n<=Parameters.N;n++){
				out_ind.write("C" + Integer.toString(n) +" : " + Double.toString(Conce[n]) + "\r\n");
			}
			out_ind.write("\r\n");
			out_ind.write("Stat distribution: \r\n");
			out_ind.write(ln + "\r\n");
			ln = "";
			for(int k=1;k<=Parameters.K;k++){
				ln += ":" +Double.toString(statpro[k]) +": stdv : " +Double.toString(Math.sqrt(sigmastatpro[k])) + "\r\n"; 
			}
			out_ind.write("profits: \r\n");
			out_ind.write(ln + "\r\n");
			System.out.println("Total Investment : " 	+ Double.toString(TotalInvestment));
			if (Parameters.Small){
				out_excel.write(Integer.toString(Parameters.N) + ":" + 
					Double.toString(TotalInvestment) + ":" + 
					Double.toString(ProducerSurplus) + ":" + 
					Double.toString(ConsumerSurplus) + ":" +					
					Double.toString(Conce[1]) + ":" + 
					Double.toString(Conce[2]) + ":" + 
					Double.toString(EntryRate) + ":" + 
					Double.toString(budslack) + ":" + 
					Double.toString(Parameters.gamma) +":" + 
					Double.toString(runningtime/1000) + ":" + 
					Integer.toString(niter) + "\r\n");
				out_excel.write(":" + 
						Double.toString(Math.sqrt(sigmaTotalInvestment)) + ":" + 
						Double.toString(Math.sqrt(sigmaProducerSurplus)) + ":" + 
						Double.toString(Math.sqrt(sigmaConsumerSurplus)) + ":" + " : " + " : " +											 
						Double.toString(Math.sqrt(sigmaEntryRate)) + "\r\n \r\n");
			}else {
				out_excel.write(Integer.toString(Parameters.N) + ":" + 
						Double.toString(TotalInvestment) + ":" + 
						Double.toString(ProducerSurplus) + ":" + 
						Double.toString(ConsumerSurplus) + ":" +					
						Double.toString(Conce[6]) + ":" + 
						Double.toString(Conce[12]) + ":" + 
						Double.toString(EntryRate) + ":" + 
						Double.toString(budslack) + ":" + 
						Double.toString(Parameters.gamma) +":" + 
						Double.toString(runningtime/1000) + ":" + 
						Integer.toString(niter) + "\r\n");
				out_excel.write(":" + 
						Double.toString(Math.sqrt(sigmaTotalInvestment)) + ":" + 
						Double.toString(Math.sqrt(sigmaProducerSurplus)) + ":" + 
						Double.toString(Math.sqrt(sigmaConsumerSurplus)) + ":" + " : " + " : " +											 
						Double.toString(Math.sqrt(sigmaEntryRate)) + "\r\n \r\n");
			}
			out_excel.close();
		} catch(Exception e){
			System.out.println(e.toString());
		}
	}
	
	
	/**
	 * Computes simulated long run indicators of interest (for log file).
	 * @param sample Sampled set of states.
	 * @param niter Iteration number.
	 * @param runningtime Algorithm's total running time.
	 * @param status CPLEX solver exit status.
	 */
	public void excelIndicator(state[] sample , int niter, double runningtime, String status){ 	
		try{
			String name = Parameters.excelname;
			if (Parameters.relaxbudget){ name +="relax";}
			FileWriter fstream_excel = new FileWriter(name + ".txt", true);	// Output file for long run indicators
			BufferedWriter out_excel = new BufferedWriter(fstream_excel);
			double[] statdist 	     = new double[Parameters.K+1];
			double[] sigmastatdist 	 = new double[Parameters.K+1];
			double[] statpro 	     = new double[Parameters.K+1];
			double[] sigmastatpro 	 = new double[Parameters.K+1];
			int[] s;								// Incumbent sampled state
			double TotalInvestment = 0;				// Mean total investment
			double sigmaTotalInvestment = 0;		// Mean total investment
			double ProducerSurplus = 0;				// Mean producer surplus
			double sigmaProducerSurplus = 0;		// Mean producer surplus
			double ConsumerSurplus = 0;				// Mean consumer surplus
			double sigmaConsumerSurplus = 0;		// Mean consumer surplus
			double EntryRate	   = 0;				// Mean entry rate
			double sigmaEntryRate  = 0;	
			double total_freq	   = 0;				// Normalizing factor
			double[] Conce		   = new double[Parameters.N+1]; // compute industry concentration			
			double[] FreqConce	   = new double[Parameters.N+1]; // compute industry concentration
			double Conceaux 	   = 0;				// compute industry concentration
			int Indaux = 0; int ncount = 0;			// compute industry concentration
			double aux_freq;  
			double[] surplus 	   = new double[2]; // Consumer and Producer surplus
			double[] mshare;						// Market shares vector for incumbent sampled state
			Profit profit = new Profit();		
					
			// Go over all sampled states
			for(int i =1 ;i<sample.length;i++){
				aux_freq = (double)sample[i].frequency;
				total_freq += aux_freq;
				s = (int[]) sample[i].index.clone();
				
				// Compute consumer surplus and market shares for incumbent state
				if(s[0]<Parameters.N){
					EntryRate +=  (Parameters.N-s[0])*aux_freq*Fdist.FdistEntry(sample[i].Strat.entry);
				}
				if(s[0]>0){
					statdist[0] +=s[0]*aux_freq;
					surplus =profit.surplus(s);
					mshare =profit.mshare(s);
					for(int k=1;k<s.length;k++){
						statdist[k] += s[k]*aux_freq;
						if(s[k]>0){
							statpro[k]+=aux_freq*sample[i].profit[k];
							// Computing mean simulated indicators
							TotalInvestment += aux_freq*s[k]*Parameters.invArray[sample[i].Strat.policy[k]];
						}
					}
					ConsumerSurplus += aux_freq*surplus[0];
					ProducerSurplus += aux_freq*surplus[1];
					
					// Compute industry concentration
					Conceaux	= 0;
					Indaux	 	= Parameters.K;
					ncount      = 1;
					while (s[0]>0){
						if (s[Indaux]>0){
							Conceaux 		  += mshare[Indaux];
							Conce[ncount] 	  += Conceaux*aux_freq;
							FreqConce[ncount] += aux_freq;
							ncount ++;
							s[0]--; s[Indaux]--;
						} else{
							Indaux--;
						}
					}
				}
			}
			for(int n=1;n<=Parameters.N;n++){
				if(FreqConce[n]>0){
					Conce[n] = Conce[n]/FreqConce[n];
				} else {
					Conce[n] = 0;
				}
			}
			TotalInvestment		= TotalInvestment/total_freq;
			EntryRate 			= EntryRate/total_freq;
			ProducerSurplus 	= ProducerSurplus/total_freq;
			ConsumerSurplus 	= ConsumerSurplus/total_freq;
			for(int k=0;k<=Parameters.K;k++){
				statdist[k] 	= statdist[k]/total_freq; 
			}
			for(int k=1;k<=Parameters.K;k++){
				statpro[k] 	= statpro[k]/total_freq; 
			}
			// compute variance of indicators
			int auxp;
			for(int i =1 ;i<sample.length;i++){
				aux_freq = (double)sample[i].frequency;				
				s = (int[]) sample[i].index.clone();				
				// Compute consumer surplus and market shares for incumbent state
				if(s[0]<Parameters.N){
					sigmaEntryRate +=  aux_freq*Math.pow((Parameters.N-s[0])*Fdist.FdistEntry(sample[i].Strat.entry)- EntryRate,2);
				}
				if(s[0]>0){
					sigmastatdist[0] +=Math.pow(s[0]-statdist[0],2)*aux_freq;
					surplus =profit.surplus(s);
					auxp=0;
					for(int k=1;k<s.length;k++){
						sigmastatdist[k] +=Math.pow(s[k]-statdist[k],2)*aux_freq;
						sigmastatpro[k] +=Math.pow(sample[i].profit[k]-statpro[k],2)*aux_freq;
						if(s[k]>0){
							// Computing mean simulated indicators
							auxp += s[k]*Parameters.invArray[sample[i].Strat.policy[k]];							
						}
					}
					sigmaTotalInvestment += aux_freq*Math.pow(auxp-TotalInvestment,2);
					sigmaConsumerSurplus += aux_freq*Math.pow(surplus[0]-ConsumerSurplus,2);
					sigmaProducerSurplus += aux_freq*Math.pow(surplus[1]-ProducerSurplus,2);					
				}
			}
			sigmaTotalInvestment	= sigmaTotalInvestment/total_freq;
			sigmaEntryRate 			= sigmaEntryRate/total_freq;
			sigmaProducerSurplus 	= sigmaProducerSurplus/total_freq;
			sigmaConsumerSurplus 	= sigmaConsumerSurplus/total_freq;
			for(int k=0;k<=Parameters.K;k++){
				sigmastatdist[k] 	= sigmastatdist[k]/total_freq; 
			}
			for(int k=1;k<=Parameters.K;k++){
				sigmastatpro[k] 	= sigmastatpro[k]/total_freq; 
			}
			
			// Write mean simulated indicators normalizing by total frequency			
			System.out.println("Total Investment : " 	+ Double.toString(TotalInvestment));
			if (Parameters.Small){
				out_excel.write(Integer.toString(Parameters.N) + ":" + 
						Double.toString(TotalInvestment) + ":" + 
						Double.toString(ProducerSurplus) + ":" + 
						Double.toString(ConsumerSurplus) + ":" + 
						Double.toString(Conce[1]) + ":" + 
						Double.toString(Conce[2]) + ":" + 
						Double.toString(EntryRate) + ":" + 
						Double.toString(Parameters.gamma) +":" + 
						Double.toString(runningtime/1000) + ":" + 
						Integer.toString(niter) + "\r\n");
				out_excel.write(":" + 
						Double.toString(Math.sqrt(sigmaTotalInvestment)) + ":" + 
						Double.toString(Math.sqrt(sigmaProducerSurplus)) + ":" + 
						Double.toString(Math.sqrt(sigmaConsumerSurplus)) + ":" + " : " + " : " +											 
						Double.toString(Math.sqrt(sigmaEntryRate)) + "\r\n");
			} else{
				out_excel.write(Integer.toString(Parameters.N) + ":" + 
						Double.toString(TotalInvestment) + ":" + 
						Double.toString(ProducerSurplus) + ":" + 
						Double.toString(ConsumerSurplus) + ":" + 
						Double.toString(Conce[6]) + ":" + 
						Double.toString(Conce[12]) + ":" + 
						Double.toString(EntryRate) + ":" + 
						Double.toString(Parameters.gamma) +":" + 
						Double.toString(runningtime/1000) + ":" + 
						status + ":" +
						Integer.toString(niter) + "\r\n");
				out_excel.write(":" + 
						Double.toString(Math.sqrt(sigmaTotalInvestment)) + ":" + 
						Double.toString(Math.sqrt(sigmaProducerSurplus)) + ":" + 
						Double.toString(Math.sqrt(sigmaConsumerSurplus)) + ":" + " : " + " : " +											 
						Double.toString(Math.sqrt(sigmaEntryRate)) + "\r\n");
			}
			out_excel.close();
		} catch(Exception e){
			System.out.println(e.toString());
		}
	}
	
	
	/**
	 * Prints out algorithm's progress.
	 * @param JBase Simulated value function for a firm using the approximate best response to the incumbent strategy when competitors use the incumbent strategy.
	 * @param JRresponse Simulated value function for a firm when industry use the incumbent strategy.
	 * @param Tolerr Tolerance error for the convergence criteria.
	 * @param niter Iteration number.
	 * @param out_sum Buffered writer for reporting algorithm evolution.
	 */
	public void ReportProgress(double JBase, double JRresponse,double Tolerr, int niter, BufferedWriter out_sum){
		try{			
			// Showing and writing convergence criteria related quantities 
			System.out.println("JBase      = "+ Double.toString(JBase));
			System.out.println("JRresponse = "+ Double.toString(JRresponse));			
			System.out.println("errTol     = "+ Double.toString(Tolerr));
			out_sum.write(Integer.toString(niter)+ ";" + Double.toString(JBase) + ";" + Double.toString(JRresponse) +"\r\n");
		} catch(Exception e){
			System.out.println(e.toString());
		}
	}
	
		
	/**
	 * Computes total number of states use in LP formulation.
	 * @param samplesize1 Size of the sampled set of states.
	 * @param samplesize2 Size of the additional set of states.
	 */
	public void TotalNumberStates(int samplesize1,int samplesize2){
		double Nstate = 1;
		DecimalFormat df = new DecimalFormat("00.00");
		
		// Computing total number of states (can not use pre-computed combinatorial numbers)
		for(int k=1;k<Parameters.K+1;k++){
			Nstate = (double)Nstate*(Parameters.N+k)/k;
		}
		Nstate--;
		String per = df.format((samplesize1)/Nstate*100.00);
		
		// Showing percentage of sampled states relative to the total number of states
		System.out.println("Sampling " + Integer.toString(samplesize1) + " states out of " + Double.toString(Math.round(Nstate)) +" = " + per + "%");
		System.out.println("Adding  " + Integer.toString(samplesize2) + " states");		
	}
	
	
	/**
	 * Computes maximum number of inner iterations in computation of approximate best response.
	 * @param niter Iteration count.
	 */
	public int maxInner(int niter){
		return Math.min(Parameters.InnermaxIt,Math.max(1,2+niter-Parameters.InnerLag));
	}	
	
	/**
	* Reallocates an array with a new size, and copies the contents.
	* of the old array to the new array.
	* @param oldArray  the old array, to be reallocated.
	* @param newSize   the new array size.
	* @return          A new array with the same contents.
	*/
	//@SuppressWarnings("unchecked")
	private static Object resizeArray (Object oldArray, int newSize) {
		int oldSize = java.lang.reflect.Array.getLength(oldArray);
	    Class elementType = oldArray.getClass().getComponentType();
	    Object newArray = java.lang.reflect.Array.newInstance(elementType,newSize);
	    int preserveLength = Math.min(oldSize,newSize);
	    if (preserveLength > 0)
	    	System.arraycopy (oldArray,0,newArray,0,preserveLength);
	    return newArray; 
	}


	/**
	 * Computes combinatorial coefficients use to index states for small instances.
	 */
	static void getCombinatorial(){
		Combinatorial = new double[Math.max(Parameters.N+1,Parameters.K+1)][Math.max(Parameters.N+1,Parameters.K+1)];
		for(int k=0;k<Math.max(Parameters.N+1,Parameters.K+1);k++){
			Combinatorial[k][0] = 1;
			for(int n=1;n<Math.max(Parameters.N+1,Parameters.K+1);n++){
				Combinatorial[k][n] = (double)Combinatorial[k][n-1]*(n+k)/n;
			}
		}
	}
	
	
	/**
	 * Returns all states in state space, and their strategies.
	 * @param q_sep  History of solutions to the approximate linear programs.
	 * @return set of possible states in the state space with their strategies.
	 */
	public state[] allStates(ApproxLP ALP, sol_sep[] q_sep){		
		try{
			int[] s;
			Profit profit   = new Profit();
			strategy Strat  = new strategy();
			double[] prof   = new double[Parameters.K+1];			
			state[] sample  = new state[(int) Combinatorial[Parameters.N][Parameters.K]];		
			for(int j = 1;j<(int) Combinatorial[Parameters.N][Parameters.K];j++){
				s      = ALP.inverseIndex(j,Combinatorial);				
				prof   = profit.profit(s);												
				Strat  = ALP.IterativeStrategy(s, q_sep, q_sep.length-1);					
				sample[j] = new state(s);
				sample[j].setStrategy(Strat);
				sample[j].setProfit(prof);
			}
			return sample;
		} catch(Exception e){
			System.out.println(e.toString());
			return null;
		}
	}
	
}