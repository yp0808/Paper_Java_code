package SALPMPE;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

import java.io.*;
import java.util.Random;

import SALPMPE.Parameters;
import SALPMPE.strategy;


/**
 * Defines an instance of SALP program, builds model and recovers coefficients in approximation.
 * @version 2.0
 */
public class ApproxLP {
	/** CPLEX LP model */IloCplex model; 
	/** Initial strategy (compact representation) */strategy compact = new strategy(); 
	/** Model decision variables (separable functions generalization)*/IloNumVar[][][] v_sep; 
	/** Model decision variable  (base value)*/IloNumVar v_const; 
	/** Slack variable for SALP implementation*/IloNumVar[][] v_slack; 
	/** Transition probability class*/prob P = new prob(); 
	/** Profit auxiliary array*/Profit profit; 
	/** Random number generator for simulating chain evolution*/Random generator = new Random(); 
	/** Auxiliary variable (keeps track of initial simulation point)*/int xini; 
	/** Auxiliary variable (keeps track of initial simulation point)*/int jini; 
	/** Coarseness of linear interpolation of separable functions*/int[][] nodes = new int[Parameters.K+1][]; 
	/** First level of differentiation in coarseness of approximation*/int[]   exeHi; 
	/** Second level of differentiation in coarseness of approximation*/int exeLow; 
	/** Third level of differentiation in coarseness of approximation*/int[] exeMid;  
	
	
	/**
	 * Defines a new IloCplex model, loads profit function, sets coarseness of approximation and loads initial strategy.
	 */
	public ApproxLP(){
		profit = new Profit();
		setExeNodes();
		loadCompact();
	}
	
	
	/**
	 * Constructs and solves the approximating linear program using the sampled set of states. 
	 * @param sample Sampled set of states used to construct the approximating linear program.
	 * @param altsample Additional set of states that complete use of basis functions.
	 * @param exit Exit strategies for states in sample set.
	 * @param altexit Exit strategies for states in additional sample set. 
	 * @param constraint indicates if a state in sample set will generate a set of constraints in the LP (limits size of LPs without affecting sampling).
	 */
	public void buildModel(state[] sample, state[] altsample, double[][] exit, double[][] altexit, boolean[] constraint){
		try{			
			double[][][] coef_sep = new double[Parameters.K+1][Parameters.K+1][Parameters.N];
			double aux_const      = 0;
			double aux_freq;
			double aux_obj;
			int[] s;			
			// Creating LP model
			model   = new IloCplex();											
			v_sep   = new IloNumVar[Parameters.K+1][Parameters.K+1][];
			v_slack = new IloNumVar[Parameters.K+1][];
			v_const = model.numVar(-Parameters.varBound, Parameters.varBound);			
			for(int x=1;x<=Parameters.K;x++){
				v_slack[x] = model.numVarArray(sample.length+altsample.length-1,0,Parameters.varBound);
				for(int k=1;k<=Parameters.K;k++){
					v_sep[x][k] = model.numVarArray((Parameters.N), -Parameters.varBound, Parameters.varBound);		
					for(int n=0;n<Parameters.N;n++){					
						coef_sep[x][k][n] = 0; 
					}
				}
			}		
			
			// reduce base functions
			if (Parameters.reducearch){reduceBaseFunctions();}   		
			
			// Objective function			
			IloLinearNumExpr objexpr = model.linearNumExpr();
			IloLinearNumExpr budexpr = model.linearNumExpr();
			
			// First we go over sampled states
			for(int i =1;i<sample.length;i++){
				if(constraint[i]){
					if(Parameters.unifCvector){aux_freq  = 1.0;}
					else{aux_freq  = (double) sample[i].frequency;}					
					s         = (int[]) sample[i].index.clone();		
					for(int x=1;x<=Parameters.K;x++){
						if(s[x]>0){
							// SALP implementation
							if (Parameters.relaxbudget){
								objexpr.addTerm(Parameters.lambda*aux_freq, v_slack[x][i]);
							} else{
								budexpr.addTerm(aux_freq, v_slack[x][i]);
							}
							aux_obj    = aux_freq;
							aux_const += aux_obj;
							s[x]--; s[0]--;
							for(int k=1;k<=Parameters.K;k++){
								coef_sep[x][k][s[k]] += aux_obj;
							}
							s[x]++;s[0]++;
						}
					}	
				}
			}
			
			// Now we go over additional states
			for(int i =1;i<altsample.length;i++){
				aux_freq  = (double) altsample[i].frequency;
				s         = (int[]) altsample[i].index.clone();		
				for(int x=1;x<=Parameters.K;x++){
					if(s[x]>0){					
						if (Parameters.relaxbudget){
							objexpr.addTerm(Parameters.lambda*aux_freq, v_slack[x][sample.length+i-1]);
						} else{
							budexpr.addTerm(aux_freq, v_slack[x][sample.length+i-1]);
						}
						aux_obj    = aux_freq;
						aux_const += aux_obj;
						s[x]--; s[0]--;
						for(int k=1;k<=Parameters.K;k++){
							coef_sep[x][k][s[k]] += aux_obj;
						}
						s[x]++;s[0]++;
					}
				}	
			}
			
			for(int x=1;x<=Parameters.K;x++){
				for(int k=1;k<=Parameters.K;k++){
					for(int n=0;n<Parameters.N;n++){
						objexpr.addTerm(coef_sep[x][k][n],v_sep[x][k][n]);						
					}
				}
			}
			objexpr.addTerm(aux_const, v_const);
			
			// Add constraints to LP for set of sampled states
			for(int i =1;i<sample.length;i++){
				if(constraint[i]){																
					s	= (int[]) sample[i].index.clone();
					for(int x=1;x<=Parameters.K;x++){
						if(s[x]>0){
							s[x]--; s[0]--;
							expectedValue(x,s, sample[i].Strat, sample[i].profit[x],exit[i][x],i);
							s[x]++;s[0]++;
						}
					}	
				}												
			}
			// Add Add constraints to LP for set of additional states
			for(int i =1;i<altsample.length;i++){															
				s	= (int[]) altsample[i].index.clone();
				for(int x=1;x<=Parameters.K;x++){
					if(s[x]>0){
						s[x]--; s[0]--;
						expectedValue(x,s, altsample[i].Strat, altsample[i].profit[x],altexit[i][x],sample.length+i-1);
						s[x]++;s[0]++;
					}	
				}												
			}							
			// SALP implementation
			if (!Parameters.relaxbudget){model.addGe(Parameters.budget*aux_const,budexpr);}
			
			// Solving approximate LP
			model.addMinimize(objexpr);
			
			// CPLEX configuration parameters
			// Barrier model recommended when theta>0, however one must limit threads to 1 for large instances due to a bug in CPLEX code
			// relax barrier running parameters for future CPLEX version, after checking the bug is gone.
			// Use default parameters for theta =0, primal algorithm usually gets much faster running times in this case.
			model.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Barrier);	 
		    model.setParam(IloCplex.IntParam.BarCrossAlg, IloCplex.Algorithm.None);
			model.setParam(IloCplex.IntParam.PreDual,1);			
			model.setParam(IloCplex.IntParam.Threads, 1);
		    model.solve();
		}catch(IloException e) {
			System.out.println(e.toString());
		}		
	}
	
	
	/**
	 * Reports memory usage. Motivated by CPLEX bug.
	 * Adjust and un-comment code according to account and system specifics.
	 * @param i iteration count.
	 * @param n inner iteration count.
	 */
	public void reportMemory(int i, int n){
		try{
			System.gc();
			//String s= null;
			Runtime runtime = Runtime.getRuntime();
			FileWriter fstream_mem = new FileWriter("mem.txt", true);
			BufferedWriter out_mem = new BufferedWriter(fstream_mem);			
			double totalmemory = (runtime.totalMemory()-runtime.freeMemory())/Math.pow(2, 20);
			out_mem.write("iteration " + n+ " inner " + i+ " - " + totalmemory + "  Mbs in JVM "  + model.getCplexStatus() + " \r\n");
			System.out.println("iteration " + n+ " inner " + i+ " - " + totalmemory + "  Mbs in JVM");
			
			//Process p = runtime.exec("top -b -u dsaure05 -n 1");
			//p.waitFor();
			//BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));		
			System.out.println("Here is the standard output of the command:\n");
            //while ((s = stdInput.readLine()) != null) {
            //    if(s.contains("java")){out_mem.write(s+ "\r\n");}
            //}
			out_mem.close();
			//p.destroy();
		} catch(Exception e){
			System.out.println(e.toString());
		}
	}
	
	
	/**
	 * Returns the solution to the approximate linear program.
	 * @return Solution to the approximate linear program.
	 */
	public sol_sep getSolution(){
		sol_sep solution = new sol_sep();
		for(int x=1;x<=Parameters.K;x++){
			for(int k=1;k<=Parameters.K;k++){
				for(int n=0;n<Parameters.N;n++){
					
					// Recovering optimal variable value if used
					try{
						solution.r_sep[x][k][n] = model.getValue(v_sep[x][k][n]);
					}catch(IloException e){
						solution.r_sep[x][k][n] = 0;
					}
				}
			}
		}
		try{
			solution.r_const = model.getValue(v_const);
		}catch(IloException e){
			solution.r_const  = 0;
		}
		return solution;
	}
	
	
	/**
	 * Returns value of slack variable in SALP implementation.
	 * @param sample Sampled set of states used to construct the approximating linear program.
	 * @param altsample Additional set of states that complete use of basis functions.
	 * @param constraint indicates if a state in sample set will generate a set of constraints in the LP (limits size of LPs without affecting sampling).
	 * @return slack variable in SALP implementation.
	 */
	public double getSlack(state[] sample, state[] altsample, boolean[] constraint){
		double budslack = 0;
		double totslack = 0;
		for(int i=1;i<sample.length;i++){
			if (constraint[i]){
				for (int x=1;x<=Parameters.K;x++){
					try{
						budslack += sample[i].frequency*model.getValue(v_slack[x][i]);
						totslack += sample[i].frequency;
					}catch(IloException e){}
				}
			}
		}
		for(int i=1;i<altsample.length;i++){
			for (int x=1;x<=Parameters.K;x++){
				try{
					budslack += altsample[i].frequency*model.getValue(v_slack[x][i+sample.length-1]);
					totslack += altsample[i].frequency;
				}catch(IloException e){}
			}
		}
		return (budslack/(totslack+0.00001));
	}
	
	
	/**
	 * Returns exit strategy implied by LP solution.
	 * @param sample Sampled set of states used to construct the approximating linear program.
	 * @param q_firm solution to LP program.
	 * @param i iteration count.
	 * @return exit strategy for states in LP formulation.
	 */
	public double[][] getExit(state[] sample, sol_sep q_firm, int i){
		int[] s;
		double[][] exit = new double[sample.length][Parameters.K+1];
		double[]   aux  = new double[2];
		for(int j=1;j<sample.length;j++){
			if (i>1){
				s = (int[]) sample[j].index.clone();
				for(int k=1;k<=Parameters.K;k++){
					if(s[k]>0){
						s[0]--;s[k]--;
						aux        = StratFromR(k, s, sample[j].Strat, q_firm);
						exit[j][k] = aux[1];
						s[0]++;s[k]++;
					}
				}
			} else{
				exit[j] = (double[]) sample[j].Strat.exit.clone();
			}
		}		
		return exit;
	}
	
	
	/**
	 * Compare exit strategies implied by successive LP solution and decide on convergence of LP heuristic.
	 * @param sample Sampled set of states used to construct the approximating linear program.
	 * @param q_firm solution to previous LP program.
	 * @param q_aux solution to current LP program.
	 * @param i iteration count.
	 * @param out_sum report file handle.
	 * @return true when exit strategies converge.
	 */
	public boolean CompareExit(state[] sample, sol_sep q_firm, sol_sep q_aux, int i, BufferedWriter out_sum){
		int[] s;
		double totalfreq   = 0;
		double Error       = 0;
		double aux         = 0;
		double[] auxArray  = new double[2];
		if (Parameters.entryexit){
			for(int j=1;j<sample.length;j++){
				totalfreq += sample[j].frequency;
				s = (int[]) sample[j].index.clone();
				for(int k=1;k<=Parameters.K;k++){
					if(s[k]>0){
						aux = 0;
						s[0]--;s[k]--;
						if(i==1){
							aux += sample[j].Strat.exit[k]; 	
						}else{
							auxArray   = StratFromR(k, s, sample[j].Strat, q_aux);
							aux       += auxArray[1];
						}
						auxArray   = StratFromR(k, s, sample[j].Strat, q_firm);
						aux -= auxArray[1];
						Error += Math.abs(aux)*sample[j].frequency/Math.abs(Math.max(1,auxArray[1]));
						s[0]++;s[k]++;
					}
				}
			}
		} else{
			totalfreq =1;
		}
		Error = Error/totalfreq;;
		if (Error > Parameters.errTol){
			try{
				System.out.println("Inner: " + Integer.toString(i)+ " -- " + Double.toString(Error) +"\r\n");
				out_sum.write("Inner: " + Integer.toString(i)+ ";" + Double.toString(Error) +"\r\n");
			}catch(Exception e){
				System.out.println(e.toString());
			}
			return false;
		} else{
			return true;
		}
		
	}
	
	
	/**
	 * Construct the constraints for the approximate linear program associated to particular industry state. 
	 * @param x Firm quality state.
	 * @param s Competitors quality state.
	 * @param strategy Incumbent strategy followed by competitors.
	 * @param pi Profit for the incumbent firm.
	 * @param exit incumbent firm's exit strategy.
	 * @param sample_index index of current state within the sampled states.
	 */
	public void expectedValue(int x, int[] s, strategy Strat, double pi, double exit, int sample_index){
		try{
			double[][] coeff = new double[Parameters.K+1][Parameters.N];
			double PnoExit   = Fdist.FdistExit(exit); 
			double Eexit 	 = Fdist.EdistExit(exit);
			// Compute constraints coefficients (Bellman operator part)
			coeff =computeCoeff(s,Strat,1);
			
			// addaptive investment grid
			//*********************************************************
			double[] thisarray = new double[Parameters.localLength];
			double up = Math.min(Parameters.maxInv,Parameters.invArray[Strat.policy[x]]+ Parameters.maxInv/(2*Parameters.reduction));
			double lo = Math.max(0.0,Parameters.invArray[Strat.policy[x]]- Parameters.maxInv/(2*Parameters.reduction));
			if(lo==0.0){
				up = Parameters.maxInv/Parameters.reduction;
			}
			if(up==Parameters.maxInv){
				lo = Parameters.maxInv*(1-1/Parameters.reduction);
			}			
			for(int h=0;h<Parameters.localLength;h++){
				thisarray[h] = (double) (lo + h*(up-lo)/(Parameters.localLength-1));
			}
			//*********************************************************
			for(int h=0;h<Parameters.localLength;h++){
				
				// Write constraint
				IloLinearNumExpr expr = model.linearNumExpr(pi - PnoExit*Parameters.invCost*thisarray[h] + Eexit);
				P.setProb(thisarray[h],x,Parameters.gamma);
				if (x == 1){
					for(int k=1;k<=Parameters.K;k++){
						for(int n=0;n<Parameters.N;n++){
							if (coeff[k][n]>0){
								expr.addTerm((PnoExit*P.Up*Parameters.beta*coeff[k][n]), v_sep[x+1][k][n]);
								expr.addTerm((PnoExit*(P.Down+P.Stay)*Parameters.beta*coeff[k][n]),v_sep[x][k][n]);
							}	
						}
						expr.addTerm(-1.0 , v_sep[x][k][s[k]]);
					}
				} else if (x==Parameters.K){
					for(int k=1;k<=Parameters.K;k++){
						for(int n=0;n<Parameters.N;n++){
							if(coeff[k][n]>0){							
								expr.addTerm((PnoExit*P.Down*Parameters.beta*coeff[k][n]),v_sep[x-1][k][n]);
								expr.addTerm((PnoExit*(P.Up+P.Stay)*Parameters.beta*coeff[k][n]),v_sep[x][k][n]);
							}				
						}
						expr.addTerm(-1.0 , v_sep[x][k][s[k]]);
					}
				} else{
					for(int k=1;k<=Parameters.K;k++){
						for(int n=0;n<Parameters.N;n++){
							if(coeff[k][n]>0){							
								expr.addTerm((PnoExit*P.Up*Parameters.beta*coeff[k][n]),v_sep[x+1][k][n]);
								expr.addTerm((PnoExit*P.Stay*Parameters.beta*coeff[k][n]),v_sep[x][k][n]);
								expr.addTerm((PnoExit*P.Down*Parameters.beta*coeff[k][n]),v_sep[x-1][k][n]);
							}
						}
						expr.addTerm(-1.0 , v_sep[x][k][s[k]]);
					}
				}
				expr.addTerm((PnoExit*Parameters.beta-1.0), v_const);
				expr.addTerm(-1.0, v_slack[x][sample_index]);
				// Add constraint to LP model
				model.addGe(0, expr);												
			}
		}catch(IloException e){
			System.out.println(e.toString());
		}		
	}

	
	/**
	 * Sampling routine. Generates the set of states to be used to formulate the approximate linear program.
	 * @param q_sep  History of solutions to the approximate linear programs.
	 * @param inGamma appreciation factor in transition dynamics.
	 * @param obli probability that a firm will select to use initial (compact) strategy.
	 * @return Set of sampled states to be used to formulate the approximate linear program.
	 */	
	public state[] sampleStates(sol_sep[] q_sep, double inGamma, double obli){		
		int[] u; 
		int[] s;
		boolean repeated;
		int sample_count;
		double rand;
		double Pexit;
		double PINIexit;
		double PINIentry;
		prob PINI = new prob();
		double Pentry;
		strategy Strat = new strategy();
		double[] prof  = new double[Parameters.K+1];
		state[] sample = new state[1];
		
		
		for(int rep =1; rep<=Parameters.nreplic;rep++){
			sample_count = 0;
			s            = (int[]) Parameters.s_0.clone();
			if(Math.IEEEremainder(rep, 100)==0){
				System.out.println("Period " + Integer.toString(rep) + "\r\n");
			}
			for(int n=0; n<Parameters.transcient + Parameters.spacing*Parameters.samplesize;n++){				
				
				prof   = profit.profit(s);
				Strat  = IterativeStrategy(s, q_sep, q_sep.length-1);
				Pentry = Fdist.FdistEntry(Strat.entry); 
				PINIentry = Fdist.FdistEntry(compact.entry);
				if(n>=Parameters.transcient & sample_count>=Parameters.spacing){
					sample_count = 0;
					repeated     = false;
										
					for(int i=1;i<sample.length;i++){
						if (sample[i].isEqual(s)){
							sample[i].increaseFrequency();
							repeated = true;
							break;
						} 
					}
					if (!repeated){
						sample =(state[])resizeArray(sample,sample.length+1);
						sample[sample.length-1] = new state(s);
						sample[sample.length-1].setStrategy(Strat);
						sample[sample.length-1].setProfit(prof);
					}
				}
				sample_count++;
				
				u = (int[]) s.clone();
				for(int k=1; k<=Parameters.K;k++){								
					if(u[k]>0){
						P.setProb(Parameters.invArray[Strat.policy[k]], k,inGamma);
						Pexit = (1.0-Fdist.FdistExit(Strat.exit[k]));
						PINI.setProb(Parameters.invArray[compact.policy[k]], k, inGamma);
						PINIexit = (1.0-Fdist.FdistExit(compact.exit[k]));
						if (k==1){					
							for(int i=1; i<=u[k];i++){
								rand = generator.nextDouble();
								if(rand>=obli){
									rand = generator.nextDouble();
									if(rand <Pexit){
										s[k]--;
										s[0]--;
									}else if(rand < Pexit +(1-Pexit)*P.Up){
										s[k]--;
										s[k+1]++;
									}
								} else{
									rand = generator.nextDouble();
									if(rand <PINIexit){
										s[k]--;
										s[0]--;
									}else if(rand < PINIexit +(1-PINIexit)*PINI.Up){
										s[k]--;
										s[k+1]++;
									}																							
								}
							}
						}else if (k==Parameters.K){
							for(int i=1; i<=u[k];i++){
								rand = generator.nextDouble();
								if(rand>=obli){
									rand = generator.nextDouble();
									if(rand <Pexit){
										s[k]--; 
										s[0]--;
									}else if(rand < Pexit + (1-Pexit)*P.Down){
										s[k]--;
										s[k-1]++;
									}
								}else{
									rand = generator.nextDouble();
									if(rand <PINIexit){
										s[k]--; 
										s[0]--;
									}else if(rand < PINIexit + (1-PINIexit)*PINI.Down){
										s[k]--;
										s[k-1]++;
									}
								}
							}
						}else{
							for(int i=1; i<=u[k];i++){
								rand = generator.nextDouble();
								if(rand>=obli){
									rand = generator.nextDouble();
									if(rand < Pexit){
										s[k]--;
										s[0]--;
									}else if(rand < Pexit + (1-Pexit)*P.Down){
										s[k]--;
										s[k-1]++;
									} else if (rand< Pexit + (1-Pexit)*(P.Down+P.Up)){
										s[k]--;
										s[k+1]++;
									}
								}else{
									rand = generator.nextDouble();
									if(rand < PINIexit){
										s[k]--;
										s[0]--;
									}else if(rand < PINIexit + (1-PINIexit)*PINI.Down){
										s[k]--;
										s[k-1]++;
									} else if (rand< PINIexit + (1-PINIexit)*(PINI.Down+PINI.Up)){
										s[k]--;
										s[k+1]++;
									}									
								}
							}
						}
					}
				}
				if(u[0]<Parameters.N){
					for (int a=1;a<=Parameters.N-u[0];a++){
						rand = generator.nextDouble();
						if(rand>=obli){
							rand = generator.nextDouble();
							if (rand < Pentry){
								s[Parameters.x_e]++;s[0]++;
							}
						}else{
							rand = generator.nextDouble();
							if (rand < PINIentry){
								s[Parameters.x_e]++;s[0]++;
							}							
						}
					}
				}
			}
		}
		System.out.println("Simulation finished!\r\n");
		System.out.println("Sampled states =" + Integer.toString(sample.length-1)+ "\r\n");	
		return sample;
	}
	
	
	/**
	 * Add states to set of sampled states so that all basis functions appear in the objective function of the LP.  
	 * @param q_sep  History of solutions to the approximate linear programs.
	 * @param sample Sampled set of states used to construct the approximating linear program.
	 * @return Set of sampled states to be used to formulate the approximate linear program.
	 * @param constraint indicates if a state in sample set will generate a set of constraints in the LP (limits size of LPs without affecting sampling).
	 */
	public state[] sampleRest(sol_sep[] q_sep, state[] sample, boolean[] constraint){	
		
		state[] altsample = new state[1];
		int[] s= new int[Parameters.K+1];
		boolean[][][] Noneed = new boolean[Parameters.K+1][Parameters.K+1][Parameters.N];
		for (int i=1;i<sample.length;i++){
			if (constraint[i]){
				s = sample[i].index;
				for(int x=1;x<=Parameters.K;x++){
					if (s[x]>0){
						s[x]--;s[0]--;
						for(int k=1;k<=Parameters.K;k++){
							Noneed[x][k][s[k]]=true;					
						}
						s[x]++;s[0]++;
					}
				}
			}
		}
		int[] newind;
		int Nnodes;
		strategy Strat = new strategy();
		double[] prof  = new double[Parameters.K+1];
		boolean taken;
		for(int x=1;x<=Parameters.K;x++){			
			for(int k=1;k<=Parameters.K;k++){
				Nnodes = nodes[k].length;
				for(int n=0;n<Nnodes;n++){
					taken = false;
					if(n>0){						
						for(int w=nodes[k][n-1]+1; w<=nodes[k][n];w++){
							if (Noneed[x][k][w] == true){taken =true; break;}						
						}
					}else {
						for(int w=0; w<=nodes[k][n];w++){
							if (Noneed[x][k][w] == true){taken =true; break;}						
						}
					}
					if (!taken){
						if (CoefTaken(altsample,x,k,n)){taken = true;}
					}
					if (!taken){
						newind 		 = new int[Parameters.K+1];
						newind[x]    = 1;
						newind[k]    += nodes[k][n];
						newind[0]    = 1+nodes[k][n];
						
						prof   		= profit.profit(newind);
						Strat  		= IterativeStrategy(newind, q_sep, q_sep.length-1);
						altsample 	= (state[])resizeArray(altsample,altsample.length+1);
						altsample[altsample.length-1] = new state(newind);								
						altsample[altsample.length-1].setStrategy(Strat);
						altsample[altsample.length-1].setProfit(prof);							
					}
				}
			}		
		}					
		return altsample;
	}
		
	/**
	 * Simulates performance of an strategy against the incumbent strategy. Used to compute the convergence criteria.
	 * @param sample Sampled set of states used to construct the approximating linear program.
	 * @param q_sep History of solutions to the approximate linear programs.
	 * @param mode "Base" to test incumbent strategy, "R" to test approximate best response.
	 * @param constraint indicates if a state in sample set will generate a set of constraints in the LP (limits size of LPs without affecting sampling). 
	 * @return Simulated expected value function for a firm using a strategy assuming everyone else is using the incumbent strategy.
	 */
	public double ComputeV(state[] sample, sol_sep[] q_sep, String  mode, boolean[] constraint){
		double Value = 0;
		double auxV  = 0;
		int[] s;
		int[] u;
		int x;
		double rand;
		double Pexit;
		double Pentry;
		double weight;
		double total;
		strategy industry_strat = new strategy();
		strategy firm_strat	    = new strategy();
		double[] auxStrat   	= new double[2];
		double[] prof = new double[Parameters.K+1];
		
		// Go over all sampled states and simulate expected continuation value
		for(int i=1;i<sample.length;i++){
			// do n replications per state to simulate expected value			
			if(constraint[i]){
			total = 0;
			for(int k=1;k<=Parameters.K;k++){
				if(sample[i].index[k]>0){total++;}
			}
			for(int k=1;k<=Parameters.K;k++){
				if(sample[i].index[k]>0){
					weight 	= (double) sample[i].frequency/total;
					auxV 	= 0;					
					for(int n=1;n<=Parameters.Vreplic;n++){
						//setting initial state						
						s = (int[]) sample[i].index.clone();
						s[k]--; s[0]--;
						x = k;
						// simulation scheme
						for(int t=1;t<=Parameters.Vperiods;t++){
							s[x]++;s[0]++;
							prof   = profit.profit(s);
							industry_strat  = IterativeStrategy(s, q_sep, q_sep.length-2);				// compute industry strategy
							s[x]--;s[0]--;
							if(mode.equals("Base")){
								firm_strat.setStrategy(industry_strat);											// compute firm strategy	
							} else{
								auxStrat = StratFromR(x, s, industry_strat, q_sep[q_sep.length-1]);
								firm_strat.policy[x] = (int) auxStrat[0];
								firm_strat.exit[x]   = auxStrat[1];
							}
							// Contribution of this periods to computation of V
							auxV += prof[x]*Math.pow(Parameters.beta, t-1);
																					
							// Determine next competitors state
							u = (int[]) s.clone();
							Pentry = (double) Fdist.FdistEntry(industry_strat.entry);
							for(int q=1; q<=Parameters.K;q++){					 
								if(u[q]>0){
									Pexit = (1 - Fdist.FdistExit(industry_strat.exit[q]));
									P.setProb(Parameters.invArray[industry_strat.policy[q]], q,Parameters.gamma);
									if (q==1){					
										for(int j=1; j<=u[q];j++){
											rand = generator.nextDouble();
											if(rand< Pexit){
												s[q]--;s[0]--;
											}else if(rand < Pexit + (1-Pexit)*P.Up){
												s[q]--;
												s[q+1]++;
											}				
										}
									}else if (q==Parameters.K){
										for(int j=1; j<=u[q];j++){
											rand = generator.nextDouble();
											if(rand < Pexit){
												s[q]--;s[0]--;
											}else if(rand < Pexit + (1-Pexit)*P.Down){
												s[q]--;
												s[q-1]++;
											}												
										}
									}else{
										for(int j=1; j<=u[q];j++){
											rand = generator.nextDouble();
											if(rand < Pexit){
												s[q]--;s[0]--;
											}else if(rand < Pexit + (1-Pexit)*P.Down){
												s[q]--;
												s[q-1]++;
											} else if (rand< Pexit + (1-Pexit)*(P.Down+P.Up)){
												s[q]--;
												s[q+1]++;
											}												
										}
									}
								}
							}
							if(u[0]<Parameters.N-1){
								for(int a=1;a<=Parameters.N-1-u[0];a++){
									rand = generator.nextDouble();
									if(rand < Pentry){
										s[Parameters.x_e]++; s[0]++;
									}
								}
							}
							rand  = generator.nextDouble();
							Pexit = (1 - Fdist.FdistExit(firm_strat.exit[x]));
							
							if (rand<Pexit){
								// Collect scrap value and exit industry
								auxV -= Math.pow(Parameters.beta, t-1)*Math.log(rand)/Parameters.epsn;
								break;
							} else {
								// pay investment and advance to next period
								auxV -= Parameters.c*Math.pow(Parameters.beta, t-1)*Parameters.invArray[firm_strat.policy[x]];
								rand  = generator.nextDouble();
								P.setProb(Parameters.invArray[firm_strat.policy[x]], x,Parameters.gamma);
								if (x==1){
									if (rand<P.Up){x++;}
								}else if (x==Parameters.K){
									if (rand<P.Down){x--;}
								}else{
									if (rand<P.Down){
										x--;
									}else if(rand > P.Down+P.Up){
										x++;
									}
								}									
							}
						}						
					}
					Value += weight*auxV/((double) Parameters.Vreplic);
				}
			}
			}
		}
		return Value;	
}	
	
	
	/**
	 * Oracle M: Computes incumbent strategy based on previous solutions to the approximate linear programs.
	 * @param s Incumbent industry state.
	 * @param q_sep History of solutions to the approximate linear programs.
	 * @param period Iteration for which one wants to reconstruct approximate best response.
	 * @return Approximate best response strategy for specified iteration.
	 */
	public strategy IterativeStrategy(int[] s, sol_sep[] q_sep, int period){ 	
		strategy Strat    = new strategy();
		strategy newStrat = new strategy();											
		int[] u           = new int[Parameters.K+1];
		u                 = (int[]) s.clone();
		double[] auxStrat = new double[2];
		Strat.setStrategy(compact);		
		for(int t=1;t<=period;t++){
			for(int x=1;x<=Parameters.K;x++){
				newStrat.policy[x]= 0;
				newStrat.exit[x]  = 0.0;
				if (u[x]>0){
					u[x]--;u[0]--;
					
					// Get best response strategy for period t
					auxStrat           = StratFromR(x, u, Strat, q_sep[t]);
					newStrat.policy[x] = (int) auxStrat[0];
					newStrat.exit[x]   = auxStrat[1];
					u[x]++;u[0]++;
				}
			}
			if(u[0]<Parameters.N){
				newStrat.entry = EntryFromR(u, Strat, q_sep[t]);	
			}else {
				newStrat.entry = 0;
			}
			Strat.setStrategy(newStrat);
		}	
		return Strat;
	}		
	
	
	/**
	 * Computes greedy response against a solution for the approximate linear program.  
	 * @param x Firm state.
	 * @param s Competitors state.
	 * @param q_sep  Solution for the approximate linear program.
	 * @return Greedy response against a solution for the approximate linear program.
	 */
	public double[] StratFromR(int x, int[] s, strategy Strat, sol_sep q_sep){
		double[][] coeff = new double[Parameters.K+1][Parameters.N];
		double[] values  = new double[Parameters.invLength];
		int pol          = 0;
		double[] PolEx   = new double[2];
		double max; 
		
		// Compute Bellman operator coefficients
		coeff            = computeCoeff(s, Strat,1);
		for(int h=0;h<Parameters.invLength;h++){
			
			// Reconstructing LP constraint
			values[h] = -Parameters.invCost*Parameters.invArray[h];							
			P.setProb(Parameters.invArray[h],x,Parameters.gamma);
			if (x == 1){
				for(int k=1;k<=Parameters.K;k++){
					for(int n=0;n<Parameters.N;n++){
						if (coeff[k][n]>0){
							values[h] += P.Up*Parameters.beta*coeff[k][n]*q_sep.r_sep[x+1][k][n];
							values[h] += (P.Down+P.Stay)*Parameters.beta*coeff[k][n]*q_sep.r_sep[x][k][n];										
						}	
					}
				}
				values[h] += Parameters.beta*q_sep.r_const;
			} else if (x==Parameters.K){
				for(int k=1;k<=Parameters.K;k++){
					for(int n=0;n<Parameters.N;n++){
						if(coeff[k][n]>0){
							values[h] += P.Down*Parameters.beta*coeff[k][n]*q_sep.r_sep[x-1][k][n];
							values[h] += (P.Up+P.Stay)*Parameters.beta*coeff[k][n]*q_sep.r_sep[x][k][n];
						}
					}
				}
				values[h] += Parameters.beta*q_sep.r_const;
			} else{
				for(int k=1;k<=Parameters.K;k++){
					for(int n=0;n<Parameters.N;n++){
						if(coeff[k][n]>0){
							values[h] += P.Up*Parameters.beta*coeff[k][n]*q_sep.r_sep[x+1][k][n];
							values[h] += P.Stay*Parameters.beta*coeff[k][n]*q_sep.r_sep[x][k][n];
							values[h] += P.Down*Parameters.beta*coeff[k][n]*q_sep.r_sep[x-1][k][n];
						}
					}
				}
				values[h] += Parameters.beta*q_sep.r_const;
			}						
		}
		
		// Looking for greedy maximizing investment decision 
		max = values[0]; 
		for(int h=1;h<Parameters.invLength;h++){
			if (values[h]>max){				
				max = values[h];		
				pol = h;			
			}
		}
		PolEx[0] = pol;
		PolEx[1] = Math.max(max, 0);	// exponential distribution scrap value
		return PolEx;
	}	
	
	
	/**
	 * Computes greedy response against a solution for the approximate linear program.  
	 * @param s Industry state.
	 * @param strategy industry strategy. 
	 * @param q_sep Solution for the approximate linear program.
	 * @return New entry threshold based on value function approximation.
	 */
	public double EntryFromR(int[] s, strategy Strat, sol_sep q_sep){
		double[][] coeff = new double[Parameters.K+1][Parameters.N];
		double en	= 0;
		
		if(s[0]<Parameters.N){
			
			// Compute Bellman operator coefficients
			coeff   = computeCoeff(s, Strat,1);// Modification August 2010 Firm assumes will enter
			en 		= Parameters.beta*q_sep.r_const;
			// Looking for implied entry strategy
			for(int k=1;k<=Parameters.K;k++){
				for(int n=0;n<Parameters.N;n++){
					if (coeff[k][n]>0){
						en += Parameters.beta*coeff[k][n]*q_sep.r_sep[Parameters.x_e][k][n];										
					}	
				}
			}
		}
		en =Math.max(en, 0); // exponential distribution setup cost		
		return en;
	}
	
	
	/**
	 * Computes coefficients necessary to calculate the Bellman operator.
	 * @param s Industry state.  
	 * @param strategy Strategy.
	 * @return Array of coefficients associated to each base function. 
	 */
	public double[][] computeCoeff(int[] s, strategy Strat, int out){
		double[][] coeff   = new double[Parameters.K+1][Parameters.N+1-out];
		double[][][] multi = new double[Parameters.K+1][][];
		double Pexit;
		double Pentry = Fdist.FdistEntry(Strat.entry); 
		double[] multiEntry = getmultiEntry(Pentry,Parameters.N-s[0]-out);
		for(int k=1;k<=Parameters.K;k++){
			multi[k] = new double[s[k]+1][2];
			P.setProb(Parameters.invArray[Strat.policy[k]],k,Parameters.gamma);
			Pexit = (1.0-Fdist.FdistExit(Strat.exit[k]));
			multi[k] = getmulti(s[k],k,Pexit); 
			for(int n=0;n<Parameters.N;n++){
				coeff[k][n]=0;
			}
		}
		
		// Computing coefficients taking advantage of transition probability structure (only jumps to adjacent states)
		for(int k=1;k<=Parameters.K;k++){
			if(k==Parameters.x_e){
				if(k==1){
					for(int n=0;n<=Parameters.N-s[0]-out;n++){
						for(int stay=0;stay<=s[k];stay++){
							for(int down=0;down<=s[k+1];down++){
								coeff[k][stay+down+n] += multi[k][stay][1]*multi[k+1][down][0]*multiEntry[n];
							}
						}
					}	
				}else if(k==Parameters.K){
					for(int n=0;n<=Parameters.N-s[0]-out;n++){
						for(int up=0;up<=s[k-1];up++){
							for(int stay=0;stay<=s[k];stay++){
								coeff[k][up+stay+n] += multi[k-1][up][2]*multi[k][stay][1]*multiEntry[n];
							}
						}
					}	
				}else {
					for(int n=0;n<=Parameters.N-s[0]-out;n++){
						for(int up=0;up<=s[k-1];up++){
							for(int stay=0;stay<=s[k];stay++){
								for(int down=0;down<=s[k+1];down++){
									coeff[k][up+stay+down+n] += multi[k-1][up][2]*multi[k][stay][1]*multi[k+1][down][0]*multiEntry[n];
								}	
							}
						}
					}	
				}
			}else{
				if(k==1){
					for(int stay=0;stay<=s[k];stay++){
						for(int down=0;down<=s[k+1];down++){
							coeff[k][stay+down] += multi[k][stay][1]*multi[k+1][down][0];
						}
					}
				}else if(k==Parameters.K){
					for(int up=0;up<=s[k-1];up++){
						for(int stay=0;stay<=s[k];stay++){
							coeff[k][up+stay] += multi[k-1][up][2]*multi[k][stay][1];
						}
					}
				}else {	
					for(int up=0;up<=s[k-1];up++){
						for(int stay=0;stay<=s[k];stay++){
							for(int down=0;down<=s[k+1];down++){
								coeff[k][up+stay+down] += multi[k-1][up][2]*multi[k][stay][1]*multi[k+1][down][0];
							}	
						}
					}
				}
			}
		}	
		return coeff;	
	}
	
	
	/**
	 * Compute multinomial probabilities for possible transitions.
	 * @param n Number of firms.
	 * @param k Quality of the firm.
	 * @param exit Exit probability.
	 * @return  Array of probabilities for possible transitions.
	 */
	public double[][] getmulti(int n,int k, double Pexit){
		double aux;
		double[][] prob = new double[n+1][3];
		if(k==1){
			if(((1-Pexit)*P.Up)==1){
				prob[n][2] = 1;
				for(int i=0;i<n;i++){
					prob[i][2] = 0;
				}
			} else if(((1-Pexit)*P.Up)==0){
				prob[0][2] = 1;
				for(int i=n;i>0;i--){
					prob[i][2] = 0;
				}
			} else{
				prob[0][2] = Math.pow((1-(1-Pexit)*P.Up),n);
				for(int i=1;i<=n;i++){
					aux = (double)(n-i+1)/i;
					prob[i][2] = prob[i-1][2]*aux*(1-Pexit)*P.Up/(1-(1-Pexit)*P.Up);
				}
			}
			if(((1-Pexit)*(P.Stay+P.Down))==1){
				prob[n][1] = 1;
				for(int i=0;i<n;i++){
					prob[i][1] = 0;
				}
			} else if(((1-Pexit)*(P.Stay+P.Down))==0){
				prob[0][1] = 1;
				for(int i=n;i>0;i--){
					prob[i][1] = 0;
				}
			} else{
				prob[0][1] = Math.pow((1-(1-Pexit)*(P.Stay+P.Down)),n);
				for(int i=1;i<=n;i++){
					aux = (double)(n-i+1)/i;
					prob[i][1] = prob[i-1][1]*aux*(1-Pexit)*(P.Stay+P.Down)/(1-(1-Pexit)*(P.Stay+P.Down));
				}
			}
		} else if (k==Parameters.K){
			if(((1-Pexit)*P.Down)==1){
				prob[n][0] = 1;
				for(int i=0;i<n;i++){
					prob[i][0] = 0;
				}
			} else if(((1-Pexit)*P.Down)==0){
				prob[0][0] = 1;
				for(int i=n;i>0;i--){
					prob[i][0] = 0;
				}
			} else{
				prob[0][0] = Math.pow((1-(1-Pexit)*P.Down),n);
				for(int i=1;i<=n;i++){
					aux = (double)(n-i+1)/i;
					prob[i][0] = prob[i-1][0]*aux*(1-Pexit)*P.Down/(1-(1-Pexit)*P.Down);
				}
			}
			if(((1-Pexit)*(P.Stay+P.Up))==1){
				prob[n][1] = 1;
				for(int i=0;i<n;i++){
					prob[i][1] = 0;
				}
			} else if(((1-Pexit)*(P.Stay+P.Up))==0){
				prob[0][1] = 1;
				for(int i=n;i>0;i--){
					prob[i][1] = 0;
				}
			} else{
				prob[0][1] = Math.pow((1-(1-Pexit)*(P.Stay+P.Up)),n);
				for(int i=1;i<=n;i++){
					aux = (double)(n-i+1)/i;
					prob[i][1] = prob[i-1][1]*aux*(1-Pexit)*(P.Stay+P.Up)/(1-(1-Pexit)*(P.Stay+P.Up));
				}
			}
		} else {
			if(((1-Pexit)*P.Down)==1){
				prob[n][0] = 1;
				for(int i=0;i<n;i++){
					prob[i][0] = 0;
				}
			} else if(((1-Pexit)*P.Down)==0){
				prob[0][0] = 1;
				for(int i=n;i>0;i--){
					prob[i][0] = 0;
				}
			} else{
				prob[0][0] = Math.pow((1-(1-Pexit)*P.Down),n);
				for(int i=1;i<=n;i++){
					aux = (double)(n-i+1)/i;
					prob[i][0] = prob[i-1][0]*aux*(1-Pexit)*P.Down/(1-(1-Pexit)*P.Down);
				}
			}
			if(((1-Pexit)*P.Up)==1){
				prob[n][2] = 1;
				for(int i=0;i<n;i++){
					prob[i][2] = 0;
				}
			} else if(((1-Pexit)*P.Up)==0){
				prob[0][2] = 1;
				for(int i=n;i>0;i--){
					prob[i][2] = 0;
				}
			} else{
				prob[0][2] = Math.pow((1-(1-Pexit)*P.Up),n);
				for(int i=1;i<=n;i++){
					aux = (double)(n-i+1)/i;
					prob[i][2] = prob[i-1][2]*aux*(1-Pexit)*P.Up/(1-(1-Pexit)*P.Up);
				}
			}
			if(((1-Pexit)*P.Stay)==1){
				prob[n][1] = 1;
				for(int i=0;i<n;i++){
					prob[i][1] = 0;
				}
			} else if(((1-Pexit)*P.Stay)==0){
				prob[0][1] = 1;
				for(int i=n;i>0;i--){
					prob[i][1] = 0;
				}
			} else{
				prob[0][1] = Math.pow((1-(1-Pexit)*P.Stay),n);
				for(int i=1;i<=n;i++){
					aux = (double)(n-i+1)/i;
					prob[i][1] = prob[i-1][1]*aux*(1-Pexit)*P.Stay/(1-(1-Pexit)*P.Stay);
				}
			}	
		}
		return prob;	
	}

	/**
	 * Compute entry probabilities as a function of number of entrants.
	 * @param Pentry Individual firm's entry probability.
	 * @param n number of potential entrants.
	 * @return distribution of possible entrants.
	 */
	public double[] getmultiEntry(double Pentry, int n){
		double[] prob = new double[n+1];
		double aux;
		if(Pentry==1){
			prob[n] = 1;
			for(int i=0;i<n;i++){
				prob[i] = 0;
			}
		} else if(Pentry==0){
			prob[0] = 1;
			for(int i=n;i>0;i--){
				prob[i] = 0;
			}
		} else{
			prob[0] = Math.pow((1-Pentry),n);
			for(int i=1;i<=n;i++){
				aux = (double)(n-i+1)/i;
				prob[i] = prob[i-1]*aux*Pentry/(1-Pentry);
			}
		}
		return prob;
	}
	
		
	/**
	 * Defines the transition dynamics for the model (subclass).    
	 * @version 2.0
	 */
	public static class prob{
		/** Probability of increasing quality level*/double Up;
		/** Probability of maintaining quality level*/double Stay;
		/** Probability of decreasing quality level*/double Down;
		
		
		/**
		 * Initialize the transition probability structure. 
		 */
		public prob(){
			Up = 0; Down = 0; Stay = 0;
		}
		
		
		/**
		 * Sets up transition probabilities depending on the investment decision and the current quality of the firm.
		 * This is the model used in Farias, Saure and Weintraub (2011).
		 * @param investment Investment decision.
		 * @param level Quality of the firm.
		 * @param inGamma Appreciation factor
		 */
		public void setProb(double investment, int level, double inGamma){
			double aux = 1 + Parameters.alpha*investment;
			if (level==1){
				Up   = Parameters.alpha*investment/aux;
				Stay = 1-Up;				
				Down = 0;				
				Up = Up*(1-inGamma)+inGamma; Stay = Stay*(1-inGamma); Down = Down*(1-inGamma);
			} else if (level == Parameters.K){
				Up   = 0;
				Stay = ((1-Parameters.delta)+Parameters.alpha*investment)/aux;
				Down = 1 - Stay;
			} else {
				Up   = ((1-Parameters.delta)*Parameters.alpha*investment)/aux;
				Stay = ((1-Parameters.delta)+Parameters.delta*Parameters.alpha*investment)/aux;
				Down = 1 - Up - Stay;
				Up = Up*(1-inGamma)+inGamma; Stay = Stay*(1-inGamma); Down = Down*(1-inGamma);
			}				
		}
	}
	
	
	/**
	* Reallocates an array with a new size, and copies the contents.
	* of the old array to the new array.
	* @param oldArray  the old array, to be reallocated.
	* @param newSize   the new array size.
	* @return          A new array with the same contents.
	* @author          Christian d'Heureuse.
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
	 * Index for states when dealing when small instances. 
	 * @param j             Index of the state.
	 * @param Combinatorial Combinatorial numbers (pre-computed for efficiency).
	 * @return              Industry state associated to index j.
	 */
	public int[] inverseIndex(int j, double[][] Combinatorial){
		int[] s = new int[Parameters.K+1];
		int[] A = new int[Parameters.K+1];
		int aux;
		for(int k=Parameters.K;k>0;k--){
			A[k] = 0;
			aux  =0;
			for(int h=k+1;h<=Parameters.K;h++){
				if(A[h]>0) aux += Combinatorial[A[h]-1][h];					
			}
			if (j+1-aux > Combinatorial[0][k]){
				for(int h=1; h <=Parameters.N; h++){
					if (Combinatorial[h][k]>=j+1-aux){
						A[k] = h;
						break;
					}
				}
			}	
		}
		s[0] = A[Parameters.K];
		for(int k=1; k<Parameters.K; k++){
			s[k] = A[Parameters.K-k+1] -A[Parameters.K-k]; 
		}
		s[Parameters.K] = A[1];	
		return s;
	}
	
	
	/**
	 * Returns index of a state when dealing with small instances.
	 * @param s             Industry state.
	 * @param Combinatorial Combinatorial numbers (pre-computed for efficiency).
	 * @return              Index associated to industry state s.
	 */
	public int returnIndex(int[] s, double[][] Combinatorial ){
		int j;
		int aux = s[0]-1;
		if (s[0] == 0){ 
			return 0;
		} else if(s[Parameters.K]==Parameters.N-1){
			return (int)Combinatorial[Parameters.N-1][Parameters.K]-1;
		} else{
			j = (int)Combinatorial[aux][Parameters.K]+1;
			for(int k=1;k<Parameters.K;k++){
				aux -= s[k];
				if (aux < 0) break;
				j += (int)Combinatorial[aux][Parameters.K-k]; 
			}
			return j-1;				
		}
	}	
	
	
	/**
	 * Returns all states in state space, and their strategies.
	 * @param q_sep  History of solutions to the approximate linear programs.
	 * @param frequency sample frequency associated (uniformly) to each state.
	 * @return set of possible states in the state space with their strategies.
	 */
	public state[] allStates(sol_sep[] q_sep, double[] frecuency){		
		try{
			int[] s;
			strategy Strat  = new strategy();
			double[] prof   = new double[Parameters.K+1];
			int state_count = 1;
			
			// Generate combinatorial numbers
			double[][] Combinatorial = new double[Math.max(Parameters.N+1,Parameters.K+1)][Math.max(Parameters.N+1,Parameters.K+1)];
			for(int k=0;k<Math.max(Parameters.N+1,Parameters.K+1);k++){
				Combinatorial[k][0] = 1;
				for(int n=1;n<Math.max(Parameters.N+1,Parameters.K+1);n++){
					Combinatorial[k][n] = (int) Combinatorial[k][n-1]*(n+k)/n;
				}
			}
			state[] sample = new state[(int) Combinatorial[Parameters.N][Parameters.K]+1];
			
			// Go over all possible values for pair firm state-competitors state				
			for(int j = 0;j<(int) Combinatorial[Parameters.N][Parameters.K];j++){					
				s      = inverseIndex(j,Combinatorial);
				
				// Compute profit and strategy vectors for incumbent industry state
				prof   = profit.profit(s);												
				Strat  = IterativeStrategy(s, q_sep, q_sep.length-1);	
				sample[state_count] = new state(s);
				sample[state_count].setStrategy(Strat);
				sample[state_count].setProfit(prof);
				if (!Parameters.unifCvector){
					sample[state_count].setFrequency(frecuency[j]);
				}
				state_count++;
			}
			return sample;
		} catch(Exception e){
			System.out.println(e.toString());
			return null;
		}
	}
	
	
	/**
	 * Limits size of sampled states used in LP formulation.
	 * @param sample Sampled set of states used to construct the approximating linear program.
	 * @param maxSize maximum number of states used in LP formulation.
	 * @return binary vector indicating if associated state is to be used when formulating the LP.
	 */
	static boolean[] getConstraints(state[] sample, int maxSize){
		int sizeSample = sample.length;
		boolean[] Constraint= new boolean[sizeSample];
		double max;
		int aux;
		if (sizeSample-1<=maxSize){
			for(int i=1;i<sizeSample;i++){
				Constraint[i] = true;
			}
		} else{
			double[] freq = new double[sizeSample];
			for(int i=1;i<sizeSample;i++){
				Constraint[i] = false;
				freq[i] 	  = sample[i].frequency;  
			}
			for (int n=1;n<=maxSize;n++){
				aux = 0;
				max = -1;	
				for(int i=1;i<sizeSample;i++){
					if (freq[i]>=max){
						max = freq[i];
						aux = i;
					}
				}
				Constraint[aux] = true;
				freq[aux] 		= -1;
			}
		}
		return Constraint;
	}
	
	
	/**
	 * Indicates whether a particular basis function is being covered by the set of sampled states.
	 * @param sample Sampled set of states used to construct the approximating linear program.
	 * @param x firm's quality state.
	 * @param k basis function's quality state.
	 * @param n number of firms in quality state k. 
	 * @return true if basis function is being covered by the sample.
	 */
	public boolean CoefTaken(state[] sample,int x,int k,int n){
		boolean auxbool = false;
		int lowbound  = -1;
		int upbound   = nodes[k][n];
		if (n>0){lowbound = nodes[k][n-1];}
		for(int i=1;i<sample.length;i++){
			if(sample[i].index[x]>0){
				if(x==k){
					if((sample[i].index[k]<=upbound+1) && (sample[i].index[k]>lowbound+1)){
						auxbool = true;
						break;
					}
				}else{
					if((sample[i].index[k]<=upbound) && (sample[i].index[k]>lowbound)){
						auxbool = true;
						break;
					}
				}
			}
		}
		return auxbool;
	}
		

	/**
	 * Reduces size of approximating architecture.
	 */
	public void reduceBaseFunctions(){
		
		System.out.println("Using new architecture!");
		IloLinearNumExpr pwlinear;
		int lobound;
		int upbound;
		for(int x=1;x<=Parameters.K;x++){
			for(int k=1;k<=Parameters.K;k++){
				for(int n=0;n<nodes[k].length;n++){
					lobound = -1;
					if (n>0){lobound = nodes[k][n-1];}
					upbound = nodes[k][n];
					if (upbound-lobound>1){
						for(int i=lobound+1;i<upbound;i++){
							try{
								pwlinear = null;
								pwlinear = model.linearNumExpr();
								pwlinear.addTerm(i-lobound, v_sep[x][k][lobound]);
								pwlinear.addTerm(lobound-upbound, v_sep[x][k][i]);
								pwlinear.addTerm(i-lobound, v_sep[x][k][upbound]);		            
								model.addEq(0.0, pwlinear);				
							}catch(IloException e) {
								System.out.println(e.toString());					
							}
						}
					}			
				}
			}
		}		
	}
	
	/**
	 * Sets coarseness of approximating architecture.
	 * customize for capacity competition model with 40 firms. Modify according to model specifics or set reducearch=false.
	 */
	public void setExeNodes(){
		int[] coarse1   = {0,1,2,3,5,7,8,10,12,15,19,39};
		int[] coarse2   = {0,1,2,3,4,5,6,12,39};	

		int[] Lexe1  = {1,2,3,4,5};
		int[] Lexe2  = {6,7,8,9,10};
		
		for(int i=0;i<Lexe1.length;i++){
			nodes[Lexe1[i]] = coarse1.clone();
		}
		for(int i=0;i<Lexe2.length;i++){
			nodes[Lexe2[i]] = coarse2.clone();
		}		
	}
	
	
	/**
	 * Loads initial strategy from compact.txt file.
	 * first line contains entry threshold.
	 * next K lines contain investment decision for each quality level.
	 * next K lines contain exit thresholds for each quality level.
	 */
	public  void loadCompact(){
		try{
			FileInputStream fstream = new FileInputStream("compact.txt");
  		    //FileInputStream fstream = new FileInputStream("C:\\workspace\\src\\SALPMPE\\compact.txt"); 		 	// uncomment and modify if running locally			
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;
			int count = 0;
			double aux;
			while ((strLine = br.readLine()) != null){
				if(count == 0){
					compact.entry = Double.parseDouble(strLine);
				}else if (count <=Parameters.K){
					aux = Double.parseDouble(strLine);
					compact.policy[count] = (int) Math.min(Math.max(Math.floor(aux/Parameters.di),0),Parameters.invLength-1);
				}else{
					compact.exit[count-Parameters.K] = Double.parseDouble(strLine);
				}
				count++;
		    }
		    in.close();
		    System.out.println("compact strategy read!\r\n");
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
	}
		
}