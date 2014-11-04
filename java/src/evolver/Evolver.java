package evolver;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import javax.swing.text.DateFormatter;

import com.sun.org.apache.xerces.internal.impl.dv.dtd.NMTOKENDatatypeValidator;

import evaluators.GRNGenomeEvaluator;
import evaluators.GRNSinusEvaluator;
import evaluators.GRNStableMediumOutputEvaluator;
import evaluators.IntertwinedSpiralsSubsequence;
import evaluators.IntertwinedSpiralsVariableTime;
import evaluators.MichalSignalProcessExp1;
import evaluators.DoublingFrequencyEvaluator;
import evaluators.MichalSignalProcessExp1;
import evaluators.MichalSignalProcessExp2;
import evaluators.MichalSignalProcessExp3;
import evaluators.NMMultirobotCCP;
import evaluators.PingPrediction;
import evaluators.RadbotEvaluator;
import evaluators.RetinaEvaluator;
import evaluators.CoverageControl;
import evaluators.IntertwinedSpirals;
import evaluators.IntertwinedSpiralsVariableTime;
import evaluators.WeatherPrediction;
import grn.GRNProtein;
import operators.GRNAddGeneMutationOperator;
import operators.GRNAligningCrossoverOperator_ParentCountProb;
import operators.GRNAligningCrossoverOperator_v1;
import operators.GRNAligningCrossoverOperator_v1b;
import operators.GRNAligningCrossoverOperator_v2;
import operators.GRNCrossoverOperator;
import operators.GRNDeleteGeneMutationOperator;
import operators.GRNGeneMutationOperator;
import operators.GRNMutationOperator;
import operators.GRNOnePointCrossoverOperator;

public class Evolver {
	protected int currentGen=0;

	ArrayList<GRNCrossoverOperator> crossovers = new ArrayList<GRNCrossoverOperator>();
	ArrayList<GRNMutationOperator> mutations = new ArrayList<GRNMutationOperator>();
	public static double mutationProbabilitySum=-1;

	ArrayList<Species> species = new ArrayList<Species>();
	

	public Random rng;
	public int randomSeed=(int)(new Date().getTime());
	public int maxNumGen=Integer.MAX_VALUE;
	public GRNGenomeEvaluator evaluator=null;
	public int populationSize = 500;
	public double speciationThreshold = 0.15;
	public double crossoverRate = 0.6;
	public double mutationRate = 0.3;//10714;
	public int speciesMinSize=15;
	public int duplicateInit = (int)(speciesMinSize*2);
	public int nMaxTries=10;
	public int tournamentSize=3;
	public boolean saveAllPopulation=false;
	public String folderName;
	public String experienceName="SignalProcessingExp2";
	public boolean greatInit=true;
	public boolean greatActivated=true;
	public int representativeMethod=1;// 0=keep the originator of the species, 1=random from prev gen, 2=best from prev gen
	public long uuid;
	public int speciesSizeAdjustingMethod=2;
	public boolean useDynamicsGenomeForDistance=true;

	public static int verbose = 0;// 0 is silent
	
	private void adjustSpeciationThresholds() {
		// Adjusting speciation threshold to regulate the sizes
		double speciesAvgSize=populationSize/species.size();
		for (Species s:species) {
			if (s.size()<speciesAvgSize) {
                s.speciationThreshold=Math.min(0.5,s.speciationThreshold+0.01);//speciationThreshold/10;
            } else {
                s.speciationThreshold=Math.max(0.03, s.speciationThreshold-0.01);//speciationThreshold/10;
            }
		}
	}
	
	public void run() {
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		uuid = System.nanoTime();
		folderName=experienceName+"/run_"+uuid;//dateFormat.format(new Date());
		try {
			Runtime.getRuntime().exec("mkdir "+experienceName+" "+folderName);
		} catch (Exception e) {
			// TODO: handle exception
		}
		System.out.println("---- Evolution parameters ----");
		System.out.println("Folder name:          "+folderName);
		System.out.println("Random seed:          "+randomSeed);
		System.out.println("GREAT activated:      "+greatActivated);
		System.out.println("Max iterations:       "+maxNumGen);
		System.out.println("Population size:      "+populationSize);
		System.out.println("GREAT initialization  "+greatInit);
		System.out.println("Speciation threshold: "+speciationThreshold);
		System.out.println("Species min size:     "+speciesMinSize);
		System.out.println("Tournament size:      "+tournamentSize);
		System.out.println("Crossover rate:       "+crossoverRate);
		System.out.println("Mutation rate:        "+mutationRate);
		System.out.println("Crossovers:");
		for (GRNCrossoverOperator cross : crossovers) {
			System.out.println("\t"+cross.name);
		}
		System.out.println("Mutations:");
		for (GRNMutationOperator mut : mutations) {
			System.out.println("\t"+mut.name);
		}
		System.out.println("------------------------------\n\n---- Starting evolution ----");
		
		// initialization
		currentGen=0;
		ArrayList<GRNGenome> offsprings=new ArrayList<GRNGenome>();
		while (offsprings.size()<populationSize-duplicateInit-1) {
			GRNGenome g = new GRNGenome();
			for (int i=0; i<evaluator.numGRNInputs; i++) {
				GRNGene gene=GRNGene.generateRandomGene(GRNProtein.INPUT_PROTEIN, i, rng);
//					System.err.print(gene+" ");
				g.addGene(gene);
			}
			for (int i=0; i<evaluator.numGRNOutputs; i++) {
				g.addGene(GRNGene.generateRandomGene(GRNProtein.OUTPUT_PROTEIN, i, rng));
			}
			if (greatInit) {
				// small GRNs
				g.addGene(GRNGene.generateRandomRegulatoryGene(rng));
			} else {
				// standard init
				int nbProt=(int)(rng.nextDouble()*(50-evaluator.numGRNInputs-evaluator.numGRNOutputs))+1;
				for (int i=0; i<nbProt; i++) {
					g.addGene(GRNGene.generateRandomRegulatoryGene(rng));
				}
			}
			g.beta=rng.nextDouble()*(g.betaMax-g.betaMin)+g.betaMin;
			g.delta=rng.nextDouble()*(g.deltaMax-g.deltaMin)+g.deltaMin;
			offsprings.add(g);
			for (int i=0; i<duplicateInit; i++) {
				GRNGenome mutG=new GRNGeneMutationOperator().cloneAndMutate(g, rng);
				offsprings.add(mutG);
			}
		}
		mutationProbabilitySum=0;
		for (GRNMutationOperator mut : mutations) {
			mutationProbabilitySum+=mut.probability;
		}

		// Building species
		//speciation(offsprings);
		//System.err.println(species.size());
		System.out.println("IT\tNanotime\t#Eval\t#Species\tSpecies size\tPopSize\tBest\tAvg\tWorst\tBest GRN Size\tAvg GRN size\tMax GRN size\tMin GRN size\tSpecies avgdist2rep");
		while (currentGen<maxNumGen) {
			evaluator.generation = currentGen;
			if ( evaluator.nonCacheable ) {
				for (GRNGenome grn : offsprings) {
					grn.hasBeenEvaluated=false;
				}
				//				for (int i=0; i<species.size(); i++) {
				//					species.get(i).recomputeFitness();
				//				}
			}

			if( verbose > 0 ) System.err.println( "Speciation" );
			// speciation
			speciation(offsprings);
			System.out.println(currentGen+"\t"+System.nanoTime()+"\t"+GRNGenomeEvaluator.numEvaluations+"\t"+statistics());
//			System.out.println("IT: "+currentGen+ " "+ System.nanoTime() + "("+GRNGenomeEvaluator.numEvaluations+") === "+ statistics()+ " ===");
			saveBestGRN();
			if (saveAllPopulation) {
				saveAllPopulation();
			}

			if( verbose > 0 ) System.err.println( "Adjusting fitness" );
			// adjusting fitness
			int numOff[]=new int[species.size()];
			int sumNumOff=0;
			if (speciesSizeAdjustingMethod==0) {
				double sumAdjFit[]=new double[species.size()];
				double sumSumAdjFit=0;
				// Adjusting according to the number of genomes in a structural perimeter
				for (int i=0; i<species.size(); i++) {
					Species s=species.get(i);
					sumAdjFit[i]=0;
					double bestFitness=s.getBestFitness();
					for (GRNGenome gi : s.getGenomes()) {
						double sh=0;
						for (GRNGenome gj : s.getGenomes()) {
							if (gj.distanceTo(gi, useDynamicsGenomeForDistance)<speciationThreshold) {
								sh++;
							}
						}
//						System.err.println(sh+"  "+species.size());
						gi.setAdjustedFitness(gi.getLastFitness()/sh);
						sumAdjFit[i]+=gi.getAdjustedFitness();
					}
					sumSumAdjFit+=sumAdjFit[i];
				}
				adjustSpeciationThresholds();
				// calculating offspring number for each species
				for (int i=0; i<species.size(); i++) {
					numOff[i]=(int)(sumAdjFit[i]/sumSumAdjFit*populationSize);
					sumNumOff+=numOff[i];
					//System.err.print(numOff[i]+" ");
				}			
			} else if (speciesSizeAdjustingMethod==1) {
				double sumAdjFit[]=new double[species.size()];
				double sumSumAdjFit=0;

				// Adjusting according to the number of genome number of genomes with comparable fitnesses
				for (int i=0; i<species.size(); i++) {
					Species s=species.get(i);
					sumAdjFit[i]=0;
					double bestFitness=s.getBestFitness();
					for (GRNGenome gi : s.getGenomes()) {
						double sh=0;
						for (GRNGenome gj : s.getGenomes()) {
							if (Math.abs(gi.getLastFitness()-gj.getLastFitness())/bestFitness<speciationThreshold) {
								sh++;
							}
						}
						gi.setAdjustedFitness(gi.getLastFitness()/sh);
						sumAdjFit[i]+=gi.getAdjustedFitness();
					}
					sumSumAdjFit+=sumAdjFit[i];
				}
				adjustSpeciationThresholds();
				// calculating offspring number for each species
				for (int i=0; i<species.size(); i++) {
					numOff[i]=(int)(sumAdjFit[i]/sumSumAdjFit*populationSize);
					sumNumOff+=numOff[i];
					//System.err.print(numOff[i]+" ");
				}			
			} else if (speciesSizeAdjustingMethod==2) {
				double sumAdjFit[]=new double[species.size()];
				double sumSumAdjFit=0;

				// Adjusting according to the species size
				for (int i=0; i<species.size(); i++) {
					Species s=species.get(i);
					double spSize=s.size();
					sumAdjFit[i]=0;
					for (GRNGenome gi : s.getGenomes()) {
						gi.setAdjustedFitness(gi.getLastFitness()/spSize);
						sumAdjFit[i]+=gi.getAdjustedFitness();
					}
					sumSumAdjFit+=sumAdjFit[i];
				}
				adjustSpeciationThresholds();
				// calculating offspring number for each species
				for (int i=0; i<species.size(); i++) {
					numOff[i]=(int)(sumAdjFit[i]/sumSumAdjFit*populationSize);
					sumNumOff+=numOff[i];
					//System.err.print(numOff[i]+" ");
				}			
			} else if (speciesSizeAdjustingMethod==3) {
				double propBad[]=new double[species.size()];
				double sumPropBad=0;
				for (int i=0; i<species.size(); i++) {
					Species s = species.get(i);
					propBad[i]=0;
					for (GRNGenome g : s.getGenomes()) {
						if (g.getLastFitness()>g.parentBestFitness) {
							propBad[i]+=1.0;
						}
					}
					propBad[i]/=s.size();
					sumPropBad+=propBad[i];
				}
				sumPropBad/=species.size();
				for (int i=0; i<species.size(); i++) {
					numOff[i]=(int)((double)(populationSize)*propBad[i]/sumPropBad);
					sumNumOff+=numOff[i];
					System.err.print(numOff[i]+"  ");
				}
			} else if (speciesSizeAdjustingMethod==4) {
				double sumAdjFit[]=new double[species.size()];
				double sumSumAdjFit=0;

				// Adjusting according to the species size
				for (int i=0; i<species.size(); i++) {
					Species s=species.get(i);
					double spSize=s.size();
					sumAdjFit[i]=0;
					for (GRNGenome gi : s.getGenomes()) {
						gi.setAdjustedFitness(gi.getLastFitness()/spSize);
						sumAdjFit[i]+=gi.getAdjustedFitness();
					}
					sumSumAdjFit+=sumAdjFit[i];
				}
				adjustSpeciationThresholds();
				// calculating offspring number for each species
				for (int i=0; i<species.size(); i++) {
					numOff[i]=5;
					sumNumOff+=5;
				}
				while (sumNumOff<populationSize) {
					double rnd=rng.nextDouble()*sumSumAdjFit;
					int i=0;
					while (i<species.size() && rnd<sumAdjFit[i]) {
						rnd-=sumAdjFit[i];
						i++;
					}
					numOff[i]++;
					sumNumOff++;
				}
/*				for (int i=0; i<species.size(); i++) {
					//numOff[i]=(int)((1.0-sumAdjFit[i]/species.get(i).size()/sumSumAdjFit)*populationSize);
					//sumNumOff+=numOff[i];
					System.err.print(numOff[i]+" ");
				}*/			
			}
//			System.err.println();
			
			// correcting approximation errors
			while (sumNumOff<populationSize) {
				//int rnd=(int)(Math.random()*species.size());
				int rnd=(int)(rng.nextDouble()*species.size());
				numOff[rnd]+=1;
				sumNumOff++;
			}

			if( verbose > 0 ) System.err.println( "Renewing species" );
			// renewing species
			offsprings.clear();
			int nGenomes;
			for (int i=0; i<species.size(); i++) {
				Species s=species.get(i);
				nGenomes=0;
				// adding crossovered genomes
				int nbFailedCross=0;
				while (nGenomes+nbFailedCross<crossoverRate*numOff[i]) {
					GRNGenome g1=s.tournamentSelect(tournamentSize,rng);
					GRNGenome g2=s.tournamentSelect(tournamentSize,rng);
					int tries=0;
					while (g1==g2 && tries++<nMaxTries) {
						g2=s.tournamentSelect(tournamentSize,rng);
					}
					if (tries<10) {
						//int indexCrossover=(int)(Math.random()*crossovers.size());
						int indexCrossover=(int)(rng.nextDouble()*crossovers.size());
						GRNGenome child = crossovers.get(indexCrossover).reproduce(g1, g2, rng);
						offsprings.add(child);
						nGenomes++;
					} else {
						nbFailedCross++;
					}
				}
				// adding mutated genomes
				while (nGenomes<(mutationRate+crossoverRate)*numOff[i]) {
					GRNGenome selectGen = s.tournamentSelect(tournamentSize,rng);
					GRNGenome mutatedGen=null;
					while (mutatedGen==null) {
						// chosing a mutation operator randomly
						//int indexMutation = (int)(Math.random()*mutations.size());
						int indexMutation = selectMutationOperator(); //(int)(rng.nextDouble()*mutations.size());
						// mutate
						mutatedGen=mutations.get(indexMutation).cloneAndMutate(selectGen, rng);
					}
					offsprings.add(mutatedGen);
					nGenomes++;
				}
				// adding best genome
				if (numOff[i]>0) {
					offsprings.add(s.getBestGenome());
					nGenomes++;
					//System.err.print(s.getBestFitness()+"  ");
				}
				// completing with existing genomes
				while (nGenomes<numOff[i]) {
					offsprings.add(s.tournamentSelect(tournamentSize, rng));
					nGenomes++;
				}
			}
			//			System.err.println("offSize "+offsprings.size());
			currentGen++;
		}
	}
	
	protected int selectMutationOperator() {
		double rnd=rng.nextDouble()*mutationProbabilitySum;
		int mutIndex=0;
		while (mutIndex<mutations.size() && rnd>mutations.get(mutIndex).probability) {
			rnd-=mutations.get(mutIndex).probability;
			mutIndex++;
		}
		if (verbose>0) System.err.println("Mutation selected: "+mutations.get(mutIndex).name+" with index "+mutIndex);
		return mutIndex;
	}

	protected String statistics() {
/*		String res="";
		for (int i=0; i<Species.speciesNextId; i++) {
			Species fs=null;
			for (Species s:species) {
				if (s.speciesId==i) {
					fs=s;
				}
			}
			if (fs==null) res+="0\t";
			else res+=fs.size()+"\t";
		}
		return res;
*/		
		int nGenomes=0;
		double fitMax=-Double.MAX_VALUE;
		double fitMin=Double.MAX_VALUE;
		double fitSum=0;
		int bestGRNSize=0;
		double avgGRNSize=0;
		int minGRNSize=Integer.MAX_VALUE;
		int maxGRNSize=0;
		String speciesSizes="(";
		String avgdist2rep="(";
		for (Species s : species) {
			nGenomes+=s.size();
			if (s.getBestFitness()>fitMax) {
				fitMax=s.getBestFitness();
				bestGRNSize=s.getBestGenome().size();
			}
			for (GRNGenome g : s.getGenomes()) {
				avgGRNSize+=g.size();
				minGRNSize=Math.min(minGRNSize, g.size());
				maxGRNSize=Math.max(maxGRNSize, g.size());
			}
			
			fitMin=Math.min(s.getWorstFitness(), fitMin);
			fitSum+=s.getFitnessSum();
			speciesSizes+=s.speciesId+":"+s.size()+"["+s.speciationThreshold+"],";
			avgdist2rep+=s.getAvgDist2Rep()+",";
		}
		speciesSizes+=")";
		avgdist2rep+=")";
//		return "Species: " +species.size()+speciesSizes+" -- PopSize: "+nGenomes+" -- FitMax/Avg/Min: "+fitMax+"/"+fitSum/nGenomes+"/"+fitMin;
		return species.size()+"\t"+speciesSizes+"\t"+nGenomes+"\t"+fitMax+"\t"+fitSum/nGenomes+"\t"+fitMin+"\t"+bestGRNSize+"\t"+avgGRNSize/populationSize+"\t"+maxGRNSize+"\t"+minGRNSize+"\t"+avgdist2rep;

	}

	protected void saveBestGRN() {
		double fitMax=-Double.MAX_VALUE;
		GRNGenome bestGenome=null;
		for (Species s : species) {
			if (s.getBestFitness()>fitMax) {
				bestGenome=s.getBestGenome();
				fitMax=s.getBestFitness();
			}
		}
		try {
			evaluator.buildGRNFromGenome(bestGenome).writeToFile(folderName+"/grn_"+currentGen+"_"+fitMax+".grn");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	protected void saveAllPopulation() {
		try {
			// making population directory
			Runtime.getRuntime().exec("mkdir "+folderName+"/population_"+currentGen);
			// making species directory
			for (int i=0; i<species.size(); i++) {
				Process p=Runtime.getRuntime().exec("mkdir "+folderName+"/population_"+currentGen+"/species_"+i);
				// saving all grns of the species
				//				System.err.println("Saving species "+i);
				Species s=species.get(i);
				s.sortGenomes();
				Collection<GRNGenome> genomes = s.getGenomes();
				Iterator<GRNGenome> g=genomes.iterator();
				boolean folderNotReady=true;
				while (folderNotReady) {
					try {
						folderNotReady=p.exitValue()!=0;
						if (folderNotReady) {
							System.err.println("Error during mkdir proccess. Retrying...");
							p=Runtime.getRuntime().exec("mkdir "+folderName+"/population_"+currentGen+"/species_"+i);
						}
					} catch (Exception e) {
						//						System.err.println("Folder not yet ready");
						Thread.sleep(1);
					}
				}
				for (int j=0; g.hasNext(); j++) {
					int tries=0;
					boolean success=false;
					while (tries<nMaxTries && !success) {
						try {
							evaluator.buildGRNFromGenome(g.next()).writeToFile(folderName+"/population_"+currentGen+"/species_"+i+"/genome_"+j+".grn");
							success=true;
						} catch (Exception e) {
							//System.err.println(folderName+"/population_"+currentGen+"/species_"+i+"/genome_"+j+".grn");
							//System.err.println("Folder not ready... waiting a second");
							Thread.sleep(1000);
							tries++;
						}
					}
					if (!success) {
						System.err.println("Couldn't write file "+folderName+"/population_"+currentGen+"/species_"+i+"/genome_"+j+".grn");
					}
				}
			}
			Runtime.getRuntime().exec("./popcompressor.sh "+folderName+" population_"+currentGen);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	protected void speciation(List<GRNGenome> population) {
		for (Species s : species) {
			if( representativeMethod == 1 ) {
				s.chooseRandomRepresentative( rng );
			} else if( representativeMethod == 2 ) {
				s.chooseBestRepresentative( rng );
			}
			s.removeAllGenomes();
		}
		
		if( verbose > 0 ) System.err.println( "Speciation: Adding to species" );
	
		for (GRNGenome g : population) {
			Species bestMatch=null;
			double minDist=Double.MAX_VALUE;
			for (Species s : species) {
				double currentDist=g.distanceTo(s.representativeGenome, useDynamicsGenomeForDistance);
				if (currentDist<minDist && currentDist<s.speciationThreshold) {
					bestMatch=s;
					minDist=currentDist;
				}
			}
			if (bestMatch==null || minDist>bestMatch.speciationThreshold) {
				// Creating a new Species
				Species s=new Species(g, speciationThreshold, true, evaluator);
				species.add(s);
			} else {
				// add the genome to the best matching species
				bestMatch.addGenome(g, false);
			}
		}
	
		if (verbose>0) System.err.println("Speciation: "+statistics());

		if( verbose > 0 ) System.err.println( "Speciation: Removing empty species" );
		// removing empty species
		int numGenomeDelete=0;//populationSize-population.size();
		int numSpeciesDeleted=0;
		for (int i=0; i<species.size(); i++) {
			Species s = species.get(i);
			//System.out.println( s.size() );
			if (s.size()<speciesMinSize && species.size() > 1) {
				numGenomeDelete+=s.size();
				species.remove(s);
				i--;
				numSpeciesDeleted++;
			}
		}
		if( verbose > 0 ) System.err.println( "Speciation: "+numSpeciesDeleted+" species deleted");
		if( verbose > 0 ) System.err.println( "Speciation: Adding " + numGenomeDelete + " genomes" );
		// creating new genomes to keep the popSize constant
		while (numGenomeDelete>0) {
			//int rnd=(int)(Math.random()*species.size());
			int rnd=(int)(rng.nextDouble()*species.size());
			Species s=species.get(rnd);
			GRNGenome g=s.tournamentSelect(tournamentSize,rng);
			boolean geneAdd=false;
			int tries=0;
			while (!geneAdd && tries++<nMaxTries) {
				// chosing a mutation operator randomly
				int indexMutation = (int)(rng.nextDouble()*mutations.size());
				// mutate
				GRNGenome mutatedGen=mutations.get(indexMutation).cloneAndMutate(g,rng);
				if (mutatedGen!=null) {
					geneAdd=s.addGenome(mutatedGen, true);
				}
			}
			if (geneAdd) {
				numGenomeDelete--;
			}
		}
		if( verbose > 0 ) System.err.println( "Speciation: complete" );
	}

	public int testPopSize() {
		int popSize=0;
		for (Species s : species) {
			popSize+=s.size();
		}
		return popSize;		
	}

	public static void main(String args[]) throws Exception {
		Evolver e=new Evolver();
		boolean customExpName=false;
		String sshTarget = "";
		
		int xover = 1;
		
		double addMutationProb=0.33;
		double delMutationProb=0.33;
		double changeMutationProb=0.33;
		int addMutationMaxSize=Integer.MAX_VALUE;
		int delMutationMinSize=0;
		
		// CL-arguments
		for (int k = 0; k < args.length; k++ ) {
//			System.err.print("\t" + args[k]);

			if ( args[k].compareTo("maxNumGen") == 0) {
				e.maxNumGen = Integer.parseInt( args[k+1] );
			} else if ( args[k].compareTo("populationSize") == 0 ) {
				e.populationSize = Integer.parseInt( args[k+1] );
			} else if ( args[k].compareTo("speciationThreshold")  == 0) {
				e.speciationThreshold = Double.parseDouble( args[k+1] );
			} else if ( args[k].compareTo("crossoverRate") == 0 ) {
				e.crossoverRate = Double.parseDouble( args[k+1] );
			} else if ( args[k].compareTo("mutationRate") == 0 ) {
				e.mutationRate = Double.parseDouble( args[k+1] );
			} else if ( args[k].compareTo("speciesMinSize") == 0 ) {
				e.speciesMinSize = Integer.parseInt( args[k+1] );
			} else if ( args[k].compareTo("nMaxTries") == 0 ) {
				e.nMaxTries= Integer.parseInt( args[k+1] );
			} else if ( args[k].compareTo("saveAllPopulation") == 0 ) {
				e.saveAllPopulation= Boolean.parseBoolean( args[k+1] );
			} else if ( args[k].compareTo("randomSeed") == 0 ) {
				e.randomSeed= Integer.parseInt( args[k+1] );
			} else if ( args[k].compareTo("experienceName") == 0 ) {
				e.experienceName = args[k+1];
				customExpName=true;
			} else if ( args[k].compareTo("sshTarget") == 0 ) {
				sshTarget= args[k+1];
			} else if ( args[k].compareTo("great")==0) {
				e.greatActivated=Boolean.parseBoolean(args[k+1]);			
			} else if ( args[k].compareTo("greatInit")==0) {
				e.greatInit=Boolean.parseBoolean(args[k+1]);			
			} else if ( args[k].compareTo("xover")==0) {
				xover=Integer.parseInt(args[k+1]);				
			} else if ( args[k].compareTo("representativeMethod")==0) {
				e.representativeMethod=Integer.parseInt(args[k+1]);
			} else if ( args[k].compareTo("speciesSizeAdjustingMethod")==0) {
				e.speciesSizeAdjustingMethod=Integer.parseInt(args[k+1]);
			} else if (args[k].compareTo("duplicateInit")==0) {
				e.duplicateInit=Integer.parseInt(args[k+1]);
			} else if (args[k].compareTo("verbose")==0) {
				e.verbose=Integer.parseInt(args[k+1]);
			} else if (args[k].compareTo("addMutationProbability")==0) {
				addMutationProb=Double.parseDouble(args[k+1]);
			} else if (args[k].compareTo("deleteMutationProbability")==0) {
				delMutationProb=Double.parseDouble(args[k+1]);
			} else if (args[k].compareTo("changeMutationProbability")==0) {
				changeMutationProb=Double.parseDouble(args[k+1]);
			} else if (args[k].compareTo("addMutationMaxSize")==0) {
				addMutationMaxSize=Integer.parseInt(args[k+1]);
			} else if (args[k].compareTo("deleteMutationMinSize")==0) {
				delMutationMinSize=Integer.parseInt(args[k+1]);
			}
		}

		e.rng = new java.util.Random( e.randomSeed );
		//e.greatActivated=false;
		
		if (!customExpName) {
			//e.evaluator=new CoverageControl( args, e.rng );//new MichalSignalProcessExp3();
			e.evaluator = new IntertwinedSpirals( args );
			e.experienceName=e.evaluator.name;
		} else {
			if( e.experienceName.compareTo("DoublingFrequencyEvaluator") == 0) e.evaluator=new DoublingFrequencyEvaluator();  
			else if( e.experienceName.compareTo("MichalSignalProcessExp1") == 0 ) e.evaluator=new MichalSignalProcessExp1();  
			else if( e.experienceName.compareTo("MichalSignalProcessExp2") == 0 ) e.evaluator=new MichalSignalProcessExp2();  
			else if( e.experienceName.compareTo("MichalSignalProcessExp3") == 0 ) e.evaluator=new MichalSignalProcessExp3();  
			else if( e.experienceName.compareTo("CoverageControl") == 0 ) e.evaluator=new CoverageControl( args, e.rng );  
			else if( e.experienceName.compareTo("IntertwinedSpirals") == 0 ) e.evaluator=new IntertwinedSpirals( args );
			else if( e.experienceName.compareTo("IntertwinedSpiralsSubsequence") == 0 ) e.evaluator=new IntertwinedSpiralsSubsequence();
			else if( e.experienceName.compareTo("IntertwinedSpiralsVariableTime") == 0 ) e.evaluator=new IntertwinedSpiralsVariableTime( args );
			else if( e.experienceName.compareTo("WeatherPrediction") == 0 ) e.evaluator=new WeatherPrediction(1);
			else if( e.experienceName.compareTo("Radbot") == 0 ) e.evaluator=new RadbotEvaluator();
			else if( e.experienceName.compareTo("PingPrediction") == 0 ) e.evaluator=new PingPrediction(0);
			//e.experienceName=e.evaluator.name;
		}
		
		if (e.greatActivated) {
			System.out.println("GREAT!");
			System.out.println("Xover: " + xover);
			if( xover==0 )
				e.crossovers.add(new GRNAligningCrossoverOperator_v1());
			else if (xover==1)
				e.crossovers.add(new GRNAligningCrossoverOperator_v1b());
			else if (xover==2)
				e.crossovers.add(new GRNAligningCrossoverOperator_v2());
			else if (xover==3) 
				e.crossovers.add(new GRNOnePointCrossoverOperator());
			else if (xover==4)
				e.crossovers.add(new GRNAligningCrossoverOperator_ParentCountProb());
			e.mutations.add(new GRNAddGeneMutationOperator(addMutationMaxSize, addMutationProb));
			e.mutations.add(new GRNDeleteGeneMutationOperator(Math.max(delMutationMinSize, e.evaluator.numGRNInputs+e.evaluator.numGRNOutputs+1), delMutationProb));
			e.mutations.add(new GRNGeneMutationOperator(changeMutationProb));
			e.greatInit=true;
		} else {
			System.out.println("GRN+AG");
			if( xover==0 )
				e.crossovers.add(new GRNAligningCrossoverOperator_v1());
			else if (xover==1)
				e.crossovers.add(new GRNAligningCrossoverOperator_v1b());
			else if (xover==2)
				e.crossovers.add(new GRNAligningCrossoverOperator_v2());
			else if (xover==3) 
				e.crossovers.add(new GRNOnePointCrossoverOperator());
			if (addMutationMaxSize<Integer.MAX_VALUE) {
				e.mutations.add(new GRNAddGeneMutationOperator(addMutationMaxSize, addMutationProb));
			} else {
				e.mutations.add(new GRNAddGeneMutationOperator(50, addMutationProb));
			}
			e.mutations.add(new GRNDeleteGeneMutationOperator(Math.max(delMutationMinSize, e.evaluator.numGRNInputs+e.evaluator.numGRNOutputs+1), delMutationProb));
			e.mutations.add(new GRNGeneMutationOperator(changeMutationProb));
			e.speciationThreshold=1000;
			e.greatInit=false;
		}



//		e.verbose=1;
		
		//e.evaluator=new NMMultirobotCCP(args,e.rng);
		//e.evaluator=new RetinaEvaluator();
		//e.evaluator=new GRNSinusEvaluator();
		e.run();

		if( sshTarget.compareTo("") != 0 ) {
			try {			
				Runtime.getRuntime().exec("tar czf "+e.experienceName+"_"+e.uuid+".tar.gz "+e.folderName);
				Runtime.getRuntime().exec("scp "+e.experienceName+"_"+e.uuid+".tar.gz "+sshTarget);
			} catch (Exception exc) {
				exc.printStackTrace();
			}
		}
	}
	
}
