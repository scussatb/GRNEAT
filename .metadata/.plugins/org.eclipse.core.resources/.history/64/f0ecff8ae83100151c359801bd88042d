package evolver;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import javax.swing.text.DateFormatter;

import com.sun.org.apache.xerces.internal.impl.dv.dtd.NMTOKENDatatypeValidator;

import evaluators.GRNGenomeEvaluator;
import evaluators.GRNSinusEvaluator;
import evaluators.GRNStableMediumOutputEvaluator;
import evaluators.MichalSignalProcessExp1;
import evaluators.MichalSignalProcessExp2;
import evaluators.MichalSignalProcessExp3;
import evaluators.NMMultirobotCCP;
import evaluators.RetinaEvaluator;
import grn.GRNProtein;

import operators.GRNAddGeneMutationOperator;
import operators.GRNAligningCrossoverOperator;
import operators.GRNCrossoverOperator;
import operators.GRNDeleteGeneMutationOperator;
import operators.GRNGeneMutationOperator;
import operators.GRNMutationOperator;

public class Evolver_neural {
	protected int currentGen=0;

	ArrayList<GRNCrossoverOperator> crossovers = new ArrayList<GRNCrossoverOperator>();
	ArrayList<GRNMutationOperator> mutations = new ArrayList<GRNMutationOperator>();

	ArrayList<Species> species = new ArrayList<Species>();

	public Random rng;
	public int randomSeed=(int)(new Date().getTime());
	public int maxNumGen=Integer.MAX_VALUE;
	public GRNGenomeEvaluator evaluator=null;
	public int populationSize = 500;
	public double speciationThreshold = 0.3;
	public double crossoverRate = 0.8;
	public double mutationRate = 0.1;
	public int speciesMinSize=15;
	public int nMaxTries=10;
	public boolean saveAllPopulation=false;
	public String folderName;
	//public String experienceName="NMMultirobotCCP";
	public String experienceName="RetinaEvaluator";
	public long uuid;

	public void run() {
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		//uuid = System.nanoTime();
		folderName=experienceName+"/run_"+dateFormat.format(new Date());
		System.out.println("Folder name: "+folderName);
		try {
			Runtime.getRuntime().exec("mkdir "+experienceName+" "+folderName);
		} catch (Exception e) {
			// TODO: handle exception
		}
		
		// initialization
		currentGen=0;
		ArrayList<GRNGenome> offsprings=new ArrayList<GRNGenome>();
		while (offsprings.size()<populationSize) {
			GRNGenome g = new GRNGenome();
			for (int i=0; i<evaluator.numGRNInputs; i++) {
				g.addGene(GRNGene.generateRandomGene(GRNProtein.INPUT_PROTEIN, i, rng));
			}
			for (int i=0; i<evaluator.numGRNOutputs; i++) {
				g.addGene(GRNGene.generateRandomGene(GRNProtein.OUTPUT_PROTEIN, i, rng));
			}
			g.addGene(GRNGene.generateRandomRegulatoryGene(rng));
			offsprings.add(g);
		}

		// Building species
		//speciation(offsprings);
		//System.err.println(species.size());

		
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
			
			// speciation
			speciation(offsprings);
			System.out.println("IT: "+currentGen+ " "+ System.nanoTime() + "("+GRNGenomeEvaluator.numEvaluations+") === "+ statistics()+ " ===");
			saveBestGRN();
			if (saveAllPopulation) {
				saveAllPopulation();
			}

			// adjusting fitness
			double sumAdjFit[]=new double[species.size()];
			double sumSumAdjFit=0;
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
			// calculating offspring number for each species
			int numOff[]=new int[species.size()];
			int sumNumOff=0;
			for (int i=0; i<species.size(); i++) {
				numOff[i]=(int)(sumAdjFit[i]/sumSumAdjFit*populationSize);
				sumNumOff+=numOff[i];
				//System.err.print(numOff[i]+" ");
			}
			while (sumNumOff<populationSize) {
				//int rnd=(int)(Math.random()*species.size());
				int rnd=(int)(rng.nextDouble()*species.size());
				numOff[rnd]+=1;
				sumNumOff++;
			}
			// renewing species
			offsprings.clear();
			int nGenomes;
			for (int i=0; i<species.size(); i++) {
				Species s=species.get(i);
				nGenomes=0;
				// adding crossovered genomes
				int nbFailedCross=0;
				while (nGenomes+nbFailedCross<crossoverRate*numOff[i]) {
					GRNGenome g1=s.tournamentSelect(7,rng);
					GRNGenome g2=s.tournamentSelect(7,rng);
					int tries=0;
					while (g1==g2 && tries++<nMaxTries) {
						g2=s.tournamentSelect(7,rng);
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
					GRNGenome selectGen = s.tournamentSelect(7,rng);
					GRNGenome mutatedGen=null;
					while (mutatedGen==null) {
						// chosing a mutation operator randomly
						//int indexMutation = (int)(Math.random()*mutations.size());
						int indexMutation = (int)(rng.nextDouble()*mutations.size());
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
					offsprings.add(s.randomSelect(rng));
					nGenomes++;
				}
			}
//			System.err.println("offSize "+offsprings.size());
			currentGen++;
		}
	}

	protected String statistics() {
		int nGenomes=0;
		double fitMax=-Double.MAX_VALUE;
		double fitMin=Double.MAX_VALUE;
		double fitSum=0;
		String speciesSizes="(";
		for (Species s : species) {
			nGenomes+=s.size();
			fitMax=Math.max(s.getBestFitness(), fitMax);
			fitMin=Math.min(s.getWorstFitness(), fitMin);
			fitSum+=s.getFitnessSum();
			speciesSizes=speciesSizes+s.size()+",";
		}
		speciesSizes=speciesSizes+")";
		return "Species: " +species.size()+speciesSizes+" -- PopSize: "+nGenomes+" -- FitMax/Avg/Min: "+fitMax+"/"+fitSum/nGenomes+"/"+fitMin;

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
			s.removeAllGenomes();
		}
		for (GRNGenome g : population) {
			Iterator<Species> speciesIt = species.iterator();
			boolean foundMySpecies=false;
			Species sp=null;
			while (speciesIt.hasNext() && !foundMySpecies) {
				sp=speciesIt.next();
				foundMySpecies=sp.addGenome(g, true);
			}
			if (!foundMySpecies) {
				// Creating a new Species
				Species s=new Species(g, speciationThreshold, true, evaluator);
				//s.addGenome(g, false);
				species.add(s);
				//System.err.println("New species "+s+" size:"+s.size()+"  "+species.size());
			} else {
				//System.err.println("Added to "+sp+" size:"+sp.size()+"  "+species.size());
			}
		}
		System.err.println("Species: "+statistics());

		// removing empty species
		int numGenomeDelete=0;
		for (int i=0; i<species.size(); i++) {
			Species s = species.get(i);
			//System.out.println( s.size() );
			if (s.size()<speciesMinSize && species.size() > 1) {
				numGenomeDelete+=s.size();
				species.remove(s);
				i--;
			}
		}
		
		// creating new genomes to keep the popSize constant
		while (numGenomeDelete>0) {
			//int rnd=(int)(Math.random()*species.size());
			int rnd=(int)(rng.nextDouble()*species.size());
			Species s=species.get(rnd);
			GRNGenome g=s.tournamentSelect(7,rng);
			boolean geneAdd=false;
			int tries=0;
			while (!geneAdd && tries++<nMaxTries) {
				// chosing a mutation operator randomly
				//int indexMutation = (int)(Math.random()*mutations.size());
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
	}

	public int testPopSize() {
		int popSize=0;
		for (Species s : species) {
			popSize+=s.size();
		}
		return popSize;		
	}

	public static void main(String args[]) {
		Evolver e=new Evolver();
		
		String sshTarget = "";
		// CL-arguments
		for (int k = 0; k < args.length; k++ ) {
			System.err.print("\t" + args[k]);
			
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
			} else if ( args[k].compareTo("sshTarget") == 0 ) {
				sshTarget= args[k+1];
			}
		}		
		
		e.rng = new java.util.Random( e.randomSeed );
		e.evaluator=new MichalSignalProcessExp3();
		//e.evaluator=new NMMultirobotCCP(args,e.rng);
		//e.evaluator=new RetinaEvaluator();
		//e.evaluator=new GRNSinusEvaluator();
		e.crossovers.add(new GRNAligningCrossoverOperator());
		e.mutations.add(new GRNAddGeneMutationOperator());
		e.mutations.add(new GRNDeleteGeneMutationOperator());
		e.mutations.add(new GRNGeneMutationOperator());
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
