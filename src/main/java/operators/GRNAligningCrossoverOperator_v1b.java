package operators;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import evolver.GRNGene;
import evolver.GRNGenome;


public class GRNAligningCrossoverOperator_v1b extends GRNCrossoverOperator {
	protected double threshold=0.15;
	
	public GRNAligningCrossoverOperator_v1b() {
		name="AligningCrossover_v1b("+threshold+")";
	}

	public GRNGenome reproduce(GRNGenome parent1, GRNGenome parent2, Random rng) {
//		System.err.println("Align");
		GRNGenome offspring = new GRNGenome();
		// Chose randomly the input genes from the parents
		Collection<GRNGene> inputGenes=parent1.getInputGenes();
		for (GRNGene gp1 : inputGenes) {
			if (rng.nextDouble()<0.5) {
				offspring.addGene(gp1);
			} else {
				offspring.addGene(parent2.getInputGeneConnectTo(gp1.getConnectTo()));
			}
		}
		// Chose randomly the output genes from the parents
		Collection<GRNGene> outputGenes=parent1.getOutputGenes();
		for (GRNGene gp1 : outputGenes) {
			if (rng.nextDouble()<0.5) {
				offspring.addGene(gp1);
			} else {
				offspring.addGene(parent2.getOutputGeneConnectTo(gp1.getConnectTo()));
			}
		}
		// Aligning regulatory genes
//		Collection<GRNGene> regGP1=parent1.getRegulatoryGenes();
		ArrayList<GRNGene> regGP1=new ArrayList<GRNGene>(parent1.getRegulatoryGenes());
		Collections.shuffle(regGP1, rng);
		ArrayList<GRNGene> regGP2=new ArrayList<GRNGene>(parent2.getRegulatoryGenes());
		Collections.shuffle(regGP2, rng);
		for (int i=0; i<regGP1.size(); i++) {
			GRNGene gp1=regGP1.get(i);
			double minDist=99999;
			GRNGene pairedGP2=null;
			int paired2Index=-1;
			for (int j=0; j<regGP2.size(); j++) {
				GRNGene gp2=regGP2.get(j);
				if (gp2!=null) {
					// not already paired
					double dist=gp1.distanceTo(gp2);
//					System.err.println(dist);
					if (dist<threshold && dist<minDist) {
						pairedGP2=gp2;
						minDist=dist;
						paired2Index=j;
					}
				}
			}
			if (pairedGP2!=null) {
				// choosing one of them
				if (rng.nextDouble()<0.5) {
					offspring.addGene(gp1);
				} else {
					offspring.addGene(pairedGP2);
				}
//				System.err.println(gp1+" and "+pairedGP2+" have been paired");
				// marking as paired
				regGP1.set(i, null);
				regGP2.set(paired2Index, null);
			}
		}
		// add some with remaining material
		// adding from both parents
		for (int i=0; i<regGP1.size(); i++) {
			if (regGP1.get(i)!=null) {
				offspring.addGene(regGP1.get(i));
			}
		}
		for (int i=0; i<regGP2.size(); i++) {
			if (regGP2.get(i)!=null) {
				offspring.addGene(regGP2.get(i));
			}
		}
		
		// crossing the dynamics coefficients
		if (rng.nextDouble()<0.5) {
			offspring.setBeta(parent1.getBeta());
		} else {
			offspring.setBeta(parent2.getBeta());
		}
		if (rng.nextDouble()<0.5) {
			offspring.setDelta(parent1.getDelta());
		} else {
			offspring.setDelta(parent2.getDelta());
		}

		offspring.parentBestFitness=Math.max(parent1.getLastFitness(), parent2.getLastFitness());

		offspring.hasBeenModified();
		return offspring;
	}

}
