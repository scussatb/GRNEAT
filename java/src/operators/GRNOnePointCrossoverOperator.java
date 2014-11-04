package operators;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;

import evolver.GRNGene;
import evolver.GRNGenome;


public class GRNOnePointCrossoverOperator extends GRNCrossoverOperator implements Comparator<GRNGene>{
	
	public GRNOnePointCrossoverOperator() {
		name="OnePointCrossover";
	}
	
	public GRNGenome reproduce(GRNGenome parent1, GRNGenome parent2, Random rng) {
//		System.err.println("One point");
		GRNGenome offspring = new GRNGenome();
		int point=(int)(rng.nextDouble()*Math.min(parent1.size(), parent2.size()));
		boolean p1First=rng.nextDouble()<0.5;
		// reordering genes
		ArrayList<GRNGene> p1Genes=new ArrayList<GRNGene>(parent1.size());
		ArrayList<GRNGene> p2Genes=new ArrayList<GRNGene>(parent2.size());
		ArrayList<GRNGene> sortInputsP1=new ArrayList<GRNGene>(parent1.getInputGenes());
		Collections.sort(sortInputsP1, this);
		ArrayList<GRNGene> sortOutputP1=new ArrayList<GRNGene>(parent1.getOutputGenes());
		Collections.sort(sortOutputP1, this);
		ArrayList<GRNGene> sortInputsP2=new ArrayList<GRNGene>(parent2.getInputGenes());
		Collections.sort(sortInputsP2, this);
		ArrayList<GRNGene> sortOutputP2=new ArrayList<GRNGene>(parent2.getOutputGenes());
		Collections.sort(sortOutputP2, this);
		ArrayList<GRNGene> shuffRegP1=new ArrayList<GRNGene>(parent1.getRegulatoryGenes());
		Collections.shuffle(shuffRegP1);
		ArrayList<GRNGene> shuffRegP2=new ArrayList<GRNGene>(parent2.getRegulatoryGenes());
		Collections.shuffle(shuffRegP2);
		p1Genes.addAll(sortInputsP1);
		p1Genes.addAll(sortOutputP1);
		p1Genes.addAll(shuffRegP1);
		p2Genes.addAll(sortInputsP2);
		p2Genes.addAll(sortOutputP2);
		p2Genes.addAll(shuffRegP2);

		int i=0;
		while (i<point) {
			if (p1First) {
				offspring.addGene(p1Genes.get(i));
			} else {
				offspring.addGene(p2Genes.get(i));
			}
			i++;
		}
		while (i<(p1First?parent2.size():parent1.size())) {
			if (p1First) {
				offspring.addGene(p2Genes.get(i));
			} else {
				offspring.addGene(p1Genes.get(i));
			}
			i++;
		}

		offspring.parentBestFitness=Math.max(parent1.getLastFitness(), parent2.getLastFitness());
	
		offspring.hasBeenModified();
		return offspring;
	}

	@Override
	public int compare(GRNGene o1, GRNGene o2) {
		return o1.getConnectTo()-o2.getConnectTo();
	}

}
