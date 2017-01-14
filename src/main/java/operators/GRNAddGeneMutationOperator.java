package operators;

import java.util.Random;

import evolver.GRNGene;
import evolver.GRNGenome;

public class GRNAddGeneMutationOperator extends GRNMutationOperator {
	public int grnMaxSize=Integer.MAX_VALUE;
	
	public GRNAddGeneMutationOperator() {
		grnMaxSize=Integer.MAX_VALUE;
		name="AddGeneMutation(maxSize="+grnMaxSize+", prob="+probability+")";
	}
	
	public GRNAddGeneMutationOperator(double prob) {
		grnMaxSize=Integer.MAX_VALUE;
		this.probability=prob;
		name="AddGeneMutation(maxSize="+grnMaxSize+", prob="+probability+")";
	}
	
	public GRNAddGeneMutationOperator(int nGRNMaxSize) {
		grnMaxSize=nGRNMaxSize;
		name="AddGeneMutation(maxSize="+grnMaxSize+", prob="+probability+")";
	}
	
	public GRNAddGeneMutationOperator(int nGRNMaxSize, double prob) {
		grnMaxSize=nGRNMaxSize;
		probability=prob;
		name="AddGeneMutation(maxSize="+grnMaxSize+", prob="+probability+")";
	}
	
	@Override
	public boolean mutateByModifying(GRNGenome aGenome, Random rng) {
		if (aGenome.size()<grnMaxSize) {
			aGenome.addGene(GRNGene.generateRandomRegulatoryGene(rng));
			aGenome.parentBestFitness=aGenome.getLastFitness();
			aGenome.hasBeenModified();
			return true;
		} else {
			return false;
		}
	}
	
}
