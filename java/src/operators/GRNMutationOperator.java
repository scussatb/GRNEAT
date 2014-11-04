package operators;

import java.util.Random;

import evolver.GRNGenome;

public abstract class GRNMutationOperator {
	public String name="SuperMutation!";
	public double probability=0.33;
	
	public GRNGenome cloneAndMutate(GRNGenome aGenome, Random rng) {
		GRNGenome offspring = aGenome.clone();
		if (mutateByModifying(offspring, rng)) {
			return offspring;
		} else {
			return null;
		}
	}

	public abstract boolean mutateByModifying(GRNGenome aGenome, Random rng);
}
