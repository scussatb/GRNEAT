package operators;

import java.util.Random;

import evolver.GRNGenome;


public abstract class GRNCrossoverOperator {
	public String name="SuperCrossover!";
	
	public abstract GRNGenome reproduce(GRNGenome parent1, GRNGenome parent2, Random rng);
}
