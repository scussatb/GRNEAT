package evaluators;

import evolver.GRNGenome;
import grn.GRNModel;

public class GRNStableMediumOutputEvaluator extends GRNGenomeEvaluator {
	
	public GRNStableMediumOutputEvaluator(int numInput) {
		numGRNInputs=numInput;
		numGRNOutputs=1;
		name="StableMediumOutputExp";
	}

	@Override
	public double evaluate(GRNGenome aGenome) {

		double fitness=1.0;
		GRNModel grn = buildGRNFromGenome(aGenome);
		
		grn.reset();
		grn.evolve(25);
		
		double lastOutputConcentration = grn.proteins.get(numGRNInputs).concentration;
		for (int nStep=0; nStep<25; nStep++) {
			for (int i=0; i<numGRNInputs; i++) {
				grn.proteins.get(i).concentration=Math.random()*0.25/numGRNInputs;
			}
			grn.evolve(25);
			fitness+=Math.abs(lastOutputConcentration-grn.proteins.get(numGRNInputs).concentration);
			lastOutputConcentration=grn.proteins.get(numGRNInputs).concentration;
		}
		
		fitness=-fitness+Math.abs(lastOutputConcentration-0.75);
//		fitness=fitness==0?Double.MAX_VALUE:1.0/fitness;
		System.err.println(grn.toString()+"  => fitness="+fitness);
		aGenome.setNewFitness(fitness);

		return fitness;
	}

}
