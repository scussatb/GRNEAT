package evaluators;

import evolver.GRNGenome;
import grn.GRNModel;

public class GRNSinusEvaluator extends GRNGenomeEvaluator {
	
	public static int numEvaluations=0;
	
	public GRNSinusEvaluator() {
		numGRNInputs=2;
		numGRNOutputs=1;
		name="SinusExperience";
	}

	@Override
	public double evaluate(GRNGenome aGenome) {

		double fitness=0.0;
		GRNModel grn = buildGRNFromGenome(aGenome);
		
		grn.reset();
		grn.evolve(25);
		

		for (int nStep=0; nStep<100; nStep++) {
			grn.proteins.get(0).concentration=0.125*(double)nStep/100.0;
			grn.proteins.get(1).concentration=0.125*(Math.sin((double)nStep)+1.0)/2.0;
			grn.evolve(1);
			fitness+=Math.abs(Math.sin((double)nStep)-(grn.proteins.get(2).concentration*8.0-1.0));
		}
		
		//System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		aGenome.setNewFitness(-fitness);

		numEvaluations++;
		return -fitness;
	}

}
