package evaluators;

import evolver.GRNGenome;
import grn.GRNModel;

public class MichalSignalProcessExp1 extends GRNGenomeEvaluator {
	double coef=3.0;
	
	public MichalSignalProcessExp1() {
		numGRNInputs=0;
		numGRNOutputs=1;
		name="SignalProcessingExp1";
	}

	@Override
	public double evaluate(GRNGenome aGenome) {

		double fitness=evolveGRN(buildGRNFromGenome(aGenome), false);
		
		//System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		aGenome.setNewFitness(fitness);

		GRNGenomeEvaluator.numEvaluations++;
		return fitness;
	}
	
	double evolveGRN(GRNModel grn, boolean printTrace) {
		double fitness=0.0;
		
		grn.reset();
		grn.evolve(25);
		double ots[]={0,0,0,0,0,0,0,0,0,0};
		int nbEvents=0;
		int lastEvent=0;
		for (int nStep=0; nStep<600; nStep++) {
			double dt;
			grn.evolve(1);
			if (nStep<100 || nStep>475) {
				dt=0;
			} else {
				dt=(Math.sin((double)(nStep-100)/75.0*2.0*Math.PI-Math.PI/2.0)+1.0)/2.0;
			}
			for (int nst=0; nst<9; nst++) {
				ots[nst]=ots[nst+1];
			}
			ots[9]=Math.min(1.0,grn.proteins.get(0).concentration*coef);
			if (nStep-lastEvent>10 && ((ots[0]<0.5 && ots[9]>=0.5) || (ots[0]>=0.5 && ots[9]<0.5))) {
				nbEvents++;
				//System.out.println("Event "+lastEvent+"  "+nStep);
				lastEvent=nStep;
			}
			fitness+=Math.abs(dt-Math.min(1.0,grn.proteins.get(0).concentration*coef))*(1+dt);
			if (printTrace) {
				System.out.println(nStep+"\t"+Math.min(1.0,grn.proteins.get(0).concentration*coef)+"\t"+dt);
			}
			if (grn.proteins.get(0).concentration!=grn.proteins.get(0).concentration) {
				System.err.println("nan");
			}
		}
		double S;
		if (nbEvents==0 || nbEvents>20) {
			S=0;
		} else {
			S=1.0-(double)(Math.abs(nbEvents-10))/10.0;
		}
		if (printTrace) {
			System.out.println(nbEvents+"  "+ S);
		}
		fitness*=1.0/(1.0+S);
		
		//System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		return -fitness;
	}
	
	public static void main(String args[]) throws Exception {
		GRNModel grn = GRNModel.loadFromFile("run_2013-03-06_10-12-08/grn_3128_-104.4891263181337.grn");
		MichalSignalProcessExp1 eval=new MichalSignalProcessExp1();
		double fit=eval.evolveGRN(grn, true);
		System.out.println("Fitness="+fit);
	}

}
