package evaluators;

import evolver.GRNGenome;
import grn.GRNModel;

public class MichalSignalProcessExp2 extends GRNGenomeEvaluator {
	public String name = "MichalSignalProcessingExp2Evaluator";
	
	double coef=3.0;
	
	public MichalSignalProcessExp2() {
		numGRNInputs=1;
		numGRNOutputs=1;
	}

	@Override
	public double evaluate(GRNGenome aGenome) {
		GRNModel grn=buildGRNFromGenome(aGenome);
		double fitness=evolveGRN(grn, false, 125.0, 1000)+
				evolveGRN(grn, false, 500, 1000)+
				evolveGRN_zero(grn, false, 1000);
		
		//System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		aGenome.setNewFitness(fitness);

		GRNGenomeEvaluator.numEvaluations++;
		return fitness;
	}
	
	double evolveGRN_zero(GRNModel grn, boolean printTrace, int nStepMax) {
		double fitness=0.0;
		
		grn.reset();
		grn.evolve(25);
		for (int nStep=0; nStep<nStepMax; nStep++) {
			double dt;
			grn.proteins.get(0).concentration=0.0;
			grn.evolve(1);
			dt=0;
			fitness+=Math.abs(dt-Math.min(1.0,grn.proteins.get(1).concentration*coef))*(1+dt);
			if (printTrace) {
				System.out.println(nStep+"\t"+Math.min(1.0,grn.proteins.get(1).concentration*coef)+"\t"+dt);
			}
			if (grn.proteins.get(1).concentration!=grn.proteins.get(1).concentration) {
				System.err.println("nan");
			}
		}
		
		//System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		return -fitness;
		
	}
	
	double evolveGRN(GRNModel grn, boolean printTrace, double halfFreq, int nStepMax) {
		double fitness=0.0;
		
		grn.reset();
		grn.evolve(25);
		double ots[]={0,0,0,0,0,0,0,0,0,0};
		int nbEvents=0;
		int lastEvent=0;
		for (int nStep=0; nStep<nStepMax; nStep++) {
			double dt;
			grn.proteins.get(0).concentration=Math.sin(nStep*Math.PI/halfFreq-Math.PI/2)/2+0.5;
			grn.evolve(1);
			dt=Math.sin(2.0*nStep*Math.PI/halfFreq-Math.PI/2)/2+0.5;;
			for (int nst=0; nst<9; nst++) {
				ots[nst]=ots[nst+1];
			}
			ots[9]=Math.min(1.0,grn.proteins.get(1).concentration*coef);
			if (nStep-lastEvent>10 && ((ots[0]<0.5 && ots[9]>=0.5) || (ots[0]>=0.5 && ots[9]<0.5))) {
				nbEvents++;
				//System.out.println("Event "+lastEvent+"  "+nStep);
				lastEvent=nStep;
			}
			fitness+=Math.abs(dt-Math.min(1.0,ots[9]))*(1+dt);
			if (printTrace) {
				System.out.println(nStep+"\t"+Math.min(1.0,ots[9])+"\t"+dt);
			}
			if (ots[9]!=ots[9]) {
				System.err.println("nan");
			}
		}
		double S;
		double nbEventDesired=(int)(2*(double)nStepMax/halfFreq);
		if (nbEvents==0 || nbEvents>nbEventDesired*2) {
			S=0;
		} else {
			S=1.0-(double)(Math.abs(nbEvents-nbEventDesired))/nbEventDesired;
		}
		if (printTrace) {
			System.out.println(nbEvents+"  "+ S+ "  " +nbEventDesired);
		}
		fitness*=1.0/(1.0+S);
		
		//System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		return -fitness;
	}
	
	double evolveGRN_expE(GRNModel grn, boolean printTrace) {
		double fitness=0.0;
		
		grn.reset();
		grn.evolve(25);
		double ots[]={0,0,0,0,0,0,0,0,0,0};
		int nbEvents=0;
		int lastEvent=0;
		for (int nStep=0; nStep<1000; nStep++) {
			double dt;
			if (nStep<75) {
				dt=(double)nStep/75.0;
			} else if (nStep<125) {
				dt=1.0;
			} else if (nStep<200) {
				dt=1.0-(double)(nStep-125)/75.0;
			} else if (nStep<500) {
				dt=0.0;
			} else if (nStep<575) {
				dt=(double)(nStep-500)/75.0;
			} else if (nStep<825) {
				dt=1.0;
			} else if (nStep<900) {
				dt=1.0-(double)(nStep-825)/75.0;
			} else {
				dt=0.0;
			}
			grn.proteins.get(0).concentration=dt;
			grn.evolve(1);
			for (int nst=0; nst<9; nst++) {
				ots[nst]=ots[nst+1];
			}
			ots[9]=Math.min(1.0,grn.proteins.get(1).concentration*coef);
			if (nStep-lastEvent>10 && ((ots[0]<0.5 && ots[9]>=0.5) || (ots[0]>=0.5 && ots[9]<0.5))) {
				nbEvents++;
				//System.out.println("Event "+lastEvent+"  "+nStep);
				lastEvent=nStep;
			}
			if (printTrace) {
				System.out.println(nStep+"\t"+Math.min(1.0,ots[9])+"\t"+dt);
			}
			if (ots[9]!=ots[9]) {
				System.err.println("nan");
			}
		}
		double S;
		double nbEventDesired=4;
		if (nbEvents==0 || nbEvents>nbEventDesired*2) {
			S=0;
		} else {
			S=1.0-(double)(Math.abs(nbEvents-nbEventDesired))/nbEventDesired;
		}
		if (printTrace) {
			System.out.println(nbEvents+"  "+ S+ "  " +nbEventDesired);
		}
		fitness*=1.0/(1.0+S);
		
		//System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		return 0;
	}
	
	public static void main(String args[]) throws Exception {
		GRNModel grn = GRNModel.loadFromFile("SignalProcessing_exp2/run_2013-03-07_12-21-10/grn_766_-39.80738103881875.grn");
		MichalSignalProcessExp2 eval=new MichalSignalProcessExp2();
		double fit=eval.evolveGRN(grn, true, 125, 1000);
//		double fit=eval.evolveGRN(grn, true, 500, 1000);
//		double fit=eval.evolveGRN_zero(grn, true, 1000);
//		double fit=eval.evolveGRN(grn, true, 7500, 15000);
//		double fit=eval.evolveGRN(grn, true, 83.33333, 1000);
//		double fit=eval.evolveGRN_expE(grn, true);
		System.out.println("Fitness="+fit);
	}

}
