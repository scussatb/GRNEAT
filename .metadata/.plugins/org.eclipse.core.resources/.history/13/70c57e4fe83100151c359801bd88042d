package evaluators;

import evolver.GRNGenome;
import grn.GRNModel;
import grn.PopGRNModel;

public class DoublingFrequencyEvaluator extends GRNGenomeEvaluator {
	double coef=3.0;
	
	public DoublingFrequencyEvaluator() {
		numGRNInputs=1;
		numGRNOutputs=1;
		name="DoublingFrequency";
	}

	@Override
	public double evaluate(GRNGenome aGenome) {
		PopGRNModel grn=buildGRNFromGenome(aGenome);
		double fitness=evolveGRN(grn, false, 125.0, 1000)+
				evolveGRN(grn, false, 500, 1000)+
				evolveGRN_zero(grn, false, 1000);
		
//		System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		aGenome.setNewFitness(fitness);

		GRNGenomeEvaluator.numEvaluations++;
		return fitness;
	}
	
	double evolveGRN_zero(PopGRNModel grn, boolean printTrace, int nStepMax) {
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
	
	double evolveGRN(PopGRNModel grn, boolean printTrace, double halfFreq, int nStepMax) {
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
				System.out.println(nStep+"\t"+Math.min(1.0,ots[9])+"\t"+dt+"\t"+(Math.sin(nStep*Math.PI/halfFreq-Math.PI/2)/2+0.5));
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
//		if (printTrace) {
//			System.out.println(nbEvents+"  "+ S+ "  " +nbEventDesired);
//		}
		fitness*=1.0/(1.0+S);
		
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
		// Generalization
		DoublingFrequencyEvaluator eval=new DoublingFrequencyEvaluator();
		double freqs[]=	{25, 	50, 	75, 	100, 	150, 	200, 	300, 	600, 	750, 	1000, 	2500, 	5000};
		int nSteps[]=	{1000, 	1000, 	1000, 	1000, 	1000, 	1000, 	1000,	2000,	3000,	5000,	10000, 	20000};
		GRNModel greatGRN[]=new GRNModel[26];
		GRNModel gaGRN[]=new GRNModel[26];
		GRNModel esGRN[]=new GRNModel[26];
		String greatFiles[]={
				"grn_400_-24.032413792966416.grn",	"grn_400_-36.18069885237963.grn",	"grn_400_-49.683578952349976.grn",
				"grn_400_-26.112305596316183.grn",	"grn_400_-38.89061830096909.grn",	"grn_400_-53.542204357043744.grn",
				"grn_400_-27.270695985811074.grn",	"grn_400_-40.6689486646976.grn",	"grn_400_-57.223087472391796.grn",
				"grn_400_-28.312808462393164.grn",	"grn_400_-42.37529203086706.grn",	"grn_400_-63.92041219782782.grn",
				"grn_400_-29.63803876716149.grn",	"grn_400_-42.447527821070686.grn",	"grn_400_-64.42734769024332.grn",
				"grn_400_-30.061168348913856.grn",	"grn_400_-43.84709896427495.grn",	"grn_400_-67.73794638177456.grn",
				"grn_400_-31.39998462760445.grn",	"grn_400_-44.49552387295133.grn",	"grn_400_-73.21817349521785.grn",
				"grn_400_-32.912628867005225.grn",	"grn_400_-46.38323716735894.grn",	"grn_400_-84.97704796945374.grn",
				"grn_400_-35.56781845622969.grn",	"grn_400_-46.77947040104973.grn"
		};
		String gaFiles[]={
				"grn_400_-102.12137449791348.grn",	"grn_400_-112.53453690768704.grn",	"grn_400_-141.27718695626322.grn",
				"grn_400_-102.82442929332575.grn",	"grn_400_-115.7039600392347.grn",	"grn_400_-154.28011585694824.grn",
				"grn_400_-103.73962545073668.grn",	"grn_400_-116.38050870294144.grn",	"grn_400_-178.05868150165992.grn",
				"grn_400_-104.68399689381289.grn",	"grn_400_-117.32167642089448.grn",	"grn_400_-208.47805762059133.grn",
				"grn_400_-105.6920548414619.grn",	"grn_400_-120.2923106290028.grn",	"grn_400_-227.40653901468474.grn",
				"grn_400_-106.92543794634534.grn",	"grn_400_-120.37455007511271.grn",	"grn_400_-284.21871700686336.grn",
				"grn_400_-107.37579208887657.grn",	"grn_400_-120.55776517651988.grn",	"grn_400_-75.1036524885922.grn",
				"grn_400_-108.17765830768083.grn",	"grn_400_-122.42478888253218.grn",	"grn_400_-99.31781410623695.grn",
				"grn_400_-109.75912476747042.grn",	"grn_400_-123.85488384879586.grn"
		};
		String esFiles[]={
				"grn_400_-115.9205504345521.grn",	"grn_400_-31.623983169294878.grn",	"grn_400_-55.68528029408594.grn",
				"grn_400_-119.78949654978919.grn",	"grn_400_-39.89252993638848.grn",	"grn_400_-55.9036173159826.grn",
				"grn_400_-140.29655880602.grn",	"grn_400_-42.12274277390592.grn",	"grn_400_-56.834260167166555.grn",
				"grn_400_-140.90696982292357.grn",	"grn_400_-42.5213474358587.grn",	"grn_400_-57.70810040844396.grn",
				"grn_400_-17.612311802917183.grn",	"grn_400_-42.75268622968957.grn",	"grn_400_-58.874799089305775.grn",
				"grn_400_-22.62268508767353.grn",	"grn_400_-44.02606416025216.grn",	"grn_400_-63.381544278393264.grn",
				"grn_400_-22.661231975204014.grn",	"grn_400_-44.226094918580294.grn",	"grn_400_-66.91868344723783.grn",
				"grn_400_-26.640982270052017.grn",	"grn_400_-48.65906720712924.grn",	"grn_400_-87.02573938628746.grn",
				"grn_400_-27.986204483660597.grn",	"grn_400_-51.23359974959173.grn",
		};
		double greatAvg[]=new double[greatGRN.length+1];
		double gaAvg[]=new double[gaGRN.length+1];
		for (int i=0; i<greatGRN.length; i++) {
			greatGRN[i]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher016-hyperion006/Generalization/DF/GREAT/"+greatFiles[i]);
			gaGRN[i]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher016-hyperion006/Generalization/DF/GA/"+gaFiles[i]);
			esGRN[i]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher016-hyperion006/Generalization/DF/ES/"+esFiles[i]);
		}
		System.out.println("GREAT");
		for (int j=0; j<greatGRN.length; j++) {
			
			System.out.print("\t"+greatFiles[j]);
		}
		System.out.println();
		for (int i=0; i<freqs.length; i++) {
			System.out.print(freqs[i]);
			for (int j=0; j<greatGRN.length; j++) {
				double greatCurFit=-eval.evolveGRN(greatGRN[j], false, freqs[i], nSteps[i]);
				System.out.print("\t"+greatCurFit);
			}
			System.out.println();
		}
		System.out.println("\n\nGA");
		for (int j=0; j<gaGRN.length; j++) {
			System.out.print("\t"+gaFiles[j]);
		}
		System.out.println();
		for (int i=0; i<freqs.length; i++) {
			System.out.print(freqs[i]);
			for (int j=0; j<gaGRN.length; j++) {
				double gaCurFit=-eval.evolveGRN(gaGRN[j], false, freqs[i], nSteps[i]);
				System.out.print("\t"+gaCurFit);
			}
			System.out.println();
		}
		System.out.println("\n\nES");
		for (int j=0; j<esGRN.length; j++) {
			System.out.print("\t"+esFiles[j]);
		}
		System.out.println();
		for (int i=0; i<freqs.length; i++) {
			System.out.print(freqs[i]);
			for (int j=0; j<esGRN.length; j++) {
				double esCurFit=-eval.evolveGRN(esGRN[j], false, freqs[i], nSteps[i]);
				System.out.print("\t"+esCurFit);
			}
			System.out.println();
		}
		System.exit(0);
		
		
		
/*		String grnFile="default.grn";
		double halfFreq=125;
		int steps=1000;
		int mode=0;
				
		for (int k=0; k<args.length; k++) {
			if (args[k].compareTo("mode") == 0) {
				k++;
				mode=Integer.parseInt(args[k]);
			} else if (args[k].compareTo("grn")==0) {
				k++;
				grnFile=args[k];
			} else if (args[k].compareTo("halfFreq")==0) {
				k++;
				halfFreq=Double.parseDouble(args[k]);
			} else if (args[k].compareTo("step")==0) {
				k++;
				steps=Integer.parseInt(args[k]);
			} else {
				System.out.println("Unrecognized option: "+args[k]);
			}
		}
		GRNModel grn = GRNModel.loadFromFile(grnFile);
		DoublingFrequencyEvaluator eval=new DoublingFrequencyEvaluator();
		double fit=0;
		switch (mode) {
		case 0:
			fit=eval.evolveGRN(grn, true, halfFreq, steps);
			break;
		case 1:
			fit=eval.evolveGRN_expE(grn, true);
			break;
		case 2:
			fit=eval.evolveGRN_zero(grn, true, steps);
		default:
			System.out.println("Unrecognized mode: "+mode);
			break;
		}
//		System.out.println(grn.toString());
		
		
		
//		=eval.evolveGRN(grn, true, 125, 1000);
//		double fit=eval.evolveGRN(grn, true, 500, 1000);
//		double fit=eval.evolveGRN_zero(grn, true, 1000);
//		double fit=eval.evolveGRN(grn, true, 7500, 15000);
//		double fit=eval.evolveGRN(grn, true, 83.33333, 1000);
//		double fit=eval.evolveGRN_expE(grn, true);
//		System.out.println("Fitness="+fit);*/
	}

}
