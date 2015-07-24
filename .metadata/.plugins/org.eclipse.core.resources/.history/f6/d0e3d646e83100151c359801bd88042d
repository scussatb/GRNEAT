package evaluators;

import java.util.Vector;

import evolver.GRNGene;
import evolver.GRNGenome;
import grn.GRNModel;
import grn.GRNProtein;
import grn.PopGRNModel;

public abstract class GRNGenomeEvaluator {
	public int generation=0;// some evaluators are dynamic problems dependent on generation number
	public boolean nonCacheable = false;// some evaluators cannot reuse fitness (i.e. dynamic problems)

	public int numGRNInputs=0;
	public int numGRNOutputs=0;
	public static int numEvaluations=0;

	public String name="SuperClass";
	
	public abstract double evaluate(GRNGenome aGenome);
	
	public PopGRNModel buildGRNFromGenome(GRNGenome aGenome) {
		Vector<GRNProtein> prots=new Vector<GRNProtein>();
		for (GRNGene gi : aGenome.getInputGenes()) {
			prots.add(gi.getProtein());
		}
		for (GRNGene go : aGenome.getOutputGenes()) {
			prots.add(go.getProtein());
		}
		for (GRNGene gr : aGenome.getRegulatoryGenes()) {
			prots.add(gr.getProtein());
		}
		PopGRNModel p=new PopGRNModel(prots, aGenome.getBeta(), aGenome.getDelta());
		return p;
	}

}
