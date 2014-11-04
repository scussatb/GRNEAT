package operators;

import evolver.GRNGene;
import evolver.GRNGenome;
import grn.GRNProtein;

import java.util.Collection;
import java.util.Random;


public class GRNGeneMutationOperator extends GRNMutationOperator {
	public GRNGeneMutationOperator() {
		name="ChangeGeneMutation(prob="+probability+")";
	}

	public GRNGeneMutationOperator(double prob) {
		this.probability=prob;
		name="ChangeGeneMutation(prob="+probability+")";
	}
	
	@Override
	public boolean mutateByModifying(GRNGenome aGenome, Random rng) {
		Collection<GRNGene> allGenes=aGenome.getAllGenes();
		int mutateIndex = (int)(rng.nextDouble()*allGenes.size()+2);
		for (GRNGene g : allGenes) {
			if (mutateIndex == 0) {
				switch ((int)(rng.nextDouble()*3)) {
				case 0: 
					g.setProtId((int)(rng.nextDouble()*GRNProtein.IDSIZE));
					break;
				case 1: 
					g.setProtEnh((int)(rng.nextDouble()*GRNProtein.IDSIZE));
					break;
				case 2: 
					g.setProtInh((int)(rng.nextDouble()*GRNProtein.IDSIZE));
					break;
				default: 
					g.setDisabled(rng.nextDouble()<0.5);
					break;
				}
			}
			mutateIndex--;
		}
		if (mutateIndex==0) {
			// mutate beta
			aGenome.setBeta(rng.nextDouble()*(aGenome.getBetaMax()-aGenome.getBetaMin())+aGenome.getBetaMin());
		}
		if (mutateIndex==1) {
			// mutate delta
			aGenome.setDelta(rng.nextDouble()*(aGenome.getDeltaMax()-aGenome.getDeltaMin())+aGenome.getDeltaMin());
		}
		aGenome.hasBeenModified();
		return true;
	}

}
