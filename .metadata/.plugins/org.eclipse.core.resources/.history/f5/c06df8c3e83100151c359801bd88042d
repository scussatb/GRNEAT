package evolver;

import grn.GRNProtein;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;

public class GRNGenome implements Comparable {	
	protected Hashtable<Integer, GRNGene> inputGenes;
	protected Hashtable<Integer, GRNGene> outputGenes;
	protected Hashtable<Long, GRNGene> regulatoryGenes;
	protected Hashtable<Long, GRNGene> allGenes;
	protected double beta;
	protected double delta;
	protected double betaMax=2.0;
	protected double betaMin=0.5;
	protected double deltaMax=2.0;
	protected double deltaMin=0.5;
//	protected int nbEval=0;
	
	protected boolean hasBeenEvaluated=false;
	protected double lastFitnessValue;
	protected double adjustedFitness;
	public double parentBestFitness;
		
	public GRNGenome() {
		inputGenes=new Hashtable<Integer, GRNGene>();
		outputGenes=new Hashtable<Integer, GRNGene>();
		regulatoryGenes=new Hashtable<Long, GRNGene>();
		allGenes=new Hashtable<Long, GRNGene>();
		beta=1.0;
		delta=1.0;
	}
	
	public GRNGenome(Enumeration<GRNGene> nGenes, double nBeta, double nDelta) {
		this();
		while (nGenes.hasMoreElements()) {
			addGene(nGenes.nextElement());
		}
		beta=nBeta;
		delta=nDelta;
	}
	
	public GRNGenome(List<GRNGene> nGenes, double nBeta, double nDelta) {
		this();
		for (GRNGene g : nGenes) {
			addGene(g);
		}
		beta=nBeta;
		delta=nDelta;
	}
	
	public void addGene(GRNGene nGene) {
		GRNGene nGeneClone=nGene.clone();
		switch (nGeneClone.protType) {
		case GRNProtein.INPUT_PROTEIN:
			inputGenes.put(nGeneClone.getConnectTo(), nGeneClone);
			break;
		case GRNProtein.OUTPUT_PROTEIN:
			outputGenes.put(nGeneClone.getConnectTo(), nGeneClone);
			break;
		case GRNProtein.REGULATORY_PROTEIN:
			regulatoryGenes.put(nGeneClone.getId(), nGeneClone);
			break;
		}
		allGenes.put(nGeneClone.getId(), nGeneClone);
		hasBeenEvaluated=false;
	}

	public int compareTo(Object o) {
		if (o instanceof GRNGenome) {
			GRNGenome genome2 = (GRNGenome) o;
			int distance=0;
			// comparing inputs
			for (GRNGene gi1 : inputGenes.values()) {
				for (GRNGene gi2 : genome2.inputGenes.values()) {
					if (gi1.compareTo(gi2)!=0) {
						distance++;
					}
				}
			}
			// comparing outputs
			for (GRNGene go1 : outputGenes.values()) {
				for (GRNGene go2 : genome2.outputGenes.values()) {
					if (go1.compareTo(go2)!=0) {
						distance++;
					}
				}
			}
			// comparing regulatories
			for (GRNGene gr1 : regulatoryGenes.values()) {
				for (GRNGene gr2 : genome2.regulatoryGenes.values()) {
					if (gr1.compareTo(gr2)!=0) {
						distance++;
					}
				}
			}
			
			return distance;
		} else {
			return Integer.MAX_VALUE;
		}
	}
	
	public double distanceTo(GRNGenome g, boolean compareDynamicsCoeff) {
		if (g==null) {
			System.err.println("(GRNGenome.distanceTo) Genome is null!");
		}
		double distance=0;
		GRNGenome gs, gl;
		if (size()<g.size()) {
			gs=this;
			gl=g;
		} else {
			gs=g;
			gl=this;
		}
		// comparing inputs
		for (GRNGene gi1 : gl.inputGenes.values()) {
			GRNGene gi2=gs.inputGenes.get(gi1.connectTo);
			if (gi2==null) {
				System.err.println("(GRNGenome.distanceTo) gi2 is null! gi1.connectTo="+gi1.connectTo);
				System.err.println(gs.toString());
			}
			distance+=gi1.distanceTo(gi2);
		}
		// comparing outputs
		for (GRNGene go1 : gl.outputGenes.values()) {
			GRNGene go2=gs.outputGenes.get(go1.connectTo);
			distance+=go1.distanceTo(go2);
		}
		// comparing regulatories
		for (GRNGene gr1 : gl.regulatoryGenes.values()) {
			double minDist=Double.MAX_VALUE;
			for (GRNGene gr2 : gs.regulatoryGenes.values()) {
				minDist=Math.min(minDist,gr1.distanceTo(gr2));
			}
			distance+=minDist;
		}
		if (compareDynamicsCoeff) {
			// take beta and delta to the distance calculation
			distance+=Math.abs(beta-g.beta)/(betaMax-betaMin);
			distance+=Math.abs(delta-g.delta)/(deltaMax-deltaMin);
			return ((double)distance)/(gl.size()+2.0);
		} else {
			return ((double)distance)/gl.size();
		}
	}

	public double getBetaMax() {
		return betaMax;
	}

	public void setBetaMax(double betaMax) {
		this.betaMax = betaMax;
	}

	public double getBetaMin() {
		return betaMin;
	}

	public void setBetaMin(double betaMin) {
		this.betaMin = betaMin;
	}

	public double getDeltaMax() {
		return deltaMax;
	}

	public void setDeltaMax(double deltaMax) {
		this.deltaMax = deltaMax;
	}

	public double getDeltaMin() {
		return deltaMin;
	}

	public void setDeltaMin(double deltaMin) {
		this.deltaMin = deltaMin;
	}

	public GRNGenome clone() {
		return new GRNGenome(allGenes.elements(), beta, delta);
	}

	public double getLastFitness() {
		return lastFitnessValue;
	}
	
	public void setNewFitness(double aFitness) {
		this.lastFitnessValue=aFitness;
		hasBeenEvaluated=true;
//		nbEval++;
//		System.err.println(this+" has been evaluated "+nbEval+" times");
	}
	
	public boolean hasBeenEvaluated() {
		return hasBeenEvaluated;
	}
	
	public GRNGene containsGene(GRNGene aGene) {
		return allGenes.get(aGene.getId());
	}
	
	public GRNGene containsGene(Integer geneId) {
		return allGenes.get(geneId);
	}
	
	public GRNGene getInputGeneConnectTo(Integer sensorId) {
		return inputGenes.get(sensorId);
	}
	
	public GRNGene getOutputGeneConnectTo(Integer outputId) {
		return outputGenes.get(outputId);
	}
	
	public Collection<GRNGene> getInputGenes() {
		return inputGenes.values();
	}
	
	public Collection<GRNGene> getOutputGenes() {
		return outputGenes.values();
	}
	
	public Collection<GRNGene> getRegulatoryGenes() {
		return regulatoryGenes.values();
	}
	
	public Collection<GRNGene> getAllGenes() {
		return allGenes.values();
	}
	
	public double getBeta() {
		return beta;
	}
	
	public void setBeta(double nBeta) {
		if (beta!=nBeta) {
			hasBeenEvaluated=false;
			beta=nBeta;
		}
	}
	
	public double getDelta() {
		return delta;
	}
	
	public void setDelta(double nDelta) {
		if (delta!=nDelta) {
			hasBeenEvaluated=false;
			delta=nDelta;
		}
	}
	
	public int size() {
		return allGenes.size();
	}
	
	public String toString() {
		String res="{";
		for (GRNGene g : allGenes.values()) {
			res+=g.toString();
		}
		res+=" ; "+beta+" ; "+delta+"}";
		return res;
	}
	
	public void removeAllGenes() {
		hasBeenEvaluated=false;
		allGenes.clear();
	}
	
	public void removeGene(GRNGene g) {
		switch (g.protType) {
		case GRNProtein.INPUT_PROTEIN:
			inputGenes.remove(g.connectTo);
			break;
		case GRNProtein.OUTPUT_PROTEIN:
			inputGenes.remove(g.connectTo);
			break;
		case GRNProtein.REGULATORY_PROTEIN:
			inputGenes.remove(g.getId());
			break;
		}
		allGenes.remove(g.getId());
		hasBeenEvaluated=false;
	}
	
	public boolean removeRandomlyRegulatoryGene(Random rng) {
		if (regulatoryGenes.size()==0) return false;
		//int deleteIndex=(int)(Math.random()*regulatoryGenes.size());
		int deleteIndex=(int)(rng.nextDouble()*regulatoryGenes.size());
		GRNGene g=(GRNGene)regulatoryGenes.values().toArray()[deleteIndex];
		removeGene(g);
		hasBeenEvaluated=false;
		return true;
	}
	
	public double getAdjustedFitness() {
		return adjustedFitness;
	}
	
	public void setAdjustedFitness(double nAdjFit) {
		adjustedFitness=nAdjFit;
	}
	
	public void hasBeenModified() {
		hasBeenEvaluated=false;
	}
}
