package evolver;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.Random;
import java.util.TreeMap;

import operators.GRNCrossoverOperator;
import sun.misc.Signal;

import evaluators.GRNGenomeEvaluator;

public class Species implements Cloneable {
	protected ArrayList<GRNGenome> genomes=null;
	protected boolean genomesAreSorted=true;
	protected double speciationThreshold=0.3;
	protected GRNGenome representativeGenome=null;
	protected GRNGenomeEvaluator evaluator=null;
	protected double fitnessSum=0;
	protected GRNGenome bestGenome=null;
	protected double bestFitness=-Double.MAX_VALUE;
	public boolean representativeIsFirst=true;		
	public boolean replacement=true;	
	public boolean useDynamicsGenomeForDistance=true;
	public int speciesId;
	static public int speciesNextId=0;
		
	public Species(GRNGenome aRepresentative, double aSpeciationThreshold, boolean addRepresentative, GRNGenomeEvaluator anEvaluator) {
		evaluator=anEvaluator;
		speciationThreshold=aSpeciationThreshold;
		representativeGenome=aRepresentative.clone();
		genomes=new ArrayList<GRNGenome>();
		if (addRepresentative) {
			addGenome(aRepresentative, false);
		}
		speciesId=speciesNextId++;
	}
	
	public void recomputeFitness() {
		for( int k = 0; k < genomes.size(); k++ ) {
			//genomes.set(k, evaluator.evaluate( genomes.get(k) ));
			evaluator.evaluate( genomes.get(k) );
		}
	}
	
	public boolean addGenome(GRNGenome aGenome, boolean testSpecies) {
//		System.err.println(representativeGenome.distanceTo(aGenome, true));
		if (aGenome==null) {
			System.err.println("(Species.addGenome) Genome is null!");
		}
		if (representativeGenome==null) {
			System.err.println("(Species.addGenome) Representative is null!");	
		}
		if (!testSpecies || representativeGenome.distanceTo(aGenome, useDynamicsGenomeForDistance)<=speciationThreshold) {
			if (!aGenome.hasBeenEvaluated()) {
				evaluator.evaluate(aGenome);
			}
			if( representativeGenome == null ) {
				representativeGenome = aGenome;
			}
			genomes.add(aGenome);
			//System.err.println(genomes.size());
			genomesAreSorted=false;
			fitnessSum+=aGenome.getLastFitness();
			if (aGenome.getLastFitness()>bestFitness) {
				bestFitness=aGenome.getLastFitness();
				bestGenome=aGenome;
				if (!representativeIsFirst) {
					representativeGenome=bestGenome.clone();
				}
			}
			return true;
		} else {
			return false;
		}
	}
	
	public Collection<GRNGenome> getGenomes() {
		return genomes;
	}
	
	public int size() {
		return genomes.size();
	}
	
	public void removeWorstGenome() {
		if (!genomesAreSorted) {
			sortGenomes();
		}
		testSorted();
		fitnessSum-=genomes.get(genomes.size()-1).getLastFitness();
		genomes.remove(genomes.size()-1);
	}
	
	public double getWorstFitness() {
		if (!genomesAreSorted) {
			sortGenomes();
		}
		testSorted();
		return genomes.get(genomes.size()-1).getLastFitness();
	}
	
	public double getBestFitness() {
		if (!genomesAreSorted) {
			sortGenomes();
		}
		testSorted();
		return genomes.get(0).getLastFitness();
	}
	
	public GRNGenome getBestGenome() {
		if (!genomesAreSorted) {
			sortGenomes();
		}
		testSorted();
		return genomes.get(0);		
	}
	
	public boolean testSorted() {
		GRNGenome last=genomes.get(0);
		boolean isSorted=true;
		boolean comparSort=true;
		for (int i=1; i<genomes.size() && isSorted; i++) {
			isSorted=genomes.get(i).getLastFitness()<=last.getLastFitness();
			comparSort=(new GRNGenomeFitnessComparator()).compare(last,genomes.get(i))>0;
			last=genomes.get(i);
		}
//		if (!isSorted) {
//			System.err.println(isSorted+"  "+comparSort);
//		}
		return isSorted;
	}
	
	public GRNGenome rouletteWheelSelect(Random rng) {
		// Roulette wheel selection

		//double rouletteValue=Math.random()*fitnessSum;
		double rouletteValue=rng.nextDouble()*fitnessSum;
		if (!genomesAreSorted) {
			sortGenomes();
		}
		testSorted();
		Iterator<GRNGenome> genome=genomes.iterator();
		GRNGenome res=genome.next();
		while (res.getLastFitness()<rouletteValue && genome.hasNext()) {
			rouletteValue-=res.getLastFitness();
			res=genome.next();
		}
		return res;
	}
	
	public GRNGenome randomSelect(Random rng) {
		//return genomes.get((int)(Math.random()*genomes.size()));
		return genomes.get((int)(rng.nextDouble()*genomes.size()));
	}
	
	public GRNGenome tournamentSelect(int nParticipant, Random rng) {
		ArrayList<GRNGenome> tournament=new ArrayList<GRNGenome>();
		while (tournament.size()<nParticipant) {
			// select randomly a participant
			GRNGenome p = randomSelect(rng);
			if (!tournament.contains(p) || replacement) {
				tournament.add(p);
			}
		}
		Collections.sort(tournament, new GRNGenomeFitnessComparator());
		return tournament.get(0);
	}
	
	public void sortGenomes() {
		Collections.sort(genomes,new GRNGenomeFitnessComparator());
		genomesAreSorted=true;
		//System.err.println(genomes.get(0).getLastFitness()>=genomes.get(1).getLastFitness());
	}
	
	private class GRNGenomeFitnessComparator implements Comparator<GRNGenome> {
		@Override
		public int compare(GRNGenome o1, GRNGenome o2) {
			double f1, f2;
			if (!o1.hasBeenEvaluated()) {
				f1=evaluator.evaluate(o1);
			}
			f1=o1.getLastFitness();
			if (!o2.hasBeenEvaluated()) {
				f2=evaluator.evaluate(o2);
			}
			f2=o2.getLastFitness();
			return f1<f2?1:f1>f2?-1:0;
		}
		
	}
	
	public GRNGenome getGenome(int index) {
		return genomes.get(index);
	}
	
	public Species clone() {
		Species s=new Species(representativeGenome, speciationThreshold, false, evaluator);
		s.genomes.addAll(this.genomes);
		s.genomesAreSorted=s.genomesAreSorted;
		return s;
	}
	
	public void chooseRandomRepresentative( Random RNG ) {		
		representativeGenome=genomes.get( RNG.nextInt( genomes.size() ) ).clone();
	}
	
	public void chooseBestRepresentative( Random RNG ) {		
		representativeGenome=bestGenome.clone();
	}
	
	public void removeAllGenomes() {
		fitnessSum=0;
		bestGenome=null;
		bestFitness = -Double.MAX_VALUE;
		//representativeGenome=null;
		genomes.clear();
	}

	public double getFitnessSum() {
		return fitnessSum;
	}
	
	public double getAvgDist2Rep() {
		double sum = 0;
		for( GRNGenome aGenome : genomes ) {
			//double d = representativeGenome.distanceTo(aGenome, true); System.err.println( d + "\t" + ( representativeGenome.distanceTo(aGenome, true)<=speciationThreshold ) + "\t" + speciationThreshold ); sum += d;						
			sum += representativeGenome.distanceTo(aGenome, useDynamicsGenomeForDistance);
		}
		return sum / genomes.size();
	}
	
	public double getFitnessAvg() {
		return fitnessSum/size();
	}
}
