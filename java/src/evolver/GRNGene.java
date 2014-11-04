package evolver;

import java.util.Random;

import grn.GRNProtein;

public class GRNGene implements Comparable {
	public double protCoef=0.75;
	public double enhCoef=0.125;
	public double inhCoef=0.125;
	
	protected static int currentInnovId=0;
	protected Long geneId;
	protected boolean isDisabled=false;

	protected int protId;
	protected int protEnh;
	protected int protInh;
	protected int protType;
	protected int connectTo;
	
	static long nextId=0;
	
	public boolean isDisabled() {
		return isDisabled;
	}
	
	public void enable() {
		isDisabled=false;
	}
	
	public void disable() {
		isDisabled=true;
	}
	
	public long getId() {
		return geneId;
	}
	
	public GRNGene() {
		this(false, 0, 0, 0, GRNProtein.REGULATORY_PROTEIN, 0);
	}
	
	public GRNGene(boolean nIsDisabled, int nProtId, int nProtEnh, int nProtInh, int nProtType, int nConnectTo) {
		geneId=nextId++;
		isDisabled=nIsDisabled;
		protId=nProtId;
		protEnh=nProtEnh;
		protInh=nProtInh;
		protType=nProtType;
		connectTo=nConnectTo;
	}
	
	public GRNGene(int nProtId, int nProtEnh, int nProtInh) {
		this(false, nProtId, nProtEnh, nProtInh, GRNProtein.REGULATORY_PROTEIN, 0);
	}
	
	public GRNGene(boolean nIsDisabled, int nProtId, int nProtEnh, int nProtInh) {
		this(nIsDisabled, nProtId, nProtEnh, nProtInh, GRNProtein.REGULATORY_PROTEIN, 0);
	}
	
	public GRNProtein getProtein() {
		return new GRNProtein(protId, protType, 0.0, protEnh, protInh);
	}
	
	public GRNGene clone() {
		return new GRNGene(isDisabled,protId, protEnh, protInh, protType, connectTo);
	}

	public int compareTo(Object o) {
		if (o instanceof GRNGene) {
			GRNGene g = (GRNGene) o;
			if (g.protType==protType) {
				return protId-g.protId;
			} else {
				return Integer.MAX_VALUE;
			}
		} else {
			return Integer.MAX_VALUE;
		}
	}
	
	public double distanceTo(GRNGene g) {
		if (g==null) {
			System.err.println("(GRNGene.distanceTo)  g is null!!!");
		}
		if (g.protType==protType) {
			return (protCoef*(double)(Math.abs(protId-g.protId))+
					enhCoef*(double)(Math.abs(protEnh-g.protEnh))+
					inhCoef*(double)(Math.abs(protInh-g.protInh)))/GRNProtein.IDSIZE;
		} else {
			return 1.0;
		}
	}
	
	public int getConnectTo() {
		return connectTo;
	}
	
	public int getProtId() {
		return protId;
	}

	public void setProtId(int protId) {
		this.protId = protId;
	}

	public int getProtEnh() {
		return protEnh;
	}

	public void setProtEnh(int protEnh) {
		this.protEnh = protEnh;
	}

	public int getProtInh() {
		return protInh;
	}

	public void setProtInh(int protInh) {
		this.protInh = protInh;
	}

	public int getProtType() {
		return protType;
	}

	public void setProtType(int protType) {
		this.protType = protType;
	}

	public void setDisabled(boolean isDisabled) {
		this.isDisabled = isDisabled;
	}
	
	public static GRNGene generateRandomGene(int nProtType, int connectTo, Random rng) {
		/*return new GRNGene(
				Math.random()<0.5, 
				(int)(Math.random()*GRNProtein.IDSIZE),
				(int)(Math.random()*GRNProtein.IDSIZE),
				(int)(Math.random()*GRNProtein.IDSIZE),
				nProtType,
				connectTo);*/
		return new GRNGene(
				rng.nextDouble()<0.5, 
				(int)(rng.nextDouble()*GRNProtein.IDSIZE),
				(int)(rng.nextDouble()*GRNProtein.IDSIZE),
				(int)(rng.nextDouble()*GRNProtein.IDSIZE),
				nProtType,
				connectTo);
	}
	
	public static GRNGene generateRandomRegulatoryGene(Random rng) {
		return generateRandomGene(GRNProtein.REGULATORY_PROTEIN, 0, rng);
	}
	
	public String toString() {
		return "["+protId+","+protEnh+","+protInh+","+protType+","+isDisabled+","+connectTo+","+geneId+"]";
	}
	
}
