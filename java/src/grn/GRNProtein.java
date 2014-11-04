package grn;

import java.io.Serializable;

import evolver.GRNGene;


public class GRNProtein implements Serializable {
	private static final long serialVersionUID = 1L;

	public static final int INPUT_PROTEIN = 0;
	public static final int REGULATORY_PROTEIN = 2;
	public static final int OUTPUT_PROTEIN = 1;
	public static final int IDSIZE = 32;
	
	public int id;
	public double concentration;
	public int enhancer;
	public int inhibiter;
	public int type;
	
	public GRNProtein() {
		id=0;
		concentration=0;
		enhancer=0;
		inhibiter=0;
		type=0;
	}
	
	public GRNProtein(int ID, int typ, double conc, int enh, int inh) {
		id=ID;
		concentration=conc;
		enhancer=enh;
		inhibiter=inh;
		type=typ;
	}
	
	public GRNProtein(GRNProtein p) {
		id=p.id;
		concentration=p.concentration;
		enhancer=p.enhancer;
		inhibiter=p.inhibiter;
		type=p.type;
	}
	
	public void setConcentration(double newConc) {
		concentration = newConc;
	}
	
	double getConcentration() {
		return concentration;
	}
	
	int getID() {
		return id;
	}
	
	public String toString() {
		return "["+id+","+enhancer+","+inhibiter+","+concentration+","+type+"]";
	}
	
	public double distanceTo(GRNProtein g, double protCoef, double enhCoef, double inhCoef ) {
		if (g==null) {
			System.err.println("(GRNProtein.distanceTo)  g is null!!!");
		}
		return (protCoef*(double)(Math.abs(id-g.id))+
				enhCoef*(double)(Math.abs(enhancer-g.enhancer))+
				inhCoef*(double)(Math.abs(inhibiter-g.inhibiter)))/GRNProtein.IDSIZE;
	}

}
