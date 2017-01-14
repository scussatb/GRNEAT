package grn;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Random;
import java.util.Vector;


public class GRNModel implements Serializable {
	private static final long serialVersionUID = 1L;

	double beta;
	double delta;
	int currentStep;
	double maxEnhance;
	double maxInhibit;
	double enhanceMatching[][];
	double inhibitMatching[][];

	public Vector<GRNProtein> proteins;

	public GRNModel(Vector<GRNProtein> p, double b, double d) {
		proteins=p;
		beta=b;
		delta=d;

		maxEnhance=0;
		maxInhibit=0;

		updateSignatures();
	}

	public GRNModel() {
		proteins=new Vector<GRNProtein>();
		beta=1;
		delta=1;
	}

	public GRNModel(GRNModel grn) {
		proteins=new Vector<GRNProtein>();
		for (int i=0; i<grn.proteins.size(); i++) {
			proteins.add(new GRNProtein(grn.proteins.get(i)));
		}

		beta=grn.beta;
		delta=grn.delta;

		maxEnhance=grn.maxEnhance;
		maxInhibit=grn.maxInhibit;

		enhanceMatching=new double[grn.enhanceMatching.length][grn.enhanceMatching[0].length];
		for (int i=0; i<grn.enhanceMatching.length; i++) {
			for (int j=0; j<grn.enhanceMatching[i].length; j++) {
				enhanceMatching[i][j]=grn.enhanceMatching[i][j];
			}
		}
		inhibitMatching=new double[grn.inhibitMatching.length][grn.inhibitMatching[0].length];
		for (int i=0; i<grn.inhibitMatching.length; i++) {
			for (int j=0; j<grn.inhibitMatching[i].length; j++) {
				inhibitMatching[i][j]=grn.inhibitMatching[i][j];
			}
		}

		currentStep=grn.currentStep;
	}

	public static GRNModel loadFromFile(String fileName) throws Exception {
		FileInputStream fis = new FileInputStream(fileName);
		ObjectInputStream ois = new ObjectInputStream(fis);

		GRNModel res = (GRNModel)ois.readObject();

		fis.close();
		ois.close();

		return res;
	}

	public Vector<GRNProtein> getProteins() {
		return proteins;
	}

	//		inline unsigned int getNbProteins() {return nbProteins;};
	public void evolve(int nbSteps) {
		double enhance, inhibit;
		int i,j,k;
		double sumConcentration;
		//		for (i=0; i<proteins.size(); i++) {
		//			cerr << proteins[i].concentration << "\t";
		//		}
		//		cerr << proteins.size() << endl;


		for (int step=0; step<nbSteps; step++) {
			Vector<GRNProtein> nextProteins=new Vector<GRNProtein>();

			// Calculating next step protein concentrations
			for (j=0; j<proteins.size(); j++) {
				// For an input protein, the concentration is not calculated
				if (proteins.get(j).type==GRNProtein.INPUT_PROTEIN) {
					nextProteins.add(proteins.get(j));
				} else {
					enhance=0;
					inhibit=0;
					// Calculating enhancing and inhibiting factor for the current protein
					for (k=0; k<proteins.size(); k++) {
						if (proteins.get(k).type!=GRNProtein.OUTPUT_PROTEIN) {
							enhance+=proteins.get(k).concentration*enhanceMatching[j][k];
							inhibit+=proteins.get(k).concentration*inhibitMatching[j][k];
						}
					}
					// if (j=5) cout << enhance << "   " << inhibit << endl;
					// Calculating the next concentration of current protein 
					nextProteins.add(new GRNProtein(proteins.get(j).id, proteins.get(j).type, Math.max(0.0,proteins.get(j).concentration+delta/proteins.size()*(enhance-inhibit)), proteins.get(j).enhancer, proteins.get(j).inhibiter));
				}
			}

			//		for (i=0; i<nextProteins.size(); i++) {
			//			cout << nextProteins[i].id << "=" << nextProteins[i].concentration << "\t";
			//		}
			//		cout << endl;
			// Reajusting proteins concentration so that their sum stay equal to 1
			sumConcentration = 0;
			for (i=0; i<proteins.size(); i++) {
				if (proteins.get(i).type!=GRNProtein.INPUT_PROTEIN) {
					sumConcentration+=nextProteins.get(i).concentration;
				}
			}
			//cout << sumConcentration << "\t" << nbDiv << "\t" << (sumConcentration-1)/nbDiv << "\t";
			if (sumConcentration!=0) {
				for (i=0; i<proteins.size(); i++) {
					if (proteins.get(i).type!=GRNProtein.INPUT_PROTEIN) {
						nextProteins.get(i).concentration/=sumConcentration;
					}
				}
			}

			// Finalizing the step by switching the vector
			proteins=nextProteins;

			currentStep++;
		}

	}

	public void duplicateProteins(int nbDup, double mutProb, Random rng) {
		//srand(time(NULL))
		int nbProteins=proteins.size();
		for (int dup=0; dup<nbDup; dup++) {
			for (int i=0; i<nbProteins; i++) {
/*				proteins.add(new GRNProtein(
						Math.random()>mutProb?proteins.get(i).id:(int)(Math.random()*GRNProtein.IDSIZE),
								proteins.get(i).type,
								proteins.get(i).concentration,
								Math.random()>mutProb?proteins.get(i).enhancer:(int)(Math.random()*GRNProtein.IDSIZE),
										Math.random()>mutProb?proteins.get(i).inhibiter:(int)(Math.random()*GRNProtein.IDSIZE)));*/
				proteins.add(new GRNProtein(
						rng.nextDouble()>mutProb?proteins.get(i).id:(int)(rng.nextDouble()*GRNProtein.IDSIZE),
								proteins.get(i).type,
								proteins.get(i).concentration,
								rng.nextDouble()>mutProb?proteins.get(i).enhancer:(int)(rng.nextDouble()*GRNProtein.IDSIZE),
										rng.nextDouble()>mutProb?proteins.get(i).inhibiter:(int)(rng.nextDouble()*GRNProtein.IDSIZE)));
				
			}
		}
		updateSignatures();

	}


	private void updateSignatures() {
		// calculating signatures
		enhanceMatching = new double[proteins.size()][proteins.size()];
		inhibitMatching = new double[proteins.size()][proteins.size()];
		for (int j=0; j<proteins.size(); j++) {
			for (int k=0; k<proteins.size(); k++) {
				enhanceMatching[j][k] = GRNProtein.IDSIZE-Math.abs(proteins.get(j).enhancer-proteins.get(k).id);
				maxEnhance=Math.max(maxEnhance, enhanceMatching[j][k]);
				inhibitMatching[j][k] = GRNProtein.IDSIZE-Math.abs(proteins.get(j).inhibiter-proteins.get(k).id);
				maxInhibit=Math.max(maxInhibit, inhibitMatching[j][k]);
			}
		}
		for (int j=0; j<proteins.size(); j++) {
			for (int k=0; k<proteins.size(); k++) {
				enhanceMatching[j][k]=Math.exp(beta*(double)enhanceMatching[j][k]-(double)maxEnhance);
				inhibitMatching[j][k]=Math.exp(beta*(double)inhibitMatching[j][k]-(double)maxInhibit);
			}
		}
	}

	public void writeToFile(String fileName) throws Exception {
		FileOutputStream fos = new FileOutputStream(fileName);
		ObjectOutputStream oos = new ObjectOutputStream(fos);

		oos.writeObject(this);

		oos.close();
		fos.close();
	}

	public GRNModel copy() {
		GRNModel grn = new GRNModel();
		for (int i=0; i<this.proteins.size(); i++) {
			grn.proteins.add(new GRNProtein(this.proteins.get(i)));
		}

		grn.beta=this.beta;
		grn.delta=this.delta;

		grn.maxEnhance=this.maxEnhance;
		grn.maxInhibit=this.maxInhibit;

		grn.enhanceMatching=new double[this.enhanceMatching.length][this.enhanceMatching[0].length];
		for (int i=0; i<this.enhanceMatching.length; i++) {
			for (int j=0; j<this.enhanceMatching[i].length; j++) {
				grn.enhanceMatching[i][j]=this.enhanceMatching[i][j];
			}
		}
		grn.inhibitMatching=new double[this.inhibitMatching.length][this.inhibitMatching[0].length];
		for (int i=0; i<this.inhibitMatching.length; i++) {
			for (int j=0; j<this.inhibitMatching[i].length; j++) {
				grn.inhibitMatching[i][j]=this.inhibitMatching[i][j];
			}
		}

		grn.currentStep=this.currentStep;

		return grn;

	}

	public void reset() {
		for (GRNProtein p : proteins) {
			p.concentration=1.0/(double)proteins.size();
		}
	}

	public String toString() {
		String res="";
		for (GRNProtein p : proteins) {
			res+=p.toString();
		}
		return res;
	}

	public long size() {
		return proteins.size();
	}
	
	public double distanceTo(GRNModel g, boolean compareDynamicsCoeff, double protCoef, double enhCoef, double inhCoef, double betaMax, double betaMin, double deltaMax, double deltaMin ) {
		if (g==null) {
			System.err.println("(GRNModel.distanceTo) Genome is null!");
		}
		double distance=0;
		GRNModel gs, gl;
		if (size()<g.size()) {
			gs=this;
			gl=g;
		} else {
			gs=g;
			gl=this;
		}
		
		// comparing all types
		for (GRNProtein gi1 : gl.proteins) {
			double minD = Double.MAX_VALUE; 
			for (GRNProtein gi2 : gs.proteins) {
				double d = gi1.distanceTo( gi2, protCoef, enhCoef, inhCoef);
				if( gi1.type == gi2.type && d < minD )
					minD = d;
			}
			distance+= minD;
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
	

	static public void main(String args[]) {
		/* For compiling a JAR that outputs the distance between 2 genomes */
		try {
			GRNModel g1 = GRNModel.loadFromFile( args[0] );
			GRNModel g2 = GRNModel.loadFromFile( args[1] );
		
			boolean compDynamicsCoef = Boolean.parseBoolean( args[2] );
			double protCoef = Double.parseDouble( args[3] );
			double enhCoef  = Double.parseDouble( args[4] );
			double inhCoef  = Double.parseDouble( args[5] );
			double betaMax  = Double.parseDouble( args[6] );
			double betaMin  = Double.parseDouble( args[7] );
			double deltaMax  = Double.parseDouble( args[8] );
			double deltaMin  = Double.parseDouble( args[9] );
			
			System.out.println( g1.distanceTo( g2, compDynamicsCoef, protCoef, enhCoef, inhCoef, betaMax, betaMin, deltaMax, deltaMin ) );
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

}
