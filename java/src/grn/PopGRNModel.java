package grn;

import evolver.GRNGene;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Random;
import java.util.Vector;


public class PopGRNModel extends GRNModel implements Serializable {
	private static final long serialVersionUID = 1L;

	double beta;
	double delta;
	int currentStep;
	double maxEnhance;
	double maxInhibit;
	double enhanceMatching[][];
	double inhibitMatching[][];

	public Vector<GRNProtein> proteins;

	public PopGRNModel(Vector<GRNProtein> p, double b, double d) {
		proteins=p;
		beta=b;
		delta=d;

		maxEnhance=0;
		maxInhibit=0;
		
		updateSignatures();
	}

	public PopGRNModel() {
		proteins=new Vector<GRNProtein>();
		beta=1;
		delta=1;
	}

	public PopGRNModel(PopGRNModel grn) {
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

	public static PopGRNModel loadFromFile(String fileName) throws Exception {
		FileInputStream fis = new FileInputStream(fileName);
		ObjectInputStream ois = new ObjectInputStream(fis);

		PopGRNModel res = (PopGRNModel)ois.readObject();

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
				enhanceMatching[j][k] = Integer.bitCount(proteins.get(j).enhancer^proteins.get(k).id);
				maxEnhance=Math.max(maxEnhance, enhanceMatching[j][k]);
				inhibitMatching[j][k] = Integer.bitCount(proteins.get(j).inhibiter^proteins.get(k).id);
				maxInhibit=Math.max(maxInhibit, inhibitMatching[j][k]);
			}
		}
		for (int j=0; j<proteins.size(); j++) {
			for (int k=0; k<proteins.size(); k++) {
				enhanceMatching[j][k]=Math.exp(beta*((double)enhanceMatching[j][k]-(double)maxEnhance));
				inhibitMatching[j][k]=Math.exp(beta*((double)inhibitMatching[j][k]-(double)maxInhibit));
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

	public PopGRNModel copy() {
		PopGRNModel grn = new PopGRNModel();
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
	
	public double distanceTo(PopGRNModel g, boolean compareDynamicsCoeff, double protCoef, double enhCoef, double inhCoef, double betaMax, double betaMin, double deltaMax, double deltaMin ) {
		if (g==null) {
			System.err.println("(GRNModel.distanceTo) Genome is null!");
		}
		double distance=0;
		PopGRNModel gs, gl;
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
		for (int i=0; i<128; i++) {
			for (int j=0; j<128; j++) {
				System.out.println(i+"\t"+j+"\t"+Integer.bitCount(i^j));
			}
		}
	}
/*		Random rnd=new Random();
		int nProts=rnd.nextInt(50);
		Vector<GRNProtein> prots=new Vector<GRNProtein>();
		for (int i=0; i<nProts; i++) {
			prots.add(new GRNProtein(rnd.nextInt(GRNProtein.IDSIZE), GRNProtein.REGULATORY_PROTEIN, 1.0/nProts, rnd.nextInt(GRNProtein.IDSIZE), rnd.nextInt(GRNProtein.IDSIZE)));
		}
		PopGRNModel grn=new PopGRNModel(prots, rnd.nextDouble()*2, rnd.nextDouble()*2);
		for (int i=0; i<10000; i++) {
			grn.evolve(1);
			System.out.print(i+"\t");
			for (int p=0; p<grn.size(); p++) {
				System.out.print(grn.proteins.get(p).concentration+"\t");
			}
			System.out.println();
		}
	}*/
}
