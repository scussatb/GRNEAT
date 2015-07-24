package evaluators;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;

import evolver.GRNGenome;
import grn.GRNModel;

public class RadbotEvaluator extends GRNGenomeEvaluator {
	public String resourceDirectory = "resources/Radbot";
	RadbotWorld rbw;
	public RadbotEvaluator() {
		numGRNInputs=4;
		numGRNOutputs=4;
		
		name="Radbot";
		
		try {
			rbw=new RadbotWorld(resourceDirectory+"/world_rnd1.rbw");
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	
	public class Radbot {
		public RadbotWorld rbWorld;
		public double levelA;
		public double levelB;
		public double levelS;
		public double levelC;
		public double levelD;

		public double deltaC;
		public double deltaD;

		public double minA, maxA;
		public double minB, maxB;
		public double minS, maxS;

		public double maxSpeed;
		public double actualSpeed;

		public double positionInWorld;

		public double deltaA_B1;
		public double deltaB_B1;
		public double deltaS_B1;
		public double deltaC_B1;
		public double deltaD_B1;
		// Behavior 2
		public double deltaA_B2;
		public double deltaB_B2;
		public double deltaS_B2;	
		public double deltaC_B2;
		public double deltaD_B2;
		// Behavior 3
		public double deltaA_B3;
		public double deltaB_B3;
		public double deltaS_B3;
		public double deltaC_B3;
		public double deltaD_B3;
		public double deltaSpeed_B3; //speed
		// Behavior 4
		public double deltaA_B4;
		public double deltaB_B4;
		public double deltaS_B4;
		public double deltaC_B4;
		public double deltaD_B4;
		public double deltaSpeed_B4; //speed

		public Radbot(RadbotWorld rbw, String fileName) throws Exception {
			rbWorld=rbw;
			String line[];
			BufferedReader reader=new BufferedReader(new FileReader(fileName));
			String readLine;
			while ((readLine=reader.readLine()) != null) {
				line=readLine.split("\t");
				if (line.length==2) {
					if (line[0].compareTo("levelA")==0) {
						levelA=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("minA")==0) {
						minA=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("maxA")==0) {
						maxA=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("levelB")==0) {
						levelB=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("minB")==0) {
						minB=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("maxB")==0) {
						maxB=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("levelS")==0) {
						levelS=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("minS")==0) {
						minS=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("maxS")==0) {
						maxS=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaA_B1")==0) {
						deltaA_B1=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaA_B2")==0) {
						deltaA_B2=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaA_B3")==0) {
						deltaA_B3=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaA_B4")==0) {
						deltaA_B4=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaB_B1")==0) {
						deltaB_B1=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaB_B2")==0) {
						deltaB_B2=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaB_B3")==0) {
						deltaS_B3=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaB_B4")==0) {
						deltaB_B4=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaS_B1")==0) {
						deltaS_B1=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaS_B2")==0) {
						deltaS_B2=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaS_B3")==0) {
						deltaS_B3=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("deltaS_B4")==0) {
						deltaS_B4=Integer.parseInt(line[1]);
					} else if (line[0].compareTo("maxSpeed")==0) {
						maxSpeed=Integer.parseInt(line[1]);
					} else {
						throw new Exception("Malformed file: Unknown token "+line[0]);
					}
				}
			}
			positionInWorld=0;
		}
		
		public boolean canContinue() {
			return rbWorld.isInTheWorld(positionInWorld) && levelS<maxS;
		}
		
		public void stepForward(int nbIt, double c1, double c2, double c3, double c4) {
			for (int i=0; i<nbIt; i++) {
				if (canContinue()) {
					applyBehavior1(c1);
					applyBehavior2(c2);
					applyBehavior3(c3);
					applyBehavior4(c4);
					positionInWorld+=actualSpeed;
				}
			}
		}
		// Stop C - Consume B - Produce A - Leave D get in
		void applyBehavior1(double percent) {
			
		  // Enable the shield => consume B, produce A and protect from C
		  if (levelB>0) {
		    levelB += deltaB_B1*percent;
		    levelA += deltaA_B1*percent;
		    if (levelA>maxA) {
		      levelA=maxA;
		    }
		    levelS += rbWorld.getRadC(positionInWorld)*(1.0-percent);
		  } else {
		    // Shield disabled => max damages
		    levelS += rbWorld.getRadC(positionInWorld);
		  }
		}

		// Stop D - Consume B - Produce A - Leave C get in
		void applyBehavior2(double percent) {
		  if (levelB>0) {
		    // Enable the shield => consume B, Produce A and protect from D
		    levelB += deltaB_B2*percent;
		    levelA += deltaA_B2*percent;
		    if (levelA>maxA) {
		      levelA=maxA;
		    }
		    levelS += rbWorld.getRadD(positionInWorld)*(1.0 - percent);
		  } else {
		    // Shield disabled => max damages
		    levelS += rbWorld.getRadD(positionInWorld);
		  }
		}

		// Consume A - Regenerate S
		void applyBehavior3(double percent) {
		  // Regenerates S only if A is avalaible
		  if (levelA>0) {
		    levelA += deltaA_B3*percent;
		    levelS += deltaS_B3*percent;
		    if (levelS < minS) {
		      levelS = minS;
		    }
		  }
		}

		// Run - Consume A - Regenerates B
		void applyBehavior4(double percent) {
		  // Regenerates B and moves if A is avalaible
		  if (levelA>0) {
		    levelA += deltaA_B4*percent;
		    levelB += deltaB_B4*percent;
		    if (levelB>maxB) {
		      levelB=maxB;
		    }
		    actualSpeed = maxSpeed*percent;
		  }	
		}
		
		double getLevelANormalized() {
			return levelA/maxA;
		}

		double getLevelBNormalized() {
			return levelB/maxB;
		}

		double getLevelSNormalized() {
			return levelS/maxS;
		}


	}

	public class RadbotWorld {
		public class WorldPosition implements Comparable<WorldPosition> {
			public double position;
			public double radC;
			public double radD;
			WorldPosition(double p, double c, double d) {
				position=p;
				radC=c;
				radD=d;
			}
			@Override
			public int compareTo(WorldPosition o) {
				// TODO Auto-generated method stub
				if (position<o.position) {
					return -1;
				} else if (position>o.position) {
					return 1;
				} else {
					return 0;
				}
			}
		}

		public ArrayList<WorldPosition> positions=new ArrayList<WorldPosition>();
		public double maxRadC=-Double.MAX_VALUE;
		public double maxRadD=-Double.MAX_VALUE;
		
		public RadbotWorld(String fileName) throws Exception {
			BufferedReader reader=new BufferedReader(new FileReader(fileName));
			String readLine;
			while ((readLine=reader.readLine())!=null) {
				String line[]=readLine.split("\t");
				if (line.length==3) {
					WorldPosition p=new WorldPosition(Double.parseDouble(line[0]),Double.parseDouble(line[1]),Double.parseDouble(line[2]));
					positions.add(p);
					maxRadC=Math.max(maxRadC,p.radC);
					maxRadD=Math.max(maxRadD,p.radD);
				} else {
					throw new Exception("Malformed file");
				}
			}
			Collections.sort(positions);
		}

		public double getRadC(double position) {
			int pP=0;
			while (pP<positions.size() && positions.get(pP).position<position) {
				pP++;
			}
			if (pP>=positions.size()) {
				return positions.get(pP-1).radC;
			} else if (pP<=0) {
				return positions.get(0).radC;
			} else {
				double radC1=positions.get(pP-1).radC;
				double radC2=positions.get(pP).radC;
				double ratio=(position-positions.get(pP-1).position)/(positions.get(pP).position-positions.get(pP-1).position);
				return radC1+ratio*(radC2-radC1);
			}
		}

		public double getRadD(double position) {
			int pP=0;
			while (pP<positions.size() && positions.get(pP).position<position) {
				pP++;
			}
			if (pP>=positions.size()) {
				return positions.get(pP-1).radD;
			} else if (pP<=0) {
				return positions.get(0).radD;
			} else {
				double radD1=positions.get(pP-1).radD;
				double radD2=positions.get(pP).radD;
				double ratio=(position-positions.get(pP-1).position)/(positions.get(pP).position-positions.get(pP-1).position);
				return radD1+ratio*(radD2-radD1);
			}
		}
		
		public double getNormalizedRadC(double position) {
			return getRadC(position)/maxRadC;
		}
		
		public double getNormalizedRadD(double position) {
			return getRadD(position)/maxRadD;
		}
		
		public boolean isInTheWorld(double position) {
			return position<=positions.get(positions.size()-1).position;
		}

	}
	

	@Override
	public double evaluate(GRNGenome aGenome) {
		double fit=0;
		try {
			GRNModel grn=buildGRNFromGenome(aGenome);
			Radbot rb=new Radbot(rbw, resourceDirectory+"/radbot.rb");
			grn.evolve(25);
			int nStep=0;
			while (rb.canContinue() && nStep<50000) {
				nStep++;
				grn.proteins.get(0).concentration=0.05*rb.getLevelANormalized();
				grn.proteins.get(1).concentration=0.05*rb.getLevelBNormalized();
				grn.proteins.get(2).concentration=0.05*rb.getLevelSNormalized();
				grn.proteins.get(3).concentration=0.05*rbw.getNormalizedRadC(rb.positionInWorld);
				grn.proteins.get(4).concentration=0.05*rbw.getNormalizedRadD(rb.positionInWorld);
				grn.evolve(1);
				double percents[]={grn.proteins.get(5).concentration,
						grn.proteins.get(6).concentration,
						grn.proteins.get(7).concentration,
						grn.proteins.get(8).concentration
				};
				double sumValues=percents[0]+percents[1]+percents[2]+percents[3];
				if (sumValues!=0) {
					percents[0]/=sumValues;
					percents[1]/=sumValues;
					percents[2]/=sumValues;
					percents[3]/=sumValues;
					double threshold=0.4;
					if (percents[0]>threshold) {
						percents[0]=1;
					} else {
						percents[0]/=threshold;
					}
					if (percents[1]>threshold) {
						percents[1]=1;
					} else {
						percents[1]/=threshold;
					}
//					System.err.println(rb.positionInWorld+"\t"+rbw.getNormalizedRadC(rb.positionInWorld)+"\t"+rbw.getNormalizedRadD(rb.positionInWorld));
				}
				rb.stepForward(1, percents[0], percents[1], percents[2], percents[3]);
//				System.err.println(rb+"\t"+rb.positionInWorld);
			}			
			fit=rb.positionInWorld;
		} catch (Exception e) {
			e.printStackTrace();
			fit=0;
		}
		aGenome.setNewFitness(fit);
		numEvaluations++;
		return fit;
	}
	
	public static void main(String args[]) {
		// generate a env file
		for (double i=0; i<100000.9; i++) {
//			System.out.println(i+"\t"+(5.0*Math.sin(i/20.0)+5.0)+"\t"+(-5.0*Math.sin(i/20.0)+5.0));
			double rnd=Math.random()*10;
			for (int j=0; j<100; j++) {
				System.out.println(i+"\t"+rnd+"\t"+(10.0-rnd));
				i++;
			}
			
		}
	}



}
