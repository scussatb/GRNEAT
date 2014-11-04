package evaluators;

import evolver.Evolver;
import evolver.GRNGenome;
import grn.GRNModel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

import operators.GRNAddGeneMutationOperator;
import operators.GRNAligningCrossoverOperator;
import operators.GRNDeleteGeneMutationOperator;
import operators.GRNGeneMutationOperator;

// Pardon the lack of style
public class CoverageControl extends GRNGenomeEvaluator {
	public Random rng;

	public double max_steps = 1000;
	public double initC = 1;//0.025;
	//	public int mapID = 1;
	public boolean mapIDs[];
	public static int maxMaps=100;

	public int width = 0;
	public int height = 0;

	public int stepsPerAction = 25;
	public double epsilon = 0.1;
	public int sensorRadius = 3;

	public int numUnobstructed = 0;
	public int[] wall;
	public double[] reward;
	public double[] dReward;

	//public String mapDirectory = "src/resources/CoverageControl";
	public String mapDirectory = "resources/CoverageControl";
	public int verbose=0; // 0 is silent

	public CoverageControl( String args[], Random evolverRNG ) {
		name="CoverageControl";

		rng = evolverRNG;

		numGRNInputs = 9;
		numGRNOutputs = 4;
		long random_seed = 1;
		boolean parameterSampling = false;
		//		mapID = rng.nextInt( 100 );// there are 100 maps
		mapIDs=new boolean[maxMaps];
		for (int k = 0; k < args.length; k++ ) {
			//			System.err.println(args[k]);
			if ( args[k].compareTo("stepsPerAction") == 0) {
				stepsPerAction = Integer.parseInt( args[k+1] );
			} else if ( args[k].compareTo("mapID") == 0 ) {
				mapIDs[Integer.parseInt( args[k+1] )-1]=true;
			} else if ( args[k].compareTo("verbose") == 0 ) {
				verbose = Integer.parseInt(args[k+1]);
			}
		}

		if (verbose>0) {
			System.out.print( "mapID\t"); 
			for (int i=0; i<maxMaps; i++) {
				if (mapIDs[i]) {
					System.out.print((i+1)+" ");
				}
			}
			System.out.println();

			System.out.println( "epsilon\t" + epsilon );
			System.out.println( "stepsPerAction\t" + stepsPerAction );
			System.out.println( "sensorRadius\t" + sensorRadius );

			System.out.println( "random_seed\t" + random_seed );
		}


	}

	public void loadMap(int mapID) {
		BufferedReader br = null;
		try {
			br = new BufferedReader( new FileReader( mapDirectory + "/greatmap_coverage_" + mapID + ".txt" ) );

			String line, allLines = "";					

			height = 0; 
			width = 0;

			while ((line= br.readLine()) != null) {
				if (verbose>0) System.out.println( line );
				allLines += line;
				width = line.length();
				height++;
			}

			wall = new int[ allLines.length() ];
			for( int k = 0; k < allLines.length(); k++ ) {
				if(  allLines.charAt( k ) == '0' ) {
					wall[k] = 0;
					numUnobstructed++;
				} else {
					wall[k] = 1;
				}
			}
			if (verbose>=2) {
				// printing wall coordinates
				System.out.println("===========\nWall coordinates:");
				for (int x=0; x<width; x++) {
					for (int y=0; y<height; y++) {
						if (wall[lidx(x,y)]==1) {
							System.out.println(x+"\t"+y);
						}
					}
				}
				System.out.println("===========");
			}
		} catch( IOException e ) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}

	// Linear index
	public int lidx(int x, int y) {
		if ( x < 0 ) x = -x;
		if ( y < 0 ) y = -y;
		return ( y + x * height );
	}	

	public void printBinaryMap(int[] map) {
		for ( int y = 0; y < height; y++ ) {
			for ( int x = 0; x < width; x++ ) {
				System.out.print( map[lidx(x,y)] );
			}
			System.out.print( "\n" );
		}
	}

	public void writeMap2File(double[] map, String fn) {
		try {
			FileWriter w = new FileWriter(fn);

			for ( int y = 0; y < height; y++ ) {
				for ( int x = 0; x < width; x++ ) {
					if ( x > 0 ) w.append( "," );
					w.append( "" + map[lidx(x,y)] );
				}
				w.append( "\n" );
			}
			w.flush();
			w.close();

		} catch( IOException e ) {
			e.printStackTrace();
		}
	}

	public void initializeReward() {
		int idx;
		reward = new double[height * width];
		for ( int y = 0; y < height; y++ ) {
			for ( int x = 0; x < width; x++ ) {
				idx = lidx(x,y);
				if ( wall[idx] == 0 ) reward[idx] = initC;				
			}
		}
	}


	public class Agent {
		public int uid;
		public double epsilon;
		public int sensor_radius;

		// Agent state variables
		public double current_reward;
		public int curY;
		public int curX;
		public int state_idx;
		public int next_state_idx;
		public int prev_action;
		public int action;
		public int nextAction;
		public GRNModel grn;

		public Agent(int uid, double epsilon, int sensor_radius, int startX, int startY, GRNModel grn) {
			this.uid = uid;
			//agent.numStates = 4 ^ ( ( 4 * sensor_radius + 1 ) ); 
			this.epsilon = epsilon;
			this.sensor_radius = sensor_radius;

			// Agent state variables
			this.current_reward = 0;
			this.curY = startY;
			this.curX = startX;
			this.state_idx = 0;
			this.prev_action = 0;
			this.action = rng.nextInt(4);
			this.grn = grn;
		}


	}


	/* Update the world and swarm after agent at index sidx is updated.
	 * (NOTE) alpha, gamma, and lambda are OPTIONAL, and they will override the
	 * GRN when they are used */	
	public Agent updateAgent( Agent ag ) {
		// For typing simplicity
		//Agent ag = agents.get( idx );

		// State transition for both world and agent
		ag.current_reward = 0;

		int y = ag.curY;
		int x = ag.curX;

		switch( ag.action ) {
		case 0:
			x = x + 1;
			break;
		case 1:
			x = x - 1;
			break;
		case 2:
			y = y + 1;
			break;
		case 3:
			y = y - 1;
			break;
		}

		int ny = y % height;
		if ( ny < 0 ) ny = -ny;
		int nx = x % width;
		if ( nx < 0 ) nx = -nx;

		int nidx = lidx(nx,ny);
		//System.out.println( "Movement " + ( wall[nidx] == 0 ) + " " + ( swarmMask[nidx] == 0 ) );
		if ( ( wall[nidx] == 0 ) ) {
			ag.curY = ny;
			ag.curX = nx;
		}

		ag.current_reward = reward[nidx];
		reward[nidx] = 0;

		//System.out.println( "Rewarding " + ag.current_reward + " " + reward[nidx] + " " + reward[rng.nextInt(1000)]);

		// Sensing/feature detection
		y = ag.curY;
		x = ag.curX;
		int r = sensorRadius;

		double[] features = new double[4 * 2];

		int oidx;
		int nUp = 0;
		for (int ox = x - r; ox < x + r; ox++) {
			for (int oy = y + 1; oy < y + r; oy++) {
				oidx = lidx(ox % width,oy % height);
				features[0] += reward[oidx];
				features[4] += wall[oidx];
				nUp++;
			}
		}
		int nNbrs = nUp;
		features[0] /= nNbrs;
		features[4] /= nNbrs;

		for (int ox = x - r; ox < x + r; ox++) {
			for (int oy = y - r; oy < y - 1; oy++) {
				oidx = lidx(ox % width,oy % height);
				features[1] += reward[oidx];
				features[5] += wall[oidx];
			}
		}
		features[1] /= nNbrs;
		features[5] /= nNbrs;

		for (int ox = x + 1; ox < x + r; ox++) {
			for (int oy = y - r; oy < y + r; oy++) {
				oidx = lidx(ox % width,oy % height);
				features[2] += reward[oidx];
				features[6] += wall[oidx];
			}
		}
		features[2] /= nNbrs;
		features[6] /= nNbrs;

		for (int ox = x - r; ox < x - 1; ox++) {
			for (int oy = y - r; oy < y + r; oy++) {
				oidx = lidx(ox % width,oy % height);
				features[3] += reward[oidx];
				features[7] += wall[oidx];
			}
		}
		features[3] /= nNbrs;
		features[7] /= nNbrs;


		for( int fid = 0; fid < features.length; fid++ ) {
			ag.grn.proteins.get(fid).concentration = features[fid];
		}
		ag.grn.proteins.get(features.length).concentration = 0.25 * Math.max( 1, ag.current_reward );

		ag.grn.evolve(stepsPerAction);

		double maxAction = 0;
		int mAid = 0;
		for( int aid = 0; aid < 4; aid++ ) {
			if( maxAction < ag.grn.proteins.get(features.length + 1 + aid).concentration ) {
				mAid = aid;
				maxAction = ag.grn.proteins.get(features.length + 1 + aid).concentration;
			}	    	
		}

		//    Reinforcement is 0. Prediction of failure given by v weight.
		double reward = ag.current_reward;

		ag.action = mAid;
		//ag.state_idx = state_idx;

		return ag;

		// Replace the old agent in the swarm
		//agents.set( idx, ag );
		//System.out.println( "Updated: " + ag.curX + "," + ag.curY );
	}

	public double evaluate(GRNModel grn) {

		double fitness=0.0;
		int nbMaps=0;

		for (int mapID=1; mapID<=maxMaps; mapID++) {
			if (mapIDs[mapID-1]) {
				//				System.out.println("map "+mapID);
				numUnobstructed=0;
				loadMap(mapID);
				double localFit=0.0;
				// Initialize map properties
				//generateHallway();
				//loadMap();
				initializeReward();
				//initializeDReward();

				//writeMap2File( reward, "initialreward.csv" );
				//writeMap2File( dReward, "initialdReward.csv" );		

				Agent ag = new Agent( 0, epsilon, sensorRadius, 0, 0, grn );

				// Warm up the GRN
				if ( grn != null ) {
					grn.reset();
					grn.evolve(25);
				}

				// Run the fitness test
				for (int nStep=0; nStep < max_steps; nStep++) {
					//System.out.println( "Time: " + nStep + " " + fitness );
					//writeMap2File( reward, "reward" + nStep + ".csv" );

					// Environment dynamics
					//updateReward();

					// We update agents sequentially, thus there will be regular bias if
					//   we do not randomize the update order.

					if (verbose>=2) System.out.println(ag.curX+"\t"+ag.curY);
					ag = updateAgent( ag );


					localFit += ag.current_reward;


				}

				localFit /= (double) (numUnobstructed);
				if (verbose>0) System.out.println("Fitness "+localFit);
				fitness+=localFit;
				nbMaps++;
			}
		}
		fitness/=nbMaps;
		return fitness;

	}

	@Override
	public double evaluate(GRNGenome aGenome) {
		GRNModel grn;
		if (aGenome!=null) {
			grn=buildGRNFromGenome(aGenome);
		} else {
			grn=null;
		}
		double fitness=0.0;

		//System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		if ( grn != null ) {
			fitness=evaluate(grn);
		}
		aGenome.setNewFitness(fitness);

		numEvaluations++;
		return fitness;
	}	

	public static void main(String args[]) throws Exception {
		/*
		String greatFiles[]={
				"grn_1.grn",	"grn_14.grn",	"grn_19.grn",	"grn_23.grn",	"grn_4.grn",	"grn_9.grn",
				"grn_10.grn",	"grn_15.grn",	"grn_2.grn",	"grn_24.grn",	"grn_5.grn",
				"grn_11.grn",	"grn_16.grn",	"grn_20.grn",	"grn_25.grn",	"grn_6.grn",
				"grn_12.grn",	"grn_17.grn",	"grn_21.grn",	"grn_7.grn",
				"grn_13.grn",	"grn_18.grn",	"grn_22.grn",	"grn_3.grn",	"grn_8.grn"
		};
		String gaFiles[]={
				"grn_1.grn",	"grn_14.grn",	"grn_19.grn",	"grn_23.grn",	"grn_4.grn",	"grn_9.grn",
				"grn_10.grn",	"grn_15.grn",	"grn_2.grn",	"grn_24.grn",	"grn_5.grn",
				"grn_11.grn",	"grn_16.grn",	"grn_20.grn",	"grn_25.grn",	"grn_6.grn",
				"grn_12.grn",	"grn_17.grn",	"grn_21.grn",	"grn_7.grn",
				"grn_13.grn",	"grn_18.grn",	"grn_22.grn",	"grn_3.grn",	"grn_8.grn"
		};
		String esFiles[]={
				"grn_1.grn",	"grn_14.grn",	"grn_19.grn",	"grn_23.grn",	"grn_4.grn",	"grn_9.grn",
				"grn_10.grn",	"grn_15.grn",	"grn_2.grn",	"grn_24.grn",	"grn_5.grn",
				"grn_11.grn",	"grn_16.grn",	"grn_20.grn",	"grn_6.grn",
				"grn_12.grn",	"grn_17.grn",	"grn_21.grn",	"grn_7.grn",
				"grn_13.grn",	"grn_18.grn",	"grn_22.grn",	"grn_3.grn",	"grn_8.grn"
		};
		GRNModel greatGRN[]=new GRNModel[greatFiles.length];
		GRNModel gaGRN[]=new GRNModel[gaFiles.length];
		GRNModel esGRN[]=new GRNModel[esFiles.length];
		for (int i=0; i<esFiles.length; i++) {
			greatGRN[i]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher017-hyperion007/launcher019/Generalization/CC/GREAT/"+greatFiles[i]);
			gaGRN[i]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher017-hyperion007/launcher019/Generalization/CC/GA/"+gaFiles[i]);
			esGRN[i]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher017-hyperion007/launcher019/Generalization/CC/ES/"+esFiles[i]);
		}
		greatGRN[24]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher017-hyperion007/launcher019/Generalization/CC/GREAT/"+greatFiles[24]);
		gaGRN[24]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher017-hyperion007/launcher019/Generalization/CC/GA/"+gaFiles[24]);
		
		System.out.println("GREAT");
		for (int j=0; j<greatGRN.length; j++) {
			
			System.out.print("\t"+greatFiles[j]);
		}
		System.out.println();
		for (int i=1; i<=100; i++) {
			System.out.print("map"+i);
			String ccArgs[]={"mapID", String.valueOf(i), "verbose", "0"};
			CoverageControl cc=new CoverageControl(ccArgs, new Random());
			for (int j=0; j<greatGRN.length; j++) {
				double greatFit=cc.evaluate(greatGRN[j]);
				System.out.print("\t"+greatFit);
			}
			System.out.println();
		}

		System.out.println("\n\nGA");
		for (int j=0; j<gaGRN.length; j++) {
			
			System.out.print("\t"+gaFiles[j]);
		}
		System.out.println();
		for (int i=1; i<=100; i++) {
			System.out.print("map"+i);
			String ccArgs[]={"mapID", String.valueOf(i)};
			CoverageControl cc=new CoverageControl(ccArgs, new Random());
			for (int j=0; j<gaGRN.length; j++) {
				double gaFit=cc.evaluate(gaGRN[j]);
				System.out.print("\t"+gaFit);
			}
			System.out.println();
		}

		System.out.println("\n\nES");
		for (int j=0; j<esGRN.length; j++) {
			
			System.out.print("\t"+esFiles[j]);
		}
		System.out.println();
		for (int i=1; i<=100; i++) {
			System.out.print("map"+i);
			String ccArgs[]={"mapID", String.valueOf(i)};
			CoverageControl cc=new CoverageControl(ccArgs, new Random());
			for (int j=0; j<esGRN.length; j++) {
				double esFit=cc.evaluate(esGRN[j]);
				System.out.print("\t"+esFit);
			}
			System.out.println();
		}
		System.exit(0);*/
/**/		
//		System.err.println(args[0]);
//		GRNModel grn=GRNModel.loadFromFile(args[0]);
		GRNModel grn=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher017-hyperion007/Generalization/CC/GA/grn_7.grn");
		//		grn.reset();
		double avg=0;
		/*		String covArgs[]=new String[12];
		for (int i=0; i<5; i++) {
			covArgs[i*2]="mapID";
			covArgs[i*2+1]=String.valueOf(i+1);
		}
		covArgs[10]="verbose";
		covArgs[11]="1";*/
		String covArgs[]={"mapID", "30"/*args[1]*/, "verbose", "0"};
		CoverageControl cov=new CoverageControl(covArgs, new Random());	

		double fit;
		fit=cov.evaluate(grn);
		System.err.println("Final fitness: "+fit);
		avg+=fit;
		System.err.println("=================\n"+(avg*100.0));

	}

}
