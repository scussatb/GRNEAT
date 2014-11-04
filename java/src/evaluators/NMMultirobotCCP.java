package evaluators;

import evolver.Evolver;
import evolver.GRNGenome;
import grn.GRNModel;

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
public class NMMultirobotCCP extends GRNGenomeEvaluator {
	public Random rng;
	
	public static int numEvaluations = 0;
	public double max_steps = 10000;
	public double initC = 0.025;
	public double C = 0.0001089;
	public double num_obstacles = 0;//3;
	
	public int accumulationType = 1; 
	
	public int hall_length = 5;// in number of rooms
	public int hall_width = 5;// in pixels
	public int room_width = 9;
	public int room_depth = 9;
	public int door_radius = 2;

	public int width = 0;
	public int height = 0;
	
	public int stepsPerAction = 25;
	public double epsilon = 0.1;
	public int numAgents = 3;
	public int numStates;
	public int sensorRadius = 9;
		
	public double globalAlpha = -1;
	public double globalGamma = -1;
	public double globalLambda = -1;		
	
	public int[] wall;
	public double[] reward;
	public double[] dReward;
	
	
	public NMMultirobotCCP( String args[], Random evolverRNG ) {
		name="NMMultirobotCCP";
		
		rng = evolverRNG;
		
		numGRNInputs = 4;
		numGRNOutputs = 4;
		long random_seed = 1;
		boolean parameterSampling = false;
		for (int k = 0; k < args.length; k++ ) {
			System.err.println(args[k]);
			if ( args[k].compareTo("stepsPerAction") == 0) {
				stepsPerAction = Integer.parseInt( args[k+1] );
			} else if ( args[k].compareTo("epsilon") == 0 ) {
				epsilon = Double.parseDouble( args[k+1] );
			} else if ( args[k].compareTo("numAgents")  == 0) {
				numAgents = Integer.parseInt( args[k+1] );
			} else if ( args[k].compareTo("sensorRadius") == 0 ) {
				sensorRadius = Integer.parseInt( args[k+1] );
			} else if ( args[k].compareTo("alpha") == 0 ) {
				globalAlpha = Double.parseDouble( args[k+1] );
			} else if ( args[k].compareTo("gamma") == 0 ) {
				globalGamma = Double.parseDouble( args[k+1] );
			} else if ( args[k].compareTo("lambda") == 0 ) {
				globalLambda = Double.parseDouble( args[k+1] );
			} else if ( args[k].compareTo("environmentType") == 0 ) {
				accumulationType = Integer.parseInt( args[k+1] );
			//} else if ( args[k].compareTo("randomSeed") == 0 ) {
			//	random_seed = Integer.parseInt( args[k+1] );
			} else if ( args[k].compareTo("parameterSampling") == 0 ) {
				//parameterSampling = ( 1 == Integer.parseInt( args[k+1] ) );
				parameterSampling = Boolean.parseBoolean( args[k+1] );
			}
		}
				
		width = ( room_width + 1 ) * hall_length + 1;
		height = hall_width + 2 * ( room_depth + 2 );
		//stepsPerAction, epsilon, num_iterations, swarm_size, environment_params, sensor_radius, alpha, gamma, lambda
		
		//rng = new java.util.Random(random_seed);
		
		System.out.println( "numAgents\t" + numAgents );
		System.out.println( "epsilon\t" + epsilon );
		System.out.println( "stepsPerAction\t" + stepsPerAction );
		System.out.println( "sensorRadius\t" + sensorRadius );
		System.out.println( "globalAlpha\t" + globalAlpha );
		System.out.println( "globalLambda\t" + globalLambda );
		System.out.println( "globalGamma\t" + globalGamma );
		System.out.println( "accumulationType\t" + accumulationType );
		System.out.println( "random_seed\t" + random_seed );
		System.out.println( "parameterSampling\t" + parameterSampling );
		
		// For fixed parameter sampling
		//if( !( globalAlpha > 0 && globalLambda > 0 && globalGamma > 0 ) ) {
		if( parameterSampling ) {
			System.err.println( "Fixed parameter sampling ");
			double fitness = 0;
			double best_fitness = 0;
			double[] best_params = new double[3];
			for( double alpha = 0.1; alpha < 1; alpha += 0.05 ) {
				for( double gamma = 0.1; gamma < 1; gamma += 0.05 ) {
					for( double lambda = 0.1; lambda < 1; lambda += 0.05 ) {
						fitness = evaluate( null );
						System.out.println( "" + alpha + "\t" + gamma + "\t" + lambda + "\t" + fitness );
						if( fitness > best_fitness ) {
							best_params[0] = alpha;
							best_params[1] = gamma;
							best_params[2] = lambda;
							best_fitness = fitness;
						}
					}
				}
			}
			System.out.println( "Best parameters:");
			System.out.println( "Fitness\t" + best_fitness );
			System.out.println( "Alpha\t" + best_params[0] );
			System.out.println( "Gamma\t" + best_params[1] );
			System.out.println( "Lambda\t" + best_params[2] );
			System.exit(0);
		}

	}
	
	// Linear index
	public int lidx(int x, int y) {
		if ( x < 0 ) x = -x;
		if ( y < 0 ) y = -y;
		return ( y + x * height );
	}

	public void generateHallway() {
		width = ( room_width + 1 ) * hall_length;
		height = hall_width + 2 * ( room_depth + 2 );

		wall = new int[ height * width ];

		for ( int r = 0; r < hall_length; r++ ) {
			int center = (int)( ( ( r + 0.5 ) * ( room_width + 1 )  ) );
			int x,y;	    
			
		    // Top row
			for ( x = 0; x <= room_width; x++ ) {
				y = 0;
				//System.out.println( "Ints " + x + " "  + y + " " + lidx(x,y) + " " + width + " " + height );
				wall[lidx( x + r * ( room_width + 1 ), y )] = 1;
				
				y = room_depth + 1;
				//wall[lidx(x,y)] = 1;
				wall[lidx( x + r * ( room_width + 1 ), y )] = 1;
			}
		    
			for ( y = 0; y < room_depth + 1; y++ ) { 
				x = r * ( room_width + 1 );				
				wall[lidx(x,y)] = 1;
				
				//x = ( r + 1 ) * ( room_width + 1 ) - 1;
				//System.out.println( "Ints " + x + " "  + y + " " + lidx(x,y) + " " + width + " " + height );
				//wall[lidx(x,y)] = 1;
			}
			    
		    // Doorway
			for ( x = ( center - door_radius + 1 ); x < ( center + door_radius ); x++ ) {
				y = ( 1 + room_depth );
				//System.out.println( "Ints " + x + " "  + y + " " + lidx(x,y) + " " + width + " " + height + " " + door_radius );
				wall[lidx(x,y)] = 0;
			}
/*		    				    
				    % Obstacles 
				    [y,x] = meshgrid( 2:( env.room_depth - 1 ), 2:( env.room_width - 1 ) );
				    valid_locations = [ y(:) x(:) ];
				    
				    idx = randperm( size(valid_locations,1) );
				    obstacles = valid_locations(idx(1:env.room_obstacles),:);
				    obs_idx = ( 1 + obstacles(:,1) ) + ( ( r - 1 ) * ( env.room_width + 1 ) + obstacles(:,2) ) * env.height;
				    env.wall(obs_idx) = 1;
*/	

			// Bottom row
			for ( x = 0; x <= room_width; x++ ) {
				y = height - 1;
				//System.out.println( "Ints " + x + " "  + y + " " + lidx(x,y) + " " + width + " " + height );
				wall[lidx( x + r * ( room_width + 1 ), y )] = 1;
				
				y = height - ( room_depth + 1 );
				//wall[lidx(x,y)] = 1;
				wall[lidx( x + r * ( room_width + 1 ), y )] = 1;
			}
		    
			for ( y = height - 1; y >= height - ( room_depth + 1 ); y-- ) { 
				x = r * ( room_width + 1 );				
				wall[lidx(x,y)] = 1;
				
				//x = ( r + 1 ) * ( room_width + 1 ) - 1;
				//System.out.println( "Ints " + x + " "  + y + " " + lidx(x,y) + " " + width + " " + height );
				//wall[lidx(x,y)] = 1;
			}
			    
		    // Doorway
			for ( x = ( center - door_radius + 1 ); x < ( center + door_radius ); x++ ) {
				y = height - ( room_depth + 1 );
				//System.out.println( "Ints " + x + " "  + y + " " + lidx(x,y) + " " + width + " " + height + " " + door_radius );
				wall[lidx(x,y)] = 0;
			}
			
/*				    
				    % Obstacles 
				    [y,x] = meshgrid( 2:( env.room_depth - 1 ), 2:( env.room_width - 1 ) );
				    valid_locations = [ y(:) x(:) ];
				    
				    idx = randperm( size(valid_locations,1) );
				    obstacles = valid_locations(idx(1:env.room_obstacles),:);
				    obs_idx = ( env.height - obstacles(:,1) ) + ( ( r - 1 ) * ( env.room_width + 1 ) + obstacles(:,2) ) * env.height;
				    env.wall(obs_idx) = 1;
*/				    
		}
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
	
	public void updateReward() {
		int idx;
		for ( int y = 0; y < height; y++ ) {
			for ( int x = 0; x < width; x++ ) {
				idx = lidx(x,y);
				reward[idx] += dReward[idx];
			}
		}
	}
	
	public void initializeDReward() {
		int idx;
		dReward = new double[height * width];
		int[] nbrWalls = new int[height*width];
		for ( int y = 0; y < height; y++ ) {
			for ( int x = 0; x < width; x++ ) {
				nbrWalls[lidx(x,y)] = 
					wall[lidx( ( x + 1 ) % width,y)] + wall[lidx( ( x - 1 ) % width,y)] +
					wall[lidx(x,( y + 1 ) % height)] + wall[lidx(x,( y + 1 ) % height)]; 
			}
		}
		for ( int y = 0; y < height; y++ ) {
			for ( int x = 0; x < width; x++ ) {				
				idx = lidx(x,y);
				if( wall[idx] == 0) {
					if ( accumulationType == 1 ) dReward[idx] = C;
					if ( accumulationType == 2 ) dReward[idx] = nbrWalls[idx] * C;
				}
			}
		}
		for ( int y = 0; y < height; y++ ) {
			for ( int x = 0; x < width; x++ ) {
				idx = lidx(x,y);
				if ( wall[idx] == 0 ) reward[idx] = initC;				
			}
		}
	}
	
	public class Agent {
		public int uid;
		public int numActions;
		public int numStates;
		public double epsilon;
		public int sensor_radius;
		
		// Learned variables
		// Q-Function being learned by SARSA
		public double[] q;
		// Eligibility trace
		public double[] e;

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
		
		public double alpha;
		public double gamma;
		public double lambda;
		
		public Agent(int uid, double epsilon, int sensor_radius, int startX, int startY, GRNModel grn) {
			this.uid = uid;
			this.numActions = 4;
			//agent.numStates = 4 ^ ( ( 4 * sensor_radius + 1 ) ); 
			this.numStates = 13824 + 1;
			this.epsilon = epsilon;
			this.sensor_radius = sensor_radius;
			
			// Learned variables
			// Q-Function being learned by SARSA
			this.q = new double[numStates * numActions];
			// Eligibility trace
			this.e = new double[numStates * numActions];

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
	
	ArrayList<Agent> agents;
	
	public class ArrayIndexComparator implements Comparator<Integer>
	{
	    private final double[] array;

	    public ArrayIndexComparator(double[] array)
	    {
	        this.array = array;
	    }

	    public Integer[] createIndexArray()
	    {
	        Integer[] indexes = new Integer[array.length];
	        for (int i = 0; i < array.length; i++)
	        {
	            indexes[i] = i; // Autoboxing
	        }
	        return indexes;
	    }

	    @Override
	    public int compare(Integer index1, Integer index2)
	    {
	         // Autounbox from Integer to int to use as array indexes
	    	if ( array[index1] > (array[index2]) ) return 1; 
	    	else if ( array[index1] < (array[index2]) ) return -1;
	    	else return 0;
	    }
	}
	
	public int perm4toidx( String inp ) {
		if ( inp.compareTo("3210") == 0 ) return 0;
		if ( inp.compareTo("3201") == 0 ) return 1;
		if ( inp.compareTo("3120") == 0 ) return 2;
		if ( inp.compareTo("3102") == 0 ) return 3;
		if ( inp.compareTo("3012") == 0 ) return 4;
		if ( inp.compareTo("3021") == 0 ) return 5;
		if ( inp.compareTo("2310") == 0 ) return 6;
		if ( inp.compareTo("2301") == 0 ) return 7;
		if ( inp.compareTo("2130") == 0 ) return 8;
		if ( inp.compareTo("2103") == 0 ) return 9;
		if ( inp.compareTo("2013") == 0 ) return 10;
		if ( inp.compareTo("2031") == 0 ) return 11;
		if ( inp.compareTo("1230") == 0 ) return 12;
		if ( inp.compareTo("1203") == 0 ) return 13;
		if ( inp.compareTo("1320") == 0 ) return 14;
		if ( inp.compareTo("1302") == 0 ) return 15;
		if ( inp.compareTo("1032") == 0 ) return 16;
		if ( inp.compareTo("1023") == 0 ) return 17;
		if ( inp.compareTo("0213") == 0 ) return 18;
		if ( inp.compareTo("0231") == 0 ) return 19;
		if ( inp.compareTo("0123") == 0 ) return 20;
		if ( inp.compareTo("0132") == 0 ) return 21;
		if ( inp.compareTo("0312") == 0 ) return 22;
		if ( inp.compareTo("0321") == 0 ) return 23;
		return -1;
	}
	
	/* Update the world and swarm after agent at index sidx is updated.
	 * (NOTE) alpha, gamma, and lambda are OPTIONAL, and they will override the
	 * GRN when they are used */	
	public void updateAgent( int idx ) {
		// For typing simplicity
		Agent ag = agents.get( idx );

		int saidx = ag.state_idx * ag.numActions + ag.action;// state-action index
		// Update eligibility trace
		//System.out.println( "update_agent: " + ag.state_idx + " " + ag.numActions + " " + ag.numStates);
		ag.e[saidx] = ag.e[saidx] + 1;

		// State transition for both world and agent
		ag.current_reward = 0;
		
		int y = ag.curY;
		int x = ag.curX;
		
		switch( ag.action ) {
		    case 1:
		        x = x + 1;
		        break;
		    case 2:
		        x = x - 1;
		        break;
		    case 3:
		        y = y + 1;
		        break;
		    case 4:
		        y = y - 1;
		        break;
		}
		
		int ny = y % height;
		if ( ny < 0 ) ny = -ny;
		int nx = x % width;
		if ( nx < 0 ) nx = -nx;
		
		// If the space is free (no obstruction) then move there
		int[] swarmMask = new int[width*height];
		for( int k = 0; k < numAgents; k++ ) {
			if ( k != idx ) {
				swarmMask[lidx(agents.get( k ).curX,agents.get( k ).curY)] = 1;
			}
		}
		
		int nidx = lidx(nx,ny);
		//System.out.println( "Movement " + ( wall[nidx] == 0 ) + " " + ( swarmMask[nidx] == 0 ) );
		if ( ( wall[nidx] == 0 ) && ( swarmMask[nidx] == 0 ) ) {
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

        double[] features = new double[4 * 3];

        int oidx;
        int nUp = 0;
        for (int ox = x - r; ox < x + r; ox++) {
        	for (int oy = y + 1; oy < y + r; oy++) {
        		oidx = lidx(ox % width,oy % height);
        		features[0] += reward[oidx];
        		features[4] += wall[oidx];
        		features[8] += swarmMask[oidx];
        		nUp++;
        	}
        }
        int nNbrs = nUp;
        features[0] /= nNbrs;
        features[4] /= nNbrs;
        features[8] /= nNbrs;
        				
        for (int ox = x - r; ox < x + r; ox++) {
        	for (int oy = y - r; oy < y - 1; oy++) {
        		oidx = lidx(ox % width,oy % height);
        		features[1] += reward[oidx];
        		features[5] += wall[oidx];
        		features[9] += swarmMask[oidx];
        	}
        }
        features[1] /= nNbrs;
        features[5] /= nNbrs;
        features[9] /= nNbrs;
        
        for (int ox = x + 1; ox < x + r; ox++) {
        	for (int oy = y - r; oy < y + r; oy++) {
        		oidx = lidx(ox % width,oy % height);
        		features[2] += reward[oidx];
        		features[6] += wall[oidx];
        		features[10] += swarmMask[oidx];
        	}
        }
        features[2] /= nNbrs;
        features[6] /= nNbrs;
        features[10] /= nNbrs;

        for (int ox = x - r; ox < x - 1; ox++) {
        	for (int oy = y - r; oy < y + r; oy++) {
        		oidx = lidx(ox % width,oy % height);
        		features[3] += reward[oidx];
        		features[7] += wall[oidx];
        		features[11] += swarmMask[oidx];
        	}
        }
        features[3] /= nNbrs;
        features[7] /= nNbrs;
        features[11] /= nNbrs;
        
        
        // Epsilon greedy action selection        
        ArrayList<Integer> rord = new ArrayList<Integer>();
        double[] tempfeats = new double[4];
		int[] feature_codes = { 1, 2, 3 };
		String thisperm;
		
		int state_idx = 0;
		for ( int fidx = 0; fidx < 3; fidx++ ) {
			rord.clear(); rord.add(0); rord.add(1); rord.add(2); rord.add(3); 
			java.util.Collections.shuffle( rord );
			for ( int k = 0; k < tempfeats.length; k++ ) {
				tempfeats[k] = features[fidx * 4 + rord.get(k)];
			}
			ArrayIndexComparator comp = new ArrayIndexComparator( tempfeats );
			Integer[] sridx =  comp.createIndexArray();
			Arrays.sort( sridx, comp );
			
			thisperm = "";
			for ( int k = 0; k < tempfeats.length; k++ ) {
				int mytrueidx;
				for ( int i = 0; i < 4; i++ ) {
					if ( sridx[k] == rord.get(i) ) thisperm += i;
				}				
			}

			feature_codes[fidx] = perm4toidx( thisperm );
			state_idx += Math.pow( 24, fidx ) * feature_codes[fidx];
			
			//System.out.println( "Features: " + fidx + " " + feature_codes[fidx] + " " + ag.q.length);
		}
		
		// Exploration
		if ( rng.nextDouble() < ag.epsilon ) { 
		    ag.nextAction = rng.nextInt( ag.numActions );
		} else {// exploitation
			rord.clear(); rord.add(0); rord.add(1); rord.add(2); rord.add(3); 
			java.util.Collections.shuffle( rord );
			double[] tmpQ = new double[4];
			for ( int k = 0; k < tmpQ.length; k++ ) {
				tmpQ[k] = ag.q[state_idx * 4 + rord.get(k)];
			}
			ArrayIndexComparator comp = new ArrayIndexComparator( tmpQ );
			Integer[] sridx =  comp.createIndexArray();
			Arrays.sort( sridx, comp );
			int[] truesort = new int[4];  
			
			thisperm = "";
			for ( int k = 0; k < tempfeats.length; k++ ) {
				int mytrueidx;
				for ( int i = 0; i < 4; i++ ) {
					if ( sridx[k] == rord.get(i) ) truesort[k] = i;
				}				
			}
			
		    ag.nextAction = truesort[0];
		}

		// We only use neuromodulation for learning params, if they are specified,
		// no neuromodulation.
		ag.alpha = globalAlpha; ag.lambda = globalLambda; ag.gamma = globalGamma;
		if( !( globalAlpha < 0 && globalLambda < 0 && globalGamma < 0 ) ) {
		    // Neuromodulation starts here
		    // senses = oarSensors( updatedWorld, swarm, ag );
		    
		    double dirt_input = 0;//sum( senses(1,:) ) / 4;
		    double obstacle_input = 0;//sum( senses(2,:) ) / 4;
		    double swarm_input = 0;//sum( senses(3,:) ) / 4;
		    
		    for( int k = 0; k < 4; k++ ) {
		    	dirt_input += features[k*4];
		    	obstacle_input += features[k*4 + 1];
		    	swarm_input += features[k*4 + 2];
		    }
		    dirt_input /= 4;
		    obstacle_input /= 4;
		    swarm_input /= 4;
		    
		    ag.grn.proteins.get(0).concentration = 0.25 * dirt_input;
		    ag.grn.proteins.get(1).concentration = 0.25 * obstacle_input;
		    ag.grn.proteins.get(2).concentration = 0.25 * swarm_input;
		    ag.grn.proteins.get(3).concentration = 0.25 * Math.max( 1, ag.current_reward );
		    
		    ag.grn.evolve(stepsPerAction);

		    double norm_factor = ag.grn.proteins.get(4).concentration;

		    // The 4th protein is used to normalize
		    if ( norm_factor == 0 ) { 
		    	ag.alpha = ag.grn.proteins.get(5).concentration;// / norm_factor;
		    	ag.gamma = ag.grn.proteins.get(6).concentration;// / norm_factor;
		    	ag.lambda = ag.grn.proteins.get(7).concentration;// / norm_factor;
		    } else {
		    	ag.alpha = Math.max( 1, ag.grn.proteins.get(5).concentration / norm_factor );
		    	ag.gamma = Math.max( 1, ag.grn.proteins.get(6).concentration / norm_factor );
		    	ag.lambda = Math.max( 1, ag.grn.proteins.get(7).concentration / norm_factor );
			}
		    
		} // Neuromodulation is done

		//    Reinforcement is 0. Prediction of failure given by v weight.
		double reward = ag.current_reward;

		// We need to compute the sarsa-delta
		double delta = reward + ag.gamma * ag.q[ag.next_state_idx * ag.numActions + ag.nextAction] - 
		    ag.q[ag.state_idx * ag.numActions + ag.action];

		for (int k = 0; k < ag.q.length; k++) {
			ag.q[k] += ( delta * ag.alpha  * ag.e[k] );
			ag.e[k] = ag.e[k] * ag.gamma * ag.lambda;
		}

		ag.action = ag.nextAction;
		ag.state_idx = state_idx;

		// Replace the old agent in the swarm
		agents.set( idx, ag );
		//System.out.println( "Updated: " + ag.curX + "," + ag.curY );
	}
	
	@Override
	public double evaluate(GRNGenome aGenome) {

		double fitness=0.0;
		GRNModel grn;
		if ( aGenome != null )
			grn = buildGRNFromGenome(aGenome);
		else
			grn = null;
		
		// Initialize map properties
		generateHallway();		
		initializeReward();
		initializeDReward();
		
		//writeMap2File( reward, "initialreward.csv" );
		//writeMap2File( dReward, "initialdReward.csv" );
		
		// Initialize RL
		int[] startPos = new int[numAgents * 2];
		for ( int k = 0; k < numAgents; k++ ) {			
			startPos[2*k] = width / numAgents * k;
			startPos[2*k+1] = height / 2;
		}
		
		agents = new ArrayList<Agent>();
		for ( int k = 0; k < numAgents; k++ ) {
			agents.add(new Agent(k, epsilon, sensorRadius, startPos[2*k], startPos[2*k+1], grn));
		}
		
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
	        updateReward();
	        
	        // We update agents sequentially, thus there will be regular bias if
	        //   we do not randomize the update order.
	        
	        ArrayList<Integer> swarm_ord = new ArrayList<Integer>();
	        for ( int k = 0; k < numAgents; k++ ) swarm_ord.add(k);
	        java.util.Collections.shuffle( swarm_ord );

	        int sidx;
	        for ( int k = 0; k < numAgents; k ++ ) {	        		        	
	            sidx = swarm_ord.get(k);// Agent index
	            
	            //System.out.println( "Agent " + k + " " + agents.get(sidx).curX + "," + agents.get(sidx).curY + " " + agents.get(sidx).current_reward );
	                        
	            /*alpha = 0.1; gamma = 0.8; lambda = 0.1;
	            if neuromodulation_flag
	                [world, swarm] = update_agent( world, swarm, sidx );
	            else
	                [world, swarm] = update_agent( world, swarm, sidx, alpha, gamma, lambda );
	            end*/
	            
	            updateAgent( sidx );
	            fitness += agents.get(sidx).current_reward;
	        }			

		}
		
		//System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		if ( grn != null ) {
			aGenome.setNewFitness(fitness);
		}

		numEvaluations++;
		return fitness;
	}	

}
