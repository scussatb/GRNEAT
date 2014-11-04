package evaluators;

import neural.FeedforwardNeuralNetwork;
import neural.NeuralPosition;
import evolver.GRNGenome;
import grn.GRNModel;

/*********
 * The retina problem. See:
 *
 * N. Kashtan and U. Alon, ÒSpontaneous evolution of modularity and network motifs.Ó Proceedings of the National Academy of Sciences, 102: 13773-13779, 2005.
 *
 * Jeff Clune, Benjamin E. Beckmann, Philip K. McKinley, and Charles Ofria. Investigating whether hyperneat produces modular neural networks. In Proc. of the Genetic and Ev. Comp. Conf., New York, NY, USA, 2010. ACM Press.
 */

public class RetinaEvaluator extends GRNGenomeEvaluator {
	public String name = "RetinaEvaluator";
	
	boolean[][] inputPatterns;
	boolean[][] targetOutputPatterns;
	int validPatternCount;
	int trialCount;
	
	double inputRangeScaling = 0.1;
	int switchRate = 20;// switch task every 20 generations
	int lastSwitch = 0;
	boolean and_mode = true;
	
	
	public int numNeuralInputs;
	public int numNeuralOutputs;
	public int numNeuralHidden1;
	public int numNeuralHidden2;
	public int numNeuralHidden3;
	
	public RetinaEvaluator() {
		nonCacheable = true;
				
		int numRetinas = 2;
		
		// All remaining parameters are autodetected
		int retinaWidth = validPatternsLeft[0].length;
		
		numGRNInputs = 4;
		numGRNOutputs = 3;
		
		numNeuralInputs = numRetinas * retinaWidth;
		numNeuralOutputs = 1;// The last one is for the true/false threshold			
		numNeuralHidden1 = 4;
		numNeuralHidden2 = 2;
		
		validPatternCount = validPatternsLeft.length;
		
		trialCount = 1 << ( numGRNInputs );
		
		inputPatterns = new boolean[trialCount][numNeuralInputs];
		targetOutputPatterns = new boolean[trialCount][numRetinas];			
		
		int[] validPatternsLeftVal = new int[validPatternCount];
		int[] validPatternsRightVal = new int[validPatternCount];
		for (int i = 0; i < validPatternCount; i++) {
			int valLeft = 0;
			int valRight = 0;
			for (int x = 0; x < retinaWidth; x++) {
				valLeft <<= 1;
				valRight <<= 1;
				valLeft |= (int) validPatternsLeft[i][x];
				valRight |= (int) validPatternsRight[i][x];
			}
			validPatternsLeftVal[i] = valLeft;
			validPatternsRightVal[i] = valRight;
		}		

		for (int t = 0; t < trialCount; t++) {
			int p = t;
			
			//System.out.print( t + "\t" );
			
			for (int x = 0; x < numGRNInputs; x++, p >>= 1) {
				inputPatterns[t][x] = ( (p & 0x1) == 1 );//? inputMaxValue : inputMinValue;
				//System.out.print( inputPatterns[t][x] + " " );
			}
			
			//System.out.print( "\t" );
			
			// Fill target outputs with default value.
			for (int x = 0; x < numRetinas; x++) {
				targetOutputPatterns[t][x] = false;
			}
			
			// If this trial corresponds to a valid pattern on the left or right then set the outputs accordingly.
			int prVal = (t >> 4) & 0x0F;
			int plVal = t & 0x0F;
			for (int i = 0; i < validPatternCount; i++) {
				targetOutputPatterns[t][0] |= (validPatternsLeftVal[i] == plVal);
				targetOutputPatterns[t][1] |= (validPatternsRightVal[i] == prVal);
				//System.out.print( " " + targetOutputPatterns[t][0] + "," + targetOutputPatterns[t][0] );
			}
			//System.out.print("\n");
			
			
		}
	}
	
	private static double[][] validPatternsLeft = new double[][]
			{ { 0, 1, 1, 1 }, // 7
			  { 1, 0, 1, 1 }, // 11
			  { 0, 0, 1, 1 }, // 3
			  { 0, 0, 0, 0 }, // 0
			  { 0, 1, 0, 0 }, // 4
			  { 1, 0, 0, 0 }, // 8
			  { 0, 0, 0, 1 }, // 1
			  { 0, 0, 1, 0 } };// 2
	// Must be same dimensions as validPatternsLeft.
	private static double[][] validPatternsRight = new double[][]
			{ { 0, 1, 0, 0 }, // 4
			  { 1, 0, 0, 0 }, // 8
			  { 0, 0, 0, 1 }, // 1
			  { 0, 0, 1, 0 }, // 2
			  { 1, 1, 0, 1 }, // 13
			  { 1, 1, 1, 0 }, // 14
			  { 1, 1, 0, 0 }, // 12
			  { 0, 0, 0, 0 } };// 0
	
	
	@Override
	public double evaluate(GRNGenome aGenome) {
		if ( ( generation - ( switchRate + lastSwitch ) ) == 0 ) {
			// Force a recompute and switch the fitness test mode
			nonCacheable = true;
			and_mode = !and_mode;
			lastSwitch = generation;
		} else if ( ( generation - ( switchRate + lastSwitch ) ) > 0 ) { 
			nonCacheable = false;			
			// We can cache fitness values again
		}		
				
		GRNModel grn = buildGRNFromGenome(aGenome);
		
		// Make the neural net
		FeedforwardNeuralNetwork nn = new FeedforwardNeuralNetwork();
		
		// Make neurons
		// Input
		for( double k = -1; k <= 1; k += ( 1.0 / (double)numNeuralInputs ) ) {// 8
			nn.addNeuron( new NeuralPosition( k, 0, 1 ) );
		}
		// Output
		nn.addNeuron( new NeuralPosition( 1, 0, 0 ) );// 1
		// Hidden 1
		for( double k = -1; k <= 1; k += ( 1.0 / (double)numNeuralHidden1) ) {// 4
			nn.addNeuron( new NeuralPosition( k, 0, 0.75 ) );
		}
		// Hidden 2
		for( double k = -1; k <= 1; k += ( 1.0 / (double)numNeuralHidden2) ) {// 2
			nn.addNeuron( new NeuralPosition( k, 0, 0.5 ) );
		}
		
		// Make connections
		nn.connectFullFeedforward();
		nn.weightsFromGRN( grn );
		
		System.out.println( nn.numNeurons + "\n" + nn.depthOrder );
		
		double fitness=0.0;
		
		double trueConc = ( 1 );
		double falseConc = ( -1 );

		for (int pattern = 0; pattern < trialCount; pattern++) {
			
			for( int k = 0; k < numNeuralInputs; k++ ) {				
				nn.setNeuron( k, ( inputPatterns[pattern][k] ? trueConc : falseConc ) );
			}
			
			nn.feedforward();
			
			double output = nn.getNeuron( numNeuralInputs );
			
			//boolean leftOutput = grn.proteins.get(numGRNInputs).concentration > booleanThreshold;
			//boolean rightOutput = grn.proteins.get(numGRNInputs+1).concentration > booleanThreshold;
			
			double target = 0;	
			
			if ( and_mode ) {
				target = ( ( targetOutputPatterns[pattern][0] && targetOutputPatterns[pattern][1] ) ? 1 : -1 );
			} else {
				target = ( ( targetOutputPatterns[pattern][0] || targetOutputPatterns[pattern][1] ) ? 1 : -1 );
			}

			fitness += Math.pow( 1, ( target - output ) );
				
			//System.err.println( "Pattern " + pattern + " " + and_mode + " " + targetOutputPatterns[pattern] + " " + output );
			
		}
		
		fitness /= trialCount;
		
		//System.err.println("fitness="+fitness+"  =>  "+grn.toString());
		aGenome.setNewFitness(fitness);

		numEvaluations++;
		return fitness;
	}

}
