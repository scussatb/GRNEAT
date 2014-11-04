package neural;

import grn.GRNModel;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeMap;
import java.util.TreeSet;

public class FeedforwardNeuralNetwork {		
	public int numNeurons = 0;
	
	public TreeMap<Long,NeuralPosition> focii;// Generally this corresponds to a point-assumption of a neuron, but this need not be true
	public HashMap<Long,Double> activation;
	public HashMap<Long[],Double> weights;
	
	public TreeMap<Double,ArrayList<Long>> depthOrder;// Stores the depth-ordering (and thus the order that neurons can be updated) of neurons, bigger number means closer to input, 1 is an input.
	public TreeMap<Long,ArrayList<Long>> inputSources;// Stores a list of input neuron IDs for each neuron
	
	public FeedforwardNeuralNetwork() {
		
		clear();
	}
	
	public void clear() {
		/*
		 * Empty out all the data from this neural net, including weights.
		 */
		focii = new TreeMap<Long,NeuralPosition>();
		activation = new HashMap<Long,Double>();
		weights = new HashMap<Long[],Double>();
		
		depthOrder = new TreeMap<Double,ArrayList<Long>>();
		inputSources = new TreeMap<Long,ArrayList<Long>>();
	}
	
	public long addNeuron( NeuralPosition pos ) {
		/* 
		 * Add a neuron at position pos and
		 * return the unique identifier of this neuron.
		 */
		long nid = focii.size();
		focii.put( nid, pos );
		activation.put( nid, 0.0 );
		if( !depthOrder.containsKey( pos.z ) ) {
			ArrayList<Long> arry = new ArrayList<Long>();
			
			depthOrder.put( pos.z, arry );
		} else {
			depthOrder.get( pos.z ).add( nid );
		}
		numNeurons++;
		return nid;
	}
	
	public void setNeuron( long nid, double val ) {
		/* 
		 * Set the activation value of a neuron based upon its unique identifier
		 */
		activation.put( nid, val );
	}
	
	public double getNeuron( long nid ) {
		/*
		 * Return the activation value of a neuron
		 */
		return activation.get( nid );
	}
	
	public void setWeight( long source_id, long dest_id, double val ) {
		/* 
		 * Set the weight between the source and dest neurons
		 */
		Long[] edge = { source_id,  dest_id };
		weights.put( edge, val );
		if( inputSources.containsKey( dest_id ) ) {
			inputSources.get( dest_id ).add( source_id );
		} else {
			ArrayList<Long> al = new ArrayList<Long>();
			al.add( source_id );
			inputSources.put( dest_id, al );
		}
		
	}
	
	public void connectFullFeedforward() {
		/*
		 * Create a 0-weight connection between all feedforward neighbors
		 */
		
		Object[] objs = depthOrder.descendingKeySet().toArray();
		Double[] keys = new Double[objs.length]; 
		for( int k = 0; k < objs.length; k++ ) {
			keys[k] = (Double) objs[k];
		}
		
		ArrayList<Long> nidInLayer;
		ArrayList<Long> nidInChildLayer = depthOrder.get( keys[0] );
		
		for ( int k = 0; k < keys.length - 1; k++ ) { // Ignore output layer
			
			nidInLayer = nidInChildLayer;
			nidInChildLayer = depthOrder.get( keys[k+1] );
			
			for ( int sid = 0; sid < nidInLayer.size(); sid++ ) {
				for ( int did = 0; did < nidInChildLayer.size(); did++ ) {
					setWeight( sid, did, 0 );
				}
			}
			
		}
		
		
	}

	public void weightsFromGRN( GRNModel grn ) {
		/*
		 * Updates all the weights of a neural network with a GRN
		 */
	
		for (Map.Entry<Long[], Double> entry : weights.entrySet()) {		
			// Initialize and warm up the GRN
			grn.reset();
			grn.evolve(25);		
			
			// Lookup information about the edge
			Long[] edge = entry.getKey();
			NeuralPosition source = focii.get( edge[0] );
			NeuralPosition dest = focii.get( edge[1] );
			
			// Load the edge into the GRN inputs
			grn.proteins.get(0).setConcentration( source.x );
			grn.proteins.get(1).setConcentration( source.y );
			grn.proteins.get(2).setConcentration( dest.x );
			grn.proteins.get(3).setConcentration( dest.y );
			
			// Run the GRN to stabilize (maybe we should run until stability criteria is satisfied, not just fixed # of iterations)			
			grn.evolve(25);
			
			double posOut = grn.proteins.get(5).concentration;
			double negOut = grn.proteins.get(6).concentration;
			double existOut = grn.proteins.get(7).concentration;
			
			double value = ( posOut > negOut ? 
							( posOut > existOut ? posOut : 0 ) :
							( negOut > existOut ? -negOut : 0 ) );
												 							
			weights.put( edge, value );
			
		}
	}
	
	public void updateNeuron( long nid ) {
		/*
		 * Update a neuron's activation		
		 */
		System.out.println( inputSources.keySet() + "\n" );
		//for( ArrayList<Long> al : inputSources.values() ) {
		for ( Map.Entry<Double,ArrayList<Long>> entry : depthOrder.entrySet() ) {
			System.out.println( entry.getKey() + "\t" + entry.getValue() );
		}
		System.out.println( inputSources.size() );
		ArrayList<Long> theseSources = inputSources.get(nid);
		
		double sum = 0;
		for( int sid = 0; sid < theseSources.size(); sid++ ) {
			Long[] edge = { theseSources.get(sid), nid }; 
			sum += activation.get( theseSources.get(sid) ) * weights.get( edge );
		}
		sum = Math.tanh( sum );
		activation.put( nid, sum );
	}
	
	public void feedforward( ) {
		/*
		 * Perform a feedforward update of the neural network
		 */		
		
		for ( Map.Entry<Double,ArrayList<Long>> entry : depthOrder.entrySet() ) {
			if( entry.getKey() != 1 ) {
				ArrayList<Long> nidInLayer = entry.getValue();
				for( int k = 0; k < nidInLayer.size(); k++ ) {
					updateNeuron( nidInLayer.get(k) );
				}
			}
		}
		
	}
}
