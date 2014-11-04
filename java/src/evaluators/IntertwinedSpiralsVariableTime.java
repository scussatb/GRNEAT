package evaluators;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import evolver.GRNGenome;
import grn.GRNModel;

/*
;;;;;;;;;;;;
;; Intertwined Spirals problem
;; Introduced in "Learning to tell two spirals apart," Lang, K.J. and Witbrock, M.J., Proceedings of the 1988 Connectionist models summer school, 1988.
*/

public class IntertwinedSpiralsVariableTime extends GRNGenomeEvaluator {
	double coef=3.0;
	int num_samples = 97;// 97 is standard, and required with this implementation, which shifts both spirals to the positive quadrant
	public static long MAX_EVAL_STEPS = 100;// The GRN specifies when to stop, but this is the upper limit
	
	public boolean fuzzyFitness = false; // false means +1 point per correct answer 0 for wrong, true means float-distance from prediction to target
	
	double[] data;
	
	public IntertwinedSpiralsVariableTime( String args[] ) {
		numGRNInputs=2;
		//numGRNOutputs=2;
		numGRNOutputs=4;
		name="IntertwinedSpiralsVariableTime";
		
		for (int k = 0; k < args.length; k++ ) {
//			System.err.print("\t" + args[k]);

			if ( args[k].compareTo("fuzzyFitness") == 0) {
				fuzzyFitness = Boolean.parseBoolean( args[k+1] );
			}
		}				
		
		
		initDataset();
	}

	public void initDataset() {
		data = new double[ 194 * 3 ];
	
		double minx = 1000, miny = 1000;
		double maxx = 0, maxy = 0;
		double xshift = 6.5, yshift = 6.5;
		
		for( int k = 0; k < num_samples; k++ ) {
			double angle = ( k * Math.PI ) / 16.0;
			double radius = (6.5 * ((104.0 - k) / 104.0));
			double x1 = ( ( radius * Math.sin( angle ) ) + xshift ) / ( 2 * xshift );
			double y1 = ( ( radius * Math.cos( angle ) ) + yshift ) / ( 2 * yshift );
		
			double x2 = ( -( radius * Math.sin( angle ) ) + xshift ) / ( 2 * xshift );
			double y2 = ( -( radius * Math.cos( angle ) ) + yshift ) / ( 2 * yshift );
		
			
			data[ 3*k ] = x1;
			data[ 3*k + 1 ] = y1;
			data[ 3*k + 2 ] = 1;
			
			data[ 3*( k + num_samples ) ] = x2;
			data[ 3*( k + num_samples ) + 1] = y2;
			data[ 3*( k + num_samples ) + 2] = 0;
			
			if( x1 < minx ) minx = x1;
			if( y1 < miny ) miny = y1;
			
			if( x2 < minx ) minx = x2;
			if( y2 < miny ) miny = y2;
			
			if( x1 > maxx ) maxx = x1;
			if( y1 > maxy ) maxy = y1;
			
			if( x2 > maxx ) maxx = x2;
			if( y2 > maxy ) maxy = y2;
			
			//System.out.println( "Angle: " + angle + " radius: " + radius + " x: " + x + " y:" + y );
			
		}		
		
		//System.out.println( "minx,miny " + minx + " " + miny );
		//System.out.println( "maxx,maxy " + maxx + " " + maxy );
				
	}
	
	@Override
	// Evaluation based on a fixed # of steps
	/*public double evaluate(GRNGenome aGenome) {
		GRNModel grn=buildGRNFromGenome(aGenome);
		
		double fitness = -num_samples * 2;
		
		for( int k = 0; k < num_samples * 2; k++ ) {
			
			grn.reset();
			grn.evolve(25);// Warmup
			
			// Set inputs
			grn.proteins.get(0).concentration = 0.5 * data[ 3 * k ];
			grn.proteins.get(1).concentration = 0.5 * data[ 3 * k + 1 ];
			    			
			grn.evolve(25);// Compute answer
			
			int targetClass = (int) data[ 3 * k + 2 ];
			int predictedClass;
			
			if( grn.proteins.get(2).concentration > grn.proteins.get(3).concentration ) {
				predictedClass = 1;
			} else {
				predictedClass = 0;
			}
			    
			double error_signal  = grn.proteins.get(2).concentration - grn.proteins.get(3).concentration * ( targetClass == 0 ? 1 : -1 ); 
			
			fitness += Math.max( -0.5, Math.min( 0.5, error_signal ) );
			
			if( predictedClass == targetClass ) {
				//fitness += 1.0;
				
			}
			
			//System.err.println( k + "\t" + data[ 3 * k ] + "\t" + data[ 3 * k + 1] + "\t" + targetClass + "\t" + predictedClass );
						
		}
				
		aGenome.setNewFitness( fitness / ( num_samples * 2.0 ) );// fitness is negative 1 for each wrong answer out of 194 answers
		
		GRNGenomeEvaluator.numEvaluations++;
		return fitness;
	}*/
	
	public double evaluate(GRNGenome aGenome) {
		GRNModel warmGrn=buildGRNFromGenome(aGenome);
		
		warmGrn.reset();
		warmGrn.evolve(25);// Warmup
		
		double fitness = -num_samples * 2;
		
		for( int k = 0; k < num_samples * 2; k++ ) {
			
			GRNModel grn = warmGrn.copy();
			
			// Set inputs
			grn.proteins.get(0).concentration = 0.25 * data[ 3 * k ];
			grn.proteins.get(1).concentration = 0.25 * data[ 3 * k + 1 ];
			grn.proteins.get(2).concentration = 0.25 * 1;
			grn.proteins.get(3).concentration = 0;
			    			
			int numSteps = 0;
			while( ( grn.proteins.get(4).concentration > grn.proteins.get(5).concentration ) && ( numSteps < 100 ) ) {
				grn.evolve(1);
				numSteps++;
			}			
			
			int targetClass = (int) data[ 3 * k + 2 ];
			int predictedClass;
			
			if( grn.proteins.get(2).concentration > grn.proteins.get(3).concentration ) {
				predictedClass = 1;
			} else {
				predictedClass = 0;
			}
			    
			double target_value = ( targetClass == 0 ? 1 : -1 );
			double predicted_value  = grn.proteins.get(2).concentration - grn.proteins.get(3).concentration; 
			// error signal is negative if signs are wrong
			
			if( fuzzyFitness ) {
				fitness += Math.abs( target_value - predicted_value ) / 2.0;// because range is [-1,1]
			} else {
				if( predictedClass == targetClass ) {
					fitness += 1.0;
				}
			}		
			
			//System.err.println( k + "\t" + data[ 3 * k ] + "\t" + data[ 3 * k + 1] + "\t" + targetClass + "\t" + predictedClass );
		}
				
		aGenome.setNewFitness( fitness / ( num_samples * 2.0 ) );// fitness is negative 1 for each wrong answer out of 194 answers

		/*System.err.println( aGenome.getLastFitness() );
		
		if( aGenome.getLastFitness() == 0 ) {
			for( int k = 0; k < num_samples * 2; k++ ) {
				
				grn.reset();
				grn.evolve(25);// Warmup
				
				// Set inputs
				grn.proteins.get(0).concentration = 0.5 * data[ 3 * k ];
				grn.proteins.get(1).concentration = 0.5 * data[ 3 * k + 1 ];
				    			
				grn.evolve(25);// Compute answer
				
				int targetClass = (int) data[ 3 * k + 2 ];
				int predictedClass;
				
				if( grn.proteins.get(2).concentration > grn.proteins.get(3).concentration ) {
					predictedClass = 1;
				} else {
					predictedClass = 0;
				}
				    				
				System.err.println( k + "\t" + data[ 3 * k ] + "\t" + data[ 3 * k + 1] + "\t" + targetClass + "\t" + predictedClass );
							
			}			
		}*/
		
		GRNGenomeEvaluator.numEvaluations++;
		return fitness;
	}	
	
	public void drawSpiralImage( GRNModel grn, String filename ) {
		double maxx = 1, maxy = 1;
		//double dx = 0.05, dy = 0.05;
		double dx = 0.0025, dy = 0.0025;
		//double dx = 0.005, dy = 0.005;
		
		BufferedImage img = new BufferedImage( (int) ( maxx / dx ), (int) ( maxy / dy ), BufferedImage.TYPE_INT_RGB );
		
		GRNModel warmGrn= grn.copy();
		
		warmGrn.reset();
		warmGrn.evolve(25);// Warmup
								
			
		
		// Plot the GRN
		for( double x = 0; x < maxx; x+= dx ) {
			for( double y = 0; y < maxy; y += dy ) {
				
				grn = warmGrn.copy();
				
				grn.proteins.get(0).concentration = 0.25 * x;
				grn.proteins.get(1).concentration = 0.25 * y;
				grn.proteins.get(2).concentration = 0.25 * 1;
				grn.proteins.get(3).concentration = 0;
				
				int numSteps = 0;
				while( ( grn.proteins.get(4).concentration > grn.proteins.get(5).concentration ) && ( numSteps < MAX_EVAL_STEPS ) ) {
					grn.evolve(1);
					numSteps++;
				}			
								
				int r = (int) Math.round( 255.0 * ( grn.proteins.get(2).concentration - grn.proteins.get(3).concentration > 0 ? 
						//( grn.proteins.get(2).concentration - grn.proteins.get(3).concentration ) : 0 ) );
						1 : 0 ) );
				int g = 0;//255 * predictedClass;
				int b = (int) Math.round( 255.0 * ( grn.proteins.get(2).concentration - grn.proteins.get(3).concentration <= 0 ? 
						//( grn.proteins.get(2).concentration - grn.proteins.get(3).concentration ) : 0 ) );
						1 : 0 ) );
				
				int pixel = ( r << 16 ) | ( g << 8 ) | b; 
				
				//System.err.println( x + "\t" + y + "\t" + (int) (x/dx) + "\t" + (int) (y/dy) + "\t" + predictedClass );
				
				img.setRGB( (int) (x/dx), (int) (y/dy), pixel );
			}
		}
		
		// Now draw the the true spirals
		int markerRadius = 6;
		
		for( int k = 0; k < num_samples * 2; k++ ) {
			

			// Set inputs
			double x = data[ 3 * k ];
			double y = data[ 3 * k + 1 ];			    			

			System.out.println( "Sample " + k + " " + x + "," + y);
			
			int targetClass = (int) data[ 3 * k + 2 ];
						
			int r = 255 * targetClass;
			int g = 127;
			int b = 255 * ( 1 - targetClass);
			
			int pixel = ( r << 16 ) | ( g << 8 ) | b; 
			
			for( int mx = -markerRadius; mx <= markerRadius; mx++ ) {
				for( int my = -markerRadius; my <= markerRadius; my++ ) {					

					int markx = Math.max( 0, 
										  Math.min( (int) (maxx / dx) - 1, 
												    (int) (x/dx) + mx ) );
					int marky = Math.max( 0, 
										  Math.min( (int) (maxy / dy) - 1, 
												    (int) (y/dy) + my ) );					
					
					//System.out.println( "Marker " + markx + "," + marky );
					
					int currpix = 255;//img.getRGB( markx, marky);					
					
					img.setRGB( markx, marky, pixel );// & currpix );
				}
			}
			
			//System.err.println( k + "\t" + data[ 3 * k ] + "\t" + data[ 3 * k + 1] + "\t" + targetClass + "\t" + predictedClass );
						
		}
		
		try {
		    // retrieve image		    
		    File outputfile = new File( filename );
		    ImageIO.write(img, "png", outputfile);
		} catch (IOException e) {
		    System.out.println( "Exception: " + e );
		}
		
		
		
	}
	
	
	public static void main(String args[]) throws Exception {		
		//GRNModel grn = GRNModel.loadFromFile("IntertwinedSpirals/run_1366656050552666000/grn_1_0.0.grn");
		
		//GRNModel grn = GRNModel.loadFromFile("IntertwinedSpirals/run_1366679988257973000/grn_9_-0.39690721649484534.grn");
		
		//GRNModel grn = GRNModel.loadFromFile("grn_161_-0.29896907216494845.grn");
		//GRNModel grn = GRNModel.loadFromFile("grn_221_-0.27835051546391754.grn");
		GRNModel grn = GRNModel.loadFromFile("grn_499_-0.26288659793814434.grn");
		
		System.out.println(grn.toString());
		
		IntertwinedSpiralsVariableTime eval=new IntertwinedSpiralsVariableTime( args );
		
		System.out.println( "Drawing image" );
		eval.drawSpiralImage( grn, "testspiral.png" );
		System.out.println( "Done drawing" );
		
//		double fit=eval.evolveGRN(grn, true, 125, 1000);
//		double fit=eval.evolveGRN(grn, true, 500, 1000);
//		double fit=eval.evolveGRN_zero(grn, true, 1000);
//		double fit=eval.evolveGRN(grn, true, 7500, 15000);
//		double fit=eval.evolveGRN(grn, true, 83.33333, 1000);
				
	}

}
