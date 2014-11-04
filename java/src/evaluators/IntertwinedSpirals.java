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

public class IntertwinedSpirals extends GRNGenomeEvaluator {
	double coef=3.0;
	int num_samples = 97;// 97 is standard, and required with this implementation, which shifts both spirals to the positive quadrant
	double inputCoef=0.5;

	double[] data;
	
	public boolean fuzzyFitness = false; // false means +1 point per correct answer 0 for wrong, true means float-distance from prediction to target
	

	public IntertwinedSpirals( String args[] ) {
		numGRNInputs=2;
		numGRNOutputs=2;
		name="IntertwinedSpirals";
		
		for (int k = 0; k < args.length; k++ ) {
//			System.err.print("\t" + args[k]);

			if ( args[k].compareTo("fuzzyFitness") == 0) {
				fuzzyFitness = Boolean.parseBoolean( args[k+1] );
			}
		}				

		System.out.println("fuzzy fitness:      "+fuzzyFitness);
		
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
	public double evaluate(GRNGenome aGenome) {
		GRNModel grn=buildGRNFromGenome(aGenome);

		double fitness = -num_samples * 2;
		grn.reset();
		grn.evolve(25);// Warmup

		for( int k = 0; k < num_samples * 2; k++ ) {

			GRNModel localGRN=grn.copy();

			// Set inputs
			localGRN.proteins.get(0).concentration = inputCoef * data[ 3 * k ];
			localGRN.proteins.get(1).concentration = inputCoef * data[ 3 * k + 1 ];
			

			localGRN.evolve(25);// Compute answer

			int targetClass = (int) data[ 3 * k + 2 ];
			int predictedClass;

			if( localGRN.proteins.get(2).concentration > localGRN.proteins.get(3).concentration ) {
				predictedClass = 1;
			} else {
				predictedClass = 0;
			}

			double target_value = ( targetClass == 0 ? 1 : -1 );
			double predicted_value  = Math.tanh( grn.proteins.get(2).concentration - grn.proteins.get(3).concentration ); 
			
			/*if( fuzzyFitness ) {
				fitness += Math.abs( target_value - predicted_value ) / 2.0;// because range is [-1,1]
			} else {
				if( predictedClass == targetClass ) {
					fitness += 1.0;
				}
			}*/
			
			
			if( predictedClass == targetClass ) {
				fitness += 1.0;
			} else if( fuzzyFitness ) {
				fitness -= Math.abs( target_value - predicted_value ) / 2.0;// because range is [-1,1]
			}

			//System.err.println( k + "\t" + data[ 3 * k ] + "\t" + data[ 3 * k + 1] + "\t" + targetClass + "\t" + predictedClass );

		}

		//aGenome.setNewFitness( fitness );
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

	public void drawSpiralImage( GRNModel grn, String filename, double d) {
		double maxx = 1, maxy = 1;
		double dx = d, dy = d;
//		double dx = 0.0025, dy = 0.0025;

		BufferedImage img = new BufferedImage( (int) ( maxx / dx ), (int) ( maxy / dy ), BufferedImage.TYPE_INT_RGB );
		// Plot the GRN
		grn.reset();
		grn.evolve(25);// Warmup

		for (int i=0; i<(int) ( maxx / dx ); i++) {
//			System.err.println(i+" lines of "+(int) ( maxx / dx )+" generated");
			for (int j=0; j<(int) ( maxy / dy ); j++) {
				double x=dx*i;
				double y=dy*j;
				GRNModel localGRN=grn.copy();

				// Set inputs
				localGRN.proteins.get(0).concentration = inputCoef * ( x );
				localGRN.proteins.get(1).concentration = inputCoef * ( y );

				localGRN.evolve(25);// Compute answer

				//int targetClass = (int) data[ 3 * k + 2 ];
				int predictedClass;				

				if( localGRN.proteins.get(2).concentration > localGRN.proteins.get(3).concentration ) {
					predictedClass = 1;
				} else {
					predictedClass = 0;
				}

				int r = 255 * predictedClass;
				int g = 0;
				int b = 255 * ( 1 - predictedClass);

				int pixel = ( r << 16 ) | ( g << 8 ) | b; 

				//System.err.println( x + "\t" + y + "\t" + (int) (x/dx) + "\t" + (int) (y/dy) + "\t" + predictedClass );

				img.setRGB( i, j, pixel );
			}
		}

		try {
			// retrieve image		    
			File outputfile = new File( filename );
			ImageIO.write(img, "png", outputfile);
		} catch (IOException e) {
			System.out.println( "Exception: " + e );
		}
	}

	public void drawSampledImage(GRNModel grn, String filename, BufferedImage img, double d) {
		double fitness=-num_samples * 2;
		int markerRadius = 6;
		double maxx = 1, maxy = 1;
		double dx = d, dy = d;
		if (img==null) {
			img = new BufferedImage( (int) ( maxx / dx ), (int) ( maxy / dy ), BufferedImage.TYPE_INT_RGB );
			for (int i=0; i<(int) ( maxx / dx ); i++) {
				for (int j=0; j<(int) ( maxy / dy ); j++) {
					img.setRGB(i, j, Integer.MAX_VALUE); 
				}
			}
		}
		
		grn.reset();
		grn.evolve(25);// Warmup

		int targetClass;
		for( int k = 0; k < num_samples * 2; k++ ) {


			// Set inputs
			double x = data[ 3 * k ];
			double y = data[ 3 * k + 1 ];			    			

//			System.out.println( "Sample " + k + " " + x + "," + y);

			// GRN marks
			if (grn!=null) {
				GRNModel localGRN=grn.copy();
				// Set inputs
				localGRN.proteins.get(0).concentration = inputCoef * data[ 3 * k ];
				localGRN.proteins.get(1).concentration = inputCoef * data[ 3 * k + 1 ];
				localGRN.evolve(25);// Compute answer
				if( localGRN.proteins.get(2).concentration > localGRN.proteins.get(3).concentration ) {
					targetClass = 1;
				} else {
					targetClass = 0;
				}
			} else {
				targetClass = (int) data[ 3 * k + 2 ];
			}
			if (targetClass==(int) data[ 3 * k + 2 ]) {
				fitness+=1.0;
			}
			
			int r = 255 * targetClass;
			int g = 0;
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

					img.setRGB( markx, marky, pixel );
				}
			}

			//System.err.println( k + "\t" + data[ 3 * k ] + "\t" + data[ 3 * k + 1] + "\t" + targetClass + "\t" + predictedClass );

		}	
		System.out.println("fitness="+fitness / ( num_samples * 2.0 ));
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

//		GRNModel grn = GRNModel.loadFromFile("grn_298_-0.28865979381443296.grn");
//		GRNModel grn = GRNModel.loadFromFile("grn_344_-0.27319587628865977.grn");
		//		GRNModel grn = GRNModel.loadFromFile("IntertwinedSpirals/run_1368526741478672000/grn_118_-0.36597938144329895.grn");

//		System.out.println(grn.toString());

		GRNModel grn = GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher018/IntertwinedSpirals/run_15117926140547400/grn_499_-0.2422680412371134.grn");
		IntertwinedSpirals eval=new IntertwinedSpirals( args );
		eval.drawSampledImage(grn, "/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher018/ep.png", null, 0.0025);
		
/*		for (int i=0; i<22; i++) {
			GRNModel grn=GRNModel.loadFromFile("launcher018/Generalization/IS/GREAT/grn_"+i+".grn");
			System.out.println( "Drawing image GREAT "+i );
			eval.drawSpiralImage( grn, "spiral_extended_GREAT_"+i+".png", 0.0025 );
			eval.drawSampledImage(grn, "spiral_grn_GREAT_"+i+".png", null, 0.0025);
			File inputFile=new File("spiral_extended_GREAT_"+i+".png");
			eval.drawSampledImage(grn, "spiral_extendedSampled_GREAT_"+i+".png", ImageIO.read(inputFile), 0.0025);
			System.out.println( "Done drawing" );
			
			grn=GRNModel.loadFromFile("launcher018/Generalization/IS/GA/grn_"+i+".grn");
			System.out.println( "Drawing image GA "+i );
			eval.drawSpiralImage( grn, "spiral_extended_GA_"+i+".png", 0.0025 );
			eval.drawSampledImage(grn, "spiral_grn_GA_"+i+".png", null, 0.0025);
			inputFile=new File("spiral_extended_GA_"+i+".png");
			eval.drawSampledImage(grn, "spiral_extendedSampled_GA_"+i+".png", ImageIO.read(inputFile), 0.0025);
			System.out.println( "Done drawing" );
		}

/**/
	}
}
