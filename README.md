
# Gene Regulatory Network Evolution Through Augmenting Topologies

# Quickstart:

This code was originally made in Eclipse, and works well there with the standard "Run" mechanism.

Go to: java/evolver/Evolver.java and press Run

# To change experiments: 

Go to the main() function in Evolver.java

Change line:
			e.evaluator = new IntertwinedSpirals( args );
			
To the evaluator of your choice.

The main function was designed for command line usage in a .jar, but by swapping the evaluator class that you initialize you can easily change the problem within the IDE of your choice.

# Citing:

Cussat-Blanc, S., Harrington, K., and Pollack, J. (2015) Gene Regulatory Network Evolution Through Augmenting Topologies. IEEE Transactions on Evolutionary Computation 19(6), pp. 823 - 837.
http://dx.doi.org/10.1109/TEVC.2015.2396199
