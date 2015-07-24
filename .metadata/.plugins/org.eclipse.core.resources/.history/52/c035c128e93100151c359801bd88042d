package evaluators;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.ObjectInputStream;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.SimpleTimeZone;

import evolver.GRNGenome;
import grn.GRNModel;

public class PingPrediction extends GRNGenomeEvaluator{
	public class Data {
		public double date;
		public double value;

		public Data(double d, double v) {
			date=d;
			value=v;
		}

		public Data() {
			date=0;
			value=0;
		}
	}

	public ArrayList<Data> pingData=new ArrayList<Data>();

	public double minPing=Double.MAX_VALUE;
	public double maxPing=-Double.MAX_VALUE;
	public double avgPing=0;
	public double nEvent=0;

	public double startDate;
	public double endDate;
	public double finalDate;

	public int mode=0;
	public int verbose=0;

	SimpleDateFormat sdf=new SimpleDateFormat("yyyy-MM-dd hh:mm:ss");

	public PingPrediction(int mode) throws Exception {

		// reading the ping file
		BufferedReader br;
		String line;
		int i;
		br=new BufferedReader(new FileReader("train_ping_s2.csv"));
		// omitting first line
		br.readLine();
		line=br.readLine();
		i=0;
		while (line!=null) {
			String data[]=line.split(";");
			double date=Double.parseDouble(data[0]);
			double value=Double.parseDouble(data[1]);
			pingData.add(new Data(date, value));
//			System.err.println("read data:\t"+date+"\t"+value);
			minPing=Math.min(minPing, value);
			maxPing=Math.max(maxPing, value);
			avgPing+=value;
			i++;
			line=br.readLine();
		}
		avgPing/=pingData.size();
		for (i=1; i<pingData.size(); i++) {
			if ((pingData.get(i-1).value<avgPing && pingData.get(i).value>avgPing) ||
				(pingData.get(i-1).value>avgPing && pingData.get(i).value<avgPing)) {
				nEvent++;
			}
		}
//		System.err.println(nEvent);
		
		// train_ping.csv
//		startDate=1367851143;
//		endDate=1367851492;
//		finalDate=1367915000;
		// train_ping_s2.txt
		startDate=1362355200;
		endDate=1365379200;
		finalDate=1367928000;
		

		// GREAT init
		numGRNInputs=1;
		numGRNOutputs=1;
		this.mode=mode;
		name="PingPrediction";

	}

	public double evaluate(GRNModel grn) {
		grn.reset();
		grn.evolve(25);
		double fitness=0;
//		double step=1; // train_ping.csv
		double step=1200; // train_ping_s2.csv
		double pingError=0;
		int nTest=0;
		double pingConcs[]=new double[(int)(((double)(endDate-startDate))/(double)step)+1];
		double maxPingConc=0;
//		System.err.println("============");
		for (double date=startDate; date<endDate; date+=step) {
			Calendar gc=Calendar.getInstance();
			gc.setTime(new Date(((long)date)*1000));
			int day=gc.get(Calendar.DAY_OF_WEEK);
			boolean weekend= (day==Calendar.SATURDAY || day==Calendar.SUNDAY);
			grn.proteins.get(0).concentration=weekend?0.2:0;
			grn.evolve(1);
			pingConcs[nTest]=grn.proteins.get(numGRNInputs).concentration;
			maxPingConc=Math.max(maxPingConc, pingConcs[nTest]);
			nTest++;
		}
		maxPingConc+=0.001;
		double maxPingErr=0;
		double date=startDate;
		double prevPing=avgPing;
		double ping=avgPing;
		double maxPingObserved=-Double.MAX_VALUE;
		nTest=0;
		double realPing=0;
		int nObservedEvent=0;
		int lastEvent=0;
		for (int i=0; i<pingConcs.length; i++) {
			//ping=(pingConcs[i]/maxPingConc)*avgPing*2; // train_ping.csv
			ping=(pingConcs[i]/maxPingConc)*(maxPing-minPing)+minPing; // train_ping_s2.csv
			realPing=getPingValue(date);
			maxPingObserved=Math.max(maxPingObserved, ping);
			if (lastEvent>0 &&
				((prevPing<avgPing && ping>avgPing) ||
				(prevPing>avgPing && ping<avgPing))) {
				nObservedEvent++;
				lastEvent=0;
			} else {
				lastEvent++;
			}
			prevPing=ping;
			fitness-=Math.abs(realPing-ping);
/*			ping=(ping-minPing)/(maxPing-minPing);
			realPing=(realPing-minPing)/(maxPing-minPing);
			if (Math.abs(ping-realPing)<0.05) {
				fitness+=1.0-Math.abs(ping-realPing)/0.05;
			}/**/
			if (verbose==2) {
				System.out.println(date+"\t"+ping+"\t"+realPing);
			}
			date+=step;
		}
		double S;
		if (nObservedEvent==0 || nObservedEvent>nEvent*2) {
			S=0;
		} else {
			S=(double)(Math.abs(nObservedEvent-nEvent))/nEvent;
		}
/**/		fitness*=1.0+S;/**/
/*		fitness*=1.0/(1.0+Math.abs(maxPingObserved-maxPing)/maxPing);/**/
		fitness*=1.0+Math.abs(maxPingObserved-maxPing)/maxPing;
		if (verbose==2) {
			System.out.println("Errors: "+fitness+" ("+nTest+" ; "+nObservedEvent+"/"+nEvent+" ; "+avgPing+")");
		}
		return fitness;
	}

	@Override
	public double evaluate(GRNGenome aGenome) {
		GRNModel grn=buildGRNFromGenome(aGenome);
		double fitness=evaluate(grn);
		aGenome.setNewFitness(fitness);
		numEvaluations++;
		return fitness;
	}

	double getPingValue(double date) {
		int i=0;

		while (i<pingData.size() && date>pingData.get(i).date) {
			i++;
		}
		if (i==0) return pingData.get(0).value;
		if (i==pingData.size()) return pingData.get(i-1).value;
		//		if (pingData.get(i).date>=endDate-5) return pingData.get(i).date;
		double a=(pingData.get(i).value-pingData.get(i-1).value)/(pingData.get(i).date-pingData.get(i-1).date);
		double b=pingData.get(i-1).value-a*pingData.get(i-1).date;
		double res=a*date+b;
		if (res!=res) System.err.println("Ping="+res+"  "+i);
//		System.err.println("real ping="+res);
		return res;		
	}

	public static void main (String args[]) throws Exception {
		PingPrediction eval=new PingPrediction(0);
		// train_ping.csv
//		GRNModel grn = GRNModel.loadFromFile("PingPrediction/run_1367928992580962000/grn_146_-5134.1913657945315.grn");
		// train_ping_s2.csv
//		GRNModel grn = GRNModel.loadFromFile("PingPrediction/run_1368020745221125000/grn_103_-2.616636248833961E9.grn");
		GRNModel grn = GRNModel.loadFromFile("PingPrediction/run_1368022198862182000/grn_119_-5.123315105599529E9.grn");
		eval.verbose=2;
		eval.endDate=eval.finalDate;
//		eval.endDate=1367852143;
		eval.evaluate(grn);
	}

}
