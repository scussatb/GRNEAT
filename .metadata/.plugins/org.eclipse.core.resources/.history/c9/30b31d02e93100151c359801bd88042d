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
import java.util.Date;

import evolver.GRNGenome;
import grn.GRNModel;

public class FinancePrediction extends GRNGenomeEvaluator{
	public class Data {
		public long date;
		public double value;

		public Data(long d, double v) {
			date=d;
			value=v;
		}

		public Data() {
			date=0;
			value=0;
		}
	}

	public ArrayList<Data> data=new ArrayList<Data>();

	public double minValue=Double.MAX_VALUE;
	public double maxValue=-Double.MAX_VALUE;
	public double avgValue=0;
	
	double threshold=0.01;
	int lastSwitch=0;
	
	public long startDate;
	public long endDate;
	public long finalDate;

	public int mode=0;
	public int verbose=0;

	SimpleDateFormat sdf=new SimpleDateFormat("yyyy-MM-dd");

	public FinancePrediction() throws Exception {

		// reading the data file
		BufferedReader br;
		String line;
		int i;
		br=new BufferedReader(new FileReader("finance.csv"));
		// omitting first line
		br.readLine();
		line=br.readLine();
		i=0;
		while (line!=null) {
			String tokens[]=line.split(";");
			long date=sdf.parse(tokens[0]).getTime();
			double value=Math.min(25, Double.parseDouble(tokens[1]));
			if (Math.abs(value)>0.001) {
				data.add(new Data(date, value));
				minValue=Math.min(minValue, value);
				maxValue=Math.max(maxValue, value);
				avgValue+=value;
				i++;
			}
			line=br.readLine();
		}
		avgValue/=data.size();

		startDate=sdf.parse("2008-10-30").getTime();
		endDate=sdf.parse("2013-10-20").getTime();
		finalDate=sdf.parse("2013-10-30").getTime();

		// GREAT init
		numGRNInputs=0;
		numGRNOutputs=1;
		name="FinancePrediction";
		//nonCacheable=true;

	}

	public double evaluate(GRNModel grn) {
		grn.reset();
		grn.evolve(25);
		double fitness=0;
		long step=24*3600*1000; // 1 day in milliseconds
		double error=0;
		int nTest=0;
		double concs[]=new double[(int)(((double)(endDate-startDate))/(double)step)+1];
		double maxConc=0;

		for (long date=startDate; date<endDate; date+=step) {
			//grn.proteins.get(0).concentration=getValue(date-step);
			grn.evolve(1);
			concs[nTest]=grn.proteins.get(numGRNInputs).concentration;
			maxConc=Math.max(maxConc, concs[nTest]);
			nTest++;
		}
		maxConc+=0.001;
		double maxErr=0;
		long date=startDate;
		double value=avgValue;
		double prevValue=avgValue;
		double prevRealValue=avgValue;
		nTest=0;
		double realValue=0;
		int nbRealOsc=0;
		int nbPredOsc=0;
/*		if (generation-lastSwitch>=5) {
			threshold/=2.0;
			lastSwitch=generation;
			System.err.println("reducing th to"+threshold);
			nonCacheable=true;
		} else {
			nonCacheable=false;
		}*/
		for (int i=0; i<concs.length; i++) {
			prevValue=value;
			prevRealValue=realValue;
			value=(concs[i]/maxConc)*4+avgValue-2;
			realValue=getValue(date);
			if (Math.abs(value-realValue)<threshold) {
				fitness+=1;
			}
			if ((prevValue<avgValue && value>avgValue)||
				(prevValue>avgValue && value<avgValue)) {
				nbPredOsc++;
			}
			if ((prevRealValue<avgValue && realValue>avgValue)||
				(prevRealValue>avgValue && realValue<avgValue)) {
				nbRealOsc++;
			}
			if (verbose==2) {
				System.out.println(sdf.format(date)+"\t"+value+"\t"+realValue);
			}
			date+=step;
		}
		if (nbPredOsc>=nbRealOsc/2 && nbPredOsc<=nbRealOsc*2) {
			fitness/=Math.abs(nbPredOsc-nbRealOsc)+1;
		} else {
			fitness/=nbRealOsc*2+1;
		}
		return fitness;
	}

	@Override
	public double evaluate(GRNGenome aGenome) {
		GRNModel grn=buildGRNFromGenome(aGenome);
		//		System.err.println(fitness+"\t=>\t"+tempError+"\t+\t"+humidityError);
		double fitness=evaluate(grn);
		aGenome.setNewFitness(fitness);
		numEvaluations++;
		return fitness;
	}

	double getValue(long date) {
		int i=0;

		while (i<data.size() && date>data.get(i).date) {
			i++;
		}
		if (i==0) return data.get(0).value;
		if (i==data.size()) return data.get(i-1).value;
		if (data.get(i).date>=endDate-5) return data.get(i).date;
		double a=(data.get(i).value-data.get(i-1).value)/(data.get(i).date-data.get(i-1).date);
		double b=data.get(i-1).value-a*data.get(i-1).date;
		double res=a*date+b;
		if (res!=res) System.err.println("temp="+res+"  "+i);
		return res;		
	}

	public static void main (String args[]) throws Exception {
		FinancePrediction eval=new FinancePrediction();
		// GREAT1 => Error: 0.06381281579816524
		//		GRNModel grn = GRNModel.loadFromFile("WeatherPrediction/run_1366963891965423000/grn_303_300.0.grn");
		GRNModel grn = GRNModel.loadFromFile("FinancePrediction/run_1383126332992865000/grn_242_34.0.grn");
		// GREAT2 => Error: 0.06381281579816524
		//		GRNModel grn = GRNModel.loadFromFile("WeatherPrediction/run_1366967649950242000/grn_88_639.0.grn");

		// Error: 0.09143223892137296
		//		GRNModel grn = GRNModel.loadFromFile("WeatherPrediction/run_1366790541127817000/grn_267_818.0.grn");
		// Error: 0.07654033568813648
		//		GRNModel grn = GRNModel.loadFromFile("WeatherPrediction/run_1366829792531593000/grn_266_537.0.grn");
		// Error: 0.07205233471544902
		//				GRNModel grn = GRNModel.loadFromFile("WeatherPrediction/run_1366920000934482000/grn_3816_373.0.grn");

		eval.verbose=2;
		eval.endDate=eval.finalDate;
		eval.evaluate(grn);
		/*		double prevTemp=eval.tempData.get(0).value;
		double prevHumi=eval.humidityData.get(0).value;
		double tempError=0;
		double humidityError=0;
		long step=10*60*1000; // 10 minutes in milliseconds
		int nTest=0;
		double curTemp=0;
		double curHumi=0;
		for (long date=eval.startDate; date<eval.endDate; date+=step) {
			double w=eval.getWeatherValue(date);
//			curTemp=prevTemp+prevTemp*0.15*(eval.getWeatherValue(date)-eval.getWeatherValue(date-step));
//			curTemp=21+prevTemp*(eval.getWeatherValue(date)-eval.getWeatherValue(date-step));
			curTemp = 22.41 + 0.009011*w*w*w + 0.002608*w*w*w*w - 0.1138*w*w - 0.0002624*w*w*w*w*w;
			curHumi=prevHumi*(1+eval.getWeatherValue(date)-eval.getWeatherValue(date-step));
			System.out.println(((double)(date)/86400000.0+25569)+"\t"+eval.getWeatherValue(date)+"\t"+curTemp+"\t"+eval.getTempValue(date)+"\t"+curHumi+"\t"+eval.getHumidityValue(date));
			double nTempError=(curTemp-eval.getTempValue(date))*(curTemp-eval.getTempValue(date))/100;
			tempError+=nTempError;
			double nHumError=(curHumi-eval.getHumidityValue(date))*(curHumi-eval.getHumidityValue(date))/10000;
			humidityError+=nHumError;
			prevTemp=curTemp;
			prevHumi=curHumi;
			nTest++;
		}
		tempError=Math.sqrt(tempError/nTest);
		humidityError=Math.sqrt(humidityError/nTest);
		System.out.println((tempError+humidityError)+"  "+tempError+"  "+humidityError);*/
	}

}
