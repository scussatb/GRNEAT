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

public class WeatherPrediction extends GRNGenomeEvaluator{
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

	public ArrayList<Data> tempData=new ArrayList<Data>();
	public ArrayList<Data> humidityData=new ArrayList<Data>();
	public ArrayList<Data> weatherData=new ArrayList<Data>();

	public double minWeather=Double.MAX_VALUE;
	public double maxWeather=-Double.MAX_VALUE;
	public double minHumidity=Double.MAX_VALUE;
	public double maxHumidity=-Double.MAX_VALUE;
	public double minTemp=Double.MAX_VALUE;
	public double maxTemp=-Double.MAX_VALUE;
	public double avgTemp=0;
	public double avgHumidity=0;
	public double avgWeather=0;

	public long startDate;
	public long endDate;
	public long finalDate;

	public int mode=0;
	public int verbose=0;

	SimpleDateFormat sdf=new SimpleDateFormat("yyyy-MM-dd hh:mm:ss");

	public WeatherPrediction(int mode) throws Exception {

		// reading the temperature file
		BufferedReader br;
		String line;
		int i;
		if (mode==0 || mode==2) {
			br=new BufferedReader(new FileReader("train_temperature.csv"));
			// omitting first line
			br.readLine();
			line=br.readLine();
			i=0;
			while (line!=null) {
				String data[]=line.split(";");
				long date=sdf.parse(data[0]).getTime();
				double value=Math.min(25, Double.parseDouble(data[1]));
				if (Math.abs(value)>0.001) {
					tempData.add(new Data(date, value));
					minTemp=Math.min(minTemp, value);
					maxTemp=Math.max(maxTemp, value);
					avgTemp+=value;
					i++;
				}
				line=br.readLine();
			}
			avgTemp/=tempData.size();
		}

		if (mode==1 || mode==2) {
			// reading the humidity file
			br=new BufferedReader(new FileReader("train_humidity.csv"));
			// omitting first line
			br.readLine();
			line=br.readLine();
			i=0;
			while (line!=null) {
				String data[]=line.split(";");
				long date=sdf.parse(data[0]).getTime();
				double value=Double.parseDouble(data[1]);
				if (Math.abs(value)>0.001) {
					humidityData.add(new Data(date, value));
					minHumidity=Math.min(minHumidity,value);
					maxHumidity=Math.max(maxHumidity, value);
					avgHumidity+=value;
					i++;
				}
				line=br.readLine();
			}
			avgHumidity/=humidityData.size();
		}
		// reading the weather file
		br=new BufferedReader(new FileReader("train_weather.csv"));
		// omitting first line
		br.readLine();
		line=br.readLine();
		i=0;
		while (line!=null) {
			String data[]=line.split(";");
			long date=sdf.parse(data[0]).getTime();
			double value=Double.parseDouble(data[1]);
			weatherData.add(new Data(date, value));
			minWeather=Math.min(minWeather, value);
			maxWeather=Math.max(maxWeather, value);
			avgWeather+=value;
			i++;
			line=br.readLine();
		}
		avgWeather/=weatherData.size();

		startDate=sdf.parse("2013-02-01 14:00:00").getTime();
		endDate=sdf.parse("2013-02-18 22:59:00").getTime();
		finalDate=sdf.parse("2013-02-22 00:00:01").getTime();

		// GREAT init
		if (mode==2) {
			numGRNInputs=3;
			numGRNOutputs=2;
		} else {
			numGRNInputs=1;
			numGRNOutputs=1;
		}
		this.mode=mode;
		name="WeatherPrediction";

	}

	public double evaluate(GRNModel grn) {
		grn.reset();
		grn.evolve(25);
		double fitness=0;
		long step=10*60*1000; // 10 minutes in milliseconds
		double tempError=0;
		double humiError=0;
		int nTest=0;
		double tempConcs[]=new double[(int)(((double)(endDate-startDate))/(double)step)+1];
		double humiConcs[]=new double[(int)(((double)(endDate-startDate))/(double)step)+1];
		double maxTempConc=0;
		double maxHumiConc=0;

		for (long date=startDate; date<endDate; date+=step) {
			grn.proteins.get(0).concentration=0.2*((getWeatherValue(date)-minWeather)/(maxWeather-minWeather));
			grn.evolve(1);
			switch (mode) {
			case 0:
				tempConcs[nTest]=grn.proteins.get(numGRNInputs).concentration;
				maxTempConc=Math.max(maxTempConc, tempConcs[nTest]);
				break;
			case 1:
				humiConcs[nTest]=grn.proteins.get(numGRNInputs).concentration;
				maxHumiConc=Math.max(maxHumiConc, humiConcs[nTest]);
				break;
			case 2:
				tempConcs[nTest]=grn.proteins.get(numGRNInputs).concentration;
				maxTempConc=Math.max(maxTempConc, tempConcs[nTest]);
				humiConcs[nTest]=grn.proteins.get(numGRNInputs+1).concentration;
				maxHumiConc=Math.max(maxHumiConc, humiConcs[nTest]);
				break;
			default:
				break;
			}
			nTest++;
		}
		maxTempConc+=0.001;
		maxHumiConc+=0.001;
		double maxTempErr=0;
		double maxHumiErr=0;
		long date=startDate;
		double prevTemp=avgTemp;
		double temp=avgTemp;
		double prevHumi=avgHumidity;
		double humi=avgHumidity;
		nTest=0;
		double realTemp=0;
		double realHumi=0;
		for (int i=0; i<tempConcs.length; i++) {
			try {
				if (mode==0||mode==2) {
					prevTemp=temp;
					temp=(tempConcs[i]/maxTempConc)*20+avgTemp-10; // concentration between 0 and 0.1
					//			temp=(tempConcs[i]/maxTempConc)*20.0+prevTemp-10.0; // concentration between 0 and 0.1
					//			temp=(tempConcs[i]/maxTempConc)*20.0+getTempValue(date-7*3600*24*1000)-10.0; // concentration between 0 and 0.1
					realTemp=getTempValue(date);
					if (date>=sdf.parse("2013-02-15 00:00:00").getTime() && 
							date<=sdf.parse("2013-02-18 00:00:00").getTime()) {
						double nTempError=(temp-realTemp)*(temp-realTemp)/100;
						tempError+=nTempError;
						maxTempErr=Math.max(maxTempErr, nTempError);
						nTest++;
					}
					if (Math.abs(temp-realTemp)<0.1) {
						fitness+=1;
					}
				}
				if (mode==1||mode==2) {
					prevHumi=humi;
					humi=avgHumidity+(humiConcs[i]/maxHumiConc)*40.0-20.0;
					realHumi=getHumidityValue(date);
					if (date>=sdf.parse("2013-02-15 00:00:00").getTime() && 
							date<=sdf.parse("2013-02-18 00:00:00").getTime()) {
						double nHumiError=(humi-realHumi)*(humi-realHumi)/100;
						humiError+=nHumiError;
						maxHumiErr=Math.max(maxHumiErr, nHumiError);
						nTest++;
					}
					if (Math.abs(humi-realHumi)<1) {
						fitness+=1;
					}
				}
				if (mode==2) {
					nTest/=2;
				}
				if (verbose==2) {
					switch (mode) {
					case 0:
						System.out.println(temp+"\t"+realTemp);
						break;
					case 1:
						System.out.println(humi+"\t"+realHumi);
						break;
					case 2:
						System.out.println(((double)(date)/86400000.0+25569)+"\t"+getWeatherValue(date)+temp+"\t"+realTemp+"\t"+humi+"\t"+realHumi);
						break;
					
					}
				}
			} catch (ParseException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			date+=step;
		}
		if (verbose==2) {
			System.out.println("Errors: "+Math.sqrt(tempError/nTest)+" ("+nTest+")"+"  "+Math.sqrt(humiError/nTest));
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

	double getTempValue(long date) {
		int i=0;

		while (i<tempData.size() && date>tempData.get(i).date) {
			i++;
		}
		if (i==0) return tempData.get(0).value;
		if (i==tempData.size()) return tempData.get(i-1).value;
		if (tempData.get(i).date>=endDate-5) return tempData.get(i).date;
		double a=(tempData.get(i).value-tempData.get(i-1).value)/(tempData.get(i).date-tempData.get(i-1).date);
		double b=tempData.get(i-1).value-a*tempData.get(i-1).date;
		double res=a*date+b;
		if (res!=res) System.err.println("temp="+res+"  "+i);
		return res;		
	}

	double getHumidityValue(long date) {
		int i=0;
		while (i<humidityData.size() && date>humidityData.get(i).date) {
			i++;
		}
		if (i==0) return humidityData.get(0).value;
		if (i==humidityData.size()) return humidityData.get(i-1).value;
		double a=(humidityData.get(i).value-humidityData.get(i-1).value)/(humidityData.get(i).date-humidityData.get(i-1).date);
		double b=humidityData.get(i-1).value-a*humidityData.get(i-1).date;
		double res=a*date+b;
		if (res!=res) System.err.println("temp="+res);
		return res;		
	}

	double getWeatherValue(long date) {
		int i=0;
		while (i<weatherData.size() && date>weatherData.get(i).date) {
			i++;
		}
		if (i==0) return weatherData.get(0).value;
		if (i==weatherData.size()) return weatherData.get(i-1).value;
		double a=(weatherData.get(i).value-weatherData.get(i-1).value)/(weatherData.get(i).date-weatherData.get(i-1).date);
		double b=weatherData.get(i-1).value-a*weatherData.get(i-1).date;
		double res=a*date+b;
		if (res!=res) System.err.println("temp="+res);
		return res;		
	}

	public static void main (String args[]) throws Exception {
		WeatherPrediction eval=new WeatherPrediction(1);
		// GREAT1 => Error: 0.06381281579816524
//		GRNModel grn = GRNModel.loadFromFile("WeatherPrediction/run_1366963891965423000/grn_303_300.0.grn");
		GRNModel grn = GRNModel.loadFromFile("WeatherPrediction/run_1366982736725785000/grn_530_622.0.grn");
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
