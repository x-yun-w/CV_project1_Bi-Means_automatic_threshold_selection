#include <fstream>
#include <iostream>
#include <string>
using namespace std;

class thresholdSelection {
public:
	int numR;
	int numC;
	int minV;
	int maxV;
	int biThrV;
	int histH;
	int maxH;
	int* histAry;
	int* gaussAry;
	int* bestFitGaussAry;
	char** graph;

public:
	thresholdSelection(ifstream& in) {
		
		
	}

	int loadHist(ifstream& in) {
		int a = 0;
		int b = 0;
		int m = 0;
		while (in >> a && in >> b) {
			histAry[a] = b;
			if(b > m)
				m = histAry[a];
		}
		return m;
	}

	void copyAry(int* a, int* b) {
		for (int i = 0; i < maxV + 1; i++) {
			b[i] = a[i];
		}
	}

	void dispHist(ofstream& out) {
		out << numR << " " << numC << " " << minV << " " << maxV << "\n";
		
		for (int i = 0; i <= maxV; i++) {
			out << i << " (" << histAry[i] << "):";
			for (int j = 1; j <= histAry[i]; j++) {
				out << "+";
			}
			out << "\n";
		}
	}

	void plotHist() {
		for (int i = 0; i < maxV + 1; i++) {
			for (int j = 0; j < histAry[i]; j++) {
				graph[i][j] = '+';
			}
			
		}
	}

	

	void setZero(int* a) {
		for (int i = 0; i < maxV + 1; i++) {
			a[i] = 0;
		}
	}

	void printGraph(ofstream& out) {
		for (int i = 0; i < maxV + 1; i++) {
			for (int j = 0; j < histH + 1; j++) {
				out<< graph[i][j];
			}
			out << endl;
		}
	}
	int biGaussian(ofstream& debug) {
		debug << "Entering biGaussian method\n";
		double sum1;
		double sum2;
		double total;
		double minSumDiff;
		int offSet = (maxV - minV) / 10;
		int dividePt = offSet;
		int bestThr = dividePt;
		minSumDiff = 999999.0;
		while (dividePt < (maxV - offSet)) {
			setZero(gaussAry);
			sum1 = fitGauss(0, dividePt, debug);
			sum2 = fitGauss(dividePt, maxV, debug);
			total = sum1 + sum2;
			if (total < minSumDiff) {
				minSumDiff = total;
				bestThr = dividePt;
				copyAry(gaussAry, bestFitGaussAry);
			}
			debug << "dividePt: " << dividePt << ", sum1: " << sum1 << ", sum2: " << sum2 << ", total: " << total << ", minSumDiff: " << minSumDiff << ", bestThr: " << bestThr << "\n";
			dividePt++;
		}
		debug << "leaving biGaussian method, minSumDiff = " << minSumDiff << "   bestThr = " << bestThr << "\n";
		return bestThr;
	}

	double fitGauss(int l, int r, ofstream& debug) {
		debug << "Entering fitGauss method\n";
	
		double mean;
		double v;
		double sum = 0.0;
		double gv;
		mean = computeMean(l, r, debug);
		v = computeVar(l, r, mean, debug);
		int i = l;
		while (i <= r) {
			gv = modifiedGauss(i, mean, v, maxH);
			sum += abs(gv - (double)(histAry[i]));
			gaussAry[i] = (int)(gv);
			i++;
		}
		debug << "leaving fitGauss method, sum is: " << sum << "\n";
		return sum;
	}

	double computeMean(int l, int r, ofstream& debug) {
		debug << "Entering computeMean method\n";
		maxH = 0;
		int sum = 0;
		int numP = 0;
		int i = l;
		while (i < r) {
			sum += (histAry[i] * i);
			numP += histAry[i];

			if (histAry[i] > maxH)
				maxH = histAry[i];

			i++;
		}
		double result = ((double)sum) / ((double)numP);
		debug << "Leaving computeMean method maxHeight is: " << maxH << "   result is: " << result << "\n";
		return result;
	}

	double computeVar(int l, int r, double mean, ofstream& debug) {
		debug << "Entering computeVar method\n";
		double sum = 0.0;
		int numP = 0;
		int i = l;
		while (i < r) {
			sum += (double)(histAry[i]) * pow((double)(i) - mean, 2);
			numP += histAry[i];
			i++;
		}
		double result = sum / ((double)numP);
		debug << "Leaving computeVar method returning result: " << result << "\n";
		return result;
	}

	double modifiedGauss(int x, double mean, double v, int maxH) {
		return (double)(maxH * exp(-1 * (pow(((double)x) - mean, 2)) / (2 * v)));
	}

	void plotGaussCurves(ofstream& debug) {
		debug << "Entering plotGaussCurves () method" << endl;

		int index = 0;

		while (index <= maxV)
		{
			int end1 = 0;
			int end2 = 0;
			if (bestFitGaussAry[index] <= histAry[index]) {
				end1 = bestFitGaussAry[index];
				end2 = histAry[index];
			}
			else {
				end1 = histAry[index];
				end2 = bestFitGaussAry[index];
			}

			int i = end1;
			while (i <= end2)
			{
				graph[index][i] = '#';
				i++;
			}
			graph[index][bestFitGaussAry[index]] = '*';
			index++;
		}

		debug << "leaving plotGaussCurves()" << endl;
	}
};

int main(int argc, char* argv[]) {
	ifstream in(argv[1]);
	ofstream out1(argv[2]);
	ofstream out2(argv[3]);
	ofstream debug(argv[4]);

	thresholdSelection t = thresholdSelection(in);
	in >> t.numR >> t.numC >> t.minV >> t.maxV;
	
	t.biThrV = 0;
	t.maxH = 0;
	t.histAry = new int[t.maxV + 1]();
	t.gaussAry = new int[t.maxV + 1]();
	t.bestFitGaussAry = new int[t.maxV + 1]();
	t.setZero(t.histAry);
	t.setZero(t.gaussAry);
	int h = t.loadHist(in);
	t.histH = h;
	t.graph = new char* [(t.maxV + 1)]();

	for (int i = 0; i < t.maxV + 1; i++) {
		t.graph[i] = new char[t.histH + 1];
	}

	for (int i = 0; i < t.maxV + 1; i++) {
		for (int j = 0; j < t.histH + 1; j++) {
			t.graph[i][j] = ' ';
		}
	}
	
	out1 <<"in main(), below is the input histogram" << endl;
	t.dispHist(out1);

	t.plotHist();
	debug << "In main(), below is the Graph after plotting the histogram onto Graph" << endl;
	t.printGraph(debug);

	t.biThrV = t.biGaussian(debug);

	out2 << "The BiGaussThrVal is " << t.biThrV << endl;
	debug<<"In main(). Below is the graph showing the histogram, the best fitted Gaussian curves and the gap" << endl;
	t.plotGaussCurves(debug);
	out2 << "In main(), Below is the final Graph" << endl;
	t.printGraph(out2);

	out1.close();
	out2.close();
	debug.close();
	in.close();
}