#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <new>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include <vector>


using namespace std;

int main ( int argc, char *argv[] ){
	char arquivo[60];
	int numarquivos = atoi(argv[1]), execution = atoi(argv[2]);
	vector<vector<double> > Lmatrix;
	vector<vector<double> > Dmatrix;
	vector<int> temp;

	for(int k = 1; k <= numarquivos; k++){
		sprintf(arquivo,"../confluence/confluence%03d/out%03d-confluence-2D.dat",execution,k);
		//printf("+++++++++++++ %s +++++++++++++\n",arquivo);
		ifstream input(arquivo);
		vector<double> Lvector;
		vector<double> Dvector;
		while (!input.eof()){
			double live_confluence,dead_confluence; int tempo;
			input >> tempo >> live_confluence >> dead_confluence;
			if(!input.eof()){
				if (k == 1) temp.push_back(tempo);
				Lvector.push_back(live_confluence);
				Dvector.push_back(dead_confluence);
			}
		}
		Lmatrix.push_back(Lvector);
		Dmatrix.push_back(Dvector);
		Lvector.clear();
		Dvector.clear();
		input.close();
	}
	double Lmedia[temp.size()],Dmedia[temp.size()];//,media2[tempmax+1];
	for(int j = 0; j < Lmatrix[0].size(); j++){
		Lmedia[j] = 0;
		Dmedia[j] = 0;
	}
	for(int j = 0; j < temp.size(); j++){
		for(int k = 1; k <= numarquivos; k++){
			Lmedia[j] += Lmatrix[k-1][j];
			Dmedia[j] += Dmatrix[k-1][j];
		}
		Lmedia[j] = Lmedia[j]/(double)numarquivos;
		Dmedia[j] = Dmedia[j]/(double)numarquivos;
	}
	sprintf(arquivo,"confluence-%04d.dat",execution);
	ofstream output(arquivo);
	for(int j = 0; j < temp.size(); j++){
		double L_SquareSum=0,D_SquareSum=0;
		for(int k = 1; k <= numarquivos; k++){
			L_SquareSum += pow(Lmatrix[k-1][j]-Lmedia[j],2);
			D_SquareSum += pow(Dmatrix[k-1][j]-Dmedia[j],2);
		}
		output << std::setw(3) << temp[j] << "\t" << scientific << std::setw(12) << Lmedia[j] << "\t" << sqrt(L_SquareSum/(double)numarquivos) << "\t" << Dmedia[j] << "\t" << sqrt(D_SquareSum/(double)numarquivos) << endl;
	}
	output.close();
	return 0;
}
