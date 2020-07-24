#include "CRS.hpp"

using namespace std;

void create_crs_ordened(int row, int column, int row_ptr[], int col_ind[], double val[], double value){
	if (row_ptr[row+1] == 0)
		row_ptr[row+1] = row_ptr[row];
	row_ptr[row+1]++;
	col_ind[row_ptr[row+1]-1] = column;
	val[row_ptr[row+1]-1]= value;
}

void add_crs(int row, int column, int row_ptr[], int col_ind[], double val[], double new_value){
	for(int i = row_ptr[row]; i< row_ptr[row+1]; i++)
		if (column == col_ind[i]) val[i] += new_value;
}

void set_crs(int row, int column, int row_ptr[], int col_ind[], double val[], double new_value){
	for(int i = row_ptr[row]; i< row_ptr[row+1]; i++)
		if (column == col_ind[i]) val[i] = new_value;
}

void print_crs(int row_ptr[], int col_ind[], double val[], int num_nos, int nz_num){
	for (int i=0;i<num_nos+1;i++)
		cout << row_ptr[i] << "\t";
	cout << endl;
	for (int i=0;i<nz_num;i++)
		cout << col_ind[i] << "\t";
	cout << endl;
	for (int i=0;i<nz_num;i++)
		cout << val[i] << "\t";
	cout << endl;
}

void zera_crs(int num_nos, int nz_num, int row_ptr[], int col_ind[], double val[]){
	for (int i =0; i <num_nos+1; i++)
		row_ptr[i] = 0;
	for (int i =0; i < nz_num; i++){
		col_ind[i] = 0;
		val[i] = 0;
	}
}
