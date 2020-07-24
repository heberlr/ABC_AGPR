#ifndef CRS
#define CRS

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

void create_crs_ordened(int row, int column, int row_ptr[], int col_ind[], double val[], double value);
void add_crs(int row, int column, int row_ptr[], int col_ind[], double val[], double new_value);
void set_crs(int row, int column, int row_ptr[], int col_ind[], double val[], double new_value);
void print_crs(int row_ptr[], int col_ind[], double val[], int num_nos, int nz_num);
void zera_crs(int num_nos, int nz_num, int row_ptr[], int col_ind[], double val[]);

#endif /* end of include guard: CRS */
