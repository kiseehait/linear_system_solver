#include <iostream>
#include <vector>
#include "fraction.h"

using namespace std;

vector<vector<sic::fraction> > input;
vector<sic::fraction> output;
size_t total_row, total_col, solution_type;

void sort_row(size_t start_index, size_t target_col)
{
	vector<sic::fraction> fraction_table;
	vector<size_t> index_table;
	size_t table_size = total_row;

	for (size_t i = 0; i < table_size; i++)
	{
		fraction_table.push_back(input[i][target_col]);
		index_table.push_back(i);
	}

	bool sorted = false;
	while (!sorted)
	{
		sorted = true;
		for (size_t i = start_index; i < table_size - 1; i++)
		{
			if ((!fraction_table[i + 1].is_zero() && fraction_table[i + 1] < fraction_table[i]) || (fraction_table[i].is_zero() && !fraction_table[i + 1].is_zero()))
			{
				swap(fraction_table[i], fraction_table[i + 1]);
				swap(index_table[i], index_table[i + 1]);
				sorted = false;
			}
		}
	}
	
	vector<vector<sic::fraction> > vector_temp(input);

	for (size_t i = 0; i < table_size; i++)
	{
		input[i] = vector_temp[index_table[i]];
	}
}

void row_operation(size_t init_row, size_t init_col, size_t target_row)
{
	sic::fraction factor(input[target_row][init_col] / input[init_row][init_col]);

	for (size_t i = 0; i < total_col; i++)
	{
		input[target_row][i] = input[target_row][i] - factor * input[init_row][i];
	}
}

void reduce_row_forward(size_t target_col, size_t start_row)
{
	size_t start_index = start_row;
	if (target_col > 0) while (!input[start_index][target_col - 1].is_zero() && start_index < total_row - 1) start_index++;
	if (start_index == total_row - 1) return;

	for (size_t i = start_index + 1; i < total_row; i++)
	{
		row_operation(start_index, target_col, i);
	}
}

void reduce_row_backward(size_t target_col)
{
	size_t start_index = total_row - 1;
	while (input[start_index][target_col].is_zero() && start_index > 0) start_index--;
	if (start_index == 0) return;

	for (size_t i = start_index; i > 0; i--)
	{
		row_operation(start_index, target_col, i - 1);
	}
}

void simplify_row(size_t target_row, size_t target_col)
{
	if (input[target_row][target_col].is_zero()) return;
	
	sic::fraction factor(input[target_row][target_col]);

	for (size_t i = 0; i < total_col; i++)
	{
		input[target_row][i] /= factor;
	}
}

bool is_non_zero_col(size_t target_col)
{
	for (size_t i = 0; i < total_row; i++)
	{
		if (!input[i][target_col].is_zero()) return true;
	}
	return false;
}

void make_reduced_echelon_form(size_t mode)
{
	size_t limit;
	if (mode == 0) limit = std::min(total_row, total_col);
	else if (mode == 1) limit = std::min(total_row, total_col - 1);

	size_t col_pointer = 0;
	for (size_t row_pointer = 0; row_pointer < limit; row_pointer++)
	{
		if (is_non_zero_col(row_pointer))
		{
			sort_row(col_pointer, row_pointer);
			reduce_row_forward(row_pointer, col_pointer);
			col_pointer++;
		}
	}

	for (size_t row_pointer = limit; row_pointer > 0; row_pointer--)
	{
		reduce_row_backward(row_pointer - 1);
	}

	col_pointer = 0;
	for (size_t row_pointer = 0; row_pointer < limit; row_pointer++)
	{
		if (is_non_zero_col(row_pointer))
		{
			simplify_row(col_pointer, row_pointer);
			col_pointer++;
		}
	}
}

void calculate_output()
{
	size_t row_pointer = 0;
	size_t col_pointer = 0;
	size_t free_val_pointer = 0;
	size_t limit = std::min(total_row, total_col - 1);

	while (col_pointer < total_col - 1)
	{
		if (is_non_zero_col(col_pointer))
		{
			if (!input[row_pointer][col_pointer].is_zero()) output[row_pointer] = input[row_pointer][total_col - 1] / input[row_pointer][col_pointer];
			row_pointer++;
		}
		col_pointer++;
	}
}

void check_solution_type()
{
	solution_type = 0;
	if (total_row < total_col - 1) solution_type = 1;
	else if (total_row > total_col - 1)
	{
		for (size_t row_pointer = total_col - 1; row_pointer < total_row; row_pointer++)
		{
			if (!input[row_pointer][total_col - 1].is_zero()) solution_type = 2;
		}
	}
	else
	{
		for (size_t pointer = 0; pointer < total_row; pointer++)
		{
			if (input[pointer][pointer].is_zero()) solution_type = 1;
		}
	}
}

void get_input()
{
	cin >> total_row >> total_col;
	int input_temp;

	input.resize(total_row);
	for (size_t i = 0; i < total_row; i++)
	{
		input[i].resize(total_col);
		for (size_t j = 0; j < total_col; j++)
		{
			cin >> input_temp;
			input[i][j] = sic::fraction(input_temp, 1);
		}
	}

	output.resize(total_col);
	solution_type = 0;
}

void print_input()
{
	for (size_t i = 0; i < total_row; i++)
	{
		for (size_t j = 0; j < total_col; j++)
		{
			input[i][j].print();
			if (j + 1 < total_col) cout << "\t";
		}
		cout << endl;
	}
}

void print_output()
{
	for (size_t i = 0; i < total_col; i++)
	{
		cout << "x" << i << " = ";
		output[i].print();
		cout << endl;
	}
}

int main()
{
	get_input();
	cout << "------" << endl;
	cout << "step1" << endl;
	make_reduced_echelon_form(1);
	print_input();
	cout << "step2" << endl;
	calculate_output();
	print_output();
	cout << "step3" << endl;
	check_solution_type();
	cout << "Solution type: " << solution_type << endl;
	return 0;
}