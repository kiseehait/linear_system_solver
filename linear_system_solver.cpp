//
// Linear System Solver version 1.0.A
// Created by Seehait Chockthanyawat
//

#include <iostream>
#include <vector>
#include <string>
#include "fraction.h"

std::vector<std::vector<sic::fraction> > input; // input matrix
std::vector<sic::fraction> output; // output (particular part)
size_t total_row, total_col; // total row, column of the input matrix

size_t solution_type; // 0 = unique solution, 1 = infinitely many solution, 2 = no solution (contradiction)

size_t calculation_mode; // 0 = associate homogeneoous system, 1 = particular system

std::vector<std::vector<sic::fraction> > free_var; // free variables (homogeneous part)
std::vector<size_t> free_var_pos; // index of each free variable
size_t total_free_var; // total free variable

// sort the row, prepare the matrix before doing next reduction
void sort_row(size_t start_row_index, size_t target_col)
{
	std::vector<sic::fraction> fraction_table;
	std::vector<size_t> index_table;

	for (size_t row_pointer = 0; row_pointer < total_row; row_pointer++)
	{
		fraction_table.push_back(input[row_pointer][target_col]);
		index_table.push_back(row_pointer);
	}

	bool sorted = false;
	while (!sorted)
	{
		sorted = true;
		for (size_t row_pointer = start_row_index; row_pointer < total_row - 1; row_pointer++)
		{
			if ((!fraction_table[row_pointer + 1].is_zero() && fraction_table[row_pointer + 1] < fraction_table[row_pointer]) || (fraction_table[row_pointer].is_zero() && !fraction_table[row_pointer + 1].is_zero()))
			{
				std::swap(fraction_table[row_pointer], fraction_table[row_pointer + 1]);
				std::swap(index_table[row_pointer], index_table[row_pointer + 1]);
				sorted = false;
			}
		}
	}
	
	std::vector<std::vector<sic::fraction> > vector_temp(input);

	for (size_t row_pointer = 0; row_pointer < total_row; row_pointer++)
	{
		input[row_pointer] = vector_temp[index_table[row_pointer]];
	}
}

// row operation: -k * row(i) + row(j)
void row_operation(size_t init_row, size_t init_col, size_t target_row)
{
	sic::fraction factor(input[target_row][init_col] / input[init_row][init_col]);

	for (size_t col_pointer = 0; col_pointer < total_col; col_pointer++)
	{
		input[target_row][col_pointer] = input[target_row][col_pointer] - factor * input[init_row][col_pointer];
	}
}

// reduction to echelon form
void reduce_row_forward(size_t target_col, size_t start_row)
{
	size_t start_row_index = start_row;
	if (target_col > 0) while (!input[start_row_index][target_col - 1].is_zero() && start_row_index < total_row - 1) start_row_index++;
	if (start_row_index == total_row - 1) return;
	if (input[start_row_index][target_col].is_zero()) return;

	for (size_t row_pointer = start_row_index + 1; row_pointer < total_row; row_pointer++)
	{
		row_operation(start_row_index, target_col, row_pointer);
	}
}

// reduction to reduced echelon form
void reduce_row_backward(size_t target_col)
{
	size_t start_row_index = total_row - 1;
	while (input[start_row_index][target_col].is_zero() && start_row_index > 0) start_row_index--;
	if (start_row_index == 0) return;

	for (size_t row_pointer = start_row_index; row_pointer > 0; row_pointer--)
	{
		row_operation(start_row_index, target_col, row_pointer - 1);
	}
}

// make each leading variable equals to 1
void simplify_row(size_t target_row, size_t target_col)
{
	if (input[target_row][target_col].is_zero()) return;
	
	sic::fraction factor(input[target_row][target_col]);

	for (size_t col_pointer = 0; col_pointer < total_col; col_pointer++)
	{
		input[target_row][col_pointer] /= factor;
	}
}

// check that target column is zero column vector or not
bool is_non_zero_col(size_t target_col)
{
	for (size_t col_pointer = 0; col_pointer < total_row; col_pointer++)
	{
		if (!input[col_pointer][target_col].is_zero()) return true;
	}
	return false;
}

// check that target row is zero row vector or not
bool is_non_zero_row(size_t target_row)
{
	for (size_t row_pointer = 0; row_pointer < total_col; row_pointer++)
	{
		if (!input[target_row][row_pointer].is_zero()) return true;
	}
	return false;
}

// make reduced echelon form matrix
void make_reduced_echelon_form()
{
	size_t limit;
	if (calculation_mode == 0) limit = std::min(total_row, total_col);
	else if (calculation_mode == 1) limit = std::min(total_row, total_col - 1);

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

// calculate the output (particular part)
void calculate_output()
{
	size_t row_pointer = 0;
	size_t col_pointer = 0;
	size_t limit = std::min(total_row, total_col - 1);

	while (col_pointer < limit)
	{
		if (is_non_zero_col(col_pointer))
		{
			if (!input[row_pointer][col_pointer].is_zero()) output[row_pointer] = input[row_pointer][total_col - 1] / input[row_pointer][col_pointer];
			row_pointer++;
		}
		col_pointer++;
	}
}

// check solution type of current linear system
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

// calculate homogeneous part
void calculate_free_var()
{
	total_free_var = 0;
	free_var_pos.clear();

	for (size_t row_pointer = 0; row_pointer < total_row; row_pointer++)
	{
		if (!is_non_zero_row(row_pointer)) total_free_var++;
	}

	total_free_var += total_col - total_row - 1;

	size_t col_pointer = 0;
	while (col_pointer < total_col - 1)
	{
		if (!is_non_zero_col(col_pointer)) free_var_pos.push_back(col_pointer);
		col_pointer++;
	}

	col_pointer = total_col - total_free_var + free_var_pos.size() - 1;
	while (free_var_pos.size() < total_free_var)
	{
		free_var_pos.push_back(col_pointer);
		col_pointer++;
	}

	free_var.resize(output.size());
	for (size_t row_pointer = 0; row_pointer < total_row; row_pointer++)
	{
		free_var[row_pointer].resize(total_free_var);
	}

	for (size_t free_var_pointer = 0; free_var_pointer < total_free_var; free_var_pointer++)
	{
		for (size_t row_pointer = 0; row_pointer < total_row; row_pointer++)
		{
			free_var[row_pointer][free_var_pointer] = sic::fraction(-1, 1) * input[row_pointer][free_var_pos[free_var_pointer]];
		}
	}
}

// check all free variables in the target row is zero or not
bool is_all_free_var_zero(size_t target_row)
{
	for (size_t free_var_pointer = 0; free_var_pointer < total_free_var; free_var_pointer++)
	{
		if (!free_var[target_row][free_var_pointer].is_zero()) return false;
	}
	return true;
}

// print the matrix
void print_input()
{
	if (calculation_mode == 0)
	{
		for (size_t row_pointer = 0; row_pointer < total_row; row_pointer++)
		{
			std::cout << "|\t";
			for (size_t col_pointer = 0; col_pointer < total_col; col_pointer++)
			{
				input[row_pointer][col_pointer].print();
				if (col_pointer + 1 < total_col) std::cout << "\t";
			}
			std::cout << "\t|\n";
		}
	}
	else
	{
		for (size_t row_pointer = 0; row_pointer < total_row; row_pointer++)
		{
			std::cout << "|\t";
			for (size_t col_pointer = 0; col_pointer < total_col; col_pointer++)
			{
				input[row_pointer][col_pointer].print();
				if (col_pointer + 2 < total_col) std::cout << "\t";
				else if (col_pointer + 1 < total_col) std::cout << "\t:\t";
			}
			std::cout << "\t|\n";
		}
	}
}

// print the solution
void print_output()
{
	std::cout << "The solution is: \n\n";
	if (solution_type == 0)
	{
		for (size_t col_pointer = 0; col_pointer < total_col - 1; col_pointer++)
		{
			std::cout << "c" << col_pointer << " = ";
			output[col_pointer].print();
			std::cout << std::endl;
		}
	}
	else if (solution_type == 1)
	{
		size_t free_var_p = free_var_pos[0];
		size_t output_pointer = 0;
		for (size_t col_pointer = 0; col_pointer < total_col - 1; col_pointer++)
		{
			if (col_pointer != free_var_p)
			{
				if (is_all_free_var_zero(output_pointer) && output[output_pointer].is_zero())
				{
					std::cout << "c" << col_pointer << " = " << 0 << std::endl;
				}
				else
				{
					std::cout << "c" << col_pointer << " = ";
					if (!output[output_pointer].is_zero())
					{
						output[output_pointer].print();
						for (size_t free_var_pointer = 0; free_var_pointer < total_free_var; free_var_pointer++)
						{
							if (!free_var[output_pointer][free_var_pointer].is_zero()) 
							{
								std::cout << " + (";
								free_var[output_pointer][free_var_pointer].print();
								std::cout << ")c" << free_var_pos[free_var_pointer];
							}
						}
						std::cout << std::endl;
					}
					else
					{
						size_t free_var_pointer = 0;
						while (free_var[output_pointer][free_var_pointer].is_zero()) free_var_pointer++;
						std::cout << "(";
						free_var[output_pointer][free_var_pointer].print();
						std::cout << ")c" << free_var_pos[free_var_pointer];
						
						free_var_pointer++;
						while (free_var_pointer < total_free_var)
						{
							if (!free_var[output_pointer][free_var_pointer].is_zero())
							{
								std::cout << " + (";
								free_var[output_pointer][free_var_pointer].print();
								std::cout << ")c" << free_var_pos[free_var_pointer];
							}
							free_var_pointer++;
						}
						std::cout << std::endl;
					}
				}
			}
			else
			{
				std::cout << "c" << col_pointer << " = any real number\n";
				free_var_p++;
			}
			output_pointer++;
		}
	}
	else if (solution_type == 2)
	{
		std::cout << "Error: no solution\n";
	}
}

// set title bar of the console
void set_title() {
    char esc_start[] = { 0x1b, ']', '0', ';', 0 };
    char esc_end[] = { 0x07, 0 };
    std::cout << esc_start << "Linear System Solver version 1.0.a" << esc_end;
}

// clear the console screen
void clear_screen() {
    for (int line_pointer = 0; line_pointer < 24; line_pointer++) std::cout << "                                                                                ";
    printf("%c[1J", 0x1B);
    printf("%c[1;1H", 0x1B);
}

// print the instruction
void print_instruction()
{
	clear_screen();
	printf("%c[1;1H", 0x1B);
	std::cout << "                       Linear System Solver version 1.0.a                       ";
	std::cout << "                            by Seehait Chockthanyawat                           ";
	std::cout << std::endl << "--------------------------------------------------------------------------------\n";
	std::cout << "Please enter: \t-h for solve associate homogeneous system\n\t\t-p for solve particular system\n\t\t-e for calculate reduced echelon form matrix\n\t\t-x or other to exit\n";
}

// print the exit instruction
void print_exit_instruction()
{
	std::cout << std::endl << "--------------------------------------------------------------------------------\n";
	std::cout << "Do you want to solve other system(s)? ('Y' or 'y' = Yes, 'N' or 'n' = No): ";
	char input_from_user;
	std::cin >> input_from_user;
	if (input_from_user == 'N' || input_from_user == 'n') exit(0);
}

// get fraction input
sic::fraction get_fraction()
{
	std::string input_from_user;
	std::cin >> input_from_user;

	int top = 0;
	int bottom = 1;
	size_t divide_pointer = 0;
	while (divide_pointer < input_from_user.size())
	{
		if (input_from_user[divide_pointer] == '/') break;
		divide_pointer++;
	}

	top = stoi(input_from_user.substr(0, divide_pointer));
	if (divide_pointer + 1 < input_from_user.size()) bottom = stoi(input_from_user.substr(divide_pointer + 1, input_from_user.size() - divide_pointer - 1));

	return sic::fraction(top, bottom);

}

// get homogeneous system input from the user
void get_homogeneous_system_input()
{
	clear_screen();
	std::cout << "Number of variable(s): ";
	std::cin >> total_col;
	total_col++;

	std::cout << "Number of equation(s): ";
	std::cin >> total_row;

	clear_screen();
	input.resize(total_row);

	std::cout << "Please enter the matrix\n\n";
	std::cout << "For example, the equations is\t 4x + 5y + 6z = 0\n\t\t\t\t 8x -3y = 0\n\n";
	std::cout << "We have an augmented matrix\t|4  5  6 : 0| or |4  5  6|\n\t\t\t\t|8 -3  0 : 0|    |8 -3  0|\n";
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "You can enter input like this:\t4 5 6\n\t\t\t\t8 -3 0\n";
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "Or you can enter input like this: 4 5 6 8 -3 0\n";
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "Input: ";
	for (size_t row_pointer = 0; row_pointer < total_row; row_pointer++)
	{
		input[row_pointer].resize(total_col);
		for (size_t col_pointer = 0; col_pointer < total_col - 1; col_pointer++)
		{
			input[row_pointer][col_pointer] = get_fraction();
		}
		input[row_pointer][total_col - 1] = sic::fraction(0, 1);
	}

	output.resize(total_col);
}

// get particular system input from the user
void get_particular_system_input()
{
	clear_screen();
	std::cout << "Number of variable(s): ";
	std::cin >> total_col;
	total_col++;

	std::cout << "Number of equation(s): ";
	std::cin >> total_row;

	clear_screen();
	input.resize(total_row);

	std::cout << "Please enter the matrix\n\n";
	std::cout << "For example, the equations is\t 4x + 5y + 6z = 7\n\t\t\t\t 8x -3y = 0\n\n";
	std::cout << "We have an augmented matrix\t|4  5  6 : 7|\n\t\t\t\t|8 -3  0 : 0|\n";
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "You can enter input like this:\t4 5 6 7\n\t\t\t\t8 -3 0 0\n";
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "Or you can enter input like this: 4 5 6 7 8 -3 0 0\n";
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "Input: ";
	for (size_t row_pointer = 0; row_pointer < total_row; row_pointer++)
	{
		input[row_pointer].resize(total_col);
		for (size_t col_pointer = 0; col_pointer < total_col; col_pointer++)
		{
			input[row_pointer][col_pointer] = get_fraction();
		}
	}

	output.resize(total_col);
}

// get matrix input from the user
void get_matrix_input()
{
	clear_screen();
	std::cout << "Number of variable(s): ";
	std::cin >> total_col;

	std::cout << "Number of equation(s): ";
	std::cin >> total_row;

	clear_screen();
	input.resize(total_row);

	std::cout << "Please enter the matrix\n\n";
	std::cout << "For example, if the matrix is\t|4  5  6|\n\t\t\t\t|8 -3  0|\n";
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "You can enter input like this:\t4 5 6\n\t\t\t\t8 -3 0\n";
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "Or you can enter input like this: 4 5 6 8 -3 0\n";
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "Input: ";
	for (size_t row_pointer = 0; row_pointer < total_row; row_pointer++)
	{
		input[row_pointer].resize(total_col);
		for (size_t col_pointer = 0; col_pointer < total_col; col_pointer++)
		{
			input[row_pointer][col_pointer] = get_fraction();
		}
	}
}

// solve the linear system
void make_solution()
{
	clear_screen();
	make_reduced_echelon_form();
	calculate_output();
	calculate_free_var();
	check_solution_type();
	std::cout << "The reduced echelon form matrix is:\n\n";
	print_input();
	std::cout << "\n--------------------------------------------------------------------------------\n";
	print_output();
}

// reduce the linear system
void make_reduced_echelon_form_matrix()
{
	clear_screen();
	make_reduced_echelon_form();
	std::cout << "The reduced echelon form matrix is:\n\n";
	print_input();
}

// get selected calculation mode from the user
void get_calculation_mode_from_user()
{
	std::string instruction;
	std::cout << std::endl << "> ";
	std::cin >> instruction;
	if (instruction == "-h")
	{
		calculation_mode = 1;
		get_homogeneous_system_input();
		clear_screen();
		make_solution();
	}
	else if (instruction == "-p")
	{
		calculation_mode = 1;
		get_particular_system_input();
		clear_screen();
		make_solution();
	}
	else if (instruction == "-e")
	{
		calculation_mode = 0;
		get_matrix_input();
		make_reduced_echelon_form_matrix();		
	}
	else exit(0);
}

// main of the program
int main()
{
	set_title();
	while (true)
	{
		print_instruction();
		get_calculation_mode_from_user();
		print_exit_instruction();
	}
	return 0;
}