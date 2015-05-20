#include <iostream>
#include "fraction.h"

using namespace std;

int main()
{
	sic::fraction f1(1, 3);
	sic::fraction f2(2, 5);
	f1.print();
	f2.print();
	sic::fraction f3;

	f3 = f1 + f2;
	f3.print();
	(f1 + f2).print();

	f3 = f1 - f2;
	f3.print();
	(f1 - f2).print();

	f3 = f1 * f2;
	f3.print();
	(f1 * f2).print();

	f3 = f1 / f2;
	f3.print();
	(f1 / f2).print();
	return 0;
}