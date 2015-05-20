#ifndef SIC_FRACTION_INCLUDED
#define SIC_FRACTION_INCLUDED

#include <iostream>
#include <algorithm>

namespace sic
{

class fraction
{
protected:
	// the representation of our fraction is top/bottom
	int top;
	int bottom;

	// simplify make the simplest fraction
	void simplify()
	{
		int factor = std::__gcd(top, bottom);
		if (factor != 0)
		{
			top /= factor;
			bottom /= factor;
		}
		if (bottom < 0)
		{
			top *= -1;
			bottom *= -1;
		}
	}



public:
	// default constructor
	fraction() : top(0), bottom(1) { }

	// custom constructor
	fraction(const int t, const int b) : top(t), bottom(b) 
	{
		simplify();
	}

	// copy constructor
	fraction(const fraction& other) : top(other.top), bottom(other.bottom) { }

	// destructor
	~fraction() { }

	// operator overload
	fraction& operator=(const fraction& other)
	{
		top = other.top;
		bottom = other.bottom;
		return *this;
	}

	bool operator==(const fraction& other)
	{
		return other.top == top && other.bottom == bottom;
	}

	bool operator!=(const fraction& other)
	{
		return other.top != top || other.bottom != bottom;
	}

	fraction& operator+=(const fraction& other)
	{
		int other_top = other.top * bottom;
		top *= other.bottom;
		top += other_top;
		bottom *= other.bottom;
		simplify();
		return *this;
	}

	fraction operator+(const fraction& other) const
	{
		int other_top = other.top * bottom;
		int this_top = top * other.bottom;
		return fraction(other_top + this_top, bottom * other.bottom);
	}

	fraction& operator-=(const fraction& other)
	{
		int other_top = other.top * bottom;
		top *= other.bottom;
		top -= other_top;
		bottom *= other.bottom;
		simplify();
		return *this;
	}

	fraction operator-(const fraction& other) const
	{
		int other_top = other.top * bottom;
		int this_top = top * other.bottom;
		return fraction(this_top - other_top, bottom * other.bottom);
	}

	fraction& operator*=(const fraction& other)
	{
		top *= other.top;
		bottom *= other.bottom;
		simplify();
		return *this;
	}

	fraction operator*(const fraction& other) const
	{
		return fraction(top * other.top, bottom * other.bottom);
	}

	fraction& operator/=(const fraction& other)
	{
		top *= other.bottom;
		bottom *= other.top;
		simplify();
		return *this;
	}

	fraction operator/(const fraction& other) const
	{
		return fraction(top * other.bottom, bottom * other.top);
	}

	// modifier
	void set_top(const int& t)
	{
		top = t;
		simplify();
	}

	void set_bottom(const int& b)
	{
		bottom = b;
		simplify();
	}

	// access
	float get_float_value()
	{
		if (bottom == 0) return 0;
		return (float) top / bottom;
	}

	void print()
	{
		std::cout << top << "/" << bottom << std::endl;
	}
};

}

#endif