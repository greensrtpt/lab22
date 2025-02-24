#include<iostream>
#include<cmath>
using namespace std;

class ComplexNumber{				
	public:
		double real;
		double imag;
		explicit ComplexNumber(double, double);
		ComplexNumber operator+(const ComplexNumber &);
		ComplexNumber operator-(const ComplexNumber &);
		ComplexNumber operator*(const ComplexNumber &);
		ComplexNumber operator/(const ComplexNumber &);
		bool operator==(const ComplexNumber &);
		double abs();
		double angle();

		ComplexNumber operator+(double) const;
		ComplexNumber operator-(double) const;
		ComplexNumber operator*(double) const;
		ComplexNumber operator/(double) const;

		friend bool operator==(const ComplexNumber &, double);
		friend bool operator==(double, const ComplexNumber &);
		friend ostream &operator<<(ostream &, const ComplexNumber &);
};

ComplexNumber::ComplexNumber(double x = 0, double y = 0){
	real = x; imag = y;
}

ComplexNumber ComplexNumber::operator+(const ComplexNumber &c){
	return ComplexNumber(real + c.real, imag + c.imag);
}

ComplexNumber ComplexNumber::operator-(const ComplexNumber &c){
	return ComplexNumber(real - c.real, imag - c.imag);
}

ComplexNumber ComplexNumber::operator*(const ComplexNumber &c){
	return ComplexNumber(real * c.real - imag * c.imag, real * c.imag + imag * c.real);
}

ComplexNumber ComplexNumber::operator/(const ComplexNumber &c){
	double denominator = c.real * c.real + c.imag * c.imag;
	return ComplexNumber((real * c.real + imag * c.imag) / denominator,
                         (imag * c.real - real * c.imag) / denominator);
}

double ComplexNumber::abs(){
	return sqrt(real * real + imag * imag);
}

double ComplexNumber::angle(){
	return atan2(imag, real) * (180.0 / M_PI); 
}

bool ComplexNumber::operator==(const ComplexNumber &c){
	return (fabs(real - c.real) < 1e-9) && (fabs(imag - c.imag) < 1e-9);
}

ComplexNumber ComplexNumber::operator+(double d) const {
	return ComplexNumber(real + d, imag);
}

ComplexNumber ComplexNumber::operator-(double d) const {
	return ComplexNumber(real - d, imag);
}

ComplexNumber ComplexNumber::operator*(double d) const {
	return ComplexNumber(real * d, imag * d);
}

ComplexNumber ComplexNumber::operator/(double d) const {
	return ComplexNumber(real / d, imag / d);
}

ComplexNumber operator+(double d, const ComplexNumber &c){
	return ComplexNumber(d + c.real, c.imag);
}

ComplexNumber operator-(double d, const ComplexNumber &c){
	return ComplexNumber(d - c.real, -c.imag);
}

ComplexNumber operator*(double d, const ComplexNumber &c){
	return ComplexNumber(d * c.real, d * c.imag);
}

ComplexNumber operator/(double d, const ComplexNumber &c){
	double denominator = c.real * c.real + c.imag * c.imag;
	return ComplexNumber((d * c.real) / denominator, (-d * c.imag) / denominator);
}

bool operator==(const ComplexNumber &c, double d){
	return (fabs(c.real - d) < 1e-9) && (fabs(c.imag) < 1e-9);
}

bool operator==(double d, const ComplexNumber &c){
	return (fabs(c.real - d) < 1e-9) && (fabs(c.imag) < 1e-9);
}

ostream &operator<<(ostream &os, const ComplexNumber &c){
	if (c.imag == 0) os << c.real;
	else if (c.real == 0) os << c.imag << "i";
	else os << c.real << (c.imag > 0 ? "+" : "") << c.imag << "i";
	return os;
}