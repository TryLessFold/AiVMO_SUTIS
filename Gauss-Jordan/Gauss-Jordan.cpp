#include <iostream>
#include <vector>
#include <iomanip>
#include <string> 
#include <fstream>

using namespace std;

class Fraction {
private:
	long long gcd(long long a, long long b) {
		while (a != b) {
			if (a > b) {
				a -= b;
			}
			else {
				b -= a;
			}
		}
		return a;
	}

public:
	long long numerator, denominator;

	Fraction() {
		numerator = 0;
		denominator = 1;
	}

	Fraction(long long n, long long d) {
		if (d == 0) {
			cerr << "Denominator may not be 0." << endl;
			exit(0);
		}
		else if (n == 0) {
			numerator = 0;
			denominator = 1;
		}
		else {
			int sign = 1;
			if (n < 0) {
				sign *= -1;
				n *= -1;
			}
			if (d < 0) {
				sign *= -1;
				d *= -1;
			}

			long long tmp = gcd(n, d);
			numerator = n / tmp * sign;
			denominator = d / tmp;
		}
	}

	operator int() { return (numerator) / denominator; }
	operator float() { return ((float)numerator) / denominator; }
	operator double() { return ((double)numerator) / denominator; }
};

Fraction operator+(const Fraction& lhs, const Fraction& rhs) {
	Fraction tmp(lhs.numerator * rhs.denominator
		+ rhs.numerator * lhs.denominator,
		lhs.denominator * rhs.denominator);
	return tmp;
}

Fraction operator+=(Fraction& lhs, const Fraction& rhs) {
	Fraction tmp(lhs.numerator * rhs.denominator
		+ rhs.numerator * lhs.denominator,
		lhs.denominator * rhs.denominator);
	lhs = tmp;
	return lhs;
}

Fraction operator-(const Fraction& lhs, const Fraction& rhs) {
	Fraction tmp(lhs.numerator * rhs.denominator
		- rhs.numerator * lhs.denominator,
		lhs.denominator * rhs.denominator);
	return tmp;
}

Fraction operator-=(Fraction& lhs, const Fraction& rhs) {
	Fraction tmp(lhs.numerator * rhs.denominator
		- rhs.numerator * lhs.denominator,
		lhs.denominator * rhs.denominator);
	lhs = tmp;
	return lhs;
}

Fraction operator*(const Fraction& lhs, const Fraction& rhs) {
	Fraction tmp(lhs.numerator * rhs.numerator,
		lhs.denominator * rhs.denominator);
	return tmp;
}

Fraction operator*=(Fraction& lhs, const Fraction& rhs) {
	Fraction tmp(lhs.numerator * rhs.numerator,
		lhs.denominator * rhs.denominator);
	lhs = tmp;
	return lhs;
}

Fraction operator*(int lhs, const Fraction& rhs) {
	Fraction tmp(lhs * rhs.numerator, rhs.denominator);
	return tmp;
}

Fraction operator*(const Fraction& rhs, int lhs) {
	Fraction tmp(lhs * rhs.numerator, rhs.denominator);
	return tmp;
}

Fraction operator/(const Fraction& lhs, const Fraction& rhs) {
	Fraction tmp(lhs.numerator * rhs.denominator,
		lhs.denominator * rhs.numerator);
	return tmp;
}

std::ostream& operator<<(std::ostream& strm, const Fraction& a) {
	string tmp;
	if (a.denominator == 1) {
		tmp = to_string(a.numerator);
		strm << tmp;
	}
	else {
		tmp = to_string(a.numerator) + "/" + to_string(a.denominator);
		strm << tmp;
	}
	return strm;
}

void PrintMatrix(vector <Fraction*> matrix, int n, int k);
void GaussJordan(vector <Fraction*> matrix, int n, int k);

void jopa_medvedya(vector <Fraction*> cp, int* a, int n, int k)
{
	// 2 4
	k++;
	vector <Fraction*> m;
	for (int i = 0; i < n * k; i++)
	{
		m.push_back(new Fraction(cp[i]->numerator, cp[i]->denominator));
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Fraction* tp = m[i* k + j];
			m[i * k + j] = m[i * k + a[j]];
			m[i * k + a[j]] = tp;
		}
	}

	GaussJordan(m, n, k);
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Fraction* tp = m[i*k + n - j - 1];
			m[i *k + n - j - 1] = m[i*k + a[n-j-1]];
			m[i * k + a[n-j-1]] = tp;
		}
	}
	
	PrintMatrix(m, n, k);
	vector <Fraction*> b;
	int jopka = 0;
	for (int i = 0; i < k - 1; i++) 
	{
		if (i == a[jopka]) {
			b.push_back(m[jopka++ * k + k - 1]);
		}
		else {
			b.push_back(new Fraction(0, 1));
		}
	}
	for (int i = 0; i < k - 1; i++)
	{
		cout << *b[i] << " ";
	}
	int flag = 0;
	for (int i = 0; i < k - 1; i++) {
		if (b[i]->numerator < 0) flag = 1;
	}
	if (!flag) {
		cout << endl << "opornoe" << endl;
		cout << "------------------------------------------------" << endl << endl;
	}
	else
	{
		cout << endl << endl;
		cout << "------------------------------------------------" << endl << endl;
	}
 	cout << endl;

	
}

void vizualka(vector <Fraction*> m, int n, int k)
{
	// 4 2
	int *A, i, p;
	A = new int[k];
	for (i = 0; i < k; i++) A[i] = i;
	while (1) {
		for (i = 0; i < k; i++) cout << A[i] + 1 << " ";
		cout << endl;
		//vector<Fraction*> cp;
		// copy(m.begin(), m.end(), cp);
		jopa_medvedya(m, A, k, n);
		// 2 4

		cout << endl;
		if (A[k - 1] < n - 1) A[k - 1]++;
		else {
			for (p = k - 1; p > 0; p--)
				if (A[p] - A[p - 1] > 1) break;
			if (p == 0) break;
			A[p - 1] ++;
			for (i = p; i < k; i++) A[i] = A[i - 1] + 1;
		}
	}
}

void GaussJordan(vector <Fraction*> matrix, int n, int k)
{
	int x = 0, y = 0, f = 0;
	Fraction null_frac(0, 1);
	for (int t = 0; t < n; t += 1)
	{
		if (matrix[y * k + x]->numerator == 0)
		{
			if (t == n - 1)
				break;
			for (int i = y + 1; i < n; i++)
			{
				if (matrix[i*k+x]->numerator>0)
				{ 
					f = 1;
					for (int j = x; j < k; j++)
					{
						*matrix[y * k + j] += *matrix[i * k + j];
					}
				}
			}
			if (f == 0)
			{
				x++;
				y++;
				continue;
			}
		}
		for (int i = 0; i < n; i++)
		{
			if (y == i) continue;
			for (int j = x + 1; j < k; j++)
			{
				*matrix[i * k + j] = ((*matrix[y*k + x])*(*matrix[i*k + j]) - (*matrix[i * k + x])*(*matrix[y * k + j])) / *matrix[y * k + x];
			}
		}
		for (int i = 0; i < n; i++)
		{
			if (y == i) continue;
			*matrix[i * k + x] = null_frac;
			//cout << x << " !" << i << " !";
		}
		x++;
		y++;
		//cout << t << ")" << endl;
		//PrintMatrix(matrix, n, k);
		//cout << endl;
	}
	x = 0, y = 0;
	for (int i = 0; i < n; i++)
	{
		if (matrix[y * k + x]->numerator == 0)
		{
			x++;
			y++;
			continue;
		}
		for (int j = k - 1; j >= x; j--)
		{
			*matrix[i * k + j] = (*matrix[i * k + j]) / (*matrix[y * k + x]);
		}
		x++;
		y++;
	}
}

void PrintMatrix(vector <Fraction*> matrix, int n, int k)
{
	int x = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < k; j++)
		{
			cout << setw(10) << *matrix[(long)i *k + j];
		}
		cout << endl;
	}
}

int main() {

	int k = 1, n = 1, tmp0;
	vector <Fraction*> matrix;
	ifstream in("matrix.txt");

	in >> n >> k;
	k++;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < k; j++)
		{
			in >> tmp0;
			matrix.push_back(new Fraction(tmp0, 1));
		}
	}
	cout << "input matrix:" << endl;
	PrintMatrix(matrix, n, k);
	
	GaussJordan(matrix, n, k);
	for (int i = 0; i < n; i++)
	{
		int f = 0;
		for (int j = 0; j < k; j++) {
			if (matrix[i * k + j]->numerator != 0) {
				f = 1;
			}
		}
		if (f == 0) {
			matrix.erase(matrix.begin() + i * k, matrix.begin() + i * k + k);
			n--;
		}
	}


	//cout << "Jordan Gauss:" << endl;
	//PrintMatrix(matrix, n, k);
	// n = 2 k = 5
	


	vizualka(matrix, k - 1, n);

	return 0;
}
