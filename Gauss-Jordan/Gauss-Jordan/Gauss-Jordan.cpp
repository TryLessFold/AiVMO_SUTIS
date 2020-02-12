#include <iostream>
#include <vector>
#include <iomanip>
#include <string> 

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
        tmp = to_string(a.numerator)+"/"+ to_string(a.denominator);
        strm << tmp;
    }
    return strm;
}

void PrintMatrix(vector <Fraction*> matrix, int n, int k);

void GaussJordan(vector <Fraction*> matrix, int n, int k)
{
    cout << "-----GAUSS-JORDAN-----" << endl;
    int x = 0, y = 0;
    Fraction null_frac(0,1);
    for (int t = 0; t < n; t += 1)
    {
        if (matrix[y * k + x]->numerator == 0)
        {
            x++;
            y++;
            continue;
        }
        for (int i = 0; i < n; i++)
        {
            if (y == i) continue;
            for (int j = x+1; j < k; j++)
            {
                //cout << x << " " << y << " " << i << " " << j << " "; 
                *matrix[i * k + j] = ((*matrix[y*k + x])*(*matrix[i*k + j]) - (*matrix[i * k + x])*(*matrix[y * k + j])) / *matrix[y * k + x];
                // cout << *matrix[i * k + j] << endl;
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
        PrintMatrix(matrix, n, k);
        cout << endl;
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
        for (int j = k-1; j >= x; j--)
        {
            *matrix[i * k + j] = (*matrix[i * k + j]) / (*matrix[y * k + x]);
        }
        x++;
        y++;
    }
    cout << "-----END_GAUSS-JORDAN-----" << endl;
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

    int k = 1, n = 1, tmp0, tmp1;
    vector <Fraction*> matrix;

    while ((n > 0)&&(k > 0))
    {
        cout << "Input equation(s) number:";
        cin >> n;
        cout << "Input coefficient(s) number: ";
        cin >> k;
        k++;

        for (int i = 0; i < n; i++)
        {
            cout << "Input equation[" << i << "] coefficient(s):" << endl;
            for (int j = 0; j < k; j++)
            {
                if (j == k-1)
                    cout << "Equally:";
                else
                    cout << "Input x[" << j << "]:";
                cin >> tmp0;
                matrix.push_back(new Fraction(tmp0, 1));
            }
        }
        PrintMatrix(matrix, n, k);
        GaussJordan(matrix, n, k);
        PrintMatrix(matrix, n, k);
        matrix.clear();
    }

    return 0;
}