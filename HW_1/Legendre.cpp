//
// Created by irabiel on 8/30/2022.
//

#include <iostream>
#include <vector>
using namespace std;

double pow(double x, int n);
double Legendre(double x, int n );
void linspace(double a, double b, int N, double ** x);
void sampledLegendre(double a, double b, int N, int n, double * p);
vector<double> linspace(double a, double b, int N);
vector<double> sampledLegendre(double a, double b, int N, int n);

int main( )
{
    double a = 0;
    double b = 1;
    int N = 11;
    int n = 2;
    vector<double> p;

    double * p1 = new double[ N ];

    sampledLegendre( a,  b,  N,  n, p1);

    p = sampledLegendre( a,  b,  N,  n);

    for (int i = 0; i < N ; i++)
    {
        cout << p[i] << endl;
        //cout << p[i] << ' ' << p1[i] << ' ' << p[i] - p1[i] << endl;
    }

    delete[] p1;
    return 0;
}


double pow (double x, int n)
{
    double prod = x;
    for (int i = 0; i<n-1 ; i++){
        prod *= x;
    }
    return prod;
}

double Legendre (double x, int n )
{
    if (n == 0){return 1;}
    else if (n == 1){return x;}
    else if (n == 2){return .5 * (3*x*x - 1);}
    else if (n == 3){return .5 * (5*pow(x,3)  - 3*x);}
    else if (n == 4){return (35*pow(x,4) - 30*x*x + 3) / 8;}
    else if (n == 5){return (63 * pow(x,5) - 70*pow(x,3) + 15*x)/8;}
    else{
        cout << "Error n is an integer between 0 and 5" << endl;
        exit(1) ;
    }
}

vector<double> linspace(double a, double b, int N) {
    vector<double> x;
    for (int i = 0; i < N; i++) {
        x.push_back(i);
        x[i] = a + (b - a) * (i / (N - 1.0));
    }
    return x;
}

vector<double> sampledLegendre(double a, double b, int N, int n)
{
    vector<double> x;
    vector<double> p;
    x = linspace(a, b,N);
    for (int i = 0 ; i < N; i++)
    {
        p.push_back(i);
        p[i] = Legendre(x[i], n );
    }

    return p;
}

void linspace(double a, double b, int N, double ** x) {
    for (int i = 0; i < N; i++) {
        x[0][i] = a + (b - a) * (i / (N - 1.0));
    }
}

void sampledLegendre(double a, double b, int N, int n, double * p)
{
    double * x = new double[ N ];
    linspace(a, b,N, &x);
    for (int i = 0 ; i < N; i++)
    {
        p[i] = Legendre (x[i], n );
    }

    delete[] x;
}

