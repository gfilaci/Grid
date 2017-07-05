#include <Grid/Grid.h>

using namespace std;

int main(int argc, char *argv[]) {
    
    iPert<double,2> v,w;
    v(0)=1.;
    v(1)=2.;
    
    w(0)=3.;
    w(1)=4.;
    cout<<v<<endl;
    cout<<w<<endl;
    cout<<v*w<<endl<<endl;


    Complex im(0,1);
    iPert<iMatrix<Complex,2>,3> S,P;
    
    S(0)(0,0) = 0;
    S(0)(0,1) = 1;
    S(0)(1,0) = -2;
    S(0)(1,1) = -im;
    
    S(1)(0,0) = -23. + im;
    S(1)(0,1) = 0;
    S(1)(1,0) = -8;
    S(1)(1,1) = 1;


    P(1) = conjugate(transpose(S(1)));
    
    S(2) = transpose(S(0)*S(1));
    
    P(0) = S(2)*S(2);
    
    P(2)(0,0) = S(0)(0,0) + 2.;
    P(2)(0,1) = S(0)(0,1);
    P(2)(1,0) = S(0)(1,0);
    P(2)(1,1) = S(0)(1,1) + 2.;
    
    cout<<S*P<<endl;
    
}
/*
compare in Mathematica with
A = {{0, 1}, {-2, -I}};
B = {{-23 + I, 0}, {-8, 1}};
CC = Transpose[A.B];
S0 = A;
S1 = B;
S2 = CC;
P0 = CC.CC;
P1 = ConjugateTranspose[B];
P2 = A + 2*IdentityMatrix[2];
S0.P0 // MatrixForm
S0.P1 + S1.P0 // MatrixForm
S0.P2 + S1.P1 + S2.P0 // MatrixForm
*/
