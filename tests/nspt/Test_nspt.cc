#include <Grid/Grid.h>

using namespace std;

int main(int argc, char *argv[]) {
    
    iPert<double,2> myv, myw;
    
    myv(0)=1.;
    myv(1)=2.;
    
    myw(0)=3.;
    myw(1)=4.;
    cout<<myv<<endl;
    cout<<myw<<endl;
    cout<<myv+myw<<endl;
    cout<<myv+myw<<endl;
    
}
