#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

const int MAXSIZE = 1000000;

vector<vector <int> > matrix;
vector<vector <string> > trace;

int s(string a, string b){
  if ( a == b)
    return 2;
  else
    return -2;
}


int main(){
  string sequenceOne;
  string sequenceTwo;

  cin >> sequenceOne;
  cin >> sequenceTwo;

  int d = 1;

  int s1 = sequenceOne.length() +1;
  int s2 = sequenceTwo.length() +1;
  cout << s1 << endl;
  cout << s2 << endl;

  matrix.resize(s1);
  trace.resize(s1);
  for(int i=0; i<s1;i++){
    matrix[i].resize(s2);
    trace[i].resize(s2);
  }
  matrix[0][0] = 0;
  vector<string> x(s1);
  vector<string>  y(s2);
  for(int i=0;i<s1;i++){
    if(i>0){
      x[i] = sequenceOne.substr(i-1,1);
    }
    matrix[i][0] = -i*d;
  }
  for(int j=0;j<s2;j++){
    if(j>0)
      y[j] = sequenceTwo.substr(j-1,1);
    matrix[0][j] = -j*d;
  }
  for(int i=1;i<s1;i++){
    for(int j=1;j<s2;j++){
      int a = matrix[i-1][j-1] + s(x[i], y[j]);
      int b = matrix[i-1][j] + d;
      int c = matrix[i][j-1] + d;
      matrix[i][j] = max(a,max(b,c));
      cout <<  matrix[i][j] << endl;
      if(matrix[i][j] == a)
        trace[i][j] = "\\";
      else if(matrix[i][j] == b)
        trace[i][j] = "|";
      else
        trace[i][j] = "-";
    }
  }
  cout << matrix[s1-1][s2-1] << endl;

  int is = s1-1;
  int js = s2-1;
  string l1;
  string l2;
  string l3;
  string space = "|";
  string space2 = "_";
  string space3 = " ";
  while(is >0 && js > 0){
    cout << trace[is][js] << endl;
    cout << x[is] <<" "<< y[js] << endl;
      if(trace[is][js] == "\\"){
        l1 = x[is] + l1;
        l2 = space + l2; 
        l3 = y[js] + l3;
        is--;
        js--;
      }else if(trace[is][js] == "-"){
        l1 =  space2 + l1;
        l2 = space3 + l2; 
        l3 = y[js] + l3;
        js--;
      }else if(trace[is][js] == "|"){
        l1 =  x[is] + l1;
        l2 = space3 + l2; 
        l3 = space2 + l3;
        is--;
      }
        

    }
  
    cout <<  l1 << endl;
    cout <<  l2 << endl;
    cout <<  l3 << endl;
      
    cout << endl;
  }

