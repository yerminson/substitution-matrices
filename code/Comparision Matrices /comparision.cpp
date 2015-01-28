#include <iostream>
#include <string.h>
#include <stdlib.h>

using namespace std;

const int COL = 20;
const int ROW = 20;

void printMatrix(int matrix[][COL]){
  for(int i=0;i<ROW;i++){
    for(int j=0;j<=i;j++){
      if(matrix[i][j] < 0)
        cout << matrix[i][j] <<  " ";
      else
        cout << " " << matrix[i][j] << " ";
    }
    cout << endl;
  }
}

int main(){
  int result[ROW][COL];
  memset(result, 0, ROW * COL * sizeof(int));
  int original[ROW][COL];
  memset(original, 0, ROW * COL * sizeof(int));

  for(int i=0;i<ROW;i++){
    for(int j=0;j<=i;j++){
      cin >> result[i][j];
    }
  }

  for(int i=0;i<ROW;i++){
    for(int j=0;j<=i;j++){
      cin >> original[i][j];
    }
  }

  cout << "Result matrix" << endl; 
  printMatrix(result);
  cout << "Original(Base) Matrix" << endl;
  printMatrix(original);

  int comparasionMatrix[ROW][COL];
  memset(comparasionMatrix, 0 , ROW * COL * sizeof(int));

  for(int i=0;i<ROW;i++){
    for(int j=0;j<=i;j++){
      comparasionMatrix[i][j] = abs(original[i][j] - result[i][j]);
    }
  }

  cout << "Comparision Matrix" << endl;
  printMatrix(comparasionMatrix);

  return 0;

}