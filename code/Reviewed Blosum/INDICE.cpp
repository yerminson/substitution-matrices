#include <iostream>
using namespace std;

int INDEX(int n,int col,int row){
 return col*n - (col*(col+3))/2 - 1 + row ; 
}
int main() {
    int n = 16600;
    for(int col = 0; col < n-1; col++){
        
        for(int row = col+1; row < n;row++)
            cout <<"INDEX("<<n<<","<<col<<","<<row<<") =" << INDEX(n,col,row) << endl;  
    }
    

    
    
    return 0;
}