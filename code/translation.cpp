#include<iostream>
#include<stdio.h>
#include<string>
using namespace std;

const int LINEAS_ENCABEZADO = 7;
const string FIN_BLOQUE = "//";

string traducirSecuencia(string secuencia, int objetivo){
  int size = secuencia.size();
  int valor = (size/objetivo)  * objetivo;
  return secuencia.substr(0,valor);
}

int main(){

  ios_base::sync_with_stdio(false);
  freopen ("blocks.dat","r", stdin);
  //freopen ("pairs.dat","w", stdout);
  freopen ("triples.dat","w", stdout);

  string linea = "";
  
  for(int i=0;i<LINEAS_ENCABEZADO;i++){
    getline(cin,linea);
    cout << linea << endl;
  }

  bool inicioBloque = true;

  while(getline(cin,linea)) {
    // Si es un fin de bloque debemos empezar de nuevo a leer uno nuevo si no es el fin del archivo.
    if(linea == FIN_BLOQUE){
      cout << linea << endl;
      inicioBloque = true;
      continue;
    }

    // Si es una linea en blanco avanzamos
    if(linea.size() == 0){
      cout << linea << endl;
      continue;
    }

    //Leemos el encabezado de un bloque
    if(inicioBloque){
      // Escribimos la primera linea del encabezado del bloque ID
      cout << linea << endl ;

      //Escribimos las lineas restantes del encabezado del bloque
      for(int i=0;i<3;i++){
        getline(cin,linea);
        cout << linea << endl;
      }
      inicioBloque = false;
      continue;
    }

    string descripcion = linea.substr(0,28);

    string temp =  linea.substr(28);
    int pos = temp.find(" ");
    string secuencia = temp.substr(0,pos);
    int objetivo = 3;
    secuencia = traducirSecuencia(secuencia,objetivo);
    string numero = temp.substr(pos);
    cout  << descripcion << secuencia << numero <<endl;
  }
}