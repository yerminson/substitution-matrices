#include <QCoreApplication>
#include <QFile>
#include <iostream>
#include <QString>
#include <QTextStream>

using namespace std;

const int LINEAS_ENCABEZADO = 8;
const QString FIN_BLOQUE("//");
QString generarEspacios(int n){
    QString espacios;
    for(int i=0;i<n;i++){
        espacios.append(" ");
    }
}
int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    // Leemos la ruta del archivo
    string rawRutaArchivo = "";
    cout << "Digite la ruta del archivo: ";
    cin >> rawRutaArchivo;

    //Convertimos el std::string a QStrign
    QString rutaArchivo = QString::fromStdString(rawRutaArchivo);

    //Cargamos el archivo que vamos a leer. Esperamos recibir un archivo blocks.dat o uno con su misma estructura.
    QFile archivo(rutaArchivo);

    //Validamos que el archivo lo podamos leer :D
    if (!archivo.open(QIODevice::ReadOnly | QIODevice::Text))
        return 0;
    //Archivo donde iremos guardando la traduccion a PARES
    QFile salidaPares("pairs.dat");
    //Archivo donde iremos guardando la traduccion a PARES
    QFile salidaTripletas("triplets.dat");

    salidaPares.open(QIODevice::WriteOnly | QIODevice::Text);
    salidaTripletas.open(QIODevice::WriteOnly | QIODevice::Text);

    QTextStream outPares(&salidaPares);
    QTextStream outTripletas(&salidaTripletas);

    QString linea;

    QTextStream in(&archivo);

    // Leemos y escribimos el encabezado que consta de 7 lineas.
    for(int i = 0; i<LINEAS_ENCABEZADO;i++){
        linea = in.readLine();
        outPares << linea ;
        outTripletas << linea ;
    }
    bool inicioBloque = true;
    while (!in.atEnd()) {
        //Leemos la linea
        cin >> rawRutaArchivo;
        linea = in.readLine();

        // Si es un fin de bloque debemos empezar de nuevo a leer uno nuevo si no es el fin del archivo.
        if(linea == FIN_BLOQUE){
            outPares << linea ;
            outTripletas << linea ;
            inicioBloque = true;
            continue;
        }

        // Si es una linea en blanco avanzamos
        if(linea.isEmpty()){
            outPares << endl ;
            outTripletas << endl ;
             continue;
        }



        //Leemos el encabezado de un bloque
        if(inicioBloque){
            // Escribimos la primera linea del encabezado del bloque ID
            outPares << linea ;
            outTripletas << linea ;

            //Escribimos las lineas restantes del encabezado del bloque
            for(int i=0;i<3;i++){
                linea = in.readLine();

                outPares << linea;
                outTripletas << linea;
            }
            inicioBloque = false;
        }

        QStringList datos = linea.split(" ");
        QStringList datosFiltrados;

        foreach (QString valor, datos)
            if(valor.toStdString() != "")
                datosFiltrados.append(valor);

        //Construnccion de la linea traduccida para pares y para tripletas
        QString par;
        QString tripleta;
        QString descripcion = datosFiltrados.at(0);
        if(descripcion.size() < 21)
            descripcion.append(generarEspacios(21 - descripcion.size()));
        par.append(descripcion);
        tripleta.append(descripcion);
        QString parentesis = datosFiltrados.at(1);
        par.append(parentesis);
        tripleta.append(parentesis);
        QString numeroUno = datosFiltrados.at(2);
        if(numeroUno.size() < 5)
            numeroUno = generarEspacios(5 - numeroUno.size()) + numeroUno;
        par.append(numeroUno+" ");
        tripleta.append(numeroUno+" ");
        QString secuencia = datosFiltrados.at(3);
        par.append(secuencia);
        tripleta.append(secuencia);
        QString numeroDos = datosFiltrados.at(4);
        if(numeroDos.size() < 5)
            numeroDos = generarEspacios(5 - numeroDos.size()) + numeroDos;
        par.append(numeroDos);
        tripleta.append(numeroDos);
        outPares << par;
        outTripletas << tripleta;

        cin >> rawRutaArchivo;




    }
    archivo.close();
    salidaPares.close();
    salidaTripletas.close();
    return a.exec();
}
