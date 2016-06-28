#ifndef QAP_H_
#define QAP_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <stdio.h>
#include <math.h>

#include <stdlib.h>
#include <float.h>
#include <limits.h>  //INT_MAX
#include <algorithm>
#include <vector>

using namespace std;

struct Matrix{
	vector<vector<int> > distancia;
	vector<vector<int> > flujo;
};

struct Greedy{
  vector<int> solucion;
  int coste;
};

struct Problema{
	vector<int> solucion;
	int coste;
};

struct Potenciales{
	vector<double> datos;
	int tam;
};

//Funcion para ordenar la matriz de flujo
bool Mayor(int i, int j);

class QAP{
  private:
  	int numD; //dimension del problema
  	int Antz; //numero de hormigas
  	vector<int> solI; //solucion inicial
  	int costeI; //coste de la solucion inicial
  	float media;
  	vector<vector<double> > feromonas; //matriz de feromonas
  
  public:
  	Matrix datosP;
    Greedy solG;
  	Problema solBL;
  	Problema solSCH;
  	Problema solSHMM;
  	Potenciales fluj;
  	Potenciales dist;
  	vector<Problema> ant;
  	vector<Problema> antMM;
  	vector<vector<int> > maskTrans;
  	double feromonaIni;
    //para actualizar feromona SHMM
    float feroMax;
    float feroMin; 

    void clear();
  	void limpiarSolIni(){solI.clear();}
  	
  	//Get y Set
  	int getNumDatos(){return numD;}
  	int getCosteBL(){return solBL.coste;}
    int getCosteG(){return solG.coste;}
  	int getCosteSCH(){return solSCH.coste;}
  	int getCosteSHMM(){return solSHMM.coste;}
  	vector<int> getSolIni(){return solI;}

  	void setNumDatos(ifstream& f){f >> numD;}
  	void setNumAnts(int n){Antz = n;}

  	//Métodos para el funcionamiento de los algoritmos
    void inicializarDatos(ifstream& f);
    int costeSol(vector<int> &sol);
    void solInit();

    void AlgGreedy();
  	void busquedaLocal(Problema solucionI);
  	int generaVecinoSec(vector<int> &solAct, vector<int> &vecino, int i, int pos2);

  	void SCH_BL();
    int TransicionSCH(int u ,int nAnt);
    void actualizacionLocal(int u, int l);
  	void actualizacionGlobalSCH();
    void SHMM_BL();
    int TransicionSHMM(int u, int nAnt);
    void actualizacionGlobalSHMM();
  	void inicioFeromonas(bool tipo, int coste = 0);
  	void TruncarFeromonas();

  	void converInd(vector<double> &vec);
    int factorizaCoste(vector<int> &sol, int &r, int &s);
  	void visualizaV(vector<int> &v, int n);
  	void visualizaM(vector<vector<int> > &m, int f, int c);
  	

  	int rdtsc();
};
#endif