#include "QAP.h"
#include "random_ppio.h"

using namespace std;

#define MAX 999999;

bool Mayor(int i, int j){
	return (i > j);
}

/*			Algoritmo Greedy			*/
void QAP::AlgGreedy(){
	vector<int> sumf(numD);
	vector<int> sumd(numD);
	vector<int> maskF(numD, 0);
	vector<int> maskD(numD, 0);
	vector<int> solGr(numD, -1);
	int cont = 0, maxf = 0, indf = 0, indd = 0, f = 0, d = 0;
	int costeGr = 0;
	int mind = MAX;

	for(int i = 0; i < numD; i++){
		for(int j = 0; j < numD; j++){
			f += datosP.flujo[i][j];
			d += datosP.distancia[i][j];
		}
		sumf[i] = f;
		sumd[i] = d;
		f = 0;
		d = 0;
	}

	while(cont < numD){
		for(int i = 0; i < numD; i++){
			if(maxf < sumf[i] && maskF[i] == 0){
				maxf = sumf[i];
				indf = i;
			}
			if(mind > sumd[i] && maskD[i] == 0){
				mind = sumd[i];
				indd = i;
			}
		}
		solGr[indf] = indd;
		maskD[indd] = maskF[indf] = 1;
		indf = indd = 0;
		mind = MAX;
		maxf = 0;
		cont++;
	}

	costeGr = costeSol(solGr);
	solG.coste = costeGr;
	solG.solucion = solGr;
}

/*					Liberacion de memoria 					*/
void QAP::clear(){
	limpiarSolIni();

	for(int i = 0; i < numD; i++){
		datosP.distancia[i].clear();
		datosP.flujo[i].clear();
	}

	datosP.distancia.clear();
	datosP.flujo.clear();
	solBL.solucion.clear();
	solSCH.solucion.clear();
	solSHMM.solucion.clear();

	numD =0;
	costeI=0;
	solBL.coste = 0;
	solSCH.coste = 0;
	solSHMM.coste = 0;
}

//Redimensiono y relleno las matrices del problema 
void QAP::inicializarDatos(ifstream& f){

	feromonaIni = pow(10,-6);
	solG.solucion.resize(numD);
	datosP.distancia.resize(numD);
	datosP.flujo.resize(numD);
	feromonas.resize(numD);
	fluj.datos.resize(numD);
	dist.datos.resize(numD);
	ant.resize(Antz);
	antMM.resize(Antz);
	maskTrans.resize(Antz);

	for(int i = 0; i< Antz; i++){
		maskTrans[i].resize(numD, 0);
		ant[i].solucion.resize(numD);
		antMM[i].solucion.resize(numD);
	}
	
	for(int i = 0; i < numD; i++){
		feromonas[i].resize(numD);
		datosP.distancia[i].resize(numD);
		datosP.flujo[i].resize(numD);
		fluj.datos[i] = i;
		dist.datos[i] = i;
	}

	for(int i = 0; i < numD; i++){
		for (int j = 0; j < numD; ++j){
			f >> datosP.distancia[i][j];
			feromonas[i][j] = feromonaIni;
		}
	}

	for(int i = 0; i < numD; i++)
		for (int j = 0; j < numD; ++j)
			f >> datosP.flujo[i][j];

	fluj.tam = numD;
	dist.tam = numD;
}

//Calcula el coste de la solucion
int QAP::costeSol(vector<int> &solA){
	int costeA = 0;
	for(int i = 0; i < numD; i++)
		for(int j = 0; j < numD; j++){
			costeA += (datosP.flujo[i][j]*datosP.distancia[solA[i]][solA[j]]);
		}

	return costeA;
}

//Generador de la solucion inicial
void QAP::solInit(){
	solI.resize(numD);
	costeI = 0;
	int aux = 0, ind = 0;
	vector<int> mask(numD,0);

	for(int i = 0; i < numD; i++)
		solI[i] = -1;

	srand(rdtsc());
	int al = (rand()%50) +1;
	Set_random(al);

	while(ind < numD){
		aux = Randint(0, numD-1);
		if(mask[aux] == 0){
			solI[ind] = aux;
			mask[aux] = 1;
			ind++;
		}
	}
	costeI = costeSol(solI);
}

//Algoritmo Busqueda local
void QAP::busquedaLocal(Problema solucionI){
	solBL.solucion.resize(numD);
	vector<int> vecinoBL(numD);
	vector<int> solActBL(solucionI.solucion);
	int costeActBL = solucionI.coste, costeVecBL = 0, pos = 0, aux;
	
	for(int i = 0; i < numD; i++){
		pos = generaVecinoSec(solActBL, vecinoBL, i, pos);
		costeVecBL = costeActBL + factorizaCoste(vecinoBL, i, pos);
		
		if(costeVecBL < costeActBL){
			costeActBL = costeVecBL;

			for(int i = 0; i < numD; i++)
				solActBL[i] = vecinoBL[i];
			break;
		}
	}	
	solBL.coste = costeActBL;

	for(int i = 0; i < numD; i++)
		solBL.solucion[i] = solActBL[i];

	vecinoBL.clear();
	solActBL.clear();
}

//Operador de vecinos secuencial   
int QAP::generaVecinoSec(vector<int> &solAct, vector<int> &vecinoSec, int i, int pos2){
	int aux = 0, pos = 0, pos1 = 0;

	for(int j = 0; j < numD; j++)
		vecinoSec[j] = solAct[j];
	pos = Randint(0, numD-1);

	if(pos != i && pos != pos2){
		aux = vecinoSec[i];
		vecinoSec[i] = vecinoSec[pos];
		vecinoSec[pos] = aux;
		return pos;
	}else generaVecinoSec(solAct, vecinoSec, i, pos2);
}

/*			Algoritmo de Colonia de Hormigas			*/
void QAP::SCH_BL(){
	solSCH.solucion.resize(numD);
	solSCH.coste = INT_MAX;
	double f = 0, d = 0;
	int loc;
	int iter = 0;// num iteraciones

	for(int i = 0; i < numD; i++){
		for(int j = 0; j < numD; j++){
			f += datosP.flujo[i][j];
			d += datosP.distancia[i][j];
		}
		fluj.datos[i] = f;
		dist.datos[i] = 1.0/d; //Heurística
		f = d = 0;
	}
	vector<double> flujosU(fluj.datos);
	sort(flujosU.begin(), flujosU.end(), Mayor);
	converInd(flujosU); 

	inicioFeromonas(false);
	/*******			Proceso constructivo			*******/
	//Bucle evaluaciones
	while(iter < 25000){
		for(int i = 0; i < numD; i++){
			for(int j = 0; j < Antz; j++){
				loc = TransicionSCH(flujosU[i], j);
				ant[j].solucion[flujosU[i]] = loc;
				actualizacionLocal(flujosU[i],loc); //evaporo en los arcos que pasa la hormiga
			}
		}
		for(int i = 0; i< Antz; i++){
			/*Calcular el coste de la solucion obtenida por las hormigas*/
			ant[i].coste = costeSol(ant[i].solucion);
			iter++;
			busquedaLocal(ant[i]);
			iter++; // incrementa en 1 las evaluaciones

			if(solBL.coste < solSCH.coste){
				solSCH.coste = solBL.coste;
				solSCH.solucion = solBL.solucion;
			}
		}
		//Reinicializo la mascara de feromonas para poder buscar mas soluciones
		for(int i = 0; i< Antz; i++){
			maskTrans[i].clear();
			maskTrans[i].resize(numD,0);
		}
		//evaporo y luego deposito feromona sobre los arcos de la mejor solucion
		actualizacionGlobalSCH();
	}
}

/*			Algoritmo de Hormigas Min-Max		*/
void QAP::SHMM_BL(){ 
	solSHMM.solucion.resize(numD);
	solSHMM.coste = INT_MAX;
	double f = 0, d = 0;
	int loc, iter = 0; // para controlar el número de iteraciones

	solInit(); //Creo solucion aleatoria para determinar la Tmax
	for(int i = 0; i < numD; i++){
		for(int j = 0; j < numD; j++){
			f += datosP.flujo[i][j];
			d += datosP.distancia[i][j];
		}
		fluj.datos[i] = f;
		dist.datos[i] = 1.0/d; //Heurística
		f = d = 0;
	}

	vector<double> flujosU(fluj.datos);
	sort(flujosU.begin(), flujosU.end(), Mayor);
	converInd(flujosU); 

	inicioFeromonas(true, costeI);
	/*			Proceso constructivo			*/
	//Bucle evaluaciones
	while(iter < 25000){
		for(int i = 0; i < numD; i++){
			for(int j = 0; j < Antz; j++){
				loc = TransicionSHMM(flujosU[i], j);
				antMM[j].solucion[flujosU[i]] = loc;
			}
		}
		for(int i = 0; i< Antz; i++){
			/*Calcular el coste de la solucion obtenida por las hormigas*/
			antMM[i].coste = costeSol(antMM[i].solucion);
			iter++;
			busquedaLocal(antMM[i]);
			iter++; // incrementa en 1 las evaluaciones
			if(solBL.coste < solSHMM.coste){
				solSHMM.coste = solBL.coste;
				solSHMM.solucion = solBL.solucion;
			}
		}
		//Reinicializo la mascara de feromonas para buscar mas soluciones 
		for(int i = 0; i< Antz; i++){
			maskTrans[i].clear();
			maskTrans[i].resize(numD,0);
		}
		//evaporo y luego deposito feromona sobre los arcos de la mejor solucion
		actualizacionGlobalSHMM();
		//actualizar matriz de feromonas
		TruncarFeromonas(); 
	} 
}

//Actializacion de los valores de la matriz de feromonas
void QAP::TruncarFeromonas(){
	feroMax = 1/(0.2*solSHMM.coste);
	feroMin = feroMax/500;
	feromonaIni = feroMax; 

	for(int i = 0; i < numD; i++)
		for(int j = 0; j < numD; j++){
			if(feromonas[i][j] > feroMax)
				feromonas[i][j] = feroMax;
			if(feromonas[i][j] < feroMin)
				feromonas[i][j] = feroMin;
		}
}

void QAP::converInd(vector<double> &v){
	for(int i=0; i<numD; i++){
		for(int j=0; j<numD; j++){
			if(v[i] == fluj.datos[j]){
				v[i] = j;
				fluj.datos[j] = -1;
			}
		}
	}
}

int QAP::rdtsc(){
	__asm__ __volatile__("rdtsc");
}

//Funcion para visualizar vectores
void QAP::visualizaV(vector<int> &v, int n){
	cout << "Vector: ";
	for(int i = 0; i < n; i++)
		cout << v[i] << " ";
	cout << endl;
}

//Funcion para visualizar matrices
void QAP::visualizaM(vector<vector<int> > &v, int f, int c){
	cout << "Matriz: ";
	for(int i = 0; i < f; i++){
		cout << "º: ";
		for(int j = 0; j < c; j++)
			cout << v[i][j] << " ";
		cout << endl;
	}
}

//Funcion para calcular la factorizacion de coste
int QAP::factorizaCoste(vector<int> &sol, int &r, int &s){
	int costef = 0;

	for(int c = 0; c < numD; c++){
			costef += (datosP.flujo[r][c]*(datosP.distancia[sol[s]][sol[c]] - datosP.distancia[sol[s]][sol[c]])) +
					(datosP.flujo[s][c]*(datosP.distancia[sol[r]][sol[c]] - datosP.distancia[sol[s]][sol[c]])) +
					(datosP.flujo[c][r]*(datosP.distancia[sol[c]][sol[s]] - datosP.distancia[sol[c]][sol[r]])) +
					(datosP.flujo[c][s]*(datosP.distancia[sol[c]][sol[r]] - datosP.distancia[sol[c]][sol[s]]));
	}
	return costef;
}

//Función para inicializar la matriz de feromonas
void QAP::inicioFeromonas(bool tipo, int c){
	
	if(!tipo){ //Algoritmo SCH
		feromonaIni = pow(10, -6);
	}else{ //Algoritmo SHMM
		feroMax = 1 /(0.2 * c);
		feroMin = feroMax / 500.0;
		feromonaIni = feroMax;
	}
	for(int i = 0; i < numD; i++)
		for(int j = 0; j < numD; j++)
			feromonas[i][j] = feromonaIni;
}

int QAP::TransicionSCH(int u, int nAnt){
	int alfa = 1, beta = 2, l, indMejor;
	double aux = 0, mejorProb = 0;
	srand(rdtsc());
	int al = Randint(0, numD-1);
	float q = rand() / (double)RAND_MAX;
	
	if(maskTrans[nAnt][al] == 0){
		indMejor = al;
		if(q <= 0.8){
			for(int i = 0; i < numD; i++){
				aux = pow(feromonas[u][i], alfa) * pow(dist.datos[i], beta);
				if(mejorProb < aux && maskTrans[nAnt][i] == 0){
					mejorProb = aux;
					indMejor = i;
				}
			}
		}
		l = indMejor;
		maskTrans[nAnt][l] = 1;
		return l;
	}else TransicionSCH(u, nAnt);
}

//Evaporo feromona de cada arco por el que pasa una hormiga
void QAP::actualizacionLocal(int u, int l){
	feromonas[u][l] = ((1 - 0.2) * feromonas[u][l]) + (0.2 * feromonaIni);
}

void QAP::actualizacionGlobalSCH(){
	int l;
	for(int u = 0; u < numD; u++){
		l = solSCH.solucion[u];
		feromonas[u][l] = ((1 - 0.2) * feromonas[u][l] + (0.2 * (1.0 / solSCH.coste)));
	}
}

int QAP::TransicionSHMM(int u, int nAnt){
	int alfa = 1, beta = 2, l, indMejor;
	double aux = 0, mejorProb = 0;
	srand(rdtsc());
	int al = Randint(0, numD-1);
	double q = Rand();
	
	if(maskTrans[nAnt][al] == 0){
		indMejor = al;
		if(q > 0.8){
			for(int i = 0; i < numD; i++){
				aux = pow(feromonas[u][i], alfa) * pow(dist.datos[i], beta);
				if(mejorProb < aux && maskTrans[nAnt][i] == 0){
					mejorProb = aux;
					indMejor = i;
				}
			}
		}
		l = indMejor;
		maskTrans[nAnt][l] = 1;
		return l;
	}else TransicionSHMM(u, nAnt);
}

void QAP::actualizacionGlobalSHMM(){
	int l;
	//Evaporo de todos los arcos
	for(int i = 0; i < Antz; i++)
		for(int j = 0; j < numD; j++)
			feromonas[i][j] = 0.2 * feromonas[i][j];
	//Actualizo los arcos de la mejor solucion
	for(int u = 0; u < numD; u++){
		l = solSHMM.solucion[u];
		feromonas[u][l] = ((1 - 0.2) * feromonas[u][l] + (1.0 / solSHMM.coste));
	}
	feroMax = 1 / (0.2 * solSHMM.coste);
	feroMin = feroMax / 500;
}