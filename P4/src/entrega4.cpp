#include "QAP.h"

using namespace std;

int main(int argc, char** argv){
	if(argc < 2){
			exit(0);
	}	
	QAP problema;
	clock_t tI;
	float cps = 1000, tF;
	int nDatos;
	

	/*					Creando flujo de lectura y objeto problema					*/
	ifstream fe;
	fe.open(argv[1]);
	if(fe){
		problema.setNumDatos(fe);
	}else 
		cerr << "No se puede acceder al archivo: " << argv[1] << endl;

	problema.setNumAnts(10);

	/*					Creacion de matrices Dis y Flu					*/
	nDatos = problema.getNumDatos();

	problema.inicializarDatos(fe);
	
	cout << "###### " << argv[1] << " ######" << endl;
	problema.solInit();

	/*			Algoritmo Greedy			*/
	/*tI = clock();
	problema.AlgGreedy();
	cout<<"\tCoste: "<<problema.getCosteG()<<"\t\t";
	tF = double(clock()- tI)/(double)cps;
	cout<<"Tiempo de ejecucion: "<<tF<<"\n\n";

	cout<<"\tSol Greedy: ";
	for(int i = 0; i < nDatos; i++)
		cout<<problema.solG.solucion[i]<<", ";
	cout<<"\n\n";

	tI = clock();
	problema.SCH_BL();
	tF = double(clock()- tI)/(double)cps;

	cout << "\nCoste SCH: " << problema.getCosteSCH() << endl;
	problema.visualizaV(problema.solSCH.solucion, problema.getNumDatos());
	cout << "Tiempo: " << tF << endl;
*/
	tI = clock();
	problema.SHMM_BL();
	tF = double(clock()- tI)/(double)cps;

	cout << "\nCoste SHMM: " << problema.getCosteSHMM() << endl;
	problema.visualizaV(problema.solSHMM.solucion, problema.getNumDatos());
	cout << "Tiempo: " << tF << endl;

	fe.close();
	problema.clear();
}