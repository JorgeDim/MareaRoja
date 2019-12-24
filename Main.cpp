#pragma once

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <stdio.h>
#include <iomanip>

#include "GlutFreeGlut.h"
#include "GL/glext.h"
#include "globales.h"

#include <objidl.h>
#include <gdiplus.h>
using namespace std;
//#endif


extern int nParticulas;
extern double ThetaMax,ThetaMin,dTheta_med;
extern  int primerdrawVelGL;
//extern int const maxpasadas; 1332
//extern int const maxpasadas;
extern int GL_immediate_mode;
extern int GL_threads;


#include "comunes.h"
////////////////////////////////////WEB

void				GuardarInstantanea();

//#define PAUSA 	                       }iEtapa++;cout<<"iEtapa="<<iEtapa<<endl;if (Etapa==iEtapa) { InicioEtapaP=InicioEtapa;InicioEtapa=Etapa;Guardar=1;
//#define PAUSAF 	                       }iEtapa++;cout<<"iEtapa="<<iEtapa<<endl;if (Etapa==iEtapa) { Guardar=1;
//#define PAUSA2 	Falta_llamar_etapaS=1; }iEtapa++;cout<<"iEtapa="<<iEtapa<<endl;if (Etapa==iEtapa) { Guardar=0;

#define BEGIN_ETAPA 	               if (Etapa==iEtapa) { InicioEtapaP=InicioEtapa;InicioEtapa=Etapa;Guardar=0;
#define PAUSA 	                       }iEtapa++;if (Etapa==iEtapa) { InicioEtapaP=InicioEtapa;InicioEtapa=Etapa;Guardar=0;
#define PAUSAF 	                       }iEtapa++;if (Etapa==iEtapa) { Guardar=0;
#define PAUSA2 	llamar_etapa_Siguiente_PAUSA2=1; }iEtapa++;if (Etapa==iEtapa) { Guardar=0;
#define END_ETAPA 	               }



#if defined _WIN32 || defined _WIN64s
string  GetFileName( const string & prompt ) {
	const int BUFSIZE = 1024;
	char buffer[BUFSIZE] = {0};
	OPENFILENAME ofns = {0};
	ofns.lStructSize = sizeof( ofns );
	ofns.lpstrFile = buffer;
	ofns.nMaxFile = BUFSIZE;
	ofns.lpstrTitle = prompt.c_str();
	GetOpenFileName( & ofns );
	return buffer;
}
#endif



void Calculo_EtapaS_Test(int inicializa);
void Calculo_EtapaS_2017(int inicializa);
void Calculo_EtapaS_2019(int inicializa);

void Calculo_EtapaS_2019_animacionSW(int inicializa);
void form_TesteDeVariablesGlobales() ;
void formulario_opciones_programa();
void DeterminaEscalaInicial(grid3D* gr);
void Fn_LeeZdeT_de_Rodrigo();

void Level0_formulario_OpcionesPrograma();

void Lee_Archivo_Poly(char * file_name);
void Lee_Malla_XML(char * file_name);
void Lee_Velocidad_VTU(char * file_name);

void MyMenu(int value);


void Lectura_Bahia_y_CalculosGeometricosMalla(char * fname) {

	Lectura_Bahia_y_CalculosGeometricosMalla_char++;
	clock_t start_t, end_t;
	double total_t;

	//Etapa de Lectura de la malla

	cout<<"Lectura_Bahia_y_CalculosGeometricosMalla(): fname="<<fname<<endl;
	start_t = clock();
	gtotal->readMSH3D(fname);

	end_t = clock();
	total_t = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;
	cout<<"Lista la lectura en "<<total_t<<" seg";
	if (total_t>60) {
		int min=total_t/60; total_t-=min*60;
		cout<<" ( "<<min<<":"<<total_t<<"min)";
	}
	cout <<endl;

	sprintf(text,"%sFin Etapa\n",text);


	cout<<__LINE__<<": voy a gtotal->GeneraCarasTriPri()...(13:24 min en 76416 Triprismas)..."<<endl;

	start_t = clock();
	gtotal->GeneraCarasTriPri();

	end_t = clock();
	total_t = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;


	cout<<"Listo gtotal->GeneraCarasTriPri() ";
	cout<<" en "<<total_t<<"seg";
	if (total_t>60) {
		int min=total_t/60; total_t-=min*60;
		cout<<" ( "<<min<<":"<<total_t<<"min)";
	}
	cout <<endl;


	gtotal->minmax();
	gtotal->CentroCarasBloques();


	cout<<__LINE__<<": voy a gtotal->Poligonos_Generar_Version3() ";
	gtotal->Poligonos_Generar_Version3();
	//	Modo_DibujaCentroBloques=true;
	//		Modo_DibujaCentroCaras=true;
	err0=1e-12;


	cout<<__LINE__<<": voy a gtotal->Poligonos_Generar_Version3() ";
	gtotal->CalculaVolumen();

}



void Etapa_Lectura_de_Malla()
{
	char fn[100],fnmsh[100],fnbin[100];


	if(1+DBG)cout<<"Etapa_Lectura_de_Malla(): CasoLectura="<<CasoLectura<<" caso="<<caso<<endl;

	if (caso>1) {
		sprintf(fn,"../../../Datos/bahia-TriPrismas%d.msh3D",caso);
		sprintf(fnmsh,"../../../Datos/malla_gtotal%d.msh",caso);
		sprintf(fnbin,"../../../Datos/malla_gtotal%d.bin",caso);
	}
	else{
		sprintf(fn,"../../../Datos/bahia-TriPrismas.msh3D");
		sprintf(fnmsh,"../../../Datos/malla_gtotal.msh");
		sprintf(fnbin,"../../../Datos/malla_gtotal.bin");
	}




	switch (CasoLectura){
	case 1: //Leeer archivo generado por matlab (leeento)
		//Lectura_Bahia_y_CalculosGeometricosMalla("bahia-TriPrismas.msh3D");
		//Lectura_Bahia_y_CalculosGeometricosMalla("bahia-TriPrismas2.msh3D");
		Lectura_Bahia_y_CalculosGeometricosMalla(fn);

		binario=0;
		gtotal->write(fnmsh);
		binario=1;
		gtotal->write(fnbin);
		break;
	case 2: //Lee archivo ASCIII y escribe Binario
		gtotal->read(fnmsh);
		//			gtotal->CentroCarasBloques();
		binario=1;gtotal->write(fnbin);
		break;
	case 3: //Lee el archivo binario
		binario=1;
		gtotal->read(fnbin);
		cout<<__LINE__<<": gtotal->nV3D="<<gtotal->nV3D<<endl;

		//			binario=0;gtotal->write("../../../Datos/malla_gtotal2.msh");
		break;
	case 4: //Lee binario calculo y escribo binario
		binario=1;
		gtotal->read(fnbin);


#if 0 //No generare las caras ya que parecen bien
		cout<<"gtotal->GeneraCarasTriPri()..."<<endl;

		start_t = clock();
		gtotal->TriPri3DAnalizados=0;
		gtotal->Cara.resize(0);
		gtotal->nCaras=0;
		gtotal->GeneraCarasTriPri();

		end_t = clock();
		total_t = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;


		cout<<"Listo gtotal->GeneraCarasTriPri() ";
		cout<<" en "<<total_t<<"seg";
		if (total_t>60) {
			int min=total_t/60; total_t-=min*60;
			cout<<" ( "<<min<<":"<<total_t<<"min)";
		}
		cout <<endl;
#endif


		for (int i=1;i<gtotal->nV3D;i++) {
			gtotal->v3D[i].z *=10;
		}
		gtotal->minmax();
		gtotal->CentroCarasBloques();


		cout<<__LINE__<<": voy a gtotal->Poligonos_Generar_Version3() ";
		gtotal->Poligonos_Generar_Version3();
		//	Modo_DibujaCentroBloques=true;
		//		Modo_DibujaCentroCaras=true;
		err0=1e-12;


		cout<<__LINE__<<": voy a gtotal->Poligonos_Generar_Version3() ";
		gtotal->CalculaVolumen();

		//binario=1;gtotal->write(fnbin);
		break;
	}

	DeterminaEscalaInicial(gtotal);

	//		cout<<"E vecUEsfera="<<vecUEsfera[0]<<","<<vecUEsfera[1]<<"."<<vecUEsfera[2]<<","<<vecUEsfera[3]<<endl;

	if (PrintTiempos) {tic(); cout<<"gtotal->CalculaNormalVertex()..."<<flush;}
	gtotal->CalculaNormalVertex();
	if (PrintTiempos) {cout<<"FIN:gtotal->CalculaNormalVertex() en ";toc();}


	if(DBG)cout<<"Etapa_Lectura_de_Malla():END"<<endl;
}

void DeterminaEscalaInicial(grid3D* gr) {
	//gtotal->xmax
	gr->minmax();
	vecXEsfera[0]=(gr->xmax+gr->xmin)/2;
	vecXEsfera[1]=(gr->ymax+gr->ymin)/2;
	vecXEsfera[2]=(gr->zmax+gr->zmin)/2;
	Escala=10/(gr->xmax-gr->xmin);
}

#if 0
void Etapa_Lectura_de_Velocidad(int step)
{
	//Lectura de campo de Velocidades
	int i;

	if (PrintLecturas) cout<<"Lectura de campo de Velocidades"<<endl;

	U.resize(gtotal->nTriPrisma3D);V.resize(gtotal->nTriPrisma3D);W.resize(gtotal->nTriPrisma3D);

	if (PrintLinea)  cout<<"532"<<endl;
	FILE * pFile;
	int nn,tmp;
	char name[100];
	sprintf(name,"UVW%d_%05d.dat",caso,step);

	if (PrintLinea)  cout<<"532"<<endl;

	pFile = fopen (name,"r");

	fscanf(pFile,"%d\n",&nn);
	if (nn != gtotal->nTriPrisma3D) {
		cout<<"nn != gtotal->nTriPrisma3D: "<<nn<<" != "<<gtotal->nTriPrisma3D<<endl;

		exit(1);
	}

	if (PrintLinea) cout<<"549"<<endl;

	for (i=0;i<nn;i++) {
		double v1,v2,v3;
		fscanf(pFile,"%d %lf %lf %lf\n",&tmp,&v1,&v2,&(W[i]));
		U[i]=v1;V[i]=v2;//W[i]=v3;
		if (i<10) {
			if (PrintLecturas) cout<<tmp<<" "<<v1<<" "<<v2<<" "<<v3<<endl;
			if (PrintLecturas) cout<<tmp<<" "<<U[i]<<" "<<V[i]<<" "<<W[i]<<endl;
		}
		if (tmp!=i) {
			cout<<"tmp!=i"<<tmp<<" "<<i<<endl;
			exit(1);
		}
	}
	fclose(pFile);

	if (glui != NULL) {
		if (PanelFlimite != NULL) PanelFlimite->enable();
		if (Checkbox_particulas != NULL)Checkbox_particulas->enable();
		if (PanelParticulas != NULL)PanelParticulas->enable();
	}
}
#endif

void Etapa_Lectura_de_Velocidades(int cuantos)
{
	//Lectura de campo de Velocidades
	int i,step;

	if (PrintLecturas) cout<<"Etapa_Lectura_de_Velocidades(int cuantos)"<<endl;

	U.resize(cuantos);V.resize(cuantos);W.resize(cuantos);
	for(step=0;step<cuantos;step++) {


		if (PrintLecturas) cout<<"Lectura de campo de Velocidades"<<endl;

		U[step].resize(gtotal->nTriPrisma3D);V[step].resize(gtotal->nTriPrisma3D);W[step].resize(gtotal->nTriPrisma3D);

		if (PrintLinea)  cout<<"532"<<endl;
		FILE * pFile;
		int nn,tmp;
		char name[100];
		sprintf(name,"UVW%d_%05d.dat",caso,step);

		if (PrintLinea)  cout<<"298: name="<<name<<endl;

		pFile = fopen (name,"r");

		fscanf(pFile,"%d\n",&nn);
		if (nn != gtotal->nTriPrisma3D) {
			cout<<"nn="<<nn<<" != gtotal->nTriPrisma3D="<<gtotal->nTriPrisma3D<<endl;
			exit(1);
		}

		if (PrintLinea) cout<<"549"<<endl;

		for (i=0;i<nn;i++) {
			double v1,v2,v3;
			fscanf(pFile,"%d %lf %lf %lf\n",&tmp,&v1,&v2,&(W[step][i]));
			U[step][i]=v1;V[step][i]=v2;//W[i]=v3;
			if (i<10) {
				if (PrintLecturas) cout<<tmp<<" "<<v1<<" "<<v2<<" "<<v3<<endl;
				if (PrintLecturas) cout<<tmp<<" "<<U[step][i]<<" "<<V[step][i]<<" "<<W[step][i]<<endl;
			}
			if (tmp!=i) {
				cout<<"tmp!=i"<<tmp<<" "<<i<<endl;
				exit(1);
			}
		}
		fclose(pFile);
	}

}


void Lectura_ZdeT(int cuantos)
{
	int i,step;

	if (PrintLecturas) cout<<"Lectura_ZdeT"<<endl;

	ZdeT.resize(cuantos);
	for(step=0;step<cuantos;step++) {

		FILE * pFile;
		int nn,tmp;
		char name[100];
		sprintf(name,"Z%d_%05d.msh3D",caso,step+1);//Matlab de 1..N
		//		cout<<"name="<<name<<endl;

		if (PrintLinea)  cout<<"532"<<endl;

		pFile = fopen (name,"r");

		fscanf(pFile,"%d\n",&nn);

		ZdeT[step].resize(nn);

		if (PrintLinea) cout<<"549"<<endl;

		for (i=0;i<nn;i++) {
			double v1,v2,v3;
			fscanf(pFile,"%d %lf\n",&tmp,&v1);
			ZdeT[step][i]=v1;
			if (tmp!=i) {
				cout<<"tmp!=i"<<tmp<<" "<<i<<endl;
				exit(1);
			}
		}
		fclose(pFile);

	}
	cout<<"fin de lectura"<<endl;

}


#if ProyectoPilas
#include "EscurrimientosSuperficiales.h"
#endif


void malla1(int &Etapa, int &iEtapa) {

#if DBG==1
	//cout<<"Malla1()"<<endl;
#endif

	int i,j,iL,malla;
	double x1,x2,x2a,y1,y2,xL,yL,R,xF,yF,RF,R1,R2,RE,th,pi=4*atan(1.0);
	double LZ=Dominio_Hsup; //Basico

#if DBG==1
	cout << "nRLC=" << nRLC << endl;
	cout << "nThLC=" << nThLC << endl;
#endif
	if (nRLC>0) nR=nRLC;
	if (nThLC>0) nTh=nThLC;

	x1=Dominio_Rint;x2a=Dominio_Rmed;x2=Dominio_Xmax-x2a;



	if (Etapa==iEtapa) {

		Escala=1/Dominio_Xmax;
		printf("Mallando \n");
		myfileSalida<<"Mallando "<<endl;
		gtotal->cubo(nR,nTh,nZ,1,1,LZ);g2.cubo(nTh,nR,nZ,1,1,LZ);


		//PAUSA
		ThetaMax=atan(1)*2;
		//		printf("Etapa: Se comprime g ");cout<<endl;
		double q=exp(0.4*(log(Dominio_Rmed/Dominio_Rint)/(nR)));
		cout<<"q="<<q<<endl;
		myfileSalida<<"q="<<q<<endl;
		double Sumq=0,qi=1;
		for (i=0;i<nR;i++) {Sumq +=qi;qi*=q;}
		//cout<<"Sumq="<<Sumq<<endl;
		double L0=(Dominio_Rmed-Dominio_Rint)/Sumq;
		for (i=0;i<gtotal->nV3D;i++) {
			xL=gtotal->v3D[i].x; yL=gtotal->v3D[i].y;
			int nRL=(xL*nR+.5);
			Sumq=0;qi=1;
			for (j=0;j<nRL;j++) {Sumq +=qi;qi*=q;}
			//cout<<"Sumqi="<<Sumq<<endl;
			double Rloc=L0*Sumq;
			th=atan(yL);
			double co=cos(th);
			if (xL<2.0/nR) co=1;
			R1=(x1+Rloc)/co;
			R2=(x1+Rloc);
			R=R2+(R1-R2)*xL;
			gtotal->v3D[i].x=R*cos(th);
			gtotal->v3D[i].y=R*sin(th);
		}

		//PAUSA

		//		printf("Etapa: Se comprime g2\n");cout<<endl;
		for (i=0;i<g2.nV3D;i++) {
			xL=g2.v3D[i].x;	yL=g2.v3D[i].y;
			int nRL=(yL*nR+.5);
			Sumq=0;qi=1;
			for (j=0;j<nRL;j++) {Sumq +=qi;qi*=q;}
			//cout<<"Sumqi="<<Sumq<<endl;
			double Rloc=L0*Sumq;
			th=atan(xL);
			double co=cos(th);
			if (yL<2.0/nR) co=1;
			R1=(x1+Rloc)/co;
			R2=(x1+Rloc);
			R=R2+(R1-R2)*yL;
			g2.v3D[i].x=R*sin(th);
			g2.v3D[i].y=R*cos(th);
		}
		g1=*gtotal;
		*gtotal=g2;
		th=1;
		gtotal->GeneraCaras(-1);

		//PAUSA

		//		printf("Etapa: junto g y g2 \n");cout<<endl;
		gtotal->Junta(g1);

		gtotal->GeneraCaras(gtotal->nH3D); //A partir de ahora se agregan las caras para los bloques actuales
		PAUSA2

		//printf("Etapa: Se agrega corona 2 lado\n");cout<<endl;
		g2.cubo(nR2,nTh,nZ,x2,x2a,LZ);g2.Traslada(x2a,0,0);gtotal->Junta(g2);

		gtotal->GeneraCaras(gtotal->nH3D); //A partir de ahora se agregan las caras para los bloques acuales
		PAUSA2

		//printf("Etapa: Se copia corona lateral 2 veces\n");cout<<endl;
		g2.cubo(nTh,nR2,nZ,x2a,x2,LZ);g2.Traslada(0,x2a,0);gtotal->Junta(g2);
		g2.cubo(nR2,nR2,nZ,x2,x2,LZ);g2.Traslada(x2a,x2a,0);gtotal->Junta(g2);


		gtotal->GeneraCaras(gtotal->nH3D); //A partir de ahora se agregan las caras para los bloques acuales
		PAUSA

		sprintf(text,"%sRetoco la malla\n",text);if (glui != NULL)glui_edittext->set_text(text);


		gtotal->minmax();
		gtotal->CentroCarasBloques();

		vector<R3> SumaCoor(gtotal->nV3D);
		vector<int> NSumaCoor(gtotal->nV3D);

		double err,errG=0,x0,y0,x1,y1;
		for (int iter=0;iter<100;iter++) {
			for (i=0;i<gtotal->nV3D; i++) {
				SumaCoor[i].x=0;
				SumaCoor[i].y=0;
				NSumaCoor[i]=0;
			}

			for (i=0;i<gtotal->nCaras; i++) {
				int iv1=gtotal->Cara[i].iv[0];
				int iv2=gtotal->Cara[i].iv[1];
				int iv3=gtotal->Cara[i].iv[2];
				int iv4=gtotal->Cara[i].iv[3];

				if (sqr(gtotal->v3D[iv1].x)+sqr(gtotal->v3D[iv1].y)+1e-5<Dominio_Rint) continue;
				if (sqr(gtotal->v3D[iv2].x)+sqr(gtotal->v3D[iv2].y)+1e-5<Dominio_Rint) continue;
				if (sqr(gtotal->v3D[iv3].x)+sqr(gtotal->v3D[iv3].y)+1e-5<Dominio_Rint) continue;
				if (sqr(gtotal->v3D[iv4].x)+sqr(gtotal->v3D[iv4].y)+1e-5<Dominio_Rint) continue;

				if (fabs(gtotal->Cara[i].normalCara.z)>0.9) {
					SumaCoor[iv1].x += gtotal->v3D[iv2].x;
					SumaCoor[iv1].y += gtotal->v3D[iv2].y;
					SumaCoor[iv1].x += gtotal->v3D[iv4].x;
					SumaCoor[iv1].y += gtotal->v3D[iv4].y;
					NSumaCoor[iv1]+=2;
					SumaCoor[iv2].x += gtotal->v3D[iv3].x;
					SumaCoor[iv2].y += gtotal->v3D[iv3].y;
					SumaCoor[iv2].x += gtotal->v3D[iv1].x;
					SumaCoor[iv2].y += gtotal->v3D[iv1].y;
					NSumaCoor[iv2]+=2;
					SumaCoor[iv3].x += gtotal->v3D[iv4].x;
					SumaCoor[iv3].y += gtotal->v3D[iv4].y;
					SumaCoor[iv3].x += gtotal->v3D[iv2].x;
					SumaCoor[iv3].y += gtotal->v3D[iv2].y;
					NSumaCoor[iv3]+=2;
					SumaCoor[iv4].x += gtotal->v3D[iv1].x;
					SumaCoor[iv4].y += gtotal->v3D[iv1].y;
					SumaCoor[iv4].x += gtotal->v3D[iv3].x;
					SumaCoor[iv4].y += gtotal->v3D[iv3].y;
					NSumaCoor[iv4]+=2;
				}
			}

			errG=0;
			for (i=0;i<gtotal->nV3D; i++) {
				if (NSumaCoor[i]>4){
					x1=SumaCoor[i].x/NSumaCoor[i];
					y1=SumaCoor[i].y/NSumaCoor[i];
					x0=gtotal->v3D[i].x;
					y0=gtotal->v3D[i].y;
					err=(fabs(x0-x1)+fabs(y0-y1))/(fabs(x0)+fabs(y0));
					if (err>errG) errG=err;
					gtotal->v3D[i].x=x1;
					gtotal->v3D[i].y=y1;
				}
#if DBG==1
				cout <<"NSumaCoor[i]="<<NSumaCoor[i]<<endl;
#endif

			}
		}


		sprintf(text,"%sErrG=%e\nerrMalla=%e\n",text,errG,errMalla);if (glui != NULL)glui_edittext->set_text(text);

		if(errG>errMalla)		Etapa--;

		PAUSA

	}
}


void malla2(int &Etapa, int &iEtapa) {

	int i,j,iL,malla;
	double x1,x2,y1,y2,xL,yL,R,xF,yF,RF,R1,R2,RE,th,pi=4*atan(1.0);
	//	int nR=5,nTh=5,nZ=10;LZ=5;
	//	int nR=7,nTh=4,nZ=1;LZ=1;
	//	int nR=12,nTh=6,nZ=1;LZ=1;
	//	int nR=12,nTh=3,nZ=1;LZ=Dominio_Hsup; //Basico
	//	int nR=34,nTh=10,nZ=1;
	double LZ=Dominio_Hsup; //Basico
	//	int nR=25,nTh=10,nZ=1;LZ=Dominio_Hsup; //Refinado
	x1=Dominio_Rint;x2=Dominio_Xmax;

	cout << "malla2(Etapa=" << Etapa << ",iEtapa="<<iEtapa<<")"<<endl;
	myfileSalida << "malla2(Etapa=" << Etapa << ",iEtapa="<<iEtapa<<")"<<endl;

	if (Etapa==iEtapa) {


		switch (CualMalla) {
		case 2:
			//			nR=34;nTh=10;
			nTh=nTh/2;    //Despues se duplicará la malla ( nth/2*2=nTh)
			ThetaMax=atan(1)*4; //pi+duplicacion
			ThetaMin=-atan(1)*0;

			break;
		case 3:
			//			nR=34;nTh=1;
			ThetaMax=atan(1)/15; //
			ThetaMin=-atan(1)/15;
			//			ThetaMax=atan(1)/2; //Test
			//			ThetaMin=-atan(1)*0;

			break;
		case 4:
			nR=34;nTh=6;
			ThetaMax=atan(1)*2; //
			ThetaMin=-atan(1)*0;

			break;
		default:
			break;
		}

		dTheta_med=(ThetaMax-ThetaMin)/nTh;
		Escala=1/x2;
		printf("Etapa: Se definen g\n");
		gtotal->cubo(nR,nTh,nZ,1,1,LZ);


		//PAUSA
		printf("Etapa: Se Transforma g ");cout<<endl;
		for (i=0;i<gtotal->nV3D;i++) {
			xL=gtotal->v3D[i].x; yL=gtotal->v3D[i].y;
			//			xL *= sqrt(xL);
			if (refinoR==1)
				xL=(1-cos(pi*xL))/2;
			th=ThetaMin+(ThetaMax-ThetaMin)*yL;
			RE=x2;
			R1=x1+(RE-x1)*xL;
			R2=x1+(x2-x1)*xL;
			R=R2+(R1-R2)*xL;
			gtotal->v3D[i].x=R*cos(th);
			gtotal->v3D[i].y=R*sin(th);
		}

		if (CualMalla==2) {
			g2=*gtotal;g2.Rota90Z();g2.Rota90Z();gtotal->Junta(g2);
			ThetaMax=atan(1)*8;
			nTh=nTh*2;  // Vuelvo a valor original (nTh/2*2)
		}

		gtotal->GeneraCaras(-1);    //Inicializa el numero de caras generadas ==0
		gtotal->GeneraCaras(gtotal->nH3D); //A partir de ahora se agregan las caras para los bloques acuales
#if DBG==1
		if (DBG) cout<<"ThetaMax="<<ThetaMax<<endl;
#endif
	}
}
void Calculo_EtapaS_REAL(int inicializa)
{

	AddMensajeRight((char *)"Procesando...");
	DrawGraphics();
	//Calculo_EtapaS_Test(inicializa);
	//Calculo_EtapaS_2017(inicializa);
	//Calculo_EtapaS_2019(inicializa);
	Calculo_EtapaS_2019_animacionSW(inicializa);

	NMensajeRight--;
	AddMensajeRight((char *)"ok...");
	DrawGraphics();

}



void Calculo_EtapaS_Test(int inicializa)
{
	int iEtapa=1;
	Etapa++;

	cout <<"Etapa:"<<Etapa<<"\t iEtapa="<<iEtapa<<endl;

	BEGIN_ETAPA
	{

	}
	PAUSA
	if(1) {
		Lee_Malla_XML("../../../Datos/bahia-TriPrismas.1.xml");
		ResetView();

		//Click Muestra BC....
		formulario_BC_addDatos();
		gluiEligeBC->show();gluiEligeBC_Visible=1;
		// Click DibujaLineas

		cout<<"-->Lee_Velocidad_VTU(\"velocity000200.vtu\");"<<endl;
		Lee_Velocidad_VTU("velocity000200.vtu");
	}
	FlagDibujaLineas=0;

	PAUSA

	cout <<"Etapa b"<<endl;


	END_ETAPA

}


void Calculo_EtapaS_2017(int inicializa)
{
	int iEtapa,i,j ;
	static int inicio=0;


	if(DBG)cout<<"Calculo_EtapaS(int "<<inicializa<<"): BEGIN"<<endl;

	if (inicializa) {
		Etapa=0;
	}
	Etapa++;
	iEtapa=1;




	if (CalculoContinuo==0 && inicio==0 && llamar_etapa_Siguiente_PAUSA2==0) {
		sprintf(text,"Etapa=%d\n%s",Etapa,text0);
		if (glui != NULL) {
			glui_edittext->set_text(text);glui_edittext->redraw_window();
			glui->refresh();
			glutSetWindow(main_window);
			glutPostRedisplay();
		}
		inicio =1;
		Etapa--;

		llamar_etapa_Siguiente_PAUSA2=1;

		return;
	}


	inicio=0;
	llamar_etapa_Siguiente_PAUSA2=0;


	caso=2; //1: malla vieja, 2:Malla nueva, 3: Malla extendida
	CasoLectura=3 ; //1:LeeMatlab, 3:Leebinario

	if (Etapa==iEtapa) {

		{

			gtotal=new grid3D;
			Etapa_Lectura_de_Malla();
			for (i=0;i<gtotal->Cara.size();i++) {
				gtotal->Cara[i].nRandom=(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			}

		}
		//		PAUSA;
		{
			int bin2=1;
			if (bin2==0) {
				cout<<"voy a leer UVWdeT"<<endl;
				Etapa_Lectura_de_Velocidades(2452);//900 ParticulasZlambda
				cout<<"voy a leer zdeT"<<endl;
				Lectura_ZdeT(2452);
				cout<<"Ya lei zdeT"<<endl;
				writeUVWZdeT_Binario("../../../Datos/UVWZdeT.bin")  ;
				cout<<"Ok"<<endl;
				for (i=0;i<10;i++){
					cout <<ZdeT[100][i]<<"\t"<<U[100][i]<<"\t"<<V[100][i]<<"\t"<<W[100][i]<<endl;
				}
				for (i=0;i<10;i++){
					cout <<ZdeT[100][i]<<"\t"<<U[100][i]<<"\t"<<V[100][i]<<"\t"<<W[100][i]<<endl;
				}
				int it=ZdeT.size()-1;
				int iit=ZdeT[it].size()-1;
				cout<<"it="<<it<<"\tiit="<<iit<<endl;
				for (i=0;i<10;i++){
					cout <<ZdeT[it][iit-i]<<"\t"<<U[it][iit-i]<<"\t"<<V[it][iit-i]<<"\t"<<W[it][iit-i]<<endl;
				}
			} else {
				cout<<"Leo archivo binario"<<endl;
				readUVWZdeT_Binario("../../../Datos/UVWZdeT.bin")  ;
				cout<<"Ok"<<endl;
				for (i=0;i<10;i++){
					cout <<ZdeT[100][i]<<"\t"<<U[100][i]<<"\t"<<V[100][i]<<"\t"<<W[100][i]<<endl;
				}
				int it=ZdeT.size()-1;
				int iit=ZdeT[it].size()-1;
				cout<<"it="<<it<<"\tiit="<<iit<<endl;
				for (i=0;i<10;i++){
					cout <<ZdeT[it][iit-i]<<"\t"<<U[it][iit-i]<<"\t"<<V[it][iit-i]<<"\t"<<W[it][iit-i]<<endl;
				}
			}

			if (glui != NULL) {
				if (PanelFlimite != NULL) PanelFlimite->enable();
				if (Checkbox_particulas != NULL)Checkbox_particulas->enable();
				if (PanelParticulas != NULL)PanelParticulas->enable();
			}
			for(i=0; i<gtotal->nV3D;i++){
				gtotal->v3D[i].z0=gtotal->v3D[i].z;
				if ((i%4)==0) {
					//					cout<<i<<endl;
					//					cout << gtotal->v3D[i].x<<","  <<gtotal->v3D[i].y  <<","  <<gtotal->v3D[i].z  <<endl;
					//					cout << gtotal->v3D[i+1].x<<","<<gtotal->v3D[i+1].y<<","  <<gtotal->v3D[i+1].z  <<endl;
					//					cout << gtotal->v3D[i+2].x<<","<<gtotal->v3D[i+2].y<<","  <<gtotal->v3D[i+2].z  <<endl;
					//					cout << gtotal->v3D[i+3].x<<","<<gtotal->v3D[i+3].y<<","  <<gtotal->v3D[i+3].z  <<endl;
				} else {
					if( gtotal->v3D[i].z0 < gtotal->v3D[i-1].z0+0.01 )
						gtotal->v3D[i].z0 =gtotal->v3D[i-1].z0+0.01;
				}
			}

			cout<<"END: gtotal->CentroCarasBloques();"<<endl;
		}
		//		PAUSA;
		{
			char fnmsh[100];
#if 0
			gtotal->TriPri3DAnalizados=0;
			gtotal->Cara.resize(0);
			gtotal->DZmin(1e-2);
			gtotal->nCaras=0;
			gtotal->GeneraCarasTriPri();
			gtotal->Soporte_Vertices();
			gtotal->CentroCarasBloques();
			gtotal->Poligonos_Generar_Version3();

			//// Verificacion
			gtotal->CalculaVolumen();
			CalculaFactorMallaCCRGJSM(*gtotal);


#endif

			if(1) {
				binario=1;
				sprintf(fnmsh,"../../../Datos/malla_gtotal_Voronoi.msh");
				gtotal->read(fnmsh);
			}

			cout <<"1187 C :"<< gtotal->TriPrisma3D[0].centro;
			cout<<"gtotal->CentroCarasBloques();"<<endl;
			gtotal->CentroCarasBloques();
			cout <<"1189 C :"<< gtotal->TriPrisma3D[0].centro;
			if(0) {
				binario=0;
				sprintf(fnmsh,"../../../Datos/malla_gtotal%d.msh",caso);
				gtotal->write(fnmsh);
			}
			if(0) {
				binario=1;
				sprintf(fnmsh,"../../../Datos/malla_gtotal_Voronoi.msh");
				gtotal->write(fnmsh);
			}
		}
		///		PAUSA;
		{
			F2Nodos.resize(gtotal->nV3D,0.0);
			F2Volumenes.resize(gtotal->nVolFinito,0.0);
			F2VolumenesP.resize(gtotal->nVolFinito,0.0);
			R3 CentroMancha;
			CentroMancha.x=(gtotal->xmax+gtotal->xmin)/2;
			CentroMancha.y=(gtotal->ymax+gtotal->ymin)/2;
			CentroMancha.z=(gtotal->zmax+gtotal->zmin)/2;
			double RadioMancha=(gtotal->xmax-gtotal->xmin)/5;
			for (i=0;i<gtotal->nV3D;i++) {
				if( sqr(gtotal->v3D[i].x-CentroMancha.x)+sqr(gtotal->v3D[i].y-CentroMancha.y)<sqr(RadioMancha)) {
					F2Nodos[i]=1;
				}
			}
			for (i=0;i<gtotal->nVolFinito;i++) {
				if( sqr(gtotal->VolFinito[i].centro.x-CentroMancha.x)+sqr(gtotal->VolFinito[i].centro.y-CentroMancha.y)<sqr(RadioMancha)) {
					F2Volumenes[i]=1;
				}
			}
			FlagCalculaEvolucion=1;
		}
		PAUSA;
		// RG 25 enero
		{

			cout <<"1224 C :"<< gtotal->VolFinito[0].centro;
			FuncionF_Nodos.resize(ZdeT.size());
			FuncionF_Volumenes.resize(ZdeT.size());
			FuncionDFDt_Volumenes.resize(ZdeT.size());
			Grad_F.resize(ZdeT.size());
			Achica.resize(ZdeT.size());
			AchicaPrevia.resize(gtotal->nVolFinito);
			Achica_t.resize(gtotal->nVolFinito);
			//			cout<< "1216" <<endl;

			for (i=0;i<ZdeT.size();i++) {
				FuncionF_Nodos[i].resize(gtotal->nV3D);
				FuncionF_Volumenes[i].resize(gtotal->nVolFinito);
				FuncionDFDt_Volumenes[i].resize(gtotal->nVolFinito);
				Grad_F[i].resize(gtotal->nVolFinito);
				Achica[i].resize(gtotal->nVolFinito);
				//				cout<< "1221" <<endl;

				for ( j=0;j< gtotal->nV3D;j++) {
					FuncionF_Nodos[i][j]= ZdeT[i][j]-gtotal->v3D[j].z0;
				}
				//				cout<< "1226" <<endl;

				//Calculo de los valores de F en cada prisma: \eqref{EqFVF} (=73) del archivo TeX
				for ( j=0;j< gtotal->nVolFinito;j++) {
					FuncionF_Volumenes[i][j]=0.0;
					for( int k=0; k < 6; k++){
						FuncionF_Volumenes[i][j] += FuncionF_Nodos[i][ gtotal->TriPrisma3D[j].iv[k] ];
					}
					FuncionF_Volumenes[i][j] /= 6;
					if (i>0) FuncionDFDt_Volumenes[i-1][j]=(FuncionF_Volumenes[i][j]-FuncionF_Volumenes[i-1][j])/3600;
					if (i+1==ZdeT.size()) FuncionDFDt_Volumenes[i][j]=FuncionF_Volumenes[i][j];
				}
				//Calculo del Gradiente de F --> Grad_F[][] y Achica[][]
				for ( j=0;j< gtotal->nVolFinito;j++) {
					//Programación del gradiente de F en el Volumen j: Grad_F[i][j].x...z
					//Ver ecuaciones \eqref{GradF} (=77) del documento TeX
					MatrizSym3x3 M_AC(0.0);
					R3 LadoDerecho(0.0);
					gtotal->VolFinito[j].centro.z=gtotal->TriPrisma3D[j].centro.z;
					for( int k=0; k < 6; k++){
						R3 a_menos_c;
						a_menos_c.x = gtotal->v3D[ gtotal->TriPrisma3D[j].iv[k] ].x  - gtotal->VolFinito[j].centro.x;
						a_menos_c.y = gtotal->v3D[ gtotal->TriPrisma3D[j].iv[k] ].y  - gtotal->VolFinito[j].centro.y;
						a_menos_c.z = gtotal->v3D[ gtotal->TriPrisma3D[j].iv[k] ].z0 - gtotal->VolFinito[j].centro.z;

						LadoDerecho += a_menos_c *
								(FuncionF_Nodos[i][ gtotal->TriPrisma3D[j].iv[k] ] - FuncionF_Volumenes[i][j] );

						M_AC.a += a_menos_c.x*a_menos_c.x;
						M_AC.b += a_menos_c.x*a_menos_c.y;
						M_AC.c += a_menos_c.x*a_menos_c.z;
						M_AC.d += a_menos_c.y*a_menos_c.y;
						M_AC.e += a_menos_c.y*a_menos_c.z;
						M_AC.f += a_menos_c.z*a_menos_c.z;
					}
					Grad_F[i][j]=M_AC.inversa()*LadoDerecho;
					Achica[i][j]=1+Grad_F[i][j].z;
					if (j<-10) {
						cout << "M="<<M_AC;
						cout << "Minv="<<M_AC.inversa();
						cout << "L="<<LadoDerecho;
						cout << Grad_F[i][j];
					}
				}
			}
			sprintf(text0,"%sCalculados Grad_F(itime,iVol)\n",text0);
			sprintf(text,"%sCalculados Grad_F(itime,iVol)\n",text);if (glui != NULL) glui_edittext->set_text(text);
		}
		// RG 25 enero FIN
		PAUSA;
		{
			Etapa--; //Ultima etapa por ahora
			Calculando=1-Calculando;
			sprintf(text,"%sCalculando=%d\n",text,Calculando);
			sprintf(text,"%sNada mas\n",text);if (glui != NULL) glui_edittext->set_text(text);

		}
		PAUSA;
		{
			gtotal->generaPoligonos2Algunos(CualesRehacer);
			err0=1e-10;

			//cout<<"***Etapa="<<Etapa<<"IniciEtapa="<<InicioEtapa <<"EtapaGlobal2Local="<<EtapaGlobal2Local<<endl;

			gtotal->CalculaVolumen();

			sprintf(text,"%sFin Etapa\n",text);if (glui != NULL) glui_edittext->set_text(text);

			TiempoCalculo=Datos_dt ;
			Reiterando=0;
			cout<<"Inicializo Temperaturas = Tambiente "<<endl;
			myfileSalida<<"Inicializo Temperaturas = Tambiente "<<endl;
			TempPilaBloquesPrevia.resize(gtotal->nH3D);
			TempPilaBloques.resize(gtotal->nH3D);
			for (i=0;i<gtotal->nH3D;i++) {
				TempPilaBloquesPrevia[i]=Datos_Tambiente;
				TempPilaBloques[i]=Datos_Tambiente;
			}
		}
		PAUSA
		{
			cout<<"Etapa: Nada mas"<<endl;

			EtapaGlobal=EtapaFIN;

			FinEtapas=1;
			CalculoContinuo=0;
			sprintf(text,"%sNada mas\n",text);
			if (glui != NULL) glui_edittext->set_text(text);
			Etapa--;
		}
		PAUSAF
		{
			cout<<"Etapa: Nada mas"<<endl;

			EtapaGlobal=EtapaFIN;

			FinEtapas=1;
			CalculoContinuo=0;
			sprintf(text,"%sNada mas\n",text);if (glui != NULL) glui_edittext->set_text(text);
			Etapa--;
		}


	}


	if(DBG)cout<<"Calculo_EtapaS(int "<<inicializa<<"): END"<<endl;

}

void Etapa_LecturaMalla_gtotal_Voronoi_MSH()
{
	char fnmsh[100];

	if(1) {
		binario=1;
		sprintf(fnmsh,"../../../Datos/malla_gtotal_Voronoi.msh");
		gtotal->read(fnmsh);
		DeterminaEscalaInicial(gtotal);
	}

	gtotal->CentroCarasBloques();
	if(0) {
		binario=0;
		sprintf(fnmsh,"../../../Datos/malla_gtotal%d.msh",caso);
		gtotal->write(fnmsh);
	}
	if(0) {
		binario=1;
		sprintf(fnmsh,"../../../Datos/malla_gtotal_Voronoi.bin");
		gtotal->write(fnmsh);
	}
	sprintf(text,"%sOk. Se ha leido malla_gtotal_Voronoi.msh\n",text);if (glui != NULL) glui_edittext->set_text(text);

}

void Etapa_Define_Mancha_circular_inicial()
{
	int i;
	F2Nodos.resize(gtotal->nV3D,0.0);
	F2Volumenes.resize(gtotal->nVolFinito,0.0);
	F2VolumenesP.resize(gtotal->nVolFinito,0.0);
	R3 CentroMancha;
	CentroMancha.x=(gtotal->xmax+gtotal->xmin)/2;
	CentroMancha.y=(gtotal->ymax+gtotal->ymin)/2;
	CentroMancha.z=(gtotal->zmax+gtotal->zmin)/2;
	double RadioMancha=(gtotal->xmax-gtotal->xmin)/5;
	for (i=0;i<gtotal->nV3D;i++) {
		if( sqr(gtotal->v3D[i].x-CentroMancha.x)+sqr(gtotal->v3D[i].y-CentroMancha.y)<sqr(RadioMancha)) {
			F2Nodos[i]=1;
		}
	}
	for (i=0;i<gtotal->nVolFinito;i++) {
		if( sqr(gtotal->VolFinito[i].centro.x-CentroMancha.x)+sqr(gtotal->VolFinito[i].centro.y-CentroMancha.y)<sqr(RadioMancha)) {
			F2Volumenes[i]=1;
		}
	}
	FlagCalculaEvolucion=1;
	sprintf(text,"%sOk. Se ha definido F2Nodos\n",text);if (glui != NULL) glui_edittext->set_text(text);
	if (DBGGtotal && gtotal) cout<<"1069:gtotal,gtotal->Cara.size()="<<gtotal<<","<<gtotal->Cara.size()<<endl;
}


void Calculo_EtapaS_2019_animacionSW(int inicializa)
{

	int iEtapa=1;

	Etapa++;

	BEGIN_ETAPA
	CasoLectura=1;//Lee ascii
	//CasoLectura=2;//Lee ascii y write binario
	//CasoLectura=3;//Lee binario
	//CasoLectura=4;//Lee binario y genera ...
	caso=1;  //malla fina original
	gtotal=new grid3D;
	Etapa_Lectura_de_Malla();
	PAUSA

	Fn_LeeZdeT_de_Rodrigo();


	PAUSA;
	{
		Etapa--; //Ultima etapa por ahora
		Calculando=1-Calculando;
		sprintf(text,"%sCalculando=%d\n",text,Calculando);
		sprintf(text,"%sNada mas\n",text);if (glui != NULL) glui_edittext->set_text(text);

	}

	END_ETAPA
}


void Fn_LeeZdeT_de_Rodrigo()
{
	int i,j,ii;
	int deltai=1800;
	int ndeltai=84600/deltai;
	char dir[100],file[200];
	int LeeTXT=0;
	double x,y,z;
	ifstream file_stream;
	strcpy(dir,	"XYZProfesorJSMAmplitud2");
	ZdeT.resize(ndeltai);
	if (LeeTXT==1) {
		ii=0;

		if (PrintTiempos) {tic(); cout<<__LINE__<<":LecturaTXT... (47 seg 1 dia cada media hora)"<<endl;}
		for (i=deltai;i<=ndeltai*deltai;i+=deltai){
			cout<<__LINE__<<": i="<<i<<endl;
			ZdeT[ii].resize(gtotal->nV3D);
			sprintf(file,"%s/Tiempo%d.txt",dir,i);
			file_stream.open (file);
			for (j=0;j<gtotal->nV3D;j++) {
				//cout <<"j="<<j<<endl;
				file_stream	>>x	>>y  >>z;
				//cout <<"j="<<j<<" z="<<z<<endl;
				ZdeT[ii][j]=z;
				//cout <<"j="<<j<<" z="<<z<<endl;
			}
			ii++;
			file_stream.close();
		}

		if (PrintTiempos) {cout<<__LINE__<<":FIN LecturaTXT...en ";toc();}
		//WriteBinario

		if (PrintTiempos) {tic(); cout<<__LINE__<<":WriteBinario..."<<endl;}
		{
			sprintf(file,"%s/Tiempos.bin",dir);



			ofstream file_ostream;
			file_ostream.open (file, std::ios::out | std::ios::binary);

			for (i=0;i<ZdeT.size();i++) {
				cout<<__LINE__<<": i="<<i<<endl;
				for (j=0;j<gtotal->nV3D;j++) {
					z=ZdeT[i][j];
					file_ostream.write((char*)&z,sizeof(z));
				}
			}

			file_ostream.close();
		}//END WriteBinario
		if (PrintTiempos) {cout<<__LINE__<<":FIN WriteBinario...en ";toc();}

	} // end caso if (LeeTXT==1)
	else {
		//LeeBinario
		if (PrintTiempos) {tic(); cout<<__LINE__<<":LeeBinario..."<<endl;}
		{
			sprintf(file,"%s/Tiempos.bin",dir);
			ifstream file_istream;
			file_istream.open (file, std::ios::in | std::ios::binary);
			for (i=0;i<ZdeT.size();i++) {
				cout<<__LINE__<<": i="<<i<<endl;
				ZdeT[i].resize(gtotal->nV3D);
				for (j=0;j<gtotal->nV3D;j++) {
					file_istream.read((char*)&z,sizeof(z));
					ZdeT[i][j]=z;
				}
			}
			file_istream.close();
		}//END LeeBinario
		if (PrintTiempos) {cout<<__LINE__<<":FIN LeeBinario...en ";toc();}
	}

}

void Calculo_EtapaS_2019(int inicializa)
{
	int iEtapa=1,i,j ;
	static int inicio=0;

	Etapa++;
	sprintf(text,"");

	if(DBG)cout<<"Calculo_EtapaS_2019(int "<<inicializa<<"): BEGIN, Etapa="<<Etapa<<endl;
	if (DBGGtotal && gtotal)		cout<<"gtotal,gtotal->Cara.size()="<<gtotal<<","<<gtotal->Cara.size()<<endl;


	BEGIN_ETAPA



	gtotal=new grid3D;

	sprintf(text,"%sOk. Solo hemos definido gtotal\n",text);if (glui != NULL) glui_edittext->set_text(text);

	PAUSA
	Etapa_LecturaMalla_gtotal_Voronoi_MSH();
	sprintf(text,"%sOk. Etapa_LecturaMalla_gtotal_Voronoi_MSH();\n",text);if (glui != NULL) glui_edittext->set_text(text);
	cout<<"gtotal->nV3D="<<gtotal->nV3D<<endl;
	for (i=0;i<=100;i+=4) {
		cout <<i+1<<" "<<gtotal->v3D[i].x<<" "<<gtotal->v3D[i].y<<" "<<gtotal->v3D[i].z<<endl;
	}

	PAUSA
	lectura_archivo_UVWZdeT_bin();

	sprintf(text,"%sOk. lectura_archivo_UVWZdeT_bin();\n",text);if (glui != NULL) glui_edittext->set_text(text);

	PAUSA;
	{
		Etapa--; //Ultima etapa por ahora
		Calculando=1-Calculando;
		sprintf(text,"%sCalculando=%d\n",text,Calculando);
		sprintf(text,"%sNada mas\n",text);if (glui != NULL) glui_edittext->set_text(text);

	}
	PAUSA;
	Etapa_Define_Mancha_circular_inicial();

	PAUSA
	{
		cout<<"Etapa: Nada mas"<<endl;

		EtapaGlobal=EtapaFIN;

		FinEtapas=1;
		CalculoContinuo=0;
		sprintf(text,"%sNada mas\n",text);
		if (glui != NULL) glui_edittext->set_text(text);
		Etapa--;
	}

	PAUSAF
	{
		cout<<"Etapa: Nada mas"<<endl;

		EtapaGlobal=EtapaFIN;

		FinEtapas=1;
		CalculoContinuo=0;
		sprintf(text,"%sNada mas\n",text);if (glui != NULL) glui_edittext->set_text(text);
		Etapa--;
	}



	END_ETAPA

	if (DBGGtotal && gtotal)		cout<<"gtotal,gtotal->Cara.size()="<<gtotal<<","<<gtotal->Cara.size()<<endl;
	if(DBG)cout<<"Calculo_EtapaS_2019(int "<<inicializa<<"): END"<<endl;

}


void NuevaLecturaDeDatos() {


#if defined _WIN32 || defined _WIN64s

	string fname = GetFileName( "Number which file: " );
	/*	ifstream ifs( fname.c_str() );
	if ( ! ifs.is_open() ) {
		cerr << "cannot open " << fname << " for input" << endl;
	}
	else {
	 */
	strcpy(file_name,fname.c_str());
	SaveOrRead(file_name,SaveReadMode);

	if (SaveReadMode==2) { //2:read
		sprintf(text,"%s\nEtapa=%d\n",file_name,Etapa);
		glui_edittext->set_text(text);glui_edittext->redraw_window();
		glui->refresh();
		//			glutSetWindow(main_window);
		//			glutPostRedisplay();
	}
	/*
		string line;
		int lineno = 1;
		while( getline( ifs, line ) ) {
			cout << setw(5) << right << lineno++ << " : " << line << "\n";
		}
	 */
	//}
#endif

}

#if defined(GLUI_GLUI_H)
void control_cb_file( int control )
{
	cout<<__FUNCTION__<<"@"<<__LINE__<<" control="<<control<<endl;
	if (control==0) {
		glui_hide=1;
		glui->hide();
	}
	if (control==7) {
		strcpy(file_name,FileBrowserGlui2->archivo);
		cout <<"\nfile_name=" <<file_name<<endl;

		glui2->hide();
		SaveOrRead(file_name,SaveReadMode);

		if (SaveReadMode==2) { //2:read
			sprintf(text,"%s\nEtapa=%d\n",file_name,Etapa);
			glui_edittext->set_text(text);glui_edittext->redraw_window();
			glui->refresh();
			glutSetWindow(main_window);
			glutPostRedisplay();
		}

		//				glui2->show();
	}
	if (control==17) {   // 2176: FileBrowserGlui2 = new GLUI_FileBrowser2(panel2, "", false, 17, control_cb_file);
		strcpy(file_name,FileBrowserGlui2->archivo);
		cout <<"file_name=" <<file_name<<endl;

	}
	if (control==18) { //lectura archivo .poly
		strcpy(file_name,FileBrowserGlui2->archivo);
		cout <<"file_name=" <<file_name<<endl;
		Lee_Archivo_Poly(file_name);

	}
	if (control==19) { //lectura archivo .xml
		strcpy(file_name,FileBrowserGlui2->archivo);
		cout <<"file_name=" <<file_name<<endl;
		Lee_Malla_XML(file_name);

	}


	if (control==2210) { //lectura archivo Velocidad.vtu
		strcpy(file_name,FileBrowserGlui2->archivo);
		cout <<"file_name=" <<file_name<<endl;
		Lee_Velocidad_VTU(file_name);

	}


	glui2->hide();
}

void Lee_Archivo_Poly(char * file_name)
{

	gPoly=new grid3D;
	gPoly->readPOLY(file_name);

}

void str_replace(char *target, const char *needle, const char *replacement)
{
	char buffer[1024] = { 0 };
	char *insert_point = &buffer[0];
	const char *tmp = target;
	size_t needle_len = strlen(needle);
	size_t repl_len = strlen(replacement);

	while (1) {
		const char *p = strstr(tmp, needle);

		// walked past last occurrence of needle; copy remaining part
		if (p == NULL) {
			strcpy(insert_point, tmp);
			break;
		}

		// copy part before needle
		memcpy(insert_point, tmp, p - tmp);
		insert_point += p - tmp;

		// copy replacement string
		memcpy(insert_point, replacement, repl_len);
		insert_point += repl_len;

		// adjust pointers, move on
		tmp = p + needle_len;
	}

	// write altered string back to target
	strcpy(target, buffer);
}


void Lee_Velocidad_VTU(char * file_name)
{
	ifstream myfile;
	int i,j,k,itmp;
	cout<<"Lee_Velocidad_VTU("<<file_name<<")"<<endl;

	XMLNode xMainNode=XMLNode::openFileHelper(file_name,"VTKFile");
	XMLNode Node2=xMainNode.getChildNode("UnstructuredGrid");
	XMLNode Node3=Node2.getChildNode("Piece");
	XMLNode Node4=Node3.getChildNode("PointData");
	XMLNode Node5=Node4.getChildNode("DataArray");

	cout<<"type="<<Node5.getAttribute("type")<<endl;
	cout<<"Name="<<Node5.getAttribute("Name")<<endl;
	cout<<"NumberOfComponents="<<Node5.getAttribute("NumberOfComponents")<<endl;
	cout<<"format="<<Node5.getAttribute("format")<<endl;

	printf("nText :'%d'\n", Node5.nText());
	//	  printf("getText :'%s'\n", Node5.getText());


	char *n = strtok((char *)Node5.getText(), " ");
	i=0;
	UU.resize(gtotal->nV3D);
	VV.resize(gtotal->nV3D);
	WW.resize(gtotal->nV3D);

	for (i=0;i<UU.size();i++) {
		UU[i] = atof(n);n = strtok(NULL, " ");
		VV[i] = atof(n);n = strtok(NULL, " ");
		WW[i] = atof(n);n = strtok(NULL, " ");
		//	      cout<<"i="<<i<<"   U="<<UU[i]<<"   V="<<VV[i]<<"   W="<<WW[i]<<endl;
	}
	VelocidadEnVertex=true;
}


void Lee_Malla_XML(char * file_name) //Lee .xml y .face
{

	char buffer[1024] = { 0 };
	char *insert_point = &buffer[0];

	gtotal=new grid3D;
	gtotal->readMallaXML(file_name);
	memcpy(insert_point, file_name, strlen(file_name));
	str_replace(buffer, ".xml", ".face");

	cout <<"New_file_name="<<buffer<<endl;
	gtotal->readFACE(buffer);
	gtotal->CentroCarasBloques();


	gtotal->minmax();
	//	gtotal->GeneraCarasTetrahedros();
	//	gtotal->GeneraCarasTriPri();
	Escala=M3_Escala/(gtotal->xmax-gtotal->xmin);

	ResetView();

}


void control_cb_BC( int control )
{
	switch(control) {
	case 1:

		gluiEligeBC->hide();gluiEligeBC_Visible=0;
	}
}




void Level0_control_cb( int control )
{
	switch(control) {
	case 100:

		glui_OpcionesProgramaLevel0->hide();
	}
}


void control_cb( int control )
{

	glui->sync_live();

#if DBG==1
	printf( "control_cb: %d", control );cout<<endl;
#endif
	if (control==0) {
		glui_hide=1;
		glui->hide();
	}
	if (control==7) {
		string file_name;
		file_name = "";
		file_name = FileBrowserGlui2->get_file();
#if DBG==1
		cout <<"\nfile_name=" <<file_name<<endl;
#endif
	}


	if (control==1002) {  // Save
		SaveReadMode=1;
		NuevaLecturaDeDatos();
		/*
		fb->entrando=1;
		fb->EditText_archivo->orig_text="";

		fb->fbreaddir(".");
		fb->execute_callback();
		glui2->show();
		 */

	}
	if (control==1003) {  // Read
		SaveReadMode=2;
		/*
		fb->EditText_archivo->orig_text="";
		fb->entrando=1;
		fb->fbreaddir(".");
		fb->execute_callback();
		glui2->show();
		 */
		NuevaLecturaDeDatos();



	}

	if (control==1102) {  // Plot Particulas
		//		MODO_Ejes=!MODO_Ejes;
	}

	if (control==1101) {  // Plot Particulas
		//		MODO_CampoVelocidades=!MODO_CampoVelocidades;
	}

	if (control==8001) { //Reset FPS
		nframes2=nframes;
		clock2=clock0;
	}


	if (DBG) cout <<"control_cb(): control="<<control<<" __LINE__="<<__LINE__<<endl;

	if (control==11030) {  //Previo a 1103
		/*
		int i,j,nt,tmp;
		nt=gtotal->Cara.size();
		AlgunosI.resize(nt);
		AlgunosRND.resize(nt);
		for(i=0;i<nt;i++) {
			AlgunosI[i]=i;
			AlgunosRND[i]=rand();
		}
		for(i=0;i<nt;i++) {
			for (j=0;j<nt;j++) {
				if (AlgunosRND[i]<AlgunosRND[j]) {
					tmp=AlgunosRND[i];
					AlgunosRND[i]=AlgunosRND[j];
					AlgunosRND[j]=tmp;
					tmp=AlgunosI[i];
					AlgunosI[i]=AlgunosI[j];
					AlgunosI[j]=tmp;
				}
			}
		}
		 */
		control=1103;
	}


	if (DBG) cout <<"control_cb(): control="<<control<<" __LINE__="<<__LINE__<<endl;

	if (control==1103) {  //

		cout<< "ModoDibujaAlgunos="<< ModoDibujaAlgunos<<endl;
		if (ModoDibujaAlgunos) {
			/*
			int i,cuantos=10,j;
			cuantos=Algunos_Porcentaje*gtotal->Cara.size()/100;
			CualesCarasDibuja.resize(cuantos);
			cuantos=gtotal->h3D.size()/3;
			cuantos=std::min(500,cuantos);
			CualesBloquesDibuja.resize(cuantos);
			//			cout <<"gtotal->nH3D="<<gtotal->nH3D<<endl;
			//			int rmax=0;
			for (i=0;i<CualesCarasDibuja.size();i++) {
				CualesCarasDibuja[i]=AlgunosI[i];
				//				rmax=max(rmax,CualesCarasDibuja[i]);
			}
			for (i=0;i<CualesBloquesDibuja.size();i++) {
				CualesBloquesDibuja[i]=rand()*gtotal->h3D.size()/RAND_MAX;
				//				rmax=max(rmax,CualesCarasDibuja[i]);
			}
			//			cout<<"rmax="<<rmax<<endl;

			 */

		}
	}


	if (DBG) cout <<"control_cb(): control="<<control<<" __LINE__="<<__LINE__<<endl;

#if ProyectoPilas
	if (control==1105) {  //Calculo la proporcion de volument sobre TLimite

		cout<<"cb: Calcular proporcion de volumen sobre TLimite, EtapaGlobal="<<EtapaGlobal<<endl;
		if(EtapaGlobal==ETAPA_CALCULO_T_ARRIBA)
			calculaFracionVolumen(Temperatura);
		if(EtapaGlobal==ETAPA_CALCULO_T_PILA)
			calculaFracionVolumen(TempPilaBloques);
	}
#endif



	if (DBG) cout <<"control_cb(): control="<<control<<" __LINE__="<<__LINE__<<endl;

	if (control==1001) {
		if (DBG) cout <<"control_cb(): control="<<control<<" __LINE__="<<__LINE__<<endl;
		Calculo_EtapaS();
	}



	if (DBG) cout <<"control_cb(): control="<<control<<" __LINE__="<<__LINE__<<endl;

	if (control==9001) { // TestDeVariables
		vecUEsfera[0]=0;
		vecUEsfera[1]=0;
		vecUEsfera[2]=0;
	}




	if (DBG) cout <<"control_cb(): control="<<control<<" __LINE__="<<__LINE__<<endl;

	if (control==9002) { // TestDeVariables


		//cout<<"M \nGlobalOldFovy="<<GlobalOldFovy<<" GlobalFovy="<<GlobalFovy<<endl;
		//cout<<"M A1 vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<","<<vecUEsfera[3]<<"]"<<endl;


#if (1==1)
		GLdouble xP,yP,zP,xW,yW,zW;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective((GLdouble)GlobalOldFovy, aspect, (GLdouble)1, (GLdouble)100.0);
		glTranslated( 0, 0, -10);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslated( 0, -1, -10);

		FuncionesOpenGL::World2Win((GLdouble)vecUEsfera[0], (GLdouble)vecUEsfera[1], (GLdouble)vecUEsfera[2],
				&xP,&yP,&zP);
		//cout<<"M xyzP=[ "<<xP<<" , "<<yP<<" , "<<zP<<" ]"<<endl;
		//cout<<"M A2 vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<","<<vecUEsfera[3]<<"]"<<endl;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective((GLdouble)GlobalFovy, aspect, (GLdouble)1, (GLdouble)100.0);
		glTranslated( 0, 0, -10);
		FuncionesOpenGL::Win2World(xP,yP,zP, &xW,&yW,&zW);
		zW /= FactorZ;
		vecUEsfera[0]=xW;
		vecUEsfera[1]=yW;
		vecUEsfera[2]=zW;
		//cout<<"M A3 vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<","<<vecUEsfera[3]<<"]"<<endl;


#else

		glGetDoublev( GL_PROJECTION_MATRIX, FuncionesOpenGL::projection );
		//Coordenadas screen= Projection*(Old)vecU
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslated( 0, -1, -10);
		glGetDoublev( GL_PROJECTION_MATRIX, FuncionesOpenGL::modelview);

		MatrizTrXvector4(FuncionesOpenGL::modelview,vecUEsfera,vecDUEsfera);  //Eye_Coord
		MatrizTrXvector4(FuncionesOpenGL::projection,vecDUEsfera,vecUEsfera); //Clip_Coord
		//		vecUEsfera[0] /= vecUEsfera[3];
		//		vecUEsfera[1] /= vecUEsfera[3];
		//		vecUEsfera[2] /= vecUEsfera[3];
		//		vecUEsfera[3] /= vecUEsfera[3];
		//cout<<"M BvecDUEsfera=["<<vecDUEsfera[0]<<","<<vecDUEsfera[1]<<","<<vecDUEsfera[2]<<","<<vecDUEsfera[3]<<"]"<<endl;
		//cout<<"M CvecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<","<<vecUEsfera[3]<<"]"<<endl;

		//Encuentra nuevas matrices

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective((GLdouble)GlobalFovy, aspect, (GLdouble)1, (GLdouble)100.0);
		glTranslated( 0, 0, -10);
		glGetDoublev( GL_PROJECTION_MATRIX, FuncionesOpenGL::projection );

		InvierteMatriz(FuncionesOpenGL::projection,MatrizINV);
		MatrizTrXvector4(MatrizINV,vecUEsfera,vecDUEsfera);
		InvierteMatriz(FuncionesOpenGL::modelview,MatrizINV);
		MatrizTrXvector4(MatrizINV,vecDUEsfera,vecUEsfera);
		//		MatrizXvector4(FuncionesOpenGL::projection,vecDUEsfera,vecUEsfera);
#endif

		//cout<<"M BvecDUEsfera=["<<vecDUEsfera[0]<<","<<vecDUEsfera[1]<<","<<vecDUEsfera[2]<<","<<vecDUEsfera[3]<<"]"<<endl;
		//cout<<"M FvecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<","<<vecUEsfera[3]<<"]"<<endl;
		Escala *= GlobalFovy/GlobalOldFovy;
		GlobalOldFovy=GlobalFovy;
	}

	if (DBG) cout <<"control_cb(): control="<<control<<" __LINE__="<<__LINE__<<endl;

	if (control==10001) {
		char *str,str2[1000],str3[1000];
		str=(char *)Glui3_EditParametros->get_text();
		cout<<"str="<<str<<endl;

		glui3->hide();
		char nombre[100];
		double valor;

		string snombre,svalor,linea;

		snombre=str;
		istringstream is1( snombre );


		while (getline(is1,linea)) {
			if(DBG) cout<<"linea="<<linea<<endl;
			istringstream trozo( linea );
			string snombre,svalor;
			if (linea[0]=='#') continue;
			if (!getline( trozo, snombre, '=' )) continue;
			replaceAll(snombre, "\t","");
			replaceAll(snombre, " ","");
			if (!getline( trozo, svalor, '\n' )) continue;
			double ff=atof(svalor.c_str());

			if 		(snombre=="malla"			) {CualMalla		=ff;}
			else if (snombre=="CalculoContinuo"	) {CalculoContinuo	=ff;}
			else if (snombre=="TipoCalculo"	    ) {TipoCalculo   	=ff;}
			else if (snombre=="nR"				) {nR				=ff;}
			else if (snombre=="nR2"				) {nR2				=ff;}
			else if (snombre=="nTh"				) {nTh				=ff;}
			else if (snombre=="NDivZ"			) {NDivZ			=ff;}
			else if (snombre=="Dominio_Xmax"	) {Dominio_Xmax		=ff;}
			else if (snombre=="Dominio_Hmax"	) {Dominio_Hmax		=ff;}
			else if (snombre=="Dominio_Hsup"	) {Dominio_Hsup		=ff;}
			else if (snombre=="Dominio_Rint"	) {Dominio_Rint		=ff;}
			else if (snombre=="rhof"			) {Datos_rhof		=ff;}
			else if (snombre=="rhos"			) {Datos_rhos		=ff;}
			else if (snombre=="phi"				) {Datos_phi		=ff;}
			else if (snombre=="cf"				) {Datos_cf			=ff;}
			else if (snombre=="Tinyeccion"		) {Datos_Tinyeccion	=ff;}
			else if (snombre=="Tambiente"		) {Datos_Tambiente	=ff;}
			else if (snombre=="TLimite" 		) {TLimite      	=ff;}
			else if (snombre=="hc_superior"		) {Datos_hc_superior=ff;}
			else if (snombre=="hc_inferior"		) {Datos_hc_inferior=ff;}
			else if (snombre=="hc_lateral"		) {Datos_hc_lateral	=ff;}
			else if (snombre=="KTermofilm"		) {Datos_KTermofilm	=ff;}
			else if (snombre=="eTermofilm"		) {Datos_eTermofilm	=ff;}
			else if (snombre=="VinyeccionLtsHr"	) {VinyeccionLtsHr	=ff;}
			else if (snombre=="kf"				) {Datos_kf			=ff;}
			else if (snombre=="ks"				) {Datos_ks				=ff;}
			else if (snombre=="Kevaporacion"	) {Datos_Kevaporacion	=ff;}
			else if (snombre=="HumedadAmbiental") {Datos_HumedadAmbiental=ff;}
			else if (snombre=="DistanciaAlBorde") {Datos_DistanciaAlBorde=ff;}
			else if (snombre=="CalorLatenteFluido") {Datos_CalorLatenteFluido=ff;}
			else {
				cout <<"Comando en Archivo de datos no considerado:"<<endl;
				cout <<"variable="<< snombre<<"\t="<<ff<<endl;
			}



		}

	}
	if (control==10002) {
		glui3->hide();
	}

	if (control==20001) {
		gluiHelp2->hide();
	}

	if (control==30001) {
		gluiOpciones->hide();
	}
	if (control==30002) {
		//cout<<"PrintMouse1="<<PrintMouse<<endl;
		//if (gluiOpciones) gluiOpciones->sync_live();
		//cout<<"PrintMouse2="<<PrintMouse<<endl;
	}


	if (control==2223) 		{
		if (DBG) cout <<"if (control==2223):"<<__LINE__<<endl;
		primerdrawVelGL=1;
	}

	if (control==25000)
	{

		switch(ListaCaso) {
		case 0:
			FlagMuestraCaraSuperior=1;
			FlagCaraSuperiorTransparente=0;
			FlagCaraSuperiorTextura=0;
			FlagCaraInferiorTextura=0;
			FlagDibujaLineas=1;
			break;
		case 1:
			FlagMuestraCaraSuperior=0;
			FlagCaraSuperiorTransparente=0;
			FlagCaraSuperiorTextura=0;
			FlagCaraInferiorTextura=1;
			FlagDibujaLineas=0;
			break;
		case 2:
			FlagMuestraCaraSuperior=1;
			FlagCaraSuperiorTransparente=1;
			FlagCaraSuperiorTextura=1;
			FlagCaraInferiorTextura=1;
			FlagDibujaLineas=0;
			break;
		case 3:
			Calculo_EtapaS();
			FlagMuestraCaraSuperior=1;
			FlagCaraSuperiorTransparente=1;
			FlagCaraSuperiorTextura=1;
			FlagCaraInferiorTextura=1;
			FlagDibujaLineas=0;

			//Reset View
			Escala=M3_Escala/(gtotal->xmax-gtotal->xmin);

			for (int i=0;i<16;i++){
				MatrizRotacionGlobal[i]   =M3_MatrizRotacionGlobal[i];
				MatrizRotacionGlobalINV[i]=M3_MatrizRotacionGlobalINV[i];
			}
			cout <<"vecXEsfera=["<<vecXEsfera[0]<<","<<vecXEsfera[1]<<","<<vecXEsfera[2]<<"]"<<endl;
			cout <<"vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<"]"<<endl;
			vecDUEsfera[0]=(vecXEsfera[0]-(gtotal->xmax+gtotal->xmin)/2)*Escala;
			vecDUEsfera[1]=(vecXEsfera[1]-(gtotal->ymax+gtotal->ymin)/2)*Escala;
			vecDUEsfera[2]=(vecXEsfera[2]-(gtotal->zmax+gtotal->zmin)/2)*Escala;

			MatrizXvector4(MatrizRotacionGlobal,vecDUEsfera,vecUEsfera);

			cout <<"vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<"]"<<endl;
			break;

		case 4:
			MODO_CampoVelocidades=1;
			nParticulas=1000;
			primerdrawVelGL=1;
			FlagMuestraCaraSuperior=0;
			FlagCaraSuperiorTransparente=1;
			FlagCaraSuperiorTextura=1;
			FlagCaraInferiorTextura=1;
			FlagDibujaLineas=0;

			break;
		case 5:
			MODO_CampoVelocidades=1;
			nParticulas=5000;
			primerdrawVelGL=1;
			FlagMuestraCaraSuperior=0;
			FlagCaraSuperiorTransparente=1;
			FlagCaraSuperiorTextura=1;
			FlagCaraInferiorTextura=1;
			FlagDibujaLineas=0;

			break;

		}
		if (glui != NULL) {
			glui->sync_live();
		}

	}



	else if (control>100 && control <10000) {
		FromCB=1;
		CB_keyboard(control-100);
	}

}
#endif

void MainCalculoContinuo() {

	time_t _tm =time(NULL );

	struct tm * curtime = localtime ( &_tm );


	sprintf(nombre0,"/Soluciones/Sal_%02d.%02d_%02d.%02d.%02d_Malla=%d_nR=%d_nTg=%d_NDivZ=%d_",
			curtime->tm_mon+1,curtime->tm_mday,
			curtime->tm_hour,curtime->tm_min,curtime->tm_sec,CualMalla,nR,nTh,NDivZ);


	sprintf(nombre0,"/Soluciones/Sal_%02d.%02d_%02d.%02d.%02d_",
			curtime->tm_mon+1,curtime->tm_mday,
			curtime->tm_hour,curtime->tm_min,curtime->tm_sec);


	char nombre[100];
	sprintf(nombre,"%s_t_porcentaje.dat",nombre0);


	mkdir("/Soluciones", 0755 );
	rmdir(nombre);


	myfileVol.open (nombre);
	sprintf(nombre,"%s_Salida.dat",nombre0);
	rmdir(nombre);
	myfileSalida.open (nombre);

	myfileSalida<<mensajes0;


	sprintf(nombre,"%s_Etapa=%02d_Parametros.dat",nombre0,0);
	rmdir(nombre);
	ofstream myfile;
	myfile.open (nombre);


	LecturaArchivoDeDatos_Post(myfile);
	myfile.close();

	Calculo_EtapaS();
	while(FinEtapas==0) {

#if DBGMain
		cout<<"\n-----------------------------------------------------------------------"
				<<"\nmain():[1080] Guardar="<<Guardar<<"\tEtapa"<<Etapa<<endl;
#endif
		Calculo_EtapaS();
#if DBGMain
		cout<<"\nmain():[1082] Guardar="<<Guardar<<"\tEtapa"<<Etapa<<endl;
#endif
		if (Guardar) {

			GuardarInstantanea();
		}
	}

}
void LecturaArchivosBMP_para_Textura ()
{
	////////////////BEGIN BMP

	ifstream mybmpfilespointer;
	///Imagen 1
	mybmpfilespointer.open ("mar.bmp", std::ios::in | std::ios::binary);
	mybmpfilespointer.read(headerinfo, 54); // read the 54-byte header size_t fread ( void * ptr, size_t size, size_t count, FILE * stream );

	// extract image height and width from header
	imagenwidth = *(int*)&headerinfo[18];
	imagenheight = *(int*)&headerinfo[22];

	int size = 3 * imagenwidth * imagenheight;
	imagendata = new char[size+1]; // allocate 3 bytes per pixel
	mybmpfilespointer.read(imagendata,  size); // read the rest of the imagesdata at once
	mybmpfilespointer.close();


	///BGR --> RGB
	for (int i = 0; i < size; i += 3)
	{
		byte dummy = imagendata[i];
		imagendata[i] = imagendata[i + 2];
		imagendata[i + 2] = dummy;
	}


	///Imagen 2
	mybmpfilespointer.open ("arena2.bmp", std::ios::in | std::ios::binary);
	mybmpfilespointer.read(headerinfo, 54); // read the 54-byte header size_t fread ( void * ptr, size_t size, size_t count, FILE * stream );

	// extract image height and width from header
	imagen2width = *(int*)&headerinfo[18];
	imagen2height = *(int*)&headerinfo[22];

	size = 3 * imagen2width * imagen2height;
	imagen2data = new char[size+1]; // allocate 3 bytes per pixel
	mybmpfilespointer.read(imagen2data,  size); // read the rest of the imagesdata at once
	mybmpfilespointer.close();

	///BGR --> RGB
	for (int i = 0; i < size; i += 3)
	{
		byte dummy = imagen2data[i];
		imagen2data[i] = imagen2data[i + 2];
		imagen2data[i + 2] = dummy;
	}

	// display image height and width from header
	cout << "image width:" << imagenwidth << endl;
	cout << "image height:" << imagenheight << endl;

#if 0
	// Guardar un archivo BMP con los datos leidos (como verificacion)
	ofstream arrayfile("bmpofstream.bmp", std::ios::out | std::ios::binary); // File Creation

	arrayfile.write(headerinfo, 54);
	arrayfile.write(imagendata, size);
	arrayfile.close();
#endif

	////TEXTURE


	////////////////////////END BMP
	///
	///
}

void LecturaLineaDeComandos(int & argc,char **argv)
{

	//------------------Linea de comandos
	strcpy(FDatos,"Datos.txt");
	int i;
	if (PrintInicio) cout<<"argc="<<argc<<endl;
	if (PrintInicio) cout<<"argv[0]="<<argv[0]<<endl;
	if (PrintInicio) cout<<"_getcwd="<<getcwd<<endl;

	char cCurrentPath[FILENAME_MAX];


	if (!getcwd(cCurrentPath, sizeof(cCurrentPath)))
	{
		exit(errno);
	}

	cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

	if (PrintInicio) printf ("The current working directory is %s\n", cCurrentPath);

	sprintf(mensajes0,"argc=%d\n",argc);
	sprintf(mensajes0,"argv[0]=%s\n",argv[0]);
	if (argc<=1) LecturaArchivoDeDatos();
	for (i=1;i<argc;i++) {
#if DBG>0
		cout<<"argv[i]="<<argv[i]<<endl;
#endif

		if (strcmp(argv[i],"-Datos")==0) {
			strcpy(FDatos,argv[i+1]);i++;
			LecturaArchivoDeDatos();
		}
		else if (strcmp(argv[i],"-malla")==0     		) {CualMalla			=atoi(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-NDivZ")==0			) {NDivZ				=atoi(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-CalculoContinuo")==0	) {CalculoContinuo=1;}
		else if (strcmp(argv[i],"-TipoCalculo")==0	    ) {TipoCalculo          =atoi(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-dt")==0				) {Datos_dt				=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-Tmax")==0				) {Datos_Tmax			=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-nR")==0				) {nR					=atoi(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-nR2")==0				) {nR2					=atoi(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-nTh")==0				) {nTh					=atoi(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-Dominio_Xmax")==0		) {Dominio_Xmax			=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-Dominio_Hmax")==0		) {Dominio_Hmax			=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-Dominio_Hsup")==0		) {Dominio_Hsup			=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-Dominio_Rint")==0		) {Dominio_Rint			=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-Dominio_Rmed")==0		) {Dominio_Rmed			=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-rhof")==0				) {Datos_rhof			=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-rhos")==0				) {Datos_rhos			=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-phi")	==0				) {Datos_phi			=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-cf")==0				) {Datos_cf				=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-cs")==0				) {Datos_cs				=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-Tinyeccion")==0		) {Datos_Tinyeccion		=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-Tambiente")==0		) {Datos_Tambiente		=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-TLimite")==0	    	) {TLimite	            =atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-hc_superior")==0		) {Datos_hc_superior	=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-hc_inferior")==0		) {Datos_hc_inferior	=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-hc_lateral")==0		) {Datos_hc_lateral		=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-KTermofilm")==0		) {Datos_KTermofilm		=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-eTermofilm")==0		) {Datos_eTermofilm		=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-VinyeccionLtsHr")==0	) {VinyeccionLtsHr		=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-kf")==0				) {Datos_kf				=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-ks")==0				) {Datos_ks				=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-Kevaporacion")==0		) {Datos_Kevaporacion	=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-HumedadAmbiental")==0	) {Datos_HumedadAmbiental=atof(argv[i+1]);i++;}
		else if (strcmp(argv[i],"-DistanciaAlBorde")==0	) {Datos_DistanciaAlBorde=atof(argv[i+1]);i++;}




	}

}


void AddMenuPrincipal() {
	MainMenu=glutCreateMenu(MyMenu);
	glutAddMenuEntry("1) Leer gtotal_voronio.msh ", 101);
	glutAddMenuEntry("2) Define Mancha circular ", 102);
	glutAddMenuEntry("3) Lanza cálculo conveccion/Difusion ", 103);
	glutAddMenuEntry("-------", 0);
	glutAddMenuEntry("Open File .poly", 2);
	glutAddMenuEntry("Open File .xml", 3);
	glutAddMenuEntry("Lee Velocidad.vtu", 6);
	glutAddMenuEntry("Select Polygono (on/Off)", 4);
	glutAddMenuEntry("Muestra BC", 5);
	glutAttachMenu(GLUT_RIGHT_BUTTON);EstadoMenu=1;

}


void MyMenu(int value)
{
	cout<<"MyMenu(int value)"<<value<<endl;
	switch (value) {
	case 101:
		Etapa_LecturaMalla_gtotal_Voronoi_MSH();
		break;
	case 102:
		Etapa_Define_Mancha_circular_inicial();
		break;
	case 103:
		Lanza_thread();
		break;
	case 2:
		cout<<"Open File .poly"<<endl;
		FileBrowserGlui2->boton_ok->user_id=18;
		FileBrowserGlui2->list2->user_id=18;
		FileBrowserGlui2->user_id=18;
		FileBrowserGlui2->fbreaddir("*.poly");
		glui2->show();
		break;
	case 3:
		cout<<"Open File .xml"<<endl;
		FileBrowserGlui2->boton_ok->user_id=19;
		FileBrowserGlui2->list2->user_id=19;
		FileBrowserGlui2->user_id=19;
		FileBrowserGlui2->fbreaddir("*.xml");
		glui2->show();
		//Se llamará a:		Lee_Malla_XML(file_name);
		break;
	case 4:
		//glutAddMenuEntry("Select Polygono (on/Off)", 4);
		Add_Poligono_Selected=1-Add_Poligono_Selected;
		cout <<"Add_Poligono_Selected="<<Add_Poligono_Selected<<endl;


		break;
	case 5:
		formulario_BC_addDatos();
		gluiEligeBC->show();gluiEligeBC_Visible=1;
		break;
	case 6:
		cout<<"Lee Velocidad.vtu"<<endl;
		FileBrowserGlui2->boton_ok->user_id=2210;
		FileBrowserGlui2->list2->user_id=2210;
		FileBrowserGlui2->user_id=2210;
		FileBrowserGlui2->fbreaddir("v*.vtu");
		glui2->show();
		break;
	}

}


int main(int argc,char **argv)
{

	omp_set_nested(1);
	omp_set_dynamic(0);
	glutInit(&argc, argv);

	//gtotal=new grid3D;

	LecturaLineaDeComandos(argc, argv);
	LecturaArchivosBMP_para_Textura();

	//cout<<"nR2="<<nR2<<endl;

	myfileDBG.open ("DBG");

	if (CalculoContinuo==1) {

		MainCalculoContinuo();


	} else {
		//---------------------------

		LecturaArchivoDeDatos_Post(cout);
		int anchoT,altoT;

		anchoT=glutGet(GLUT_SCREEN_WIDTH);if (anchoT>2000) anchoT=1980;

		altoT=glutGet(GLUT_SCREEN_HEIGHT);

		//anchoT=1024;
		//altoT=768;

		if (PrintInicio) { printf("ancho=%d, alto=%d\n",anchoT,altoT);cout<<endl; }

		width=anchoT-10*0;
		height=altoT-40*0;
		glutInitWindowSize(width, height);

		glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );

		main_window =glutCreateWindow("VersionGLUT-cubo");

		glutPositionWindow(glutGet(GLUT_SCREEN_WIDTH)-width-5, glutGet(GLUT_SCREEN_HEIGHT)-height-5);

#if DBG==1
		cout<<"DBG: inicializacion"<<endl;
#endif
		inicializacion(); ///Se leen los datos físicos del problema


		if (PrintInicio) { cout<<"ok inicializacion"<<endl; }
#if FREEGLUT
		glutMouseWheelFunc(CB_mouse_whell);
#endif
		glutMouseFunc(CB_mouse);
		glutMotionFunc(CB_motion);
		glutSpecialFunc( CB_keyboardSpecial );



		glutVisibilityFunc(visible);
		//  glutKeyboardUpFunc(keyup);
		//  glutSpecialUpFunc(specialup);

#if DBG==1
		cout<<"DBG: AddMenuPrincipal"<<endl;
#endif
		AddMenuPrincipal();

		if (PrintInicio) { cout<<"ok AddMenuPrincipal"<<endl; }
#if defined(GLUI_GLUI_H)
		GLUI_Master.set_glutReshapeFunc( ResizeGraphics );
		//GLUI_Master.set_glutKeyboardFunc(CB_keyboard);
		glutKeyboardFunc(CB_keyboard);

		if (DBG) { cout<<"voy a formulario_glui()"<<endl; }
		formulario_glui();

		if (DBG) { cout<<"voy a formulario_BC()"<<endl; }
		formulario_BC();
		//		Level0_formulario_OpcionesPrograma();
		if (DBG) { cout<<"ok formulario_BC"<<endl; }

#else
		glutReshapeFunc(ResizeGraphics);
		glutKeyboardFunc(CB_keyboard);
#endif

		//	glutSetCursor(GLUT_CURSOR_CYCLE);

		glutDisplayFunc(DrawGraphics);
#if DBG==1
		cout<<"Inicio4"<<endl;
#endif


		glutMainLoop();
		printf("hola");

	}


	exit(0);
}

void GuardarInstantanea() {
	char nombre[100];
	int tmp=InicioEtapa;
	InicioEtapa=EtapaGlobal2Local;
	if (TipoCalculo==CalculoEstacionario) {
		sprintf(nombre,"%s_Etapa=%02d_Size=%d.dat",nombre0,EtapaGlobal,gtotal->nH3D);
	} else {
		sprintf(nombre,"%s_Etapa=%02d_Size=%d_t=%.0f.dat",nombre0,EtapaGlobal,gtotal->nH3D,TiempoCalculo);
	}
	rmdir(nombre);


	cout<<"Save:"<<nombre<<endl;
	myfileSalida<<"Save:"<<nombre<<endl;

#if DBGMain
	cout<<"EtapaGlobal2Local="<<EtapaGlobal2Local<<endl;

#endif
	SaveOrRead(nombre,1);


	InicioEtapa=tmp;
}



#if defined(GLUI_GLUI_H)

void Formulario_OpenFile ()
{


	if (DBG) cout<<"Formulario_OpenFile ()"<<endl;
	glui2 = GLUI_Master.create_glui("GLUI Window",GLUI_SUBWINDOW_RIGHT);

	cout<<"glui2->glui_id="<<glui2->get_glut_window_id()<<endl;



	GLUI_Panel *panel2 = glui2->add_panel( "" );
	//		panel2->set_w(300);
	//		panel2->set_h(500);

	new GLUI_StaticText(panel2,"Open File:");
	FileBrowserGlui2 = new GLUI_FileBrowser2(panel2, "", false, 17, control_cb_file);
	//FileBrowserGlui2 = new GLUI_FileBrowser(glui2, "", GLUI_PANEL_EMBOSSED, 7, control_cb);
	FileBrowserGlui2->set_w(300);
	FileBrowserGlui2->set_h(500);
	FileBrowserGlui2->set_allow_change_dir(1);
	glui2->add_column_to_panel(panel2,false);
	//		fb2 = new GLUI_List(panel2, "", GLUI_PANEL_EMBOSSED, 7, control_cb);
	glui2->add_column(false);

	FileBrowserGlui2->fbreaddir("*.txt");
	//		fbreaddir2(FileBrowserGlui2, "*.txt");
	//		glutPositionWindow(glutGet(GLUT_SCREEN_WIDTH)-10, glutGet(GLUT_SCREEN_HEIGHT)-10);
	glui2->hide();

	if (DBG) cout<<"END:Formulario_OpenFile ()"<<endl;
}

/************************************************ GLUI::GLUI() **********/

void   formulario_glui()
{
	GLUI_Spinner *tsp;
	/////GLUI:
	//	glui = GLUI_Master.create_glui( "GLUI", 0, 90, 90 ); // name, flags, x, and y

	glui = GLUI_Master.create_glui_subwindow(main_window,GLUI_SUBWINDOW_RIGHT);
	glui->set_main_gfx_window( main_window);


	GLUI_Panel *mipanel=glui->add_panel("Clipping");

	gluiClipping = glui->add_checkbox_to_panel(mipanel,"'P'= Clipping",&ClippingON, 'P'+100 ,control_cb );

	tsp=glui->add_spinner_to_panel( mipanel, "Ax:" ,GLUI_SPINNER_FLOAT, &Ax );
	tsp->set_alignment(GLUI_ALIGN_LEFT);	tsp->w=10;
	if (PrintInicio) cout<<"Antes set_w ....."<<endl;
	tsp=glui->add_spinner_to_panel( mipanel, "+By:",GLUI_SPINNER_FLOAT, &By );
	tsp->set_alignment(GLUI_ALIGN_LEFT);	tsp->w=10;
	tsp=glui->add_spinner_to_panel( mipanel, "+Cz:",GLUI_SPINNER_FLOAT, &Cz );
	tsp->set_alignment(GLUI_ALIGN_LEFT);	tsp->w=10;
	tsp=glui->add_spinner_to_panel( mipanel, " =D:" ,GLUI_SPINNER_FLOAT, &DD );
	tsp->set_alignment(GLUI_ALIGN_LEFT);	tsp->w=10;
	if (PrintInicio)  cout<<"Despues set_w"<<endl;
	GLUI_Panel *ptmp2;
#if 0
	PanelFlimite = glui->add_panel("",GLUI_PANEL_EMBOSSED);
	ptmp2 = glui->add_panel_to_panel(PanelFlimite,"",GLUI_PANEL_NONE);
	//	ptmp->set_alignment(GLUI_ALIGN_LEFT);
	glui->add_checkbox_to_panel(ptmp2,"TLimite",&TLimite_if);
	glui->add_column_to_panel(ptmp2,false);
	tsp=glui->add_spinner_to_panel(ptmp2,"" ,GLUI_SPINNER_FLOAT, &TLimite , 1105, control_cb );


	glui_porcentaje=glui->add_statictext_to_panel(PanelFlimite,"");

	PanelFlimite->disable();
#endif


#if 0
	ptmp2 = glui->add_panel("",GLUI_PANEL_NONE);
	glui->add_button_to_panel(ptmp2,"Save", 1002 ,control_cb );
	glui->add_column_to_panel(ptmp2,false);
	glui->add_column_to_panel(ptmp2,false);
	glui->add_button_to_panel(ptmp2,"Read", 1003 ,control_cb );
#endif

	//	gluiMueve		=glui->add_checkbox("[M] Mueve centro",&MueveCentro);
	//	gluiNumera		=glui->add_checkbox("'N' Numera Vertices(on/off)",NULL, 'N'+100 ,control_cb );
	//	glui->add_checkbox("[ ] Ejes ",&MODO_Ejes, 1102 ,control_cb);
	glui->add_checkbox("Vel ",&MODO_CampoVelocidades2);
	glui->add_checkbox("(Centro) ",&DibujaCentroVel);

	//	gluiNormales		=glui->add_checkbox("[N] Normales ",&ModoDibujaNormales);
	//	gluiInterior		=glui->add_checkbox("[I] Interior ",&ModoDibujaInterior);
	//	gluiBordes		=glui->add_checkbox("[B] Bordes ",&ModoDibujaFrontera);
	//glui->add_checkbox("[ ] centros ",&Modo_DibujaCentroBloques);

	ptmp2 = glui->add_panel("",GLUI_PANEL_RAISED);
	glui->add_checkbox_to_panel(ptmp2, "",&ModoDibujaAlgunos, 11030 ,control_cb);
	glui->add_column_to_panel(ptmp2,false);
	Spinner_Algunos =tsp=glui->add_spinner_to_panel(ptmp2, "%",GLUI_SPINNER_INT, &Algunos_Porcentaje,1103,control_cb);

	tsp->set_int_limits(0,100);

	if (DBG) cout<<"linea 2007"<<endl;

	glui->add_statictext( "" );

	if (DBG) cout<<"linea 2011"<<endl;

	form_TesteDeVariablesGlobales();

	if (DBG) cout<<"linea 2015"<<endl;

#if 0
	glui->add_statictext( "" );
	glui->add_button("Etapa Siguiente [Space]", 1001 ,control_cb );
	glui->add_statictext( "" );
#endif

	PanelParticulas = glui->add_panel("",GLUI_PANEL_EMBOSSED);
	Checkbox_particulas=glui->add_checkbox_to_panel(PanelParticulas,"Particulas",&MODO_CampoVelocidades, 1101 ,control_cb );
	glui->add_checkbox_to_panel(PanelParticulas,"Origen",&MODO_Origen);
	glui->add_checkbox_to_panel(PanelParticulas,"Pausa",&MODO_Pausa);


	if (DBG) cout<<"linea 2029"<<endl;

	Spinner_particulas =tsp=glui->add_spinner_to_panel(PanelParticulas,"N Particulas",GLUI_SPINNER_INT, &nParticulas,2223,control_cb);
	tsp->set_int_limits(0,20000);

	//	tsp=glui->add_spinner("Factor[0,1]",GLUI_SPINNER_FLOAT, &FactorCercania );
	//	tsp->set_int_limits(0,1);
	tsp=glui->add_spinner_to_panel(PanelParticulas,"Largo=",GLUI_SPINNER_INT, &largoCola );
	tsp->set_int_limits( 0, maxpasadas);
	//	tsp->edittext->set_w(0);




	glui->add_statictext( "" );
	//	glui->add_statictext( "");
	glui_FPS=glui->add_statictext( "FPS" );
	glui_FPS->set_alignment(GLUI_ALIGN_CENTER);
	glui->add_button("Reset", 8001 ,control_cb );
	//	glui_MSG=glui->add_statictext( "" );
	//	glui_MSG->set_alignment(GLUI_ALIGN_LEFT);


	glui_edittext= new GLUI_TextBox(glui,text);
	glui_edittext->set_w(180);
	glui_edittext->set_h(120);


	mipanel=glui->add_panel("Uso del Boton 1 del Mouse");
	glui_GrupoModoDelMouse=glui->add_radiogroup_to_panel(mipanel,&MODO_de_mover);
	glui->add_radiobutton_to_group(glui_GrupoModoDelMouse,"[t]=Trasladar");
	glui->add_radiobutton_to_group(glui_GrupoModoDelMouse,"[r]=Rotar");
	glui->add_radiobutton_to_group(glui_GrupoModoDelMouse,"[s]=Spin Z");
	glui->add_radiobutton_to_group(glui_GrupoModoDelMouse,"[e]=Escala");


	//	glui->add_button("'R'=Reset View", 'R'+100 ,control_cb );  //Eliminado en Version1
	glui->add_button("'F'=Full(on/off)", 'F'+100 ,control_cb );
	//	glui->add_button("'G'=Game", 'g'+100 ,control_cb );
	glui->add_button("'C'=Color", 'C'+100 ,control_cb );


	glui->add_button("[Z]=Hide", 'Z'+100,control_cb );


	gluiHelp		=glui->add_checkbox("[V] Verbose ",&MODO_MenuMENSAJES); //Eliminado en Version 1
	glui->add_checkbox("[ ] NumH ",&MODO_NumeraH);
	//	glui->add_checkbox("[ ] NumFF ",&MODO_NumeraFF); //Eliminado en Version 1

#if 1==1
	if (1==1) { // Formulario para abrir archivo


		if (DBG) cout<<"voy a Formulario_OpenFile"<<endl;
		Formulario_OpenFile ();


	}
#endif

	if (1==1) { // Formulario para abrir archivo


		if (DBG) cout<<"BEGIN:Formulario para abrir archivos (glui3)"<<endl;
		glui3 = GLUI_Master.create_glui("GLUI Open File",GLUI_SUBWINDOW_RIGHT,1,1);
		/*
		/// Codigo equivalente:
//		GLUI *GLUI_Master_Object::create_glui( const char *name, long flags,int x,int y )
//		{
		  GLUI *new_glui = new GLUI;
			if (DBG) cout<<"BEGIN:Formulario para abrir archivos (glui3)+1"<<endl;
		  new_glui->init("GLUI Open File",GLUI_SUBWINDOW_RIGHT, 10, 10, -1 );
			if (DBG) cout<<"BEGIN:Formulario para abrir archivos (glui3)+2"<<endl;
		  new_glui->link_this_to_parent_last( glui2 );
			if (DBG) cout<<"BEGIN:Formulario para abrir archivos (glui3)+3"<<endl;
		  glui3=new_glui;
			if (DBG) cout<<"BEGIN:Formulario para abrir archivos (glui3)+4"<<endl;
//		  return new_glui;
//		}

		 */
		if (DBG) cout<<"Formulario para abrir archivo: 2562"<<endl;

		GLUI_Panel *panel3 = glui3->add_panel( "P" );


		if (DBG) cout<<"Formulario para abrir archivo: 2571"<<endl;

		//		glui3->add_edittext_to_panel(panel3, "CalculoContinuo", GLUI_EDITTEXT_INT, &E_CalculoContinuo);
		Glui3_EditParametros= new GLUI_TextBox(panel3, StringParametros);

		Glui3_EditParametros->set_w(320);
		Glui3_EditParametros->set_h(520);
		//		Glui3_EditParametros->set_font( GLUT_BITMAP_8_BY_13);


		GLUI_Panel *panel3b = glui3->add_panel( "" );
		panel3b->set_w(320);
		glui3->add_button_to_panel(panel3b,"Cambia", 10001 ,control_cb )->set_alignment(GLUI_ALIGN_LEFT);
		glui3->add_column_to_panel(panel3b,false);
		glui3->add_statictext_to_panel(panel3b,"       ");
		glui3->add_column_to_panel(panel3b,false);
		glui3->add_button_to_panel(panel3b,"Hide", 10002 ,control_cb )->set_alignment(GLUI_ALIGN_RIGHT);

		glui3 ->hide();
		if (DBG) cout<<"END:Formulario para abrir archivo"<<endl;
	}


	if (DBG) cout<<"Voy a formulario_opciones_programa()"<<endl;
	formulario_opciones_programa();

	if (1==1) { // Formulario para mostrar Help
		gluiHelp2 = GLUI_Master.create_glui("GLUI Help",GLUI_SUBWINDOW_RIGHT);



		//		glui = GLUI_Master.create_glui_subwindow(main_window,GLUI_SUBWINDOW_RIGHT);
		//		glui->set_main_gfx_window( main_window);




		GLUI_Panel *panelH1 = gluiHelp2->add_panel( "" );

		TextBoxHelp= new GLUI_TextBox(panelH1,textHelp);
		TextBoxHelp->set_w(320);
		TextBoxHelp->set_h(520);




		GLUI_Panel *panelH2 = gluiHelp2->add_panel( "" );
		panelH2->set_w(320);
		gluiHelp2->add_statictext_to_panel(panelH2,"       ");
		gluiHelp2->add_column_to_panel(panelH2,false);
		gluiHelp2->add_button_to_panel(panelH2,"Hide", 20001 ,control_cb )->set_alignment(GLUI_ALIGN_RIGHT);

		gluiHelp2 ->hide();
	}
	//	glutIdleFunc( idleevent );
	GLUI_Master.set_glutIdleFunc( idleevent );
	//	glui_hide=false;
	//	master->hide(); //Borrar pantalla principal para test


	/*
	GLUI_Control *gluiNode=(GLUI_Control *)glui->main_panel->first_child();
	while (gluiNode!=NULL) {
		cout << "gluiNode->name="<<gluiNode->name<<endl;;
		gluiNode->set_h( gluiNode->h/2);
		gluiNode=(GLUI_Control *)gluiNode->next();
	}
	 */



	//Checkbox_particulas->deactivate();
	//Checkbox_particulas->disable();

	//PanelParticulas->disable();
	if (PrintFunciones) cout<<"formulario_glui()end"<<endl;
}

void Level0_formulario_OpcionesPrograma()
{

	glui_OpcionesProgramaLevel0 = GLUI_Master.create_glui_subwindow(main_window,GLUI_SUBWINDOW_LEFT);

	//	glui_OpcionesProgramaLevel0 = GLUI_Master.create_glui_subwindow(main_window,GLUI_SUBWINDOW_TOP);



	GLUI_Panel *panelH2 = glui_OpcionesProgramaLevel0->add_panel( "" );
	panelH2->set_w(320);
	//	gluiEligeBC->add_statictext_to_panel(panelH2,"       ");
	//	gluiEligeBC->add_column_to_panel(panelH2,false);
	GLUI_Rollout * RO_Level0 = glui_OpcionesProgramaLevel0->add_rollout_to_panel(panelH2,"Level0", 1, GLUI_PANEL_EMBOSSED);

	glui_OpcionesProgramaLevel0->add_button_to_panel(RO_Level0, "LeeXML", 100, Level0_control_cb);


	glui_OpcionesProgramaLevel0 ->hide();

}
void formulario_BC_addDatos()
{
	int i,nMax;
	char s[100];
	GLUI_Node *primero,*otro;
	GLUI *GLUI_TMP;

	if (gPoly != NULL) {
		if (gluiEligeBCpanel1!=NULL) {
			primero=gluiEligeBCpanel1->first_child() ;
			while (primero!=NULL) {
				otro=primero->next();
				primero->unlink();
				primero=otro;
			}
		}

		//		cout<<"gluiEligeBCpanel1->first_child() ="<<gluiEligeBCpanel1->first_child() <<endl;;
		//		cout<<"gPoly->nBC="<<gPoly->nBC<<endl;

		nMax=50;
		for (i=0;i<gPoly->nBC;i++) {
			sprintf(s,"[%d]",gPoly->BC[i]);
			cout<<"i="<<i<<" s="<<s<<endl;
			if (i==0)
				gluiEligeBC->add_checkbox_to_panel(gluiEligeBCpanel1,s,&(gPoly->DrawBC[i]), 30002,control_cb_BC);
			else
				gluiEligeBC->add_checkbox_to_panel(gluiEligeBCpanel1,s,&(gPoly->DrawBC[i]), 30002,control_cb_BC);

			cout<<"gluiEligeBCpanel1->first_child() ="<<gluiEligeBCpanel1->first_child() <<endl;
			if (i==nMax)
				gluiEligeBC->add_column_to_panel(gluiEligeBCpanel1,false);
			if (i==2*nMax)
				gluiEligeBC->add_column_to_panel(gluiEligeBCpanel1,false);
		}

	}

	if (gtotal != NULL) {
		if (gluiEligeBCpanel1!=NULL) {
			primero=gluiEligeBCpanel1->first_child() ;
			while (primero!=NULL) {
				otro=primero->next();
				primero->unlink();
				primero=otro;
			}
		}

		//		cout<<"gluiEligeBCpanel1->first_child() ="<<gluiEligeBCpanel1->first_child() <<endl;;
		//		cout<<"gPoly->nBC="<<gPoly->nBC<<endl;

		nMax=50;
		for (i=0;i<gtotal->nBC;i++) {
			sprintf(s,"[%d]",gtotal->BC[i]);
			cout<<"i="<<i<<" s="<<s<<endl;
			if (i==0)
				gluiEligeBC->add_checkbox_to_panel(gluiEligeBCpanel1,s,&(gtotal->DrawBC[i]), 30002,control_cb_BC);
			else
				gluiEligeBC->add_checkbox_to_panel(gluiEligeBCpanel1,s,&(gtotal->DrawBC[i]), 30002,control_cb_BC);

			cout<<"gluiEligeBCpanel1->first_child() ="<<gluiEligeBCpanel1->first_child() <<endl;
			if (i==nMax)
				gluiEligeBC->add_column_to_panel(gluiEligeBCpanel1,false);
			if (i==2*nMax)
				gluiEligeBC->add_column_to_panel(gluiEligeBCpanel1,false);
		}

	}

}

void formulario_BC()
{

	// Formulario para mostrar Help
	//	gluiEligeBC = GLUI_Master.create_glui("Elija la BC a mostrar",GLUI_SUBWINDOW_RIGHT);
	gluiEligeBC = GLUI_Master.create_glui_subwindow(main_window,GLUI_SUBWINDOW_LEFT);




	GLUI_Panel *panelH2 = gluiEligeBC->add_panel( "" );
	//	gluiEligeBC->add_statictext_to_panel(panelH2,"       ");
	//	gluiEligeBC->add_column_to_panel(panelH2,false);
	gluiEligeBC->add_button_to_panel(panelH2,"Hide", 1 ,control_cb_BC )->set_alignment(GLUI_ALIGN_RIGHT);

	cout<<"gluiEligeBC->next()="<<gluiEligeBC->next()<<endl;;
	cout<<"panelH2->next()="<<panelH2->next()<<endl;;

	panelH2->set_w(130);
	gluiEligeBC_w=130;

	gluiEligeBCpanel1 = gluiEligeBC->add_panel( "BC..." );

	gluiEligeBC->hide();gluiEligeBC_Visible=0;


}


void formulario_opciones_programa()
{

	if (DBG) cout<<"formulario_opciones_programa()"<<endl;

	// Formulario para mostrar Help
	gluiOpciones = GLUI_Master.create_glui("GLUI Opciones",GLUI_SUBWINDOW_RIGHT);


	GLUI_Panel *panel1 = gluiOpciones->add_panel( "Opciones" );


	gluiOpciones->add_checkbox_to_panel(panel1,"Print Mouse",&PrintMouse, 30002,control_cb);
	gluiOpciones->add_checkbox_to_panel(panel1,"Print Linea",&PrintLinea, 30002,control_cb);
	gluiOpciones->add_checkbox_to_panel(panel1,"Print Tiempos",&PrintTiempos, 30002,control_cb);
	gluiOpciones->add_checkbox_to_panel(panel1,"Print Caras",&PrintCaras, 30002,control_cb);
	gluiOpciones->add_checkbox_to_panel(panel1,"Print Invert",&PrintInvert, 30002,control_cb);



	GLUI_Panel *panelH2 = gluiOpciones->add_panel( "" );
	panelH2->set_w(320);
	gluiOpciones->add_statictext_to_panel(panelH2,"       ");
	gluiOpciones->add_column_to_panel(panelH2,false);
	gluiOpciones->add_button_to_panel(panelH2,"Hide", 30001 ,control_cb )->set_alignment(GLUI_ALIGN_RIGHT);

	gluiOpciones ->hide();
	if (DBG) cout<<"END:formulario_opciones_programa()"<<endl;

}

void form_TesteDeVariablesGlobales() {

	GLUI_Spinner *tsp;
	GLUI_Panel *ptmp2;



	if (DBG) cout<<"TesteDeVariablesGlobales()"<<endl;
	ptmp2 = glui->add_panel("",GLUI_PANEL_RAISED);
	//glui->add_checkbox_to_panel(ptmp2,"TSup",&FlagMuestraCaraSuperior);
	//glui->add_checkbox_to_panel(ptmp2,"TInf",&FlagMuestraCaraInferior);
	glui->add_checkbox_to_panel(ptmp2,"Grid",&FlagDibujaLineas);

#if 0
	if (DBG) cout<<"linea 2242"<<endl;

	glui->add_column_to_panel(ptmp2,false);
	GLUI_RadioGroup *RG=glui->add_radiogroup_to_panel(ptmp2, &ListaCaso,25000,control_cb);
	glui->add_radiobutton_to_group(RG, "0");
	glui->add_radiobutton_to_group(RG, "1");
	glui->add_radiobutton_to_group(RG, "2");
	glui->add_radiobutton_to_group(RG, "3");
	glui->add_radiobutton_to_group(RG, "4");
	glui->add_radiobutton_to_group(RG, "5");
	control_cb(25000);
#endif

#if 1
	///// Test de algunas Variables
	tsp=glui->add_spinner("FactorNormales",GLUI_SPINNER_FLOAT, &FactorNormales);
	tsp->set_float_limits(0,1);
	tsp=glui->add_spinner("FactorSuavidad",GLUI_SPINNER_FLOAT, &FactorSuavidad);
	tsp->set_float_limits(0,1);
	tsp=glui->add_spinner("FactorAchica",GLUI_SPINNER_FLOAT, &FactorAchica);
	tsp->set_float_limits(0,1);
	//	tsp=glui->add_spinner("FactorAchicaV",GLUI_SPINNER_FLOAT, &FactorAchicaV);
	//	tsp->set_float_limits(0,1);
	tsp=glui->add_spinner("FactorZ",GLUI_SPINNER_FLOAT, &FactorZ);
	tsp->set_float_limits(1e-5,1e5);
	tsp->set_speed(0.1);
	tsp=glui->add_spinner("FactorV",GLUI_SPINNER_FLOAT, &factorV);
	tsp->set_float_limits(0,1000000);
	tsp->set_speed(0.1);
	tsp=glui->add_spinner("FactorVmin",GLUI_SPINNER_FLOAT, &factorVmin);
	tsp->set_float_limits(0,1);
	tsp=glui->add_spinner("AnimaV",GLUI_SPINNER_INT, &Global_nia);
	tsp->set_float_limits(1,10);
	//	glui->add_spinner("Factor Vh",GLUI_SPINNER_FLOAT, &factorDespegaLineas );

#endif


#if 0
	tsp=glui->add_spinner("BGColorR",GLUI_SPINNER_INT, &BGColorR);
	tsp->set_float_limits(0,255);
	tsp=glui->add_spinner("BGColorG",GLUI_SPINNER_INT, &BGColorG);
	tsp->set_float_limits(0,255);
	tsp=glui->add_spinner("BGColorB",GLUI_SPINNER_INT, &BGColorB);
	tsp->set_float_limits(0,255);
#endif


#if 0
	//Test Materiales
	tsp=glui->add_spinner("Ambient",GLUI_SPINNER_FLOAT, &FactorAmbient);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.1);

	tsp=glui->add_spinner("Difusse",GLUI_SPINNER_FLOAT, &FactorDifusse);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.1);

	tsp=glui->add_spinner("Specular",GLUI_SPINNER_FLOAT, &FactorSpecular);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.1);

	tsp=glui->add_spinner("Emission",GLUI_SPINNER_FLOAT, &FactorEmission);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.1);

#endif


	tsp=glui->add_spinner("Claridad",GLUI_SPINNER_FLOAT, &FactorClaridad);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.1);

#if 0
	//Test Luces
	tsp=glui->add_spinner("Luz1 Ambient",GLUI_SPINNER_FLOAT, &FactorAmbientL1);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.2);

	tsp=glui->add_spinner("Luz1 Difusse",GLUI_SPINNER_FLOAT, &FactorDifusseL1);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.2);

	tsp=glui->add_spinner("Luz1 Specular",GLUI_SPINNER_FLOAT, &FactorSpecularL1);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.2);


	tsp=glui->add_spinner("Luz2 Ambient",GLUI_SPINNER_FLOAT, &FactorAmbientL2);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.2);

	tsp=glui->add_spinner("Luz2 Difusse",GLUI_SPINNER_FLOAT, &FactorDifusseL2);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.2);

	tsp=glui->add_spinner("Luz2 Specular",GLUI_SPINNER_FLOAT, &FactorSpecularL2);
	tsp->set_float_limits(0,3);
	tsp->set_speed(0.2);


#endif


#if 0
	ptmp2 = glui->add_panel("",GLUI_PANEL_RAISED);
	glui->add_button_to_panel(ptmp2,"Reset centro", 9001 ,control_cb );
	glui->add_column_to_panel(ptmp2,false);
	glui->add_column_to_panel(ptmp2,false);
	glui->add_button_to_panel(ptmp2,"Reset centro2", 9001 ,control_cb );
#endif

	tsp=glui->add_spinner("Fovy",GLUI_SPINNER_FLOAT, &GlobalFovy, 9002 ,control_cb);
	tsp->set_float_limits(0.05,90);
	tsp->set_speed(1);
#if 0
	tsp=glui->add_spinner("SizeCentros",GLUI_SPINNER_FLOAT, &GlobalCentros);
	tsp->set_float_limits(0.001,1);
	tsp->set_speed(1);

	tsp=glui->add_spinner("GL_immediate_mode",GLUI_SPINNER_INT, &GL_immediate_mode);
	tsp->set_int_limits(0,4);
#endif

	tsp=glui->add_spinner("GL_threads",GLUI_SPINNER_INT, &GL_threads);
	tsp->set_int_limits(1,20);

#if 1
	tsp=glui->add_spinner("FactorTime",GLUI_SPINNER_FLOAT, &FactorClockTime);
	tsp->set_speed(1);
#endif

	if (DBG) cout<<"END: TesteDeVariablesGlobales()"<<endl;

	return;
}

#endif




#ifdef _WIN32
#include <windows.h>
#endif
void fbreaddir2(GLUI_FileBrowser *fb, const char *d) {
	GLUI_String item;
	int i = 0;


#ifdef _WIN32

	WIN32_FIND_DATAW FN;
	HANDLE hFind;
	//char search_arg[MAX_PATH], new_file_path[MAX_PATH];
	//sprintf(search_arg, "%s\\*.*", path_name);

	hFind = FindFirstFileW(L"*.*", &FN);
	if (fb->list) {
		fb->list->delete_all();
		if (hFind != INVALID_HANDLE_VALUE) {
			do {

				size_t origsize = wcslen(FN.cFileName) + 1;
				const size_t newsize = 100;
				size_t convertedChars = 0;
				char nstring[newsize];
				//			wcstombs_s(&convertedChars, nstring, origsize, FN.cFileName, _TRUNCATE);
				wcstombs(nstring,  FN.cFileName, origsize);



				int len = wcslen(FN.cFileName);
				if (FN.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
					item = '\\';
					item += nstring;
					fb->list->add_item(i,item.c_str());
				} else {
					item = nstring;
					//		  FileBrowserGlui2->list->add_item(i,item.c_str());
				}
				i++;
			} while (FindNextFileW(hFind, &FN) != 0);

			if (GetLastError() == ERROR_NO_MORE_FILES)
				FindClose(&FN);
			else
				perror("fbreaddir");
		}
	}

	hFind = FindFirstFileW(L"*.txt", &FN);
	if (hFind != INVALID_HANDLE_VALUE) {
		do {

			size_t origsize = wcslen(FN.cFileName) + 1;
			const size_t newsize = 100;
			size_t convertedChars = 0;
			char nstring[newsize];
			//			wcstombs_s(&convertedChars, nstring, origsize, FN.cFileName, _TRUNCATE);
			wcstombs(nstring,  FN.cFileName, origsize);



			int len = wcslen(FN.cFileName);
			if (FN.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
				item = '\\';
				item += nstring;
				//		  FileBrowserGlui2->list->add_item(i,item.c_str());
			} else {
				item = nstring;
				fb->list->add_item(i,item.c_str());
			}
			i++;
		} while (FindNextFileW(hFind, &FN) != 0);

		if (GetLastError() == ERROR_NO_MORE_FILES)
			FindClose(&FN);
		else
			perror("fbreaddir");
	}

#elif defined(__GNUC__)

	DIR *dir;
	struct dirent *dirp;
	struct stat dr;

	if (FileBrowserGlui2->list) {
		FileBrowserGlui2->list->delete_all();
		if ((dir = opendir(d)) == NULL)
			perror("fbreaddir:");
		else {
			while ((dirp = readdir(dir)) != NULL)   /* open directory     */
			{
				if (!lstat(dirp->d_name,&dr) && S_ISDIR(dr.st_mode)) /* dir is directory   */
					item = dirp->d_name + GLUI_String("/");
				else
					item = dirp->d_name;

				FileBrowserGlui2->list->add_item(i,item.c_str());
				i++;
			}
			closedir(dir);
		}
	}
#endif
}
