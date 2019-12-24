#pragma once

#if defined _WIN32 || defined _WIN64s
#define mkdir(x,y) _mkdir((x))
#define rmdir(x) _rmdir((x))

#endif

#include <stdio.h>
#include <sstream>
#include <signal.h>
#include <omp.h>
#include <map>
//#include <utility> // make_pair

//#include <set>
//#include <algorithm>

#include <functional>
#include <iostream>

std::map<std::string, int>  CuentaLlamadas;

std::multimap<int,std::string> dst ;


std::map<std::string, int>::iterator itCuentaLlamadas;
std::multimap<int,std::string>::iterator it_dst;

extern double ThetaMax;
extern int step_hora;
extern vector<float> ParticulasZ;
extern vector<float> ParticulasZlambda;

extern int imagenwidth;
extern int imagenheight;

void DrawMensajesRight();
void DrawMensajesCenter();
void DrawMenuOnOff();
int Iteracion_Conveccion_Difusion();

void Lee_Velocidad_VTU(char * file_name);
void SaveOrRead(char *ifile_name,int iSaveReadMode);

void control_cb( int control );
void CalculaCoeficientes_DominioVariableZ();
double CalculaConveccionDifusionDt(vector<double> &Concentracion, vector<double> &U1,
		vector<double> &U2,vector<double> &U3,grid3D *g,int niteraciones,double err0,vector<double> &ConcentracionPrevia);
double CalculaConveccionDifusionDt_DominioVariableZ(vector<double> &Concentracion, vector<double> &U1,
		vector<double> &U2,vector<double> &U3,grid3D *g,int niteraciones,double err0,vector<double> &ConcentracionPrevia);



void AddMensaje(char* newMensaje)
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	int i;
	for (i=0;i<NMensajes-1;i++) {
		strcpy(Mensaje[i],Mensaje[i+1]);
	}
	strcpy(Mensaje[NMensajes-1],newMensaje);
}

void AddMensajeRight(char* newMensaje)
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	//	cout<<"M AddMensajeRight(char* "<<newMensaje<<")"<<endl;
	int i;
	if (NMensajeRight>=MaxMensajeRight)
		for (i=0;i<NMensajeRight-1;i++) {
			strcpy(MensajeRight[i],MensajeRight[i+1]);
		}
	else
		NMensajeRight++;
	cout<<"M NMensajeRight="<<NMensajeRight<<"  newMensaje="<<newMensaje<<endl;
	strcpy(MensajeRight[NMensajeRight-1],newMensaje);
}

void AddMensajeCenter(char* newMensaje)
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	//	cout<<"M AddMensajeRight(char* "<<newMensaje<<")"<<endl;
	int i;
	if (NMensajeCenter==MaxMensajeCenter)
		for (i=0;i<NMensajeCenter-1;i++) {
			strcpy(MensajeCenter[i],MensajeCenter[i+1]);
		}
	else
		NMensajeCenter++;
	strcpy(MensajeCenter[NMensajeCenter-1],newMensaje);
}



void rotacion(float wAngleX,float wAngleY) {

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	glLoadIdentity();
	glRotatef(wAngleX, 1.0f, 0.0f, 0.0f);
	glRotatef(wAngleY, 0.0f, 1.0f, 0.0f);
	glMultMatrixf((GLfloat *)MatrizRotacionGlobal);
	glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat*)MatrizRotacionGlobal);
	glLoadIdentity();
	glMultMatrixf((GLfloat *)MatrizRotacionGlobalINV);
	glRotatef(-wAngleY, 0.0f, 1.0f, 0.0f);
	glRotatef(-wAngleX, 1.0f, 0.0f, 0.0f);
	glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat*)MatrizRotacionGlobalINV);
	FuncionesOpenGL::modelview_calculado=false;
}


void MyDBG() {

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	int i,dif=0;
	for (i=0;i<16;i++) {
		if (fabs(MatrizRotacionGlobal[i]-MatrizRotacionGlobal0[i])>1e-5) dif=1;
	}
	if (dif>0) {
		myfileDBG << "MatrizRotacionGlobal=";
		for (i=0;i<16;i++) {
			myfileDBG <<MatrizRotacionGlobal[i]<<" ";
			MatrizRotacionGlobal0[i]=MatrizRotacionGlobal[i];
		}
		myfileDBG <<endl;
	}
}

void print_text(int x, int y, char* s) 
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	//	cout<<"M print_text(int "<<x<<", int "<<y<<", char* "<<s<<") "<<endl;
	int lines;
	char* p;
	glDisable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, width, 
			0, height, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	//	glColor3ub(128, 0, 255);
	glRasterPos2i(x, y);
	for(p = s; *p; p++) {
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *p);
	}
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glEnable(GL_DEPTH_TEST);
}

void TrasponeMatriz() {

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	MatrizRotacionGlobalINV[0]=MatrizRotacionGlobal[0];
	MatrizRotacionGlobalINV[1]=MatrizRotacionGlobal[4];
	MatrizRotacionGlobalINV[2]=MatrizRotacionGlobal[8];
	MatrizRotacionGlobalINV[3]=MatrizRotacionGlobal[12];

	MatrizRotacionGlobalINV[4]=MatrizRotacionGlobal[1];
	MatrizRotacionGlobalINV[5]=MatrizRotacionGlobal[5];
	MatrizRotacionGlobalINV[6]=MatrizRotacionGlobal[9];
	MatrizRotacionGlobalINV[7]=MatrizRotacionGlobal[13];

	MatrizRotacionGlobalINV[8]=MatrizRotacionGlobal[2];
	MatrizRotacionGlobalINV[9]=MatrizRotacionGlobal[6];
	MatrizRotacionGlobalINV[10]=MatrizRotacionGlobal[10];
	MatrizRotacionGlobalINV[11]=MatrizRotacionGlobal[14];

	MatrizRotacionGlobalINV[12]=MatrizRotacionGlobal[3];
	MatrizRotacionGlobalINV[13]=MatrizRotacionGlobal[7];
	MatrizRotacionGlobalINV[14]=MatrizRotacionGlobal[11];
	MatrizRotacionGlobalINV[15]=MatrizRotacionGlobal[15];
}

void CopiaBloquesAVertices(vector<double> FFB,vector<double> &FV,double *minF,double *maxF,int iBC) {

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	int i,j;
	vector<int> cuantos(FV.size());
	for (i=0;i<gtotal->nV3D;i++) {
		cuantos[i]=0;
		FV[i]=0;
	}
	for (i=0; i<gtotal->nH3D ; i++) {
		for( j=0 ; j<8 ; j++) {
			int iv=gtotal->h3D[i].iv[j];
			FV[iv] += FFB[i];
			cuantos[iv] ++;
		}
	}

	if (iBC>0) {

		for (i=0; i<gtotal->nCaras ; i++) {
			if (gtotal->Cara[i].iBC ==2 || gtotal->Cara[i].iBC ==3) {
				for (j=0; j<4 ; j++) {
					FV[gtotal->Cara[i].iv[j]]=0;
					cuantos[gtotal->Cara[i].iv[j]]=0;
				}
			}
		}
		for (i=0; i<gtotal->nCaras ; i++) {
			if (gtotal->Cara[i].iBC ==2 || gtotal->Cara[i].iBC ==3) {
				for (j=0; j<4 ; j++) {
					switch (iBC)
					{
					case 1:
						FV[gtotal->Cara[i].iv[j]] += gtotal->Cara[i].BC;
						cuantos[gtotal->Cara[i].iv[j]]++;
						break;
					case 2:
						FV[gtotal->Cara[i].iv[j]] += gtotal->Cara[i].BC2;
						cuantos[gtotal->Cara[i].iv[j]]++;
						if(1==0)
							cout<<"g.Cara[i].BC2="<<gtotal->Cara[i].BC2
							<<"\tcuantos[g.Cara[i].iv[j]]="<<cuantos[gtotal->Cara[i].iv[j]]
																	 <<"g.Cara[i].iv[j]="<<gtotal->Cara[i].iv[j]
																											  <<endl;
						break;
					};
				}
			}
		}
	}

	*maxF=*minF=FV[0]/cuantos[0];
	for (i=0;i<gtotal->nV3D;i++) {
		FV[i] /= cuantos[i];
		if (*minF>FV[i]) {
			*minF=FV[i];
			if(1==0) cout<<"FV[i]="<<FV[i]<<"\tcuantos[i]="<<cuantos[i]<<"\ti="<<i<<endl;
		}
		if (*maxF<FV[i]) *maxF=FV[i];
	}

}

void ProyectaVolumenesAVertices(vector<double> & FVol,vector<double> &FVert,double *minF,double *maxF,int FlagTrataBC)
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	int i,j;
	if(DBG) cout<<"ProyectaVolumenesAVertices(..): FVert.size()==gtotal->nV3D: "<<FVert.size()<<"<-->"<<gtotal->nV3D<<endl;
	vector<int> cuantos(FVert.size());
	vector<double> SumaVert(FVert.size());
	for (i=0;i<gtotal->nV3D;i++) {
		cuantos[i]=0;
		SumaVert[i]=0;
	}

	if(DBG) cout<<"220"<<endl;
	for (i=0; i<gtotal->nVolFinito ; i++) { //Triprismas
		for( j=0 ; j<6 ; j++) {
			int iv=gtotal->TriPrisma3D[i].iv[j];
			SumaVert[iv] += FVol[i];
			cuantos[iv] ++;
		}
	}

	if(DBG) cout<<"229: FlagTrataBC="<<FlagTrataBC<<endl;

	if (FlagTrataBC>0) {

		for (i=0; i<gtotal->nCaras ; i++) {
			if (gtotal->Cara[i].iBC ==2 || gtotal->Cara[i].iBC ==3) {
				for (j=0; j<gtotal->Cara[i].nvCara ; j++) {
					switch (FlagTrataBC)
					{
					case 1:
						SumaVert[gtotal->Cara[i].iv[j]] += gtotal->Cara[i].BC;
						cuantos[gtotal->Cara[i].iv[j]]++;
						break;
					case 2:
						SumaVert[gtotal->Cara[i].iv[j]] += gtotal->Cara[i].BC2;
						cuantos[gtotal->Cara[i].iv[j]]++;
						break;
					};
				}
			}
		}
	}

	*maxF=*minF=SumaVert[0]/cuantos[0];
	for (i=0;i<gtotal->nV3D;i++) {
		FVert[i] = SumaVert[i]/cuantos[i];
		if (*minF>FVert[i]) {
			*minF=FVert[i];
		}
		if (*maxF<FVert[i]) *maxF=FVert[i];
	}
	if(DBG) cout<<"ProyectaVolumenesAVertices(..):END"<<endl;

}


void tic() {

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;


	start_t = clock();

}

void toc() {

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	end_t = clock();
	total_t = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;


	cout <<total_t<<" seg";
	if (total_t>60) {
		int min=total_t/60; total_t-=min*60;
		cout<<" ( "<<min<<":"<<total_t<<"min)";
	}
	cout<<endl;
}

/* this function is run by the second thread */
void *Thread_Resolucion_Numerica_Ecuaciones(void *x_void_ptr) {

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	/* increment x to 100 */
	int *x_ptr = (int *)x_void_ptr;
	double a;
	while(++(*x_ptr) ){  //Aqui  está el cálculo, pero solo cuando se lo piden

		pthread_mutex_lock(&mutex_calculo);
		while (Inicia_calculo == 0)
			pthread_cond_wait(&cond_calculo, &mutex_calculo);  //aqui espero que me llamen
		Termine_calculo=0;
		Inicia_calculo=0;
		pthread_mutex_unlock(&mutex_calculo);

		cout<<"pthread: lanzado..."<<endl;
		//lineas de test
		a=0;
		for (int i=0;i<1000000;i++) {
			a+=cos(i)+ sin (i) * tan (i);
		}
		cout<<"pthread: a="<<a<<endl;
		cout<<"pthread: gtotal->Cara.size()="<<gtotal->Cara.size()<<endl;

		//Verdadero Codigo
		Iteracion_Conveccion_Difusion();

		//Variable de control
		Termine_calculo=1;

	}

	printf("x increment finished\n");

	/* the function must return something - NULL will do */
	return NULL;

}

void Lanza_thread() {

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	if (Termine_calculo==1) {
		pthread_mutex_lock(&mutex_calculo);
		Inicia_calculo=1;
		pthread_cond_signal(&cond_calculo);
		pthread_mutex_unlock(&mutex_calculo);

	}
}


void Convierte_Tsimulacion_en_StepHora() {

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	step_hora=Tsimulacion/deltat_hora; //simulacion cada hora????  Mensaje
	stepf=(Tsimulacion-step_hora*deltat_hora)/deltat_hora;

	//cout<<"DrawGraphics: ZdeT.size()"<<ZdeT.size()<<" step="<<step<<endl;
	if (step_hora+1>=(ZdeT.size())) {
		cout<<"DrawGraphics: step_hora+1>=(ZdeT.size()"<<ZdeT.size()<<" step_hora="<<step_hora<<endl;
		clock00=clock_Time;
		clock_TimeP=clock_Time;
		Tsimulacion=0;
		TsimulacionP=0;
		TColaP=-1000;
		step_hora=0;
		stepf=0;
	}
	cout <<"Tsimulacion="<<Tsimulacion<<" step_hora="<<step_hora<<" stepf="<<stepf<<endl;
}

void DrawGraphics()
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	if (DBG) cout<<"DrawGraphics()"<<endl;
	static int contador=0;
	char linea[100];
	int i,j;

	if (!DrawYES) return; 
	DrawYES=0;


	pthread_mutex_lock(&mutex_DrawRead);
	//Test1: pthreads

	/* this variable is our reference to the second thread */
	//pthread_t inc_x_thread;

	/* create a second thread which executes Thread_Resolucion_Numerica_Ecuaciones(&x) */
	if (thread_lanzado==0) {
		Inicia_calculo=0;Termine_calculo=1;
		if(pthread_create(&inc_x_thread, NULL, Thread_Resolucion_Numerica_Ecuaciones, &x_thread)) {

			fprintf(stderr, "Error creating thread\n");
			return;

		}
		thread_lanzado=1;
	}
	if (x_thread != x_thread_p) {
		cout <<"x_thread="<<x_thread<<endl;
		cout<<"main: gtotal="<<gtotal<<endl;
		if (gtotal)
			cout<<"main: gtotal->Cara.size()="<<gtotal->Cara.size()<<endl;
		x_thread_p=x_thread;
	}
	if (0&& gtotal)
		cout<<"main2: gtotal->Cara.size()="<<gtotal<<","<<gtotal->Cara.size()<<endl;

	if (x_thread>2 & gtotal!=NULL && gtotal->Cara.size()==0) exit(1);



	// actualizacion del t actual;
	clock_Time=clock();
	if (Calculando ==1 && ZdeT.size()>0 ) {
		Tsimulacion=TsimulacionP+(clock_Time-clock_TimeP)*FactorClockTime;
		//if (Tsimulacion>TsimulacionP+100) Tsimulacion=TsimulacionP+100; //Datos_dt <= 100 !!!
		cout<<"Tsimulacion="<<Tsimulacion<<endl;

		//Todo , JSM test de ZdeT
		if (gtotal!=NULL && ZdeT.size()>0 ) {
			for(i=0;i<gtotal->nV3D;i++) {
				gtotal->v3D[i].z=ZdeT[step_hora][i]*(1-stepf)+ZdeT[step_hora+1][i]*stepf;
			}
		}
		Convierte_Tsimulacion_en_StepHora();

	}

#if 0

	if (ZdeT.size()>0 && step_hora>240) {
		for(i=0;i<gtotal->nV3D;i++) {
			gtotal->v3D[i].vz=0;
			gtotal->v3D[i].nvz=0;
		}
		for(i=0;i<gtotal->nTriPrisma3D;i++) {
			for(j=0;j<6;j++) {
				int ip=i;
				int jj=gtotal->TriPrisma3D[i].iv[j];
				gtotal->v3D[jj].vz += (W[step_hora][ip]*(1-stepf)+W[step_hora+1][ip]*stepf);
				gtotal->v3D[jj].nvz ++;
			}
		}
		for(i=0;i<gtotal->nV3D;i++) {
			if( gtotal->v3D[i].nvz>3 && i%4 != 0) {
				gtotal->v3D[i].vz /= gtotal->v3D[i].nvz;
				gtotal->v3D[i].z += gtotal->v3D[i].vz* (Tsimulacion-TsimulacionP);
			}
		}
	}
#endif

	char  mensajeCenter[1000];
	double Tdia,Thora;

	Tdia=Tsimulacion/24/3600,Thora=(Tdia-floor(Tdia))*24;

	if (writeTofileBMP==1) {
		sprintf(mensajeCenter,"T=%dd %.1fh (frame=%03d)",(int)Tdia,Thora,iframe);

	}
	else {
		sprintf(mensajeCenter,"T=%dd %.1fh",(int)Tdia,Thora);
	}
	AddMensajeCenter(mensajeCenter);
	sprintf(mensajeCenter,"");
	AddMensajeCenter(mensajeCenter);

	if (DBG) cout<<"DrawGraphics()138"<<endl;

	glClearColor(BGColorR/255. *FactorClaridad/1.7, BGColorG/255. *FactorClaridad/1.7, BGColorB/255. *FactorClaridad/1.7, 1.0f);

	glClear(GL_COLOR_BUFFER_BIT );
	glClear(GL_DEPTH_BUFFER_BIT);


#if DBGMouse
	MyDBG();
#endif

	glDisable(GL_CLIP_PLANE0);
	if (DBG) cout<<"DrawGraphics()143"<<endl;

	glClear(GL_DEPTH_BUFFER_BIT);
	if (DBG) cout<<"DrawGraphics()166"<<endl;



	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective((GLdouble)GlobalFovy, aspect, (GLdouble)1, (GLdouble)100.0);

	//	glViewport(0, 0, width, height);
	glTranslated( 0, 0, -10);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslated( 0, -1, -10);
	glShadeModel(GL_SMOOTH);
	FuncionesOpenGL::ActivaLuz0();


	glTranslatef(vecUEsfera[0], vecUEsfera[1],vecUEsfera[2]);
	FuncionesOpenGL::material(0);
	FuncionesOpenGL::esfera(0.002*GlobalFovy,20);

	glScalef(Escala, Escala, Escala);
	if (FlagPrintEscala) {
		cout <<"Escala="<<Escala<<endl;
	}
	glMultMatrixf((GLfloat *)MatrizRotacionGlobal);

	if (0*FlagPrintMatrizRotacionGlobal) {
		cout <<"float 	  M3_Escala="<<Escala<<";"<<endl;
		cout <<"GLfloat   M3_MatrizRotacionGlobal[16]= {";
		for (i=0;i<16;i++){
			cout<<MatrizRotacionGlobal[i];
			if (i<15) cout<<" , ";
		}
		cout <<"};"<<endl;
		cout <<"GLfloat   M3_MatrizRotacionGlobalINV[16]= {";
		for (i=0;i<16;i++){
			cout<<MatrizRotacionGlobalINV[i];
			if (i<15) cout<<" , ";
		}
		cout <<"};"<<endl;
	}

	glTranslatef(-vecXEsfera[0],-vecXEsfera[1],-vecXEsfera[2]);

	if (MODO_Ejes) {
		FuncionesOpenGL::ejes();
	}


	if (0 && ColorON==2) {

		AddMensajeCenter("ColorON==2");
		sprintf(mensajeCenter,"FlagEvoluciona==%d",FlagEvoluciona);
		AddMensajeCenter(mensajeCenter);
		double minF,maxF;
		if (FlagDifusion && FlagConveccion)
			sprintf(text,"Calculo con Convección y Difusión\n");
		if (FlagDifusion==0 && FlagConveccion)
			sprintf(text,"Calculo de Convección pura\n");
		if (FlagDifusion && FlagConveccion==0)
			sprintf(text,"Calculo de Difusión pura\n");
		if (FlagDifusion==0 && FlagConveccion==0)
			sprintf(text,"Calculo de Estacionario\n");
		if (FlagEvoluciona) {
			Datos_dt=(Tsimulacion-TsimulacionP); //Todo //Corregir
			cout<<"Tsimulacion="<<Tsimulacion<<"step_hora="<<step_hora<<endl;
			if (gtotal!=NULL && Datos_dt>0) {
				for (i=0;i<gtotal->nVolFinito;i++) {
					F2VolumenesP[i]=F2Volumenes[i];
				}

#if DominioVariable
				double AchicaPrevia_min=1e20,Achicat_min=1e20;
				for (i=0;i<gtotal->nVolFinito;i++) {

					if (Achica_t_PrimeraVez==1) {
						AchicaPrevia[i]=Achica[step_hora][i];
					} else {
						AchicaPrevia[i]=Achica_t[i];
					}
					if (AchicaPrevia[i]<AchicaPrevia_min) AchicaPrevia_min=AchicaPrevia[i];
					Achica_t[i]=Achica[step_hora][i]*(1-stepf)+Achica[step_hora+1][i]*stepf;
					if (Achica_t[i]<Achicat_min) Achicat_min=Achica_t[i];
					if (i<-20)
						cout <<"0) AchicaPrevia[3]="<<AchicaPrevia[3]<<",  i="<<i<<endl;
				}

				Achica_t_PrimeraVez=0;
				cout << "Achicat_min=     "<<Achicat_min<<endl;
				cout << "AchicaPrevia_min="<<AchicaPrevia_min<<endl;
				cout <<"1) AchicaPrevia[3]="<<AchicaPrevia[3]<<endl;

#endif
				TipoCalculo=CalculoEvolucion;
				//				CalculaConveccionDifusionDt(F2Volumenes, U[step_hora],V[step_hora],W[step_hora],gtotal,30,1e-4,F2VolumenesP);

				CalculaConveccionDifusionDt(F2Volumenes, U[step_hora],V[step_hora],W[step_hora],gtotal,30,1e-7,F2VolumenesP);
				int FlagTrataBC=0;
				double minF,maxF;
				ProyectaVolumenesAVertices(F2Volumenes,F2Nodos,&minF,&maxF,FlagTrataBC) ;

				sprintf(text,"%s errG=%g ",text,Global_error);

				sprintf(text,"%s minF=%f, maxF=%f\n",text,minF,maxF);


			}
		}
		else
			Tsimulacion=TsimulacionP;

		int FlagTrataBC=0;
		ProyectaVolumenesAVertices(F2Volumenes,F2Nodos,&minF,&maxF,FlagTrataBC) ;
		sprintf(text,"%sConv/Diff=%f\n",text,Peclet);


		if (glui != NULL) glui_edittext->set_text(text);



		if (DBG) cout<<"378"<<endl;
		//FlagCalculaEvolucion=0;
	}

	//	FuncionesOpenGL::ActivaLuz0();

	if (ClippingON) {
		double eq[4];
		eq[0]=-Ax;	eq[1]=-By;	eq[2]=-Cz/FactorZ;	eq[3]=DD;
		glClipPlane(GL_CLIP_PLANE0, eq);
		glEnable(GL_CLIP_PLANE0);
	} 
	if (DBG) cout<<"DrawGraphics()380" <<" ColorON="<<ColorON <<" gtotal="<<gtotal <<endl;


	//glEnable(GL_CULL_FACE);

	factorDespegaLineas=3e-2/sqrt(Escala); ///Formula deducida con Matlab....

	if (0 && ColorON==0) {
		////////////TEXTURE 1
		glGenTextures( 1, &textureID );
		if (DBG) cout<<" 383 "<<endl;
		glBindTexture( GL_TEXTURE_2D, textureID );


		//gluBuild2DMipmaps( GL_TEXTURE_2D, 3, imagenwidth, imagenheight,GL_RGB, GL_UNSIGNED_BYTE, imagesdata );

		//cout<<"imagenwidth="<<imagenwidth<<"  imagenheight="<<imagenheight<<endl;

		glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, imagenwidth, imagenheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagendata);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

		////////////TEXTURE 2
		glGenTextures( 1, &textureID2 );
		if (DBG) cout<<"398"<<endl;
		glBindTexture( GL_TEXTURE_2D, textureID2 );


		glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, imagen2width, imagen2height, 0, GL_RGB, GL_UNSIGNED_BYTE, imagen2data);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);


		//////////ENDTEXTURE
	}



	if (gtotal!=NULL && MODO_CampoVelocidades) {
		//		glDisable(GL_CLIP_PLANE0);
		//		FuncionesOpenGL::material(0);	gtotal->drawVelGL(U,V,W);
		FuncionesOpenGL::material(0);	gtotal->drawParticulas_TriPrisma();
	}

	if(DBG) cout<<"419 ColorON="<<ColorON<<endl;

	switch (ColorON)
	{
	case 0:  //Caso Normal
		AddMensajeCenter("ColorON==0: drawGL()");
		if (gtotal!=NULL ) {
			FuncionesOpenGL::material(1);gtotal->drawGL();
		}
		if (gPoly!=NULL) {

			FuncionesOpenGL::material(1);gPoly->drawGL();
		}
		break;
	case 1:  //Dibuja Campo 1

		AddMensajeCenter("ColorON==1: drawGL(F)");
		if (gtotal!=NULL ) gtotal->drawGL(F);
		break;
	case 2:  //Dibuja Campo 2

		AddMensajeCenter("ColorON==2: drawGL(F2Nodos,0,1)");
		if (gtotal!=NULL ) gtotal->drawGL(F2Nodos,0,1); //Ya dibujado
		break;
	case 3:
		if (gtotal!=NULL ) {
			AddMensajeCenter("ColorON==3: drawGL(), drawVoronoi()");
			FuncionesOpenGL::material(1);gtotal->drawGL();
			FuncionesOpenGL::material(2);gtotal->drawVoronoi();
		}
		break;
	};

	if (DBG) cout<<"DrawGraphics()249"<<endl;
	if (gtotal!=NULL && MODO_CampoVelocidades2) {
		FuncionesOpenGL::material(1);	gtotal->drawVelGL2();
	}
	if (DBG) cout<<"DrawGraphics()253"<<endl;



	DrawMenuOnOff();

	{/////////////////bloque que imprime texto......
		glDisable(GL_CLIP_PLANE0);
		glPushAttrib( GL_LIGHTING_BIT );
		glDisable( GL_LIGHTING );

		if (DBG) cout<<"DrawGraphics()260"
				<<"\n MODO_MenuMENSAJES="<<MODO_MenuMENSAJES
				<<endl;
		if (MODO_MenuMENSAJES)
			DrawMensajes();
		DrawMensajesRight();
		DrawMensajesCenter();
		if(1==1) {
			clock2F=clockF=clock ();
			nframes++;nframes2++;
			if (DBG) cout<<"DrawGraphics()267"<<endl;

			FPS=(nframes*1.0)*CLOCKS_PER_SEC/(clockF-clock0);
			FPS2=(nframes2*1.0)*CLOCKS_PER_SEC/(clock2F-clock2);
			if (clockF-clock0>CLOCKS_PER_SEC) {
				glColor3d(.8,.9, 1);
				glRasterPos3f(-2.07*aspect,1.95,5);
				///modificado GLUT
				sprintf(s,"FPS=%.2f (%.2f)\r",FPS2,FPS);//fflush(stdout);
				//sprintf(s,"%.0f F/S",FPS);FuncionesOpenGL::Print(s,1);
				clock0=clockF;
				nframes=0;
				glui_FPS->set_name(s);
				//					printf("%s",s);
			}
			//			print_text(width-150,5,s);

		}
		glPopAttrib();
	}

	/// print BC
	///
	if (gPoly != NULL ) {
		int j=0;
		for (i=0;i<gPoly->nBC;i++) {
			if (gPoly->DrawBC[i]>0) {
				char s[100];
				sprintf(s,"[%i] -",gPoly->BC[i]);
				print_text(3,height-17*(++j+1),s);
			}
		}
	}

	if (writeTofileBMP==1) {
#if 1
		if (iframe<0) {
			char nombre[100];
			iframe=0;
			//sprintf(nombre,"D:/Temuco2019/A01/velocity%06d.vtu",iframe);
			//Lee_Velocidad_VTU(nombre);
			//DeboDibujarSiguiente=2;

			TimeSaved=step_hora+stepf-DT_toSave;
			cout<<"writeTofileBPM:(iframe<0==TRUE)"<<iframe<<" TimeSaved="<<TimeSaved<<endl;
		}
		else {
			TimeToSave=step_hora+stepf;
			if (TimeToSave>=TimeSaved+DT_toSave-1e-8) {
				TimeSaved+=DT_toSave;
				cout<<"writeTofileBPM:"<<iframe<<" TimeSaved="<<TimeSaved<<endl;
				int size = 3 * width * height;
				writedata = new char[size+1]; // allocate 3 bytes per pixel

				glPixelStorei(GL_PACK_ALIGNMENT, 1);
				glReadPixels(0, 0, width, height, GL_BGR, GL_UNSIGNED_BYTE, writedata);


				// Guardar un archivo BMP con los datos leidos (como verificacion)
				char nombre[200];
				sprintf(nombre,"anim/BMP/bmpofstream_%08d.bmp",iframe);
				ofstream arrayfile(nombre, std::ios::out | std::ios::binary); // File Creation

				headerinfo[18]=(char) width;
				headerinfo[22]=(char) height;


				arrayfile.write(headerinfo, 18);
				arrayfile.write((char *)&width, 4);
				arrayfile.write((char *)&height, 4);
				arrayfile.write((char *)(headerinfo+26), 54-18-4-4);
				arrayfile.write(writedata, size);
				delete writedata;
				arrayfile.close();

				cout<<"writeTofileBPM:END OK"<<endl;//3115
				iframe++;
				if (iframe>MaxFrames) {
					writeTofileBMP=0;
					iframe=-1;
				}
				else {
					//sprintf(nombre,"D:/Temuco2019/A01/velocity%06d.vtu",iframe);
					//Lee_Velocidad_VTU(nombre);
					DeboDibujarSiguiente=2;
				}
			}
		}
#endif

	}
	glutSwapBuffers();

	clock_TimeP=clock_Time;
	TsimulacionP=Tsimulacion;



	glDeleteTextures( 1, &textureID2 );
	glDeleteTextures( 1, &textureID );

	DrawYES=1;


	pthread_mutex_unlock(&mutex_DrawRead);
	if (0 && gtotal)
		cout<<"main2: gtotal->Cara.size()="<<gtotal<<","<<gtotal->Cara.size()<<endl;
	if (DBG) cout<<"DrawGraphics():END"<<endl;

}

void writeUVWZdeT_Binario(char * file_name)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	ofstream myfile;
	int i,j;
	int dim1,dim2,dim3;

	myfile.open (file_name, std::ios::out | std::ios::binary);
	int irw=0;
	dim1=ZdeT.size();
	dim2=ZdeT[0].size();
	dim3=U[0].size();

	double d;

	myfile.write((char*)&dim1,sizeof(dim1));cout<<dim1<<" ";
	myfile.write((char*)&dim2,sizeof(dim2));cout<<dim2<<" ";
	myfile.write((char*)&dim3,sizeof(dim3));cout<<dim3<<" ";
	for (i=0;i<dim1;i++) {
		for (j=0;j<dim2;j++) {
			d=ZdeT[i][j]; myfile.write((char*)&d,sizeof(d));	if(irw<20){cout<<d<<" "; irw++;}
		}
		for (j=0;j<dim3;j++) {
			d=U[i][j];    myfile.write((char*)&d,sizeof(d));   	if(irw<20){cout<<d<<" "; irw++;}
			d=V[i][j];    myfile.write((char*)&d,sizeof(d));   	if(irw<20){cout<<d<<" "; irw++;}
			d=W[i][j];    myfile.write((char*)&d,sizeof(d));	if(irw<20){cout<<d<<" "; irw++;}
		}
	}

	myfile.close();
}

void readUVWZdeT_Binario(char * file_name)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	ifstream myfile;
	int i,j;
	int dim1,dim2,dim3;


	cout<<__FUNCTION__<<"@"<<__LINE__<<" "<<endl;

	myfile.open (file_name, std::ios::in | std::ios::binary);
	int irw=0;


	myfile.read((char*)&dim1,sizeof(dim1));cout<<dim1<<" ";
	myfile.read((char*)&dim2,sizeof(dim2));cout<<dim2<<" ";
	myfile.read((char*)&dim3,sizeof(dim3));cout<<dim3<<" ";


	dim1=dim1/3;
	ZdeT.resize(dim1);
	U.resize(dim1);
	V.resize(dim1);
	W.resize(dim1);

	double d,d1,d2,d3;

	cout<<__FUNCTION__<<"@"<<__LINE__<<" "<<endl;

	for (i=0;i<dim1;i++) {
		//cout<<__FUNCTION__<<"@"<<__LINE__<<" "<<"i="<<i<<endl;
		if (dim3 != gtotal->nTriPrisma3D) {
			cout << "Error en archivo binario : dim3="<<dim3<<endl;
			cout <<"gtotal->nTriPrisma3D="<<gtotal->nTriPrisma3D<<endl;
		}
		ZdeT[i].resize(dim2);
		U[i].resize(dim3);
		V[i].resize(dim3);
		W[i].resize(dim3);
		//cout<<__FUNCTION__<<"@"<<__LINE__<<" "<<endl;
		for (j=0;j<dim2;j++) {
			myfile.read((char*)&d,sizeof(d)); ZdeT[i][j]=d;		if(irw<20){cout<<d<<" "; irw++;}
		}
		//cout<<__FUNCTION__<<"@"<<__LINE__<<" "<<endl;
		for (j=0;j<dim3;j++) {
			myfile.read((char*)&d1,sizeof(d1)); U[i][j]=d1;		if(irw<20){cout<<d1<<" "; irw++;}
			myfile.read((char*)&d2,sizeof(d2)); V[i][j]=d2;		if(irw<20){cout<<d2<<" "; irw++;}
			myfile.read((char*)&d3,sizeof(d3)); W[i][j]=d3;		if(irw<20){cout<<d3<<" "; irw++;}
			d=sqrt(sqr(d1)+sqr(d2)+sqr(d3));
			if (d>FactorVel) FactorVel=d;
		}
		//cout<<__FUNCTION__<<"@"<<__LINE__<<" "<<endl;
	}

	myfile.close();
}

void save(vector<double> &F,ofstream &myfile) {


	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	int i;
	myfile  <<F.size()<<endl;
	for (i=0;i<     F.size();i++) { myfile<<     F[i]<<" "; } myfile<<endl;

}

void read(vector<double> &F,ifstream &myfile) {

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	int i,sizeL;
	myfile  >> sizeL; F.resize(sizeL);
	for (i=0;i<F.size();i++) { myfile>>F[i]; }

}

void DrawMensajesDatos()
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	char str[100];
	int i;


	for (i=0;i<NMensajes;i++) {
		strcpy(Mensaje[i],"");
	}


}


void DrawMensajesDatosGui3()
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	// Nueva version....
	char str[100];
	int i;

	sprintf(StringParametros,"#Parametros usados:\n");


	sprintf(StringParametros,"%smalla           =%d\n",StringParametros,CualMalla           );
	sprintf(StringParametros,"%sCalculoContinuo =%d\n",StringParametros,CalculoContinuo	    );
	sprintf(StringParametros,"%sTipoCalculo     =%d\n",StringParametros,TipoCalculo	        );
	sprintf(StringParametros,"%snR              =%d\n",StringParametros,nR				    );
	sprintf(StringParametros,"%snR2             =%d\n",StringParametros,nR2				    );
	sprintf(StringParametros,"%snTh             =%d\n",StringParametros,nTh				    );
	sprintf(StringParametros,"%sNDivZ           =%d\n",StringParametros,NDivZ			    );
	sprintf(StringParametros,"%sDominio_Xmax    =%f\n",StringParametros,Dominio_Xmax	    );
	sprintf(StringParametros,"%sDominio_Rint    =%f\n",StringParametros,Dominio_Rint	    );
	sprintf(StringParametros,"%sDominio_Rmed    =%f\n",StringParametros,Dominio_Rmed	    );
	sprintf(StringParametros,"%sDominio_Hmax    =%f\n",StringParametros,Dominio_Hmax	    );
	sprintf(StringParametros,"%sDominio_Hsup    =%f\n",StringParametros,Dominio_Hsup		);
	sprintf(StringParametros,"%srhof            =%f\n",StringParametros,Datos_rhof			);
	sprintf(StringParametros,"%srhos            =%f\n",StringParametros,Datos_rhos			);
	sprintf(StringParametros,"%sphi             =%f\n",StringParametros,Datos_phi			);
	sprintf(StringParametros,"%scf              =%f\n",StringParametros,Datos_cf			);
	sprintf(StringParametros,"%sTinyeccion      =%f\n",StringParametros,Datos_Tinyeccion	);
	sprintf(StringParametros,"%sTambiente       =%f\n",StringParametros,Datos_Tambiente	    );
	sprintf(StringParametros,"%sTLimite        =%f\n",StringParametros,TLimite 	    );
	sprintf(StringParametros,"%shc_superior     =%f\n",StringParametros,Datos_hc_superior	);
	sprintf(StringParametros,"%shc_inferior     =%f\n",StringParametros,Datos_hc_inferior	);
	sprintf(StringParametros,"%shc_lateral      =%f\n",StringParametros,Datos_hc_lateral	);
	sprintf(StringParametros,"%sKTermofilm      =%e\n",StringParametros,Datos_KTermofilm	);
	sprintf(StringParametros,"%seTermofilm      =%e\n",StringParametros,Datos_eTermofilm	);
	sprintf(StringParametros,"%sVinyeccionLtsHr =%f\n",StringParametros,VinyeccionLtsHr	    );
	sprintf(StringParametros,"%skf              =%f\n",StringParametros,Datos_kf			);
	sprintf(StringParametros,"%sks              =%f\n",StringParametros,Datos_ks				     );
	sprintf(StringParametros,"%sKevaporacion       =%e\n",StringParametros,Datos_Kevaporacion	     );
	sprintf(StringParametros,"%sHumedadAmbiental   =%f\n",StringParametros,Datos_HumedadAmbiental    );
	sprintf(StringParametros,"%sDistanciaAlBorde   =%f\n",StringParametros,Datos_DistanciaAlBorde    );
	sprintf(StringParametros,"%sCalorLatenteFluido =%e\n",StringParametros,Datos_CalorLatenteFluido);


	if (glui3 != NULL) Glui3_EditParametros->set_text(StringParametros);

}


void DrawMenuOnOff()
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	//	cout<<"M DrawMensajesRight"<<endl;
	glPushAttrib( GL_LIGHTING_BIT );
	glDisable( GL_LIGHTING );
	glColor3d(1,0, 0);

	x0Menu=0;
	if (gluiEligeBC_Visible) x0Menu+=gluiEligeBC_w;
	y0Menu=0;
	x1Menu=x0Menu+120;
	y1Menu=y0Menu+30;
	glDisable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, width, height, 0, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	if (EstadoMenu)
		glBegin(GL_QUADS);
	else {
		glBegin(GL_LINE_LOOP);
	}
	glVertex2i(x0Menu, y0Menu);
	glVertex2i(x1Menu, y0Menu);
	glVertex2i(x1Menu, y1Menu);
	glVertex2i(x0Menu , y1Menu);

	glEnd();



	char str[100];
	int i;

	strcpy(str,"RButton=Menu");

	glColor3d(1,1, 0);
	glRasterPos2i((x0Menu+x1Menu)/2-9*strlen(str)/2,(y0Menu+y1Menu)/2+7);
	char* p;

	for(p=str; *p; p++) {
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *p);
	}



	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glEnable(GL_DEPTH_TEST);




}


void DrawMensajesCenter()
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	//	cout<<"M DrawMensajesRight"<<endl;
	int i,x0,y0,x1,y1;
	glPushAttrib( GL_LIGHTING_BIT );
	glDisable( GL_LIGHTING );
	glColor3d(1,1, 0);

	for (i=0;i<NMensajeCenter;i++){
		x0=width/2-500;
		y0=height-17*(i+1);
		print_text(x0,y0,MensajeCenter[i]);






	}
}
void DrawMensajesRight()
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	//	cout<<"M DrawMensajesRight"<<endl;
	int i;
	glPushAttrib( GL_LIGHTING_BIT );

	glDisable( GL_LIGHTING );


	glColor3d(1,1, 0.1);


	for (i=0;i<NMensajeRight;i++){
		print_text(width-640,height-17*(i+1),MensajeRight[i]);
	}
}

void DrawMensajes()
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	int i;
	glPushAttrib( GL_LIGHTING_BIT );

	glDisable( GL_LIGHTING );


	glColor3d(1,1, 0.7);


	for (i=0;i<NMensajes;i++){
		//		print_text(3,height-17*(NMensajes-i),Mensaje[i]);
		print_text(3,height-17*(i+1),Mensaje[i]);
	}




	switch (mode)
	{
	case modeT : 
		glPushMatrix();
		glColor3d(.8,.9, 1);

		glTranslatef(-2.07*aspect,-2,5);
		glRasterPos2f(0,0);

		glListBase(baseBIT);
		glCallLists(12, GL_UNSIGNED_BYTE, "(T)raslación"); 

		glPopMatrix();
		break;
	case modeR:
		glPushMatrix();
		glColor3d(.8,.9, 1);

		glTranslatef(-2.07*aspect,-2,5);
		glRasterPos2f(0,0);

		glListBase(baseBIT);
		glCallLists(10, GL_UNSIGNED_BYTE, "(R)otación"); 

		glPopMatrix();
		break;
	}
	glPopAttrib();

}


void MatrizXvector4(GLfloat *mat,GLfloat *vec,GLfloat *resultado)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	resultado[0]=mat[0]*vec[0] + mat[4]*vec[1] + mat[ 8]*vec[2] + mat[12]*vec[3];
	resultado[1]=mat[1]*vec[0] + mat[5]*vec[1] + mat[ 9]*vec[2] + mat[13]*vec[3];
	resultado[2]=mat[2]*vec[0] + mat[6]*vec[1] + mat[10]*vec[2] + mat[14]*vec[3];
	resultado[3]=mat[3]*vec[0] + mat[7]*vec[1] + mat[11]*vec[2] + mat[15]*vec[3];
}


void MatrizXvector4(GLdouble *mat,GLfloat *vec,GLfloat *resultado)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	resultado[0]=mat[0]*vec[0] + mat[4]*vec[1] + mat[ 8]*vec[2] + mat[12]*vec[3];
	resultado[1]=mat[1]*vec[0] + mat[5]*vec[1] + mat[ 9]*vec[2] + mat[13]*vec[3];
	resultado[2]=mat[2]*vec[0] + mat[6]*vec[1] + mat[10]*vec[2] + mat[14]*vec[3];
	resultado[3]=mat[3]*vec[0] + mat[7]*vec[1] + mat[11]*vec[2] + mat[15]*vec[3];
}

void MatrizTrXvector4(GLdouble *mat,GLfloat *vec,GLfloat *resultado)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	resultado[0]=mat[ 0]*vec[0] + mat[ 1]*vec[1] + mat[ 2]*vec[2] + mat[ 3]*vec[3];
	resultado[1]=mat[ 4]*vec[0] + mat[ 5]*vec[1] + mat[ 6]*vec[2] + mat[ 7]*vec[3];
	resultado[2]=mat[ 8]*vec[0] + mat[ 9]*vec[1] + mat[10]*vec[2] + mat[11]*vec[3];
	resultado[3]=mat[12]*vec[0] + mat[13]*vec[1] + mat[14]*vec[2] + mat[15]*vec[3];
}


void ZGlobal(double*v)
{
	if (Cuenta) CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	v[0]=MatrizRotacionGlobalINV[8];
	v[1]=MatrizRotacionGlobalINV[9];
	v[2]=MatrizRotacionGlobalINV[10];
}

void inicializacionGL()
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	if (DBG) cout<<"inicializacionGL()"<<endl;
	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective( 0.0 , 640.0/480.0, 0.001, 100.0);
	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
	//	glRotatef(-60,1,0,0);
	//	glRotatef(-135,0,0,1);
	glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat*)MatrizRotacionGlobal);
	glLoadIdentity();
	//	glRotatef(135,0,0,1);
	//	glRotatef(60,1,0,0);
	glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat*)MatrizRotacionGlobalINV);
	FuncionesOpenGL::modelview_calculado=false;
	if (DBG) cout<<"inicializacionGL():END"<<endl;
}
void inicializacion()
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	int i,j,iL,malla;
	double x1,x2,y1,y2,xL,yL,R,xF,yF,RF,th,pi=4*atan(1.0);

	for (i=0;i<NMensajes;i++) {
		Mensaje[i]=(char *)(new char[100]);
		strcpy(Mensaje[i],"");
	}

	for (i=0;i<MaxMensajeRight;i++) {
		MensajeRight[i]=(char *)(new char[100]);
		strcpy(MensajeRight[i],"");
	}

	for (i=0;i<MaxMensajeCenter;i++) {
		MensajeCenter[i]=(char *)(new char[100]);
		strcpy(MensajeCenter[i],"");
	}

	clock2=clock0=clock00=clock ();
	glClearColor(.9, .9, 1.0, 1.0f);
	glClearColor(BGColorR/255., BGColorG/255., BGColorB/255., 1.0f);
	glClearDepth(1.0);

	malla=0;
	if (malla>0) {
		gtotal->cubo(5,2,2);
		g1.cubo(5,2,2);
		g2.cubo(5,2,2);
		switch (malla)
		{
		case 1:
			gtotal->cubo(2,2,2);
			g2.cubo(1,2,2);
			x1=0.1;x2=1;
			y1=0.1;y2=1;
			for (i=0;i<gtotal->nV3D ;i++){
				xL=gtotal->v3D[i].x ;
				yL=gtotal->v3D[i].y ;
				gtotal->v3D[i].x = x1+xL*(x2-x1);
				gtotal->v3D[i].y = yL*(y1+xL*(y2-y1));
			}
			x1=x2;x2=1.5;
			for (i=0;i<g2.nV3D ;i++){
				xL=g2.v3D[i].x ;
				g2.v3D[i].x = x1+xL*(x2-x1);
			}
			gtotal->Junta(g2);
			g2=*gtotal;

			for (i=0;i<g2.nV3D ;i++){
				g2.v3D[i].y *= -1;
			}
			for (i=0;i<g2.nH3D ;i++){
				for (j=0;j<4;j++) {
					iL=g2.h3D[i].iv[j];
					g2.h3D[i].iv[j]=g2.h3D[i].iv[7-j];
					g2.h3D[i].iv[7-j]=iL;
				}
			}
			gtotal->Junta(g2);

			g2=*gtotal;
			g2.Rota90Z();	gtotal->Junta(g2);
			g2.Rota90Z();	gtotal->Junta(g2);
			g2.Rota90Z();	gtotal->Junta(g2);
			g2.cubo(1,1,2);
			x1=1;x2=1.5;
			for (i=0;i<g2.nV3D ;i++){
				xL=g2.v3D[i].x ;
				yL=g2.v3D[i].y ;
				g2.v3D[i].x = x1+xL*(x2-x1);
				g2.v3D[i].y = x1+yL*(x2-x1);
			}
			gtotal->Junta(g2);
			g2.Rota90Z();	gtotal->Junta(g2);
			g2.Rota90Z();	gtotal->Junta(g2);
			g2.Rota90Z();	gtotal->Junta(g2);
			break;
		case 2:
			gtotal->cubo(2,2,2);g2=*gtotal;
			x1=0.05;x2=1;y1=0.1;y2=1;
			for (i=0;i<gtotal->nV3D;i++) {
				xL=gtotal->v3D[i].x; yL=gtotal->v3D[i].y;
				R=x1+(x2-x1)*xL;
				th=pi/4*yL;
				xF=R; yF=R*yL;
				xF=xF*xL+(1-xL)*R*cos(th);
				yF=yF*xL+(1-xL)*R*sin(th);
				gtotal->v3D[i].x=xF;
				gtotal->v3D[i].y=yF;
			}
			for (i=0;i<g2.nV3D;i++) {
				xL=g2.v3D[i].x;	yL=g2.v3D[i].y;
				R=x1+(x2-x1)*yL; th=pi/4*xL;
				yF=R; xF=R*xL;
				xF=xF*yL+(1-yL)*R*sin(th);
				yF=yF*yL+(1-yL)*R*cos(th);
				g2.v3D[i].x=xF;
				g2.v3D[i].y=yF;
			}
			gtotal->Junta(g2);

			g2=*gtotal;g2.Rota90Z();gtotal->Junta(g2);g2.Rota90Z();gtotal->Junta(g2);g2.Rota90Z();gtotal->Junta(g2);
			g2=*gtotal;
			g2.Traslada(2,0,0);gtotal->Junta(g2);
			g2.Traslada(0,2,0);gtotal->Junta(g2);
			g2.Traslada(-2,0,0);gtotal->Junta(g2);


		}

		F.resize(gtotal->nV3D);
		U.resize(gtotal->nV3D);
		V.resize(gtotal->nV3D);
		W.resize(gtotal->nV3D);

		for (i=0;i<gtotal->nV3D;i++) {
			float x,y,z;
			x=gtotal->v3D[i].x;
			y=gtotal->v3D[i].y;
			z=gtotal->v3D[i].z;
			F[i]=sqrt(x*x+z*z);
		}
	}
	inicializacionGL();
}

///Cambiado para GLUT de () a void ResizeGraphics(int lwidth, int lheight)
void ResizeGraphics(int lwidth, int lheight)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	// Get new window size

	if (DBG) cout<<"ResizeGraphics()"<<endl;
	int tx, ty, tw, th;
	tx=0; ty=0;	tw= lwidth;	th = lheight;

#if defined(GLUI_GLUI_H)
	if (!glui_hide) {
		GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );
	}
#endif
	width = tw;
	height = th;

	aspect = (GLfloat) width / height;

	// Adjust graphics to window size
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective((GLdouble)GlobalFovy, aspect, (GLdouble)1, (GLdouble)100.0);
	glViewport(tx, ty, width, height);
	glTranslated( 0, 0, -10);
	glMatrixMode(GL_MODELVIEW);
	FuncionesOpenGL::projection_calculado=false;
	FuncionesOpenGL::viewport_calculado=false;
	if (DBG) cout<<"ResizeGraphics():END"<<endl;
}

void CB_mouse_whell(int button, int state, int x, int y)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	//	printf("CB_mouse_whell: button=%d, state=%d, x=%d, y=%d\n",button,state,x,y);	cout<<endl;
	switch (state)
	{
	case 1:

		Escala *= 1.2;
		FuncionesOpenGL::modelview_calculado=false;
		break;
	case -1:

		Escala /= 1.2;
		FuncionesOpenGL::modelview_calculado=false;
		break;
	}
	return;
}

bool gluInvertMatrix(const GLfloat  m[16], GLfloat  invOut[16])
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	double inv[16], det;
	int i;

	if (PrintInvert) cout<<"I gluInvertMatrix:1"<<endl;
	inv[0] = m[5]  * m[10] * m[15] -
			m[5]  * m[11] * m[14] -
			m[9]  * m[6]  * m[15] +
			m[9]  * m[7]  * m[14] +
			m[13] * m[6]  * m[11] -
			m[13] * m[7]  * m[10];

	inv[4] = -m[4]  * m[10] * m[15] +
			m[4]  * m[11] * m[14] +
			m[8]  * m[6]  * m[15] -
			m[8]  * m[7]  * m[14] -
			m[12] * m[6]  * m[11] +
			m[12] * m[7]  * m[10];

	inv[8] = m[4]  * m[9] * m[15] -
			m[4]  * m[11] * m[13] -
			m[8]  * m[5] * m[15] +
			m[8]  * m[7] * m[13] +
			m[12] * m[5] * m[11] -
			m[12] * m[7] * m[9];

	inv[12] = -m[4]  * m[9] * m[14] +
			m[4]  * m[10] * m[13] +
			m[8]  * m[5] * m[14] -
			m[8]  * m[6] * m[13] -
			m[12] * m[5] * m[10] +
			m[12] * m[6] * m[9];

	inv[1] = -m[1]  * m[10] * m[15] +
			m[1]  * m[11] * m[14] +
			m[9]  * m[2] * m[15] -
			m[9]  * m[3] * m[14] -
			m[13] * m[2] * m[11] +
			m[13] * m[3] * m[10];

	inv[5] = m[0]  * m[10] * m[15] -
			m[0]  * m[11] * m[14] -
			m[8]  * m[2] * m[15] +
			m[8]  * m[3] * m[14] +
			m[12] * m[2] * m[11] -
			m[12] * m[3] * m[10];

	if (PrintInvert) cout<<"I gluInvertMatrix:10"<<endl;
	inv[9] = -m[0]  * m[9] * m[15] +
			m[0]  * m[11] * m[13] +
			m[8]  * m[1] * m[15] -
			m[8]  * m[3] * m[13] -
			m[12] * m[1] * m[11] +
			m[12] * m[3] * m[9];

	inv[13] = m[0]  * m[9] * m[14] -
			m[0]  * m[10] * m[13] -
			m[8]  * m[1] * m[14] +
			m[8]  * m[2] * m[13] +
			m[12] * m[1] * m[10] -
			m[12] * m[2] * m[9];

	inv[2] = m[1]  * m[6] * m[15] -
			m[1]  * m[7] * m[14] -
			m[5]  * m[2] * m[15] +
			m[5]  * m[3] * m[14] +
			m[13] * m[2] * m[7] -
			m[13] * m[3] * m[6];

	inv[6] = -m[0]  * m[6] * m[15] +
			m[0]  * m[7] * m[14] +
			m[4]  * m[2] * m[15] -
			m[4]  * m[3] * m[14] -
			m[12] * m[2] * m[7] +
			m[12] * m[3] * m[6];

	inv[10] = m[0]  * m[5] * m[15] -
			m[0]  * m[7] * m[13] -
			m[4]  * m[1] * m[15] +
			m[4]  * m[3] * m[13] +
			m[12] * m[1] * m[7] -
			m[12] * m[3] * m[5];

	inv[14] = -m[0]  * m[5] * m[14] +
			m[0]  * m[6] * m[13] +
			m[4]  * m[1] * m[14] -
			m[4]  * m[2] * m[13] -
			m[12] * m[1] * m[6] +
			m[12] * m[2] * m[5];

	if (PrintInvert) cout<<"I gluInvertMatrix:20"<<endl;
	inv[3] = -m[1] * m[6] * m[11] +
			m[1] * m[7] * m[10] +
			m[5] * m[2] * m[11] -
			m[5] * m[3] * m[10] -
			m[9] * m[2] * m[7] +
			m[9] * m[3] * m[6];

	inv[7] = m[0] * m[6] * m[11] -
			m[0] * m[7] * m[10] -
			m[4] * m[2] * m[11] +
			m[4] * m[3] * m[10] +
			m[8] * m[2] * m[7] -
			m[8] * m[3] * m[6];

	inv[11] = -m[0] * m[5] * m[11] +
			m[0] * m[7] * m[9] +
			m[4] * m[1] * m[11] -
			m[4] * m[3] * m[9] -
			m[8] * m[1] * m[7] +
			m[8] * m[3] * m[5];

	inv[15] = m[0] * m[5] * m[10] -
			m[0] * m[6] * m[9] -
			m[4] * m[1] * m[10] +
			m[4] * m[2] * m[9] +
			m[8] * m[1] * m[6] -
			m[8] * m[2] * m[5];

	if (PrintInvert) cout<<"I gluInvertMatrix:30"<<endl;
	det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

	if (det == 0) {
		if (PrintInvert) cout <<"Matriz No invertible"<<endl;
		return false;
	}

	if (PrintInvert) cout<<"I gluInvertMatrix:40: det="<<det<<endl;
	det = 1.0 / det;


	if (PrintInvert) cout<<"I gluInvertMatrix:50"<<endl;
	for (i = 0; i < 16; i++)
		invOut[i] = inv[i] * det;

	if (PrintInvert) cout<<"I gluInvertMatrix:60"<<endl;
	return true;
}


bool InvierteMatriz(const GLdouble  m[16], GLdouble  invOut[16])
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	double inv[16], det;
	int i;

	inv[0] = m[5]  * m[10] * m[15] -
			m[5]  * m[11] * m[14] -
			m[9]  * m[6]  * m[15] +
			m[9]  * m[7]  * m[14] +
			m[13] * m[6]  * m[11] -
			m[13] * m[7]  * m[10];

	inv[4] = -m[4]  * m[10] * m[15] +
			m[4]  * m[11] * m[14] +
			m[8]  * m[6]  * m[15] -
			m[8]  * m[7]  * m[14] -
			m[12] * m[6]  * m[11] +
			m[12] * m[7]  * m[10];

	inv[8] = m[4]  * m[9] * m[15] -
			m[4]  * m[11] * m[13] -
			m[8]  * m[5] * m[15] +
			m[8]  * m[7] * m[13] +
			m[12] * m[5] * m[11] -
			m[12] * m[7] * m[9];

	inv[12] = -m[4]  * m[9] * m[14] +
			m[4]  * m[10] * m[13] +
			m[8]  * m[5] * m[14] -
			m[8]  * m[6] * m[13] -
			m[12] * m[5] * m[10] +
			m[12] * m[6] * m[9];

	inv[1] = -m[1]  * m[10] * m[15] +
			m[1]  * m[11] * m[14] +
			m[9]  * m[2] * m[15] -
			m[9]  * m[3] * m[14] -
			m[13] * m[2] * m[11] +
			m[13] * m[3] * m[10];

	inv[5] = m[0]  * m[10] * m[15] -
			m[0]  * m[11] * m[14] -
			m[8]  * m[2] * m[15] +
			m[8]  * m[3] * m[14] +
			m[12] * m[2] * m[11] -
			m[12] * m[3] * m[10];

	inv[9] = -m[0]  * m[9] * m[15] +
			m[0]  * m[11] * m[13] +
			m[8]  * m[1] * m[15] -
			m[8]  * m[3] * m[13] -
			m[12] * m[1] * m[11] +
			m[12] * m[3] * m[9];

	inv[13] = m[0]  * m[9] * m[14] -
			m[0]  * m[10] * m[13] -
			m[8]  * m[1] * m[14] +
			m[8]  * m[2] * m[13] +
			m[12] * m[1] * m[10] -
			m[12] * m[2] * m[9];

	inv[2] = m[1]  * m[6] * m[15] -
			m[1]  * m[7] * m[14] -
			m[5]  * m[2] * m[15] +
			m[5]  * m[3] * m[14] +
			m[13] * m[2] * m[7] -
			m[13] * m[3] * m[6];

	inv[6] = -m[0]  * m[6] * m[15] +
			m[0]  * m[7] * m[14] +
			m[4]  * m[2] * m[15] -
			m[4]  * m[3] * m[14] -
			m[12] * m[2] * m[7] +
			m[12] * m[3] * m[6];

	inv[10] = m[0]  * m[5] * m[15] -
			m[0]  * m[7] * m[13] -
			m[4]  * m[1] * m[15] +
			m[4]  * m[3] * m[13] +
			m[12] * m[1] * m[7] -
			m[12] * m[3] * m[5];

	inv[14] = -m[0]  * m[5] * m[14] +
			m[0]  * m[6] * m[13] +
			m[4]  * m[1] * m[14] -
			m[4]  * m[2] * m[13] -
			m[12] * m[1] * m[6] +
			m[12] * m[2] * m[5];

	inv[3] = -m[1] * m[6] * m[11] +
			m[1] * m[7] * m[10] +
			m[5] * m[2] * m[11] -
			m[5] * m[3] * m[10] -
			m[9] * m[2] * m[7] +
			m[9] * m[3] * m[6];

	inv[7] = m[0] * m[6] * m[11] -
			m[0] * m[7] * m[10] -
			m[4] * m[2] * m[11] +
			m[4] * m[3] * m[10] +
			m[8] * m[2] * m[7] -
			m[8] * m[3] * m[6];

	inv[11] = -m[0] * m[5] * m[11] +
			m[0] * m[7] * m[9] +
			m[4] * m[1] * m[11] -
			m[4] * m[3] * m[9] -
			m[8] * m[1] * m[7] +
			m[8] * m[3] * m[5];

	inv[15] = m[0] * m[5] * m[10] -
			m[0] * m[6] * m[9] -
			m[4] * m[1] * m[10] +
			m[4] * m[2] * m[9] +
			m[8] * m[1] * m[6] -
			m[8] * m[2] * m[5];

	det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

	if (det == 0)
		return false;

	det = 1.0 / det;

	for (i = 0; i < 16; i++)
		invOut[i] = inv[i] * det;

	return true;
}



void CB_mouse(int button, int state, int x, int y)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	xPos0 = x;
	yPos0 = y;

#if DBGMouse
	myfileDBG<<"CB_mouse: B="<<button<<" state="<<state<<" x="<<x<<" y="<<y<<endl;
#endif

	//cout<<"M CB_mouse(int "<<button<<", int "<<state<<", int "<<x<<", int "<<y<<")"<<endl;
	//	printf("CB_mouse: button=%d, state=%d, x=%d, y=%d\n",button,state,x,y);
	KeyControlAltShift = glutGetModifiers();
	//printf("glutGetModifiers=%d",KeyControlAltShift);	cout<<endl;


	//cout<<"CB_mouse: B="<<button<<" state="<<state<<" x="<<x<<" y="<<y<<endl;
	if (state==0 & x>x0Menu & x< x1Menu & y>y0Menu & y < y1Menu) {
		EstadoMenu=1-EstadoMenu;
		if (EstadoMenu) {
			glutAttachMenu(GLUT_RIGHT_BUTTON);
			//cout<<"CB_mouse: B=glutAttachMenu(GLUT_RIGHT_BUTTON);   EstadoMenu="<<EstadoMenu<<endl;
		}
		else {
			glutDetachMenu(GLUT_RIGHT_BUTTON);
			//cout<<"CB_mouse: B=glutDetachMenu(GLUT_RIGHT_BUTTON);  EstadoMenu="<<EstadoMenu<<endl;
		}
		//Lanza_thread();
	}
	if (state==GLUT_UP) {iPush=0;}
	else {
		switch (button)
		{
		case GLUT_LEFT_BUTTON:

			if (MODO_Rotacion2==0) iPush=1;
			else {
				//cout<<"MODO_Rotacion2=="<<MODO_Rotacion2<<endl;
				iPush=3;
			}

			if (KeyControlAltShift==GLUT_ACTIVE_CTRL) iPush=2;


			break;

		case GLUT_RIGHT_BUTTON:
			iPush=3;
			break;
		case GLUT_MIDDLE_BUTTON:

			iPush=2;
			break;

		case 3:

			Escala *=1.1;
			break;
		case 4:

			Escala /=1.1;
			break;


		}

		if (iPush==2) {
			if (Add_Particulas) {
				iPush=0;
				//Convertir coordenadas de la pantalla a OpenGL.....


				//float x=xPos0,y=yPos0;
				GLint viewport[4];
				GLdouble modelview[16];
				GLdouble projection[16];
				GLfloat winX, winY, winZ;
				GLdouble posX, posY, posZ;

				FuncionesOpenGL::ObtieneMatrices();

				winX = (float)x;
				winY = (float)FuncionesOpenGL::viewport[3] - (float)y;
				glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );

				FuncionesOpenGL::Win2World(winX, winY, winZ, &posX, &posY, &posZ);

				posZ /=FactorZ;

				printf("%f %f %f\n",posX,posY,posZ);


				gtotal->Particulas_Init_3(posX,posY,posZ);




			} else if (Add_Voronoi||Add_VolumenINI||Flag_MuestraCara||Add_Poligono_Selected) {
				iPush=0;
				//Convertir coordenadas de la pantalla a OpenGL.....


				//float x=xPos0,y=yPos0;
				GLint viewport[4];
				GLdouble modelview[16];
				GLdouble projection[16];
				GLfloat winX, winY, winZ;
				GLdouble posX, posY, posZ;

				FuncionesOpenGL::ObtieneMatrices();

				winX = (float)x;
				winY = (float)FuncionesOpenGL::viewport[3] - (float)y;
				glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );

				FuncionesOpenGL::Win2World(winX, winY, winZ, &posX, &posY, &posZ);
				posZ /=FactorZ;

				printf("Wind=[%.2f %.2f %.4f]\n",winX,winY,winZ);
				printf("Worl=[%.2f %.2f %.4f]\n",posX,posY,posZ);

				if (Add_Voronoi)
					gtotal->SeleccionaVolFinito(posX,posY,posZ);
				else if(Flag_MuestraCara)
					gtotal->SeleccionaYMuestraCara(posX, posY, posZ);
				else if(Add_Poligono_Selected) {
					if (gPoly != NULL)
						gPoly->SeleccionaYMuestraPolygono(posX, posY, posZ);
					if(gtotal!=NULL && gtotal->nCaras>0)
						gtotal->SeleccionaYMuestraCara(posX, posY, posZ);
				}
				else
					gtotal->SeleccionaVolFinito(posX,posY,posZ);





			} else if (MueveCentro) {
				iPush=0;
				//Convertir coordenadas de la pantalla a OpenGL.....


				//float x=xPos0,y=yPos0;
				GLfloat winX, winY, winZ;
				GLdouble posX, posY, posZ;


				FuncionesOpenGL::ObtieneMatrices();

				winX = (float)x;
				winY = (float)FuncionesOpenGL::viewport[3] - (float)y;
				//				winY = -(float)y;
				glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );

				FuncionesOpenGL::Win2World(winX, winY, winZ, &posX, &posY, &posZ);
				posZ ;
				//gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);

				//cout<<"M posX,Y,Z="<<posX<<" "<<posY<<" "<<posZ<<endl;
				vecDXEsfera[0]=posX-vecXEsfera[0];
				vecDXEsfera[1]=posY-vecXEsfera[1];
				vecDXEsfera[2]=posZ-vecXEsfera[2];
				vecXEsfera[0]=posX;
				vecXEsfera[1]=posY;
				vecXEsfera[2]=posZ;

				//cout<<"M vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<"]"<<endl;
				gluInvertMatrix(MatrizRotacionGlobal,MatrizRotacionGlobalINV);
				MatrizXvector4(MatrizRotacionGlobal,vecDXEsfera,vecDUEsfera);
				vecUEsfera[0]+=vecDUEsfera[0]*Escala;
				vecUEsfera[1]+=vecDUEsfera[1]*Escala;
				vecUEsfera[2]+=vecDUEsfera[2]*Escala;
				//cout<<"M vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<"]"<<endl;
			}
		}

	}

	//cout<<"CB_mouse()END: iPush="<<iPush<<endl;

}

void CB_motion(int x, int y)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;
	//	printf("CB_motion: iPush=%d, x=%d, y=%d",iPush,x,y);cout<<endl;
	float dx,dy;
	GLdouble winX,winY,winZ;
	xPos = x;
	yPos = y;

	if (iPush==0) return;
#if DBGMouse
	myfileDBG<<"CB_motion: x="<<x<<" y="<<y<<endl;
#endif

	int lMODO_de_mover=MODO_de_mover;



	//cout<<"iPush="<<iPush<<" lMODO_de_mover="<<lMODO_de_mover<<" KeyControlAltShift="<<KeyControlAltShift<<endl;

	if (KeyControlAltShift==GLUT_ACTIVE_CTRL)
		lMODO_de_mover=0;
	if (iPush==2) lMODO_de_mover=3;
	if (iPush==3) lMODO_de_mover=0;
	switch (lMODO_de_mover)
	{
	case 0: //Trasladar
		dx=(xPos-xPos0+0.0);
		dy=-(yPos-yPos0+0.0);
		double dxE,dyE,dzE,xP,yP,zP;

		if (PrintMouse) cout <<"M \ndxy=["<<dx<<" , "<<dy<<endl;
		//cout<<"M A1 vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<","<<vecUEsfera[3]<<"]"<<endl;


		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective((GLdouble)GlobalFovy, aspect, (GLdouble)1, (GLdouble)100.0);
		glTranslated( 0, 0, -10);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslated( 0, -1, -10);

		FuncionesOpenGL::World2Win((GLdouble)vecUEsfera[0], (GLdouble)vecUEsfera[1], (GLdouble)vecUEsfera[2],
				&xP,&yP,&zP);
		//cout<<"M xyzP=[ "<<xP<<" , "<<yP<<" , "<<zP<<" ]"<<endl;

		FuncionesOpenGL::Win2World(xP+dx, yP+dy, zP, &dxE, &dyE, &dzE);
		dzE /=FactorZ;
		if (PrintMouse) cout <<"M dxyzE=["<<dxE<<" , "<<dyE<<" , "<<dzE<<endl;

		vecUEsfera[0] =dxE;
		vecUEsfera[1] =dyE;
		vecUEsfera[2] =dzE;
		//cout<<"M A2 vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<","<<vecUEsfera[3]<<"]"<<endl;
		FuncionesOpenGL::modelview_calculado=false;
		break;
	case 1: //Rotar
		rotacion(yPos-yPos0, xPos-xPos0);


#if DBGMouse
		myfileDBG<<"Rotacion Mouse"<<endl;
		MyDBG();
#endif
		break;
	case 2: //Spin en torno a Z
		int dx1,dy1,x1,y1;
		dx1=xPos-xPos0;dy1=yPos-yPos0;
		x1=xPos0-width/2;y1=yPos0-height/2;

		FuncionesOpenGL::World2Win(vecXEsfera[0],vecXEsfera[1],vecXEsfera[2], &winX, &winY, &winZ);
		printf("xPos0=%d,winX=%f,yPos0=%d,winY=%f\n",xPos0,winX,yPos0,winY);
		x1=xPos0-winX;y1=yPos0-(height-winY);

		printf("x1=%d,y1=%d\n",x1,y1);

		wAngleZ=-(dx1*(-y1)+dy1*x1+0.0)/(x1*x1+y1*y1+1.0)*180/3.14;
		glLoadIdentity();
		glRotatef(wAngleZ, 0.0f, 0.0f, 1.0f);
		glMultMatrixf((GLfloat *)MatrizRotacionGlobal);
		glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat*)MatrizRotacionGlobal);
		glLoadIdentity();
		glMultMatrixf((GLfloat *)MatrizRotacionGlobalINV);
		glRotatef(-wAngleZ, 0.0f, 0.0f, 1.0f);
		glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat*)MatrizRotacionGlobalINV);
		FuncionesOpenGL::modelview_calculado=false;
		break;

	case 3:
		Escala *= 1-2*(yPos-yPos0+0.0)/height;
		FuncionesOpenGL::modelview_calculado=false;
		break;

	}
	yPos0=yPos;
	xPos0=xPos;
	return;
}

void CB_keyboardSpecial( int key, int x, int y )
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	// manejo de teclas especiales
	printf("CB_keyboardSpecial: char=%d, x=%d, y=%d\n",key,x,y);cout<<endl;
	switch (key) {
	case 100:
		rotacion(0,-1);
		break;
	case 101:
		rotacion(-1,0);
		break;
	case 102:
		rotacion(0,1);
		break;
	case 103:
		rotacion(1,0);
		break;

	}
	glutPostRedisplay();
}


void idleevent()
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	static int generador=0;
	if(DBG)cout<<"idleevent()"<<endl;

	if ( glutGetWindow() != main_window ) 
		glutSetWindow(main_window);  


	if(DBG)cout<<"idleevent(): Falta_llamar_etapaS="<<llamar_etapa_Siguiente_PAUSA2<<endl;
	if(llamar_etapa_Siguiente_PAUSA2){
		cout<<"idleevent(): 1938"<<endl;
		Calculo_EtapaS();
	}


	if(DBG)cout<<"idleevent(): CalculoContinuo="<<CalculoContinuo<<endl;
	if (CalculoContinuo) {
		if (Etapa==0){
			cout<<"idleevent(): 1946"<<endl;
			Calculo_EtapaS(); //Lanzar la primera etapa

			if (nRLC>0) nR=nRLC;
			if (nThLC>0) nTh=nThLC;

			time_t _tm =time(NULL );

			struct tm * curtime = localtime ( &_tm );
			sprintf(nombre0,"Soluciones/Sal_%02d.%02d_%02d.%02d.%02d_Malla=%d_nR=%d_nTg=%d_NDivZ=%d_",
					curtime->tm_mon,curtime->tm_mday,
					curtime->tm_hour,curtime->tm_min,curtime->tm_sec,CualMalla,nR,nTh,NDivZ);
		}


		cout<<"idleevent(): 1961"<<endl;
		Calculo_EtapaS();

		char nombre[100];
		sprintf(nombre,"%s_Size=%d_Etapa=%d.dat",nombre0,gtotal->nH3D,Etapa);
		SaveOrRead(nombre,1);
	} else	if (Etapa==0){
		cout<<"idleevent(): 1968"<<endl;
		Calculo_EtapaS(); //Lanzar la primera etapa
	}


	generador++;
	if (generador >0) {
		generador=0;
		//		g.GeneraCaras(g.nH3D);
		//		g.GeneraCaras();

	}
	glutPostRedisplay();

	if(DBG)cout<<"idleevent():END"<<endl;
}

void replaceAll(std::string& str, const std::string& from, const std::string& to)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	if(from.empty())
		return;
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
	}
}

void LecturaArchivoDeDatos()
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	ifstream myfile;
	ofstream myfile2;
	int numlinea=0;

	char mensajeE[10000];
	myfile2.open("dondeestoy.txt");
	myfile2<<"hola"<<endl;
	myfile2.close();
#if DBG == 1
	myfile2.open("dondeestoy.txt");
	myfile2<<"hola"<<endl;
	myfile2.close();
#endif
	if (PrintFunciones) cout<<"LecturaArchivoDeDatos():"<<FDatos<<endl;
	sprintf(mensajes0,"%sLecturaArchivoDeDatos():%s\n",mensajes0,FDatos);

	sprintf(mensajeE,"Archivo: \"%s\" \n",FDatos);

	myfile.open (FDatos);
	string linea;
	int errordatos=0;

	while (getline(myfile,linea)) {
		numlinea++;
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
		else if (snombre=="TipoCalculo"	    ) {TipoCalculo	    =ff;}
		else if (snombre=="dt"				) {Datos_dt			=ff;}
		else if (snombre=="Tmax"			) {Datos_Tmax		=ff;}
		else if (snombre=="nR"				) {nR				=ff;}
		else if (snombre=="nR2"				) {nR2				=ff;}
		else if (snombre=="nTh"				) {nTh				=ff;}
		else if (snombre=="NDivZ"			) {NDivZ			=ff;}
		else if (snombre=="Dominio_Xmax"	) {Dominio_Xmax		=ff;}
		else if (snombre=="Dominio_Hmax"	) {Dominio_Hmax		=ff;}
		else if (snombre=="Dominio_Hsup"	) {Dominio_Hsup		=ff;}
		else if (snombre=="Dominio_Rint"	) {Dominio_Rint		=ff;}
		else if (snombre=="Dominio_Rmed"	) {Dominio_Rmed		=ff;}
		else if (snombre=="rhof"			) {Datos_rhof		=ff;}
		else if (snombre=="rhos"			) {Datos_rhos		=ff;}
		else if (snombre=="phi"				) {Datos_phi		=ff;}
		else if (snombre=="cf"				) {Datos_cf			=ff;}
		else if (snombre=="cs"				) {Datos_cs			=ff;}
		else if (snombre=="Tinyeccion"		) {Datos_Tinyeccion	=ff;}
		else if (snombre=="Tambiente"		) {Datos_Tambiente	=ff;}
		else if (snombre=="TLimite"		    ) {TLimite       	=ff;}
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
		else if (snombre=="Solucion_error") {ErrorMax_Cada_t=ff;}
		else {
			cout <<"Comando en Archivo de datos no considerado:"<<endl;
			cout <<"variable="<< snombre<<"\t="<<ff<<endl;
			sprintf(mensajeE,"%sLinea=%d, Variable \"%s\" no tiene sentido (= %f)\n",mensajeE,numlinea,snombre.c_str(),ff);
			cout <<"variable="<< snombre<<"\t="<<ff<<endl;
			errordatos=1;
		}

	}


	myfile.close();
	if (errordatos==1) {

#if defined _WIN32 || defined _WIN64s
		MessageBox( 0, mensajeE, "Error Grave en Archivo de Datos", 0 );
#endif
		exit(1);
	}

}

void LecturaArchivoDeDatos_Post(ostream& os)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	if (PrintParametros) {
		os<<"#Parametros usados:\n"<<endl;


		os<<"malla              ="<<CualMalla       <<endl;
		os<<"CalculoContinuo    ="<<CalculoContinuo	<<endl;
		os<<"TipoCalculo        ="<<TipoCalculo	    <<endl;
		os<<"dt                 ="<<Datos_dt		<<endl;
		os<<"Tmax               ="<<Datos_Tmax		<<endl;
		os<<"nR                 ="<<nR				<<endl;
		os<<"nR2                ="<<nR2				<<endl;
		os<<"nTh                ="<<nTh				<<endl;
		os<<"NDivZ              ="<<NDivZ			<<endl;
		os<<"Dominio_Xmax       ="<<Dominio_Xmax	<<endl;
		os<<"Dominio_Rint       ="<<Dominio_Rint	<<endl;
		os<<"Dominio_Rmed       ="<<Dominio_Rmed	<<endl;
		os<<"Dominio_Hmax       ="<<Dominio_Hmax	<<endl;
		os<<"Dominio_Hsup       ="<<Dominio_Hsup	<<endl;
		os<<"rhof               ="<<Datos_rhof		<<endl;
		os<<"rhos               ="<<Datos_rhos		<<endl;
		os<<"phi                ="<<Datos_phi		<<endl;
		os<<"cf                 ="<<Datos_cf		<<endl;
		os<<"cs                 ="<<Datos_cs		<<endl;
		os<<"Tinyeccion         ="<<Datos_Tinyeccion<<endl;
		os<<"Tambiente          ="<<Datos_Tambiente	<<endl;
		os<<"TLimite            ="<<TLimite     	<<endl;
		os<<"hc_superior        ="<<Datos_hc_superior<<endl;
		os<<"hc_inferior        ="<<Datos_hc_inferior<<endl;
		os<<"hc_lateral         ="<<Datos_hc_lateral	<<endl;
		os<<"KTermofilm         ="<<Datos_KTermofilm	<<endl;
		os<<"eTermofilm         ="<<Datos_eTermofilm	<<endl;
		os<<"VinyeccionLtsHr    ="<<VinyeccionLtsHr	<<endl;
		os<<"kf                 ="<<Datos_kf			<<endl;
		os<<"ks                 ="<<Datos_ks				<<endl;
		os<<"Kevaporacion       ="<<Datos_Kevaporacion	<<endl;
		os<<"HumedadAmbiental   ="<<Datos_HumedadAmbiental<<endl;
		os<<"DistanciaAlBorde   ="<<Datos_DistanciaAlBorde<<endl;
		os<<"CalorLatenteFluido ="<<Datos_CalorLatenteFluido<<endl;
		os<<"Solucion_error     ="<<ErrorMax_Cada_t<<endl;
	}

	Datos_KtEt=Datos_KTermofilm/Datos_eTermofilm;
	Datos_htilde=(Datos_hc_superior*Datos_KtEt)/(Datos_hc_superior+Datos_KtEt);
	Datos_Vinyeccion=VinyeccionLtsHr/1000/3600; //% m^3/m^2/s
	Datos_km_mu=Datos_phi*Datos_kf+(1-Datos_phi)*Datos_ks;
	Datos_km2=(1-Datos_phi)*Datos_ks;
	Datos_rhom=Datos_phi*Datos_rhof+(1-Datos_phi)*Datos_rhos;
	Datos_cm=Datos_phi*Datos_cf+(1-Datos_phi)*Datos_cs;

	if (PrintParametros) {
		os<<"\n#Parametros derivados:\n"<<endl;
		os<<"# KtEt      =KTermofilm/eTermofilm                    = "<<Datos_KtEt<<endl;
		os<<"# htilde    =(hc_superior*KtEt)/(hc_superior+KtEt)    = "<<Datos_htilde<<endl;
		os<<"# Vinyeccion=VinyeccionLtsHr/1000/3600; //% m^3/m^2/s = "<<Datos_Vinyeccion<<endl;
		os<<"# km        =phi*kf+(1-phi)*ks                        = "<<Datos_km_mu<<endl;
		os<<"# km2       =(1-phi)*ks                               = "<<Datos_km2<<endl;
		os<<"# rhom      =phi*rhof+(1-phi)*rhos                    = "<<Datos_rhom<<endl;
		os<<"# cm        =phi*cf  +(1-phi)*cs                      = "<<Datos_cm<<endl;
	}


}
void SaveOrRead(char *ifile_name,int iSaveReadMode)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	char leyendo[100];
	char leyendo2[100]={'aa'};
	char ArchivoParametros[200];
	int i=0;

	if (iSaveReadMode==1) { //Save


		//Busca Nombre pila del archivo y genera el nombre de los parametros
		strcpy(leyendo2,ifile_name);
		char *st=leyendo2;

		while(strncmp(st,"Etapa",5)!=0 && i<100) {st++;i++;}
		cout<<"i="<<i<<endl;
		if (i==100) {
			sprintf(leyendo2,"%s_Etapa=%02d.dat",leyendo2,EtapaGlobal);
			sprintf(ifile_name,"%s_Etapa=%02d.dat",ifile_name,EtapaGlobal);
			st=leyendo2;i=0;
			while(strncmp(st,"Etapa",5)!=0 && i<100) {st++;i++;}
		}
		cout<<"i="<<i<<endl;
		*st=0;
		cout<<"leyendo2="<<leyendo2<<endl;
		sprintf(FDatos,"%sEtapa=00_Parametros.dat",leyendo2);
		cout<<FDatos<<endl;

		rmdir(FDatos);
		ofstream myfile;
		myfile.open (FDatos);


		LecturaArchivoDeDatos_Post(myfile);
		myfile.close();

		myfile.open (ifile_name);

		myfile.precision(17);
		VersionDatos=0.20;
		myfile<<"#Version "<<VersionDatos<<" "<<CualMalla<<endl;
		cout<<"#Version "<<VersionDatos<<" "<<CualMalla<<endl;
		if(VersionDatos>=0.1999)	{
			myfile<<EtapaGlobal<<"  "<<TiempoCalculo<<endl;
		}
		gtotal->write(myfile);
		myfile<<Etapa<<endl;
		myfile<<InicioEtapa<<endl;
		save(F,myfile);
		save(PotencialPsi,myfile);
		save(Temperatura,myfile);
		save(F2Nodos,myfile);
		save(PotencialVInfiltracion,myfile);
		save(TempPilaBloques,myfile);
		//		save(U,myfile);
		//		save(V,myfile);
		//		save(W,myfile);
		myfile.close();
	} else {   //Read
		int num1,num2;
		ifstream myfile;


#if DBGMouse
		myfileDBG << "Voy a leer"<<endl;
		MyDBG();
#endif

		//Busca Nombre pila del archivo y genera el nombre de los parametros
		strcpy(leyendo2,ifile_name);
		char *st=leyendo2;

		while(strncmp(st,"Etapa",5)!=0 && i<100) {st++;i++;}

		cout<<"i="<<i<<endl;
		*st=0;
		cout<<"leyendo2="<<leyendo2<<endl;
		sprintf(FDatos,"%sEtapa=00_Parametros.dat",leyendo2);
		cout<<FDatos<<endl;

		//Leo archivo de Parametros

		cout<<"Leo archivo de Parametros..."<<endl;

		LecturaArchivoDeDatos();
		LecturaArchivoDeDatos_Post(cout);
		CalculoContinuo=0;

		// Leo archivo de resultados0
		cout<<"Leo archivo de resultados..."<<endl;
		time_t _tm0 =time(NULL );

		cout<<"1:"<<ifile_name<<endl;
		myfile.open (ifile_name);
		cout<<"2"<<endl;
		myfile.ignore(100,' ');
		cout<<"3"<<endl;
		myfile >>VersionDatos>>CualMalla;
		cout<<"4"<<endl;
		cout<<"VersionDatos="<<VersionDatos<<"\tCualMalla="<<CualMalla<<endl;

		if(VersionDatos>=0.1999)	{
			myfile>>EtapaGlobal>>TiempoCalculo;
			cout<<"EtapaGlobal(VersionDatos>=0.20)="<<EtapaGlobal<<endl;
			cout<<"TiempoCalculo(VersionDatos>=0.20)="<<TiempoCalculo<<endl;
		}


		gtotal->read(myfile);
		if (DBG) cout<<"Malla Ok"<<endl;
		myfile>>Etapa;
		if (PrintInvert) cout<<"InicioEtapa(antes)="<<InicioEtapa<<endl;
		if(VersionDatos>0.109)	{
			myfile>>InicioEtapa;
			if (PrintInvert) cout<<"InicioEtapa(VersionDatos>0.109)="<<InicioEtapa<<endl;
		}
		else
		{
			InicioEtapa=Etapa;
			if (PrintInvert) cout<<"InicioEtapa(=Etapa)="<<InicioEtapa<<endl;
		}

		read(F,myfile);
		if (DBG) cout<<"F Ok"<<endl;
		read(PotencialPsi,myfile);
		if (DBG) cout<<"FF Ok"<<endl;
		read(Temperatura,myfile);
		if (DBG) cout<<"FF2 Ok"<<endl;
		read(F2Nodos,myfile);
		if (DBG) cout<<"F2Nodos Ok"<<endl;
		read(PotencialVInfiltracion,myfile);
		if (DBG) cout<<"PotencialVInfiltracion Ok"<<endl;
		read(TempPilaBloques,myfile);
		if (DBG) cout<<"TempPilaBloques Ok"<<endl;
		//		read(U,myfile);
		//		if (DBG) cout<<"U Ok"<<endl;
		//		read(V,myfile);
		//		if (DBG) cout<<"V Ok"<<endl;
		//		read(W,myfile);
		//		if (DBG) cout<<"W Ok"<<endl;
		myfile.close();

		err0=1e-8;
		switch (CualMalla) {
		case 1:
			ThetaMax=atan(1)*2;
			if (Etapa==7) InicioEtapa=Etapa-1;
			break;
		}

		time_t _tm1 =time(NULL );

		cout<< "lectura en "<<_tm1-_tm0<<" seg"<<endl;
		if (PanelFlimite != NULL) {
			control_cb(1105);
			PanelFlimite->execute_callback();
			PanelFlimite->enable();
		}


#if DBGMouse
		myfileDBG << "Ya Lei"<<endl;
		MyDBG();
#endif

	}
}


void visible(int vis)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	if (PrintFunciones) cout<<"visible()"<<endl;

	if (vis == GLUT_VISIBLE) {

		glutIdleFunc(idleevent);

	} else {
		glutIdleFunc(NULL);
	}
	if (PrintFunciones) cout<<"visible():END"<<endl;
}

double HumedadSaturacion(double T)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	double TVapor[25],PVapor[25];

	TVapor[ 0]=0  ; PVapor[ 0]=0.0060;
	TVapor[ 1]=5  ; PVapor[ 1]=0.0086;
	TVapor[ 2]=10 ; PVapor[ 2]=0.0121;
	TVapor[ 3]=15 ; PVapor[ 3]=0.0168;
	TVapor[ 4]=20 ; PVapor[ 4]=0.0231;
	TVapor[ 5]=25 ; PVapor[ 5]=0.0313;
	TVapor[ 6]=30 ; PVapor[ 6]=0.0419;
	TVapor[ 7]=35 ; PVapor[ 7]=0.0555;
	TVapor[ 8]=40 ; PVapor[ 8]=0.0728;
	TVapor[ 9]=45 ; PVapor[ 9]=0.0946;
	TVapor[10]=50 ; PVapor[10]=0.1217;
	TVapor[11]=55 ; PVapor[11]=0.1553;
	TVapor[12]=60 ; PVapor[12]=0.1966;
	TVapor[13]=65 ; PVapor[13]=0.2468;
	TVapor[14]=70 ; PVapor[14]=0.3075;
	TVapor[15]=75 ; PVapor[15]=0.3804;
	TVapor[16]=80 ; PVapor[16]=0.4673;
	TVapor[17]=85 ; PVapor[17]=0.5706;
	TVapor[18]=90 ; PVapor[18]=0.6918;
	TVapor[19]=95 ; PVapor[19]=0.8341;
	TVapor[20]=100; PVapor[20]=0.9999;

	int i;
	double Pt,Psi;
	if (T< TVapor[0]) Pt=PVapor[ 0];
	else if (T>= TVapor[20]) Pt=PVapor[ 20];
	else
		for (i=0;i<20;i++) {
			if (TVapor[i]<=T && T< TVapor[i+1] ) {
				double lambda=(T-TVapor[i])/(TVapor[i+1]-TVapor[i]);
				Pt=PVapor[i]+lambda*(PVapor[i+1]-PVapor[i]);
			}
		}

	Psi=0.6242*Pt / (1-Pt);
	return(Psi);


}


double Psi_ev(double Tr)
{

	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	double Psi=Datos_Kevaporacion*(HumedadSaturacion(Tr)-Datos_HumedadAmbiental*HumedadSaturacion(Datos_Tambiente));
	return(Psi);
}

double CalculaMpunto(vector<double> &F,vector<double> &Temp,vector<double> &U,vector<double> &V,
		vector<double> &W,grid3D g,int niteraciones,double err0)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	int iiter,i,j,iC;
	double SumaCoeffDif_F,SumaCoeffDif,Coeff_Difusion,err,errG,nF;
	if(DBG) cout<<"CalculaMpunto"<<endl;
	if (F.size() != g.nH3D) F.resize(g.nH3D);
	if (U.size() != g.nH3D) U.resize(g.nH3D);
	if (V.size() != g.nH3D) V.resize(g.nH3D);
	if (W.size() != g.nH3D) W.resize(g.nH3D);
	//Resuelvo la ecuación varias veces
	niteraciones=g.nH3D*10;
	if (niteraciones<30000) niteraciones=30000;
	cout <<"niteraciones="<<niteraciones<<endl;
	myfileSalida <<"niteraciones="<<niteraciones<<endl;

	if(DBG) cout<<"CalculaMpunto(): Potencial"<<endl;
	for (iiter=0 ; iiter<niteraciones ; iiter++) {
		errG=0;
		int i0=g.nH3D-1,istep=-1;
		for (i=i0 ; i<g.nH3D && i>= 0 ; i+=istep) {
			//Escribo la ecuación para cada hexahedro [i]
			SumaCoeffDif_F=0;SumaCoeffDif=0;
			for (j=0 ; j< g.h3D[i].vecino.size() ;j++) {
				//recorro cada vecino==cara [i]--[j]

				Coeff_Difusion = g.h3D[i].Poligono[j].AreaPoligono /Dominio_Hsup/ g.h3D[i].Poligono[j].Dab;
				if (g.h3D[i].tipo_vecino[j] == ES_BLOQUE) {
					SumaCoeffDif_F += Coeff_Difusion * F[ g.h3D[i].vecino[j] ] ;
					SumaCoeffDif   += Coeff_Difusion;
				} else {
					iC=g.h3D[i].vecino[j];
					//if (g.Cara[iC].iBC ==1) // Nada que hacer: Derivada normal igual a cero
					if (g.Cara[iC].iBC ==2 ||g.Cara[iC].iBC ==3 ) {
						//						if (iiter==0 ) cout<<"g.Cara[iC].iBC="<<g.Cara[iC].iBC<<endl;
						SumaCoeffDif_F += Coeff_Difusion * g.Cara[iC].BC ;
						SumaCoeffDif   += Coeff_Difusion;
					}
					if (g.Cara[iC].iBC ==10 ) { // Cara de la infiltración
						SumaCoeffDif_F -= Datos_rhof*Datos_Vinyeccion*g.h3D[i].Poligono[j].AreaPoligono;
						//Aqui ponemos la evaporación
						SumaCoeffDif_F -= Psi_ev(Temp[i])*g.h3D[i].Poligono[j].AreaPoligono;
					}

				}
			}
			nF = SumaCoeffDif_F/SumaCoeffDif;
			err=fabs(nF-F[i])/(fabs(nF)+1e-20);
			if (err > errG) errG=err;
			F[i]=nF;
			//			if (F[i]==0) F[i]=0.5;

		}
		if (errG<err0) {
			cout<<"break: errG<"<<err0<<", iiter="<<iiter<<endl;
			myfileSalida<<"break: errG<"<<err0<<", iiter="<<iiter<<endl;
			break;
		}
	}

	if(DBG) cout<<"CalculaMpunto(): Velocidades"<<endl;
	for (i=0 ; i<g.nH3D ; i++) {
		SumaCoeffDif_F=0;SumaCoeffDif=0;
		double A,AA,AB,AC,RA;
		double BA,BB,BC,RB;
		double CA,CB,CC,RC;
		double nx,ny,nz,DET;
		AA=0;AB=0;AC=0;RA=0;
		BA=0;BB=0;BC=0;RB=0;
		CA=0;CB=0;CC=0;RC=0;
		for (j=0 ; j< g.h3D[i].vecino.size() ;j++) {
			if(DBG) cout<<"CalculaMpunto(): ...1115"<<endl;
			if (g.h3D[i].tipo_vecino[j] == ES_BLOQUE) {
				Coeff_Difusion = ( F[ g.h3D[i].vecino[j] ] - F[i] ) / g.h3D[i].Poligono[j].Dab;
			} else {
				iC=g.h3D[i].vecino[j];
				Coeff_Difusion=0;
				if (g.Cara[iC].iBC ==2 ||g.Cara[iC].iBC ==3) {
					Coeff_Difusion = ( g.Cara[iC].BC - F[i] ) / g.h3D[i].Poligono[j].Dab;

				}
				if (g.Cara[iC].iBC ==10) {
					Coeff_Difusion = -Datos_rhof*Datos_Vinyeccion*Dominio_Hsup;

				}
				if (g.Cara[iC].iBC ==110) {
					Coeff_Difusion = Datos_rhof*Datos_Vinyeccion*Dominio_Hsup;

				}
			}
			A=g.h3D[i].Poligono[j].AreaPoligono;
			if(DBG) cout<<"CalculaMpunto(): ...1138: A="<<A<<endl;
			nx=g.h3D[i].Poligono[j].normal.x ;
			ny=g.h3D[i].Poligono[j].normal.y ;
			nz=g.h3D[i].Poligono[j].normal.z ;
			AA += nx*nx*A; AB+= nx*ny*A; AC+= nx*nz*A; RA += -nx*Coeff_Difusion*A;
			BA += ny*nx*A; BB+= ny*ny*A; BC+= ny*nz*A; RB += -ny*Coeff_Difusion*A;
			CA += nz*nx*A; CB+= nz*ny*A; CC+= nz*nz*A; RC += -nz*Coeff_Difusion*A;
		}
		// AA AB AC --> RA
		// BA BB BC --> RB
		// CA CB CC --> RC
		DET= AA*(BB*CC-BC*CB)+BA*(CB*AC-AB*CC)+CA*(AB*BC-BB*AC) ;
		if(DBG) cout<<"CalculaMpunto(): ...1150: i,U,det ="<<i<<","<<U.size()<<","<<DET<<endl;
		U[i]= ( RA*(BB*CC-BC*CB)+RB*(CB*AC-AB*CC)+RC*(AB*BC-BB*AC) )/DET;
		V[i]=-( RA*(BA*CC-BC*CA)+RB*(CA*AC-AA*CC)+RC*(AA*BC-BA*AC) )/DET;
		W[i]= ( RA*(BA*CB-BB*CA)+RB*(CA*AB-AA*CB)+RC*(AA*BB-BA*AB) )/DET;
	}

	if(DBG) cout<<"CalculaMpunto():END"<<endl;
	return(errG);
}


double CalculaTempSuperficie(vector<double> &Temp,vector<double> &PotencialV,vector<double> &U,
		vector<double> &V,vector<double> &W,grid3D g,int niteraciones,double err0)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	int iiter,i,j,iC;
	double SumaCoeffDif_F2,SumaCoeffDif,Coeff_Difusion,err,errG,nTemp;
	if (Temp.size() !=g.nH3D) Temp.resize(g.nH3D);
	if (U.size() !=g.nH3D) U.resize(g.nH3D);
	if (V.size() !=g.nH3D) V.resize(g.nH3D);
	if (W.size() !=g.nH3D) W.resize(g.nH3D);
	//Resuelvo la ecuación varias veces
	niteraciones=g.nH3D*50;
	if (niteraciones<30000) niteraciones=30000;
	cout <<"niteraciones="<<niteraciones<<endl;
	myfileSalida <<"niteraciones="<<niteraciones<<endl;
	for (iiter=0 ; iiter<niteraciones ; iiter++) {
		errG=0;
		for (i=0 ; i<g.nH3D ; i++) {
			if(DBG) {
				if (i==20 && iiter<40) {
					cout<<"Temp[20]="<<Temp[20]<<endl;
				}
				if (i==44) {
					cout<<"Temp[44]="<<Temp[44]<<endl;
				}
			}
			//Escribo la ecuación para cada hexahedro [i]
			SumaCoeffDif_F2=0;SumaCoeffDif=0;
			double PotencialV_i=PotencialV[i],PotencialV_j;
			for (j=0 ; j< g.h3D[i].vecino.size() ;j++) {
				if(1==0) cout <<"CalculaLaplacianoCero2: iiter,i,j="<<iiter<<" "<<i<<" "<<j
						<< "\tPotencialV.size=" << PotencialV.size()
						<<"\tg.h3D[i].vecino[j]="<<g.h3D[i].vecino[j]
																   <<"\tg.h3D[i].tipo_vecino[j]="<<g.h3D[i].tipo_vecino[j]
																														<<endl;
				//recorro cada vecino==cara [i]--[j]

				if (g.h3D[i].tipo_vecino[j] == ES_BLOQUE) {
					PotencialV_j=PotencialV[ g.h3D[i].vecino[j] ];
					if (PotencialV_j>PotencialV_i) Coeff_Difusion = g.h3D[i].Poligono[j].AreaPoligono/Dominio_Hsup / g.h3D[i].Poligono[j].Dab
							*(PotencialV_j-PotencialV_i)*Datos_cf;
					else Coeff_Difusion=0;
					SumaCoeffDif_F2 += Coeff_Difusion * Temp[ g.h3D[i].vecino[j] ] ;
					SumaCoeffDif    += Coeff_Difusion;
				} else {
					iC=g.h3D[i].vecino[j];
					//1== Derivada normal igual a cero
					if (g.Cara[iC].iBC ==2 ||g.Cara[iC].iBC ==3 ) { //2,3 son condiciones DIrichlet
						PotencialV_j=g.Cara[iC].BC;
						if (PotencialV_j>PotencialV_i) Coeff_Difusion = g.h3D[i].Poligono[j].AreaPoligono/Dominio_Hsup / g.h3D[i].Poligono[j].Dab
								*(PotencialV_j-PotencialV_i)*Datos_cf;

						//						if (iiter==0 ) cout<<"g.Cara[iC].iBC="<<g.Cara[iC].iBC<<endl;
						if (1==0) cout<<"g.Cara[iC].BC2="<<g.Cara[iC].BC2
								<<"\t.BC="<<g.Cara[iC].BC
								<<"\tg.Cara[iC].iBC="<<g.Cara[iC].iBC
								<<"\tiC="<<iC
								<<"\tCoeff_Difusion="<<Coeff_Difusion
								<<"\tPotencialV_i="<<PotencialV_i<<"PotencialV_j="<<PotencialV_j
								<<endl;
						SumaCoeffDif_F2 += g.Cara[iC].BC2 * Coeff_Difusion;
						SumaCoeffDif += Coeff_Difusion;
					}
					if (g.Cara[iC].iBC ==11 ) {
						Coeff_Difusion=Datos_htilde*g.h3D[i].Poligono[j].AreaPoligono;;
						SumaCoeffDif_F2 += Coeff_Difusion*Datos_Tambiente;
						SumaCoeffDif    += Coeff_Difusion;
						//Evaporacion
						SumaCoeffDif_F2 -= Datos_CalorLatenteFluido*Psi_ev(Temp[i])*g.h3D[i].Poligono[j].AreaPoligono;
					}

				}
			}
			nTemp = SumaCoeffDif_F2/SumaCoeffDif;
			err=fabs(nTemp-Temp[i])/(fabs(nTemp)+1e-20);
			if (err>1) {
				cout <<"i="<<i
						<<"\tiiter="<<iiter
						<< "\tTemperatura previa="<<Temp[i]
														 <<"\t Nueva="<<nTemp;
				nTemp=Temp[i]+0.3*(nTemp-Temp[i]);
				cout <<"\t Nueva2="<<nTemp<<endl;

			} else {
				nTemp=Temp[i]+0.7*(nTemp-Temp[i]);
			}
			err=fabs(nTemp-Temp[i])/(fabs(nTemp)+1e-20);
			if (err > errG) errG=err;
			Temp[i]=nTemp;
			//			if (F[i]==0) F[i]=0.5;
		}
		if (errG<err0) {
			cout<<"break: errG<"<<err0<<", iiter="<<iiter<<endl;
			myfileSalida<<"break: errG<"<<err0<<", iiter="<<iiter<<endl;
			break;
		}
	}
	for (i=0 ; i<g.nH3D ; i++) {
		double A,AA,AB,AC,RA;
		double BA,BB,BC,RB;
		double CA,CB,CC,RC;
		double nx,ny,nz,DET;
		AA=0;AB=0;AC=0;RA=0;
		BA=0;BB=0;BC=0;RB=0;
		CA=0;CB=0;CC=0;RC=0;
		for (j=0 ; j< g.h3D[i].vecino.size() ;j++) {
			if (g.h3D[i].tipo_vecino[j] == ES_BLOQUE) {
				Coeff_Difusion = ( Temp[ g.h3D[i].vecino[j] ] - Temp[i] ) / g.h3D[i].Poligono[j].Dab;
			} else {
				iC=g.h3D[i].vecino[j];
				Coeff_Difusion=0;
				if (g.Cara[iC].iBC ==2 ||g.Cara[iC].iBC ==3) {
					Coeff_Difusion = ( g.Cara[iC].BC2 - Temp[i] ) / g.h3D[i].Poligono[j].Dab;

				}
			}
			A=g.h3D[i].Poligono[j].AreaPoligono;
			nx=g.h3D[i].Poligono[j].normal.x ;
			ny=g.h3D[i].Poligono[j].normal.y ;
			nz=g.h3D[i].Poligono[j].normal.z ;
			AA += nx*nx*A; AB+= nx*ny*A; AC+= nx*nz*A; RA += -nx*Coeff_Difusion*A;
			BA += ny*nx*A; BB+= ny*ny*A; BC+= ny*nz*A; RB += -ny*Coeff_Difusion*A;
			CA += nz*nx*A; CB+= nz*ny*A; CC+= nz*nz*A; RC += -nz*Coeff_Difusion*A;
		}
		// AA AB AC --> RA
		// BA BB BC --> RB
		// CA CB CC --> RC
		DET= AA*(BB*CC-BC*CB)+BA*(CB*AC-AB*CC)+CA*(AB*BC-BB*AC) ;
		U[i]= ( RA*(BB*CC-BC*CB)+RB*(CB*AC-AB*CC)+RC*(AB*BC-BB*AC) )/DET;
		V[i]=-( RA*(BA*CC-BC*CA)+RB*(CA*AC-AA*CC)+RC*(AA*BC-BA*AC) )/DET;
		W[i]= ( RA*(BA*CB-BB*CA)+RB*(CA*AB-AA*CB)+RC*(AA*BB-BA*AB) )/DET;
	}
	return(errG);
}

double MediaArmonica(double h1,double h2)
{
	return ((h1*h2) / (h1+h2));
}

double CalculaTemperaturaPilaEnTmasDt(vector<double> &Temp,vector<double> &PotencialV,vector<double> &U,
		vector<double> &V,vector<double> &W,grid3D g,int niteraciones,double err0,vector<double> &TempPrevia)
{
	int iiter,i,j,iC;
	double SumaCoeffDif_F2,SumaCoeffDif,Coeff_Difusion,err,errG,nF,SumaVolumen;
	if (Temp.size() ==0) Temp.resize(g.nH3D);
	//Resuelvo la ecuación varias veces
	niteraciones=g.nH3D*50;
	if (niteraciones<30000) niteraciones=30000;
	cout <<"niteraciones="<<niteraciones<<endl;
	myfileSalida <<"niteraciones="<<niteraciones<<endl;
	double SumCoeffEstacionarios=0,SumCoeffTemporales=0;
	int    NCoeffEstacionarios=0;
	for (iiter=0 ; iiter<niteraciones ; iiter++) {
		errG=0;
		for (i=0 ; i<g.nH3D ; i++) {
			//Escribo la ecuación para cada hexahedro [i]
			SumaCoeffDif_F2=0;SumaCoeffDif=0;SumaVolumen=0;
			double PotencialV_i=PotencialV[i],PotencialV_j;
			for (j=0 ; j< g.h3D[i].vecino.size() ;j++) {
				if(1==0) cout <<"CalculaLaplacianoCero2: iiter,i,j="<<iiter<<" "<<i<<" "<<j
						<< "\tPotencialV_.size=" << PotencialV.size()
						<<"\tg.h3D[i].vecino[j]="<<g.h3D[i].vecino[j]
																   <<"\tg.h3D[i].tipo_vecino[j]="<<g.h3D[i].tipo_vecino[j]
																														<<endl;
				//recorro cada vecino==cara [i]--[j]

				if (g.h3D[i].tipo_vecino[j] == ES_BLOQUE) {
					PotencialV_j=PotencialV[ g.h3D[i].vecino[j] ];
					if (PotencialV_j>PotencialV_i)
						Coeff_Difusion = Datos_rhof * Datos_cf *(PotencialV_j-PotencialV_i) *
						g.h3D[i].Poligono[j].AreaPoligono / g.h3D[i].Poligono[j].Dab;
					else Coeff_Difusion=0;

					Coeff_Difusion += Datos_km_mu * g.h3D[i].Poligono[j].AreaPoligono / g.h3D[i].Poligono[j].Dab;

					SumaCoeffDif_F2 += Coeff_Difusion * Temp[ g.h3D[i].vecino[j] ] ;
					SumaCoeffDif    += Coeff_Difusion;
				} else {
					iC=g.h3D[i].vecino[j];
					double htilde_local;
					//if (g.Cara[iC].iBC==1) Se trata de Derivada normal igual a cero (nada que hacer)
					switch (g.Cara[iC].iBC) {
					case 2:
					case 3:   //Condicion Dirichlet (como un bloque con T conocida == BC2)
						PotencialV_j=g.Cara[iC].BC;
						if (PotencialV_j>PotencialV_i)
							Coeff_Difusion = Datos_rhof * Datos_cf *(PotencialV_j-PotencialV_i) *
							g.h3D[i].Poligono[j].AreaPoligono / g.h3D[i].Poligono[j].Dab;
						else Coeff_Difusion=0;

						Coeff_Difusion += Datos_km_mu * g.h3D[i].Poligono[j].AreaPoligono / g.h3D[i].Poligono[j].Dab;

						if (1==0) cout<<"g.Cara[iC].BC2="<<g.Cara[iC].BC2
								<<"\t.BC="<<g.Cara[iC].BC
								<<"\tg.Cara[iC].iBC="<<g.Cara[iC].iBC
								<<"\tiC="<<iC
								<<"\tCoeff_Difusion="<<Coeff_Difusion
								<<"\tPotencialV_i="<<PotencialV_i<<"PotencialV_j="<<PotencialV_j
								<<endl;
						SumaCoeffDif_F2 += Coeff_Difusion * g.Cara[iC].BC2;
						SumaCoeffDif    += Coeff_Difusion;
						break;

					case 11:   // Cara convectiva abajo
						htilde_local=MediaArmonica(Datos_km_mu/(g.h3D[i].Poligono[j].Dab/2),Datos_hc_inferior);

						Coeff_Difusion=htilde_local*g.h3D[i].Poligono[j].AreaPoligono;
						SumaCoeffDif_F2 += Coeff_Difusion * Datos_Tambiente;
						SumaCoeffDif    += Coeff_Difusion;
						break;
					case 12:   // Cara convectiva lateral
						htilde_local
						=MediaArmonica(Datos_km_mu/(g.h3D[i].Poligono[j].Dab/2),Datos_km2/Datos_DistanciaAlBorde);
						htilde_local=MediaArmonica(htilde_local,Datos_hc_lateral);
						Coeff_Difusion=htilde_local*g.h3D[i].Poligono[j].AreaPoligono;
						SumaCoeffDif_F2 += Coeff_Difusion * Datos_Tambiente;
						SumaCoeffDif    += Coeff_Difusion;
						break;
					}

				}
			}

			if (TipoCalculo==CalculoEvolucion) {
				///////////////////Aporte del termino temporal.....
				double CoeffT=Datos_rhom*Datos_cm*g.h3D[i].volumen/Datos_dt;

				//Estadistica (cual es mas importante)
				SumCoeffEstacionarios += SumaCoeffDif;
				SumCoeffTemporales    += CoeffT;
				NCoeffEstacionarios   ++;


				SumaCoeffDif_F2 += CoeffT*TempPrevia[i];
				SumaCoeffDif    += CoeffT;
			}
			nF = SumaCoeffDif_F2/SumaCoeffDif;
			err=fabs(nF-Temp[i])/(fabs(nF)+1e-20);
			if (err > errG) errG=err;
			Temp[i]=nF;
			//			if (F[i]==0) F[i]=0.5;
		}
		if (errG<err0) {
			cout<<"break: errG<"<<err0<<endl;
			myfileSalida<<"break: errG<"<<err0<<endl;
			break;
		}
	}
	cout<<"iiter_final="<<iiter<<endl;
	myfileSalida<<"iiter_final="<<iiter<<endl;
	if (TipoCalculo==CalculoEvolucion) {
		//Estadistica (cual es mas importante)
		cout << "<CoeffEstacionarios>="<<SumCoeffEstacionarios /NCoeffEstacionarios;
		cout << "\t<SumCoeffTemporales>="<<SumCoeffTemporales / NCoeffEstacionarios <<endl;
		myfileSalida << "<CoeffEstacionarios>="<<SumCoeffEstacionarios /NCoeffEstacionarios;
		myfileSalida << "\t<SumCoeffTemporales>="<<SumCoeffTemporales / NCoeffEstacionarios <<endl;
	}
	for (i=0 ; i<g.nH3D ; i++) {
		double A,AA,AB,AC,RA;
		double BA,BB,BC,RB;
		double CA,CB,CC,RC;
		double nx,ny,nz,DET;
		AA=0;AB=0;AC=0;RA=0;
		BA=0;BB=0;BC=0;RB=0;
		CA=0;CB=0;CC=0;RC=0;
		for (j=0 ; j< g.h3D[i].vecino.size() ;j++) {
			if (g.h3D[i].tipo_vecino[j] == ES_BLOQUE) {
				Coeff_Difusion = ( Temp[ g.h3D[i].vecino[j] ] - Temp[i] ) / g.h3D[i].Poligono[j].Dab;
			} else {
				iC=g.h3D[i].vecino[j];
				Coeff_Difusion=0;
				if (g.Cara[iC].iBC ==2 ||g.Cara[iC].iBC ==3) {
					Coeff_Difusion = ( g.Cara[iC].BC2 - Temp[i] ) / g.h3D[i].Poligono[j].Dab;

				}
			}
			A=g.h3D[i].Poligono[j].AreaPoligono;
			nx=g.h3D[i].Poligono[j].normal.x ;
			ny=g.h3D[i].Poligono[j].normal.y ;
			nz=g.h3D[i].Poligono[j].normal.z ;
			AA += nx*nx*A; AB+= nx*ny*A; AC+= nx*nz*A; RA += -nx*Coeff_Difusion*A;
			BA += ny*nx*A; BB+= ny*ny*A; BC+= ny*nz*A; RB += -ny*Coeff_Difusion*A;
			CA += nz*nx*A; CB+= nz*ny*A; CC+= nz*nz*A; RC += -nz*Coeff_Difusion*A;
		}
		// AA AB AC --> RA
		// BA BB BC --> RB
		// CA CB CC --> RC
		DET= AA*(BB*CC-BC*CB)+BA*(CB*AC-AB*CC)+CA*(AB*BC-BB*AC) ;
		U[i]= ( RA*(BB*CC-BC*CB)+RB*(CB*AC-AB*CC)+RC*(AB*BC-BB*AC) )/DET;
		V[i]=-( RA*(BA*CC-BC*CA)+RB*(CA*AC-AA*CC)+RC*(AA*BC-BA*AC) )/DET;
		W[i]= ( RA*(BA*CB-BB*CA)+RB*(CA*AB-AA*CB)+RC*(AA*BB-BA*AB) )/DET;
	}
	return(errG);
}

double CalculaConveccionDifusionDt(vector<double> &Concentracion,grid3D* g,int niteraciones,double err0,
		vector<double> &ConcentracionPrevia)
{

	if(DBG) cout<<"CalculaConveccionDifusionDt(...)"<<endl;
	if(DBG) cout <<"niteraciones="<<niteraciones<<endl;


	FlagConveccion=1;TipoCalculo=CalculoEvolucion;

	if (DBG>0) {
		cout<<"CalculaConveccionDifusionDt(";int ip=0;
		if (TipoCalculo==CalculoEvolucion) {cout <<"Dt="<<Datos_dt;ip=1;}

		if (FlagDifusion) {if (ip==1) cout <<" , ";cout <<"Difusión: mu="<<Datos_km_mu;ip=1;}
		if (FlagConveccion) {if (ip==1) cout <<" , ";cout <<"Convección";ip=1;}
		cout <<");"<<endl;
	}
	int iiter,i,j,iC;
	double err,errG,SumaVolumen,CFL_max,CFL_media;
	int nCFL_media;
	double SumCoeffEstacionarios=0,SumCoeffTemporales=0,SumaDifusion=0,SumaConveccion=0;
	int    NCoeffEstacionarios=0;

	if (Concentracion.size() ==0) Concentracion.resize(g->nVolFinito);


	//Resuelvo la ecuación varias veces
	for (iiter=0 ; iiter<niteraciones ; iiter++) {

		if(DBG) cout<<"iiter="<<iiter<<endl;
		errG=0;
		CFL_max=0;
		CFL_media=0;
		nCFL_media=0;
		int ierrG=0;

		//pthread_mutex_lock(&mutex_DrawRead);
#pragma omp parallel for num_threads(GL_threads)
		for (i=0 ; i<g->nVolFinito ; i++) {
			//Escribo la ecuación para el VolFinito [i]

			double Suma_Numerador=0,Suma_Denominador=0,SumaVolumen=0,ConcentracionNueva;
			double Suma_Dt=0,Suma_Conveccion=0,Suma_Difusion=0;
			if(DBG>1) cout<<"2456: i="<<i<<"    g->nVolFinito"<<g->nVolFinito<<"=="<<g->VolFinito.size()<<endl;

			double U1i=U[step_hora][i]*(1-stepf)+U[step_hora+1][i]*stepf;
			double U2i=V[step_hora][i]*(1-stepf)+V[step_hora+1][i]*stepf;
			double U3i=W[step_hora][i]*(1-stepf)+W[step_hora+1][i]*stepf;


			for (int j=0 ; j< g->VolFinito[i].vecino.size() ;j++) {
				//recorro cada cara que conecta al VolFinito [i] con su vecino [j]
				double Coeff_Cj=0,Cij;

				if (0 && DBG ) {
					cout <<"g->VolFinito[i].tipo_vecino[j]="<<g->VolFinito[i].tipo_vecino[j]<<endl;
				}
				if (g->VolFinito[i].tipo_vecino[j] == ES_BLOQUE) {
					int vj=g->VolFinito[i].vecino[j];

					if(FlagDifusion) {
						Cij=Datos_km_mu * g->VolFinito[i].Poligono[j].AreaPoligono / g->VolFinito[i].Poligono[j].Dab;
						Coeff_Cj += Cij;
						Suma_Difusion +=Cij;
					}
					if(FlagConveccion) {

						//cout <<"2852: U1.size()="<<U1.size()<<endl;
						//cout <<"2852: g->VolFinito.size()="<<g->VolFinito.size()<<endl;
						//cout <<"2852: g->VolFinito[i].vecino.size()="<<g->VolFinito[i].vecino.size()<<endl;

						double U1j=U[step_hora][vj]*(1-stepf)+U[step_hora+1][vj]*stepf;
						double U2j=V[step_hora][vj]*(1-stepf)+V[step_hora+1][vj]*stepf;
						double U3j=W[step_hora][vj]*(1-stepf)+W[step_hora+1][vj]*stepf;




						R3 VelocidadFrontera;
						VelocidadFrontera.x=(U1i+U1j)/2;
						VelocidadFrontera.y=(U2i+U2j)/2;
						VelocidadFrontera.z=(U3i+U3j)/2;
						//cout <<"2861"<<endl;


						double VelocidadNormal=ppunto(VelocidadFrontera,g->VolFinito[i].vecino_normal[j]);

						//cout <<"2861"<<endl;
						if(VelocidadNormal<0) {
							Cij = -VelocidadNormal*g->VolFinito[i].Poligono[j].AreaPoligono ;
							Coeff_Cj += Cij;
							Suma_Conveccion +=Cij;
						} else {
							//Coeff_Cj=0;
						}
					}

					Suma_Numerador 		+= Coeff_Cj * Concentracion[ vj ] ;
					Suma_Denominador    += Coeff_Cj;
				} //end: if (g->VolFinito[i].tipo_vecino[j] == ES_BLOQUE)
				else {
					//Esta superficie Sj es parte de la frontera.
					iC=g->VolFinito[i].vecino[j];
					double htilde_local;
					//if (g->Cara[iC].iBC==1) Se trata de Derivada normal igual a cero (nada que hacer)
					switch (g->Cara[iC].iBC) {
					case 112:
					case 113:   //Condicion Dirichlet (como un bloque con Concentracion conocida == BC2)

						if (FlagConveccion) {
							R3 VelocidadFrontera;
							VelocidadFrontera.x=U1i;
							VelocidadFrontera.y=U2i;
							VelocidadFrontera.z=U3i;

							double VelocidadNormal=ppunto(VelocidadFrontera,g->VolFinito[i].vecino_normal[j]);

							if(VelocidadNormal<0) {
								Cij      = -VelocidadNormal*g->VolFinito[i].Poligono[j].AreaPoligono ;
								Coeff_Cj += Cij;
								Suma_Conveccion +=Cij;
							}

						}
						if (FlagDifusion) {
							Cij       = Datos_km_mu* g->VolFinito[i].Poligono[j].AreaPoligono / g->VolFinito[i].Poligono[j].Dab;
							Coeff_Cj += Cij;
							Suma_Difusion +=Cij;
						}

						Suma_Numerador += Coeff_Cj * g->Cara[iC].BC2;
						Suma_Denominador    += Coeff_Cj;
						break;
					}//end switch (g->Cara[iC].iBC)

				} //end: if (g->VolFinito[i].tipo_vecino[j] == ES_CARA)

			} //end: for (int j=0 ; j< g->VolFinito[i].vecino.size() ;j++)

			if (TipoCalculo==CalculoEvolucion) {
				///////////////////Aporte del termino temporal.....
				double CoeffT= g->VolFinito[i].volumen/Datos_dt;
				Suma_Dt +=CoeffT;

				//Estadistica (cual es mas importante)
				SumCoeffEstacionarios += Suma_Denominador;
				SumCoeffTemporales    += CoeffT;
				NCoeffEstacionarios   ++;


				Suma_Numerador   += CoeffT*ConcentracionPrevia[i];
				Suma_Denominador += CoeffT;
			}


			ConcentracionNueva = Suma_Numerador/Suma_Denominador;

			double CFLi=(Suma_Conveccion+Suma_Difusion)/Suma_Dt;

			if (DBG>0 && i==7055 && (iiter%100==0)) { //Suma_Denominador<1e-5
				cout<<"Para i="<<i<<": C="<<Suma_Numerador<<
						"/"<<Suma_Denominador<<"="<<ConcentracionNueva;
				cout<<" CFLi="<<CFLi;
				cout<<" Suma_Dt="<<Suma_Dt<<" Suma_Conveccion="<<Suma_Conveccion<<" Suma_Difusion="<<Suma_Difusion;
			}
			//			err=fabs(ConcentracionNueva-Concentracion[i])/(fabs(ConcentracionNueva)+1e-20);
			err=fabs(ConcentracionNueva-Concentracion[i]);

			CFL_media+=CFLi;
			nCFL_media++;
			if (CFLi>CFL_max) {
				CFL_max=CFLi;
			}
			if (err > errG) {
				if(DBG) cout<<"err="<<err<<">errG="<<errG<<"   ConcentracionNueva="<<ConcentracionNueva<<"   Concentracion["<<i<<"]="<<Concentracion[i]<<endl;
				if(DBG) cout<<"2545: Suma_Numerador,Suma_Denominador: "<<Suma_Numerador<<" , "<<Suma_Denominador
						<< "   CP="<<ConcentracionPrevia[i]
														 << "   Ci="<<Concentracion[i]<<endl;
				errG=err;
				ierrG=i;

			}
			Concentracion[i]=ConcentracionNueva;
		}//end for (i=0 ; i<g->nVolFinito ; i++)
		//pthread_mutex_unlock(&mutex_DrawRead);


		CFL_media /= nCFL_media;

		Global_error=errG;
		if ((iiter%100==0)) { //Suma_Denominador<1e-5
			cout<<" errG="<<Global_error<<" ierrG="<<ierrG
					<<" CFLG="<<CFL_max<<" iiter="<<iiter<<endl;
			int FlagTrataBC=0;
			double minF,maxF;
			ProyectaVolumenesAVertices(F2Volumenes,F2Nodos,&minF,&maxF,FlagTrataBC) ;

		}

		if (errG<err0 && iiter>3) {
			if (DBG) cout<<"break: errG<"<<err0<<endl;
			break;
		}
		if (DBG) cout<<"errG[i]="<<errG<<endl;
	}

	if(DBG){
		cout<<"iiter_final="<<iiter;
		cout<<"  CFL_max="<<CFL_max;
		cout<<"  CFL_media="<<CFL_media;
		cout<<"  errG="<<errG<<endl;
	}
	if (TipoCalculo==CalculoEvolucion) {
		//Estadistica (cual es mas importante)
		if(DBG) {
			cout << "<CoeffDifusion>="<<SumaDifusion /NCoeffEstacionarios;
			cout << "<SumaConveccion>="<<SumaConveccion /NCoeffEstacionarios;
			cout << "<CoeffEstacionarios>="<<SumCoeffEstacionarios /NCoeffEstacionarios;
			cout << "\t<SumCoeffTemporales>="<<SumCoeffTemporales / NCoeffEstacionarios <<endl;
		}

		Peclet=SumaConveccion/SumaDifusion;

		if (SumCoeffTemporales<2*SumCoeffEstacionarios) {
			FactorClockTime *= 0.5;

			if (glui != NULL) glui->sync_live();
		}
		if (SumCoeffTemporales>10*SumCoeffEstacionarios) {
			FactorClockTime *= 2;

			if (glui != NULL) glui->sync_live();
		}
	}
#if 0
	//Calculo del Gradiente
	for (i=0 ; i<g->nVolFinito ; i++) {
		double A,AA,AB,AC,RA;
		double BA,BB,BC,RB;
		double CA,CB,CC,RC;
		double nx,ny,nz,DET;
		AA=0;AB=0;AC=0;RA=0;
		BA=0;BB=0;BC=0;RB=0;
		CA=0;CB=0;CC=0;RC=0;
		for (j=0 ; j< g->VolFinito[i].vecino.size() ;j++) {
			if (g->VolFinito[i].tipo_vecino[j] == ES_BLOQUE) {
				Coeff_Cj = ( Concentracion[ g->VolFinito[i].vecino[j] ] - Concentracion[i] ) / g->VolFinito[i].Poligono[j].Dab;
			} else {
				iC=g->VolFinito[i].vecino[j];
				Coeff_Cj=0;
				if (g->Cara[iC].iBC ==2 ||g->Cara[iC].iBC ==3) {
					Coeff_Cj = ( g->Cara[iC].BC2 - Concentracion[i] ) / g->VolFinito[i].Poligono[j].Dab;

				}
			}
			A=g->VolFinito[i].Poligono[j].AreaPoligono;
			nx=g->VolFinito[i].Poligono[j].normal.x ;
			ny=g->VolFinito[i].Poligono[j].normal.y ;
			nz=g->VolFinito[i].Poligono[j].normal.z ;
			AA += nx*nx*A; AB+= nx*ny*A; AC+= nx*nz*A; RA += -nx*Coeff_Cj*A;
			BA += ny*nx*A; BB+= ny*ny*A; BC+= ny*nz*A; RB += -ny*Coeff_Cj*A;
			CA += nz*nx*A; CB+= nz*ny*A; CC+= nz*nz*A; RC += -nz*Coeff_Cj*A;
		}
		// AA AB AC --> RA
		// BA BB BC --> RB
		// CA CB CC --> RC
		DET= AA*(BB*CC-BC*CB)+BA*(CB*AC-AB*CC)+CA*(AB*BC-BB*AC) ;
		U1[i]= ( RA*(BB*CC-BC*CB)+RB*(CB*AC-AB*CC)+RC*(AB*BC-BB*AC) )/DET;
		V[i]=-( RA*(BA*CC-BC*CA)+RB*(CA*AC-AA*CC)+RC*(AA*BC-BA*AC) )/DET;
		U3[i]= ( RA*(BA*CB-BB*CA)+RB*(CA*AB-AA*CB)+RC*(AA*BB-BA*AB) )/DET;
	}
#endif
	if(DBG) {
		double minF=1e20,maxF=-1e20;
		for (i=0;i<g->nVolFinito;i++) {
			if (minF>Concentracion[i]) minF=Concentracion[i];
			if (maxF<Concentracion[i]) maxF=Concentracion[i];

		}
		if(DBG) cout<<"2703: minF,maxF="<<minF<<" , "<<maxF<<endl;

		if(DBG) cout<<"CalculaConveccionDifusionDt(...):END"<<endl;
	}
	return(errG);
}

void CalculaCoeficientes_DominioVariableZ(grid3D* g)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	if (Coeff_H.size()==0) {
		Coeff_H.resize(g->nVolFinito);
		Coeff_Q.resize(g->nVolFinito);
		Coeff_Kn.resize(g->nVolFinito);
		Coeff_Vtilde.resize(g->nVolFinito);
		Coeff_1omega.resize(g->nVolFinito,0);
		Coeff_2omega.resize(g->nVolFinito,0);
	}
	for (int i=0 ; i<g->nVolFinito ; i++) {
		if (Coeff_H[i].size()==0) {
			Coeff_H[i].resize(g->VolFinito[i].vecino.size(),0);
			Coeff_Q[i].resize(g->VolFinito[i].vecino.size(),0);
			Coeff_Kn[i].resize(g->VolFinito[i].vecino.size(),0);
			Coeff_Vtilde[i].resize(g->VolFinito[i].vecino.size(),0);
		}
		Coeff_1omega[i]=0;
		for (int j=0 ; j< g->VolFinito[i].vecino.size() ;j++) {
			R3 Normalj=g->VolFinito[i].vecino_normal[j];
			double Areaj=g->VolFinito[i].Poligono[j].AreaPoligono ;
			double Achicaj = Achica[step_hora][i];

			Coeff_Kn[i][j] =0;
			if (g->VolFinito[i].tipo_vecino[j] == ES_BLOQUE) { //Caso: Sigma_j\subseteq \int\Omega
				int vf_j=g->VolFinito[i].vecino[j];
				double Achicaj=(Achica[step_hora][i] + Achica[step_hora][vf_j]  )/2;

				Coeff_H[i][j] = ( Grad_F[step_hora][i] + Grad_F[step_hora][vf_j] )/2;
				Coeff_Q[i][j] = ( Grad_F[step_hora][i].NormaSQR() /Achica[step_hora][i] + Grad_F[step_hora][vf_j].NormaSQR() /Achica[step_hora][vf_j]  )/2;
				Coeff_Kn[i][j] =ppunto(Coeff_H[i][j] , Normalj)-Coeff_Q[i][j]*Normalj.z;

				R3 VelocidadFrontera;

				VelocidadFrontera.x=(U[step_hora][i]+U[step_hora][vf_j])/2;
				VelocidadFrontera.y=(V[step_hora][i]+V[step_hora][vf_j])/2;
				VelocidadFrontera.z=(W[step_hora][i]+W[step_hora][vf_j])/2;

				double Parentesis=ppunto(VelocidadFrontera, Coeff_H[i][j])+(FuncionDFDt_Volumenes[step_hora][i]+FuncionDFDt_Volumenes[step_hora][vf_j])/2;
				Coeff_Vtilde[i][j] =Achicaj*ppunto(VelocidadFrontera, Normalj)-Parentesis*Normalj.z;

			} else {
				//Esta superficie Sigma_j es parte de la frontera.
				int iC=g->VolFinito[i].vecino[j];
				switch (g->Cara[iC].iBC) {
				case 112:
				case 113:   //Condicion Dirichlet (como un bloque con Concentracion conocida == BC2)

					Coeff_H[i][j] = Grad_F[step_hora][i] * 2;

					Coeff_Q[i][j] = Grad_F[step_hora][i].NormaSQR() /Achica[step_hora][i] * 2;

					Coeff_Kn[i][j] =ppunto(Coeff_H[i][j] , Normalj)-Coeff_Q[i][j]*Normalj.z;

					R3 VelocidadFrontera;
					R3 Normalj=g->VolFinito[i].vecino_normal[j];

					VelocidadFrontera.x = U[step_hora][i];
					VelocidadFrontera.y = V[step_hora][i];
					VelocidadFrontera.z = W[step_hora][i];

					double Parentesis=ppunto(VelocidadFrontera, Coeff_H[i][j])+FuncionDFDt_Volumenes[step_hora][i];
					Coeff_Vtilde[i][j] =Achicaj*ppunto(VelocidadFrontera, Normalj)-Parentesis*Normalj.z;

					break;
				}
			}
			Coeff_1omega[i]+= Coeff_Kn[i][j]* Areaj;
			if (DBGVARZ2) {
				cout <<"2823: i="<<i<<endl;
				cout<<"Coeff_2omega[i]="<<Coeff_2omega[i]<<endl;
				cout<<"Normalj="<<Normalj<<endl;
				nprintDBGVARZ++;
				if (nprintDBGVARZ>10) exit(0);
			}
			Coeff_2omega[i]+= Coeff_H[i][j]* Normalj.z* Areaj ;
		}
	}

}
double CalculaConveccionDifusionDt_DominioVariableZ(vector<double> &Concentracion, vector<double> &U1,
		vector<double> &U2,vector<double> &U3,grid3D* g,int niteraciones,double err0,vector<double> &ConcentracionPrevia)
{
	int iiter,i,j,iC;
	double Suma_Numerador,Suma_Denominador,Coeff_Cj,err,errG,ConcentracionNueva,SumaVolumen;

	if(DBG) cout<<"CalculaConveccionDifusionDt_DominioVariableZ(...)"<<endl;
	if (Concentracion.size() ==0) Concentracion.resize(g->nVolFinito);
	//Resuelvo la ecuación varias veces

	if(DBG) cout <<"niteraciones="<<niteraciones<<endl;
	double SumCoeffEstacionarios=0,SumCoeffTemporales=0,SumaDifusion=0,SumaConveccion=0,SumaTemporal=0;
	int    NCoeffEstacionarios=0;

	Datos_km_mu=0.001;

	CalculaCoeficientes_DominioVariableZ(g);

	for (iiter=0 ; iiter<niteraciones ; iiter++) {

		if(DBG) cout<<"iiter="<<iiter<<endl;
		errG=0;



		//cout <<"2) AchicaPrevia[3]="<<AchicaPrevia[3]<<endl;
#pragma omp parallel for num_threads(GL_threads)
		for (i=0 ; i<g->nVolFinito ; i++) {
			//Escribo la ecuación para el VolFinito [i]
			double volumen_i=g->VolFinito[i].volumen;
			double Coeff_Cj=0,Coeff_1D,Coeff_2D,Coeff_3D,Coeff_4D,Coeff_5D,Coeff_6D,Coeff_7D,Suma_67j;
			double LadoDerecho, Coeff_Uw;

			LadoDerecho=0; Coeff_Uw=0;
			char sCoeff_Uw[100000];
			if (DBGVARZ) sprintf(sCoeff_Uw,"[2869]Coeff_Uw=%g\n",Coeff_Uw);

			double Suma_Numerador=0,Suma_Denominador=0,SumaVolumen=0;
			if(DBG>1) cout<<"2456: i="<<i<<"    g->nVolFinito"<<g->nVolFinito<<"=="<<g->VolFinito.size()<<endl;

			//cout <<"3) AchicaPrevia[3]="<<AchicaPrevia[3]<<",   i="<<i<<endl;

			for (int j=0 ; j< g->VolFinito[i].vecino.size() ;j++) {	//recorro cada cara que conecta al VolFinito [i] con su vecino [j]
				R3 Normalj=g->VolFinito[i].vecino_normal[j];
				double Areaj =g->VolFinito[i].Poligono[j].AreaPoligono ;
				double Dabj  =g->VolFinito[i].Poligono[j].Dab ;

				if (DBG & 0) {
					cout <<"2405"<<endl;
					cout <<"g->VolFinito[i].tipo_vecino[j]="<<g->VolFinito[i].tipo_vecino[j]<<endl;
				}

#if DIFUSION==0
				if(FlagDifusion) { // Laplaciano....
					Suma_67j=0;
					if (g->VolFinito[i].tipo_vecino[j] == ES_BLOQUE) { //\Sigma_j \subseteq \int\Omega
						int vf_j=g->VolFinito[i].vecino[j];
						double volumen_j=g->VolFinito[vf_j].volumen;
						double Achicaj=(Achica[step_hora][i] + Achica[step_hora][vf_j]  )/2;

						//Parte difusion
						Coeff_1D = Datos_km_mu*Achicaj*Areaj/Dabj;
						Coeff_3D = -Datos_km_mu/4/volumen_i*Normalj.z*Areaj*Coeff_1omega[i];
						Coeff_4D = -Datos_km_mu/4/volumen_j*Normalj.z*Areaj*Areaj*(Coeff_Kn[i][j]+ppunto(Coeff_H[i][j],Normalj));
						Coeff_5D = -Datos_km_mu/4/volumen_i*ppunto(Coeff_2omega[i],Normalj)*Areaj;


						if (DBGVARZ2) {
							cout <<"i="<<i<<endl;
							cout<<"Coeff_1D="<<Coeff_1D<<endl;
							cout<<"Coeff_3D="<<Coeff_3D<<endl;
							cout<<"Coeff_4D="<<Coeff_4D<<endl;
							cout<<"Coeff_5D="<<Coeff_5D<<endl;
							cout<<"Coeff_2omega[i]="<<Coeff_2omega[i]<<endl;
							cout<<"Normalj="<<Normalj<<endl;
							nprintDBGVARZ++;
							if (nprintDBGVARZ>10) exit(0);
						}
						for (int k=0 ; k< g->VolFinito[vf_j].vecino.size() ;k++) {
							R3 Normalk=g->VolFinito[vf_j].vecino_normal[k];
							double Areak =g->VolFinito[vf_j].Poligono[k].AreaPoligono ;
							Coeff_6D= (Datos_km_mu/4)*(Areaj/volumen_j) *Coeff_Kn[i][j]*Normalk.z*Areak;
							Coeff_7D= (Datos_km_mu/4)*(Areaj/volumen_j) *Normalj.z *ppunto(Coeff_H[i][j],Normalk)*Areak;
							if (g->VolFinito[vf_j].tipo_vecino[k] == ES_BLOQUE) { //\Sigma_k \subseteq \int\Omega
								int vf_k=g->VolFinito[vf_j].vecino[k];
								if (vf_k == i) continue;
								Suma_67j+=(Coeff_6D+Coeff_7D)*(Concentracion[ vf_j ]- Concentracion[ vf_k ]);
							} else {
								int iCk=g->VolFinito[vf_j].vecino[k];
								switch (g->Cara[iCk].iBC) {
								case 112:
								case 113:   //Condicion Dirichlet (como un bloque con Concentracion conocida == BC2)
									Suma_67j+=(Coeff_6D+Coeff_7D)*(Concentracion[ vf_j ]- g->Cara[iCk].BC2);
									break;
								}
							}
						}
						Coeff_Uw += Coeff_1D+Coeff_3D+Coeff_4D+Coeff_5D;

						if (DBGVARZ) {
							sprintf(sCoeff_Uw,"%s[2933]Coeff_Uw=%g\n",sCoeff_Uw,Coeff_Uw);
							sprintf(sCoeff_Uw,"%s[2933]Coeff_1D=%g\n",sCoeff_Uw,Coeff_1D);
							sprintf(sCoeff_Uw,"%s[2933]Coeff_3D=%g\n",sCoeff_Uw,Coeff_3D);
							sprintf(sCoeff_Uw,"%s[2933]Coeff_4D=%g\n",sCoeff_Uw,Coeff_4D);
							sprintf(sCoeff_Uw,"%s[2933]Coeff_5D=%g\n",sCoeff_Uw,Coeff_5D);
						}

						LadoDerecho += (Coeff_1D+Coeff_3D+Coeff_4D+Coeff_5D)* Concentracion[ vf_j ] + Suma_67j;
						SumaDifusion += Coeff_1D+Coeff_3D+Coeff_4D+Coeff_5D;
					} else { //Sigma_j es una cara
						iC=g->VolFinito[i].vecino[j];
						double htilde_local;
						switch (g->Cara[iC].iBC) {
						case 112:
						case 113:   //Condicion Dirichlet (como un bloque con Concentracion conocida == BC2)
							double Achicaj = Achica[step_hora][i];
							if (FlagDifusion) {//OK Dominio Variable
								Coeff_2D = Datos_km_mu*Achicaj*Areaj/Dabj;
								Coeff_3D = -Datos_km_mu/4/volumen_i*Normalj.z*Areaj*Coeff_1omega[i];
								Coeff_5D = -Datos_km_mu/4/volumen_i*ppunto(Coeff_2omega[i],Normalj)*Areaj;
								Coeff_Cj 	 += Coeff_2D+Coeff_3D+Coeff_5D;
							}
							Coeff_Uw += Coeff_2D+Coeff_3D+Coeff_5D;

							if (DBGVARZ) sprintf(sCoeff_Uw,"%s[2952]Coeff_Uw=%g\n",sCoeff_Uw,Coeff_Uw);

							LadoDerecho += (Coeff_2D+Coeff_3D+Coeff_5D)* g->Cara[iC].BC2 ;
							SumaDifusion+= Coeff_2D+Coeff_3D+Coeff_5D;
							break;
						}
					} //End Sigma_j es una cara
					if(DBG>1) cout<<"2510"<<endl;

				} //Fin Difusion (Laplaciano)
#endif
#if CONVECCION==0
				if(FlagConveccion) {//FALTA Dominio Variable

					double Coeff_Vj = Coeff_Vtilde[i][j]*Areaj;

					if (g->VolFinito[i].tipo_vecino[j] == ES_BLOQUE) { //\Sigma_j \subseteq \int\Omega
						int vf_j=g->VolFinito[i].vecino[j];
						if (Coeff_Vj<0 )
							LadoDerecho += (-Coeff_Vj)* Concentracion[ vf_j ];
						else {
							Coeff_Uw += Coeff_Vj;

							if (DBGVARZ) sprintf(sCoeff_Uw,"%s[2975]Coeff_Uw=%g\n",sCoeff_Uw,Coeff_Uw);


							SumaConveccion += Coeff_Vj;
						}
					} else {
						//Esta superficie Sj es parte de la frontera.
						iC=g->VolFinito[i].vecino[j];
						double htilde_local;
						//if (g->Cara[iC].iBC==1) Se trata de Derivada normal igual a cero (nada que hacer)
						switch (g->Cara[iC].iBC) {
						case 112:
						case 113:   //Condicion Dirichlet (como un bloque con Concentracion conocida == BC2)
							if (Coeff_Vj<0 )
								LadoDerecho += (-Coeff_Vj)* g->Cara[iC].BC2;
							else {
								Coeff_Uw += Coeff_Vj;

								if (DBGVARZ) sprintf(sCoeff_Uw,"%s[2993]Coeff_Uw=%g\n",sCoeff_Uw,Coeff_Uw);


								SumaConveccion += Coeff_Vj;
							}
							break;
						}
					}
				} //Fin Conveccion
#endif
			} //fin for(j....

			if (TipoCalculo==CalculoEvolucion) {
				double CoeffT  = g->VolFinito[i].volumen/Datos_dt;
				LadoDerecho   += CoeffT*AchicaPrevia[i]*ConcentracionPrevia[i];
				Coeff_Uw      += CoeffT*Achica_t[i];

				if (DBGVARZ) {
					sprintf(sCoeff_Uw,"%s[3010]Coeff_Uw=%g\n",sCoeff_Uw,Coeff_Uw);
					sprintf(sCoeff_Uw,"%s[3010]Achica_t[i]=%g\n",sCoeff_Uw,Achica_t[i]);
					sprintf(sCoeff_Uw,"%s[3010]g->VolFinito[i].volumen=%g\n",sCoeff_Uw,g->VolFinito[i].volumen);
					sprintf(sCoeff_Uw,"%s[3010]Datos_dt=%g\n",sCoeff_Uw,Datos_dt);
					sprintf(sCoeff_Uw,"%s[3010]CoeffTemporales=%g\n",sCoeff_Uw,CoeffT*Achica_t[i]);
				}


				SumCoeffTemporales += CoeffT*Achica_t[i];

				if (ConcentracionPrevia[i]>0.5) {


					if (DBGVARZ2) {
						double AchicaPreviai=AchicaPrevia[i];
						cout <<"i="<<i<<endl;
						cout<<"Achica_t[i]="<<Achica_t[i]<<endl;
						cout<<"AchicaPrevia[i]="<<AchicaPrevia[i]<<endl;
						cout<<"CoeffT*Achica_t[i]="<<CoeffT*Achica_t[i]<<endl;
						cout<<"CoeffT*AchicaPrevia[i]*ConcentracionPrevia[i]="<<CoeffT*AchicaPrevia[i]*ConcentracionPrevia[i]<<endl;
						cout<<"Cuociente: ="<<AchicaPrevia[i]*ConcentracionPrevia[i]/Achica_t[i]<<endl;
						nprintDBGVARZ++;
						if (nprintDBGVARZ>10) exit(0);
					}
				}
			}
			NCoeffEstacionarios++;
			if (Coeff_Uw<1e-5) {
				cout<<"3014: Coeff_Uw<1e-5 ="<<Coeff_Uw<<endl;


				if (DBGVARZ) {
					cout <<sCoeff_Uw<<endl;

					nprintDBGVARZ++;
					if (nprintDBGVARZ>10) exit(0);
				}

			}
			ConcentracionNueva = LadoDerecho/Coeff_Uw;
			err=fabs(ConcentracionNueva-Concentracion[i]);
			if (err > errG) {
				errG=err;
			}
			Concentracion[i]=ConcentracionNueva;
		}
		if (errG<err0 && iiter>3) {
			if (DBGVARZ) cout<<"break: errG<"<<err0<<endl;
			break;
		}
		if (DBGVARZ) cout<<"errG[i]="<<errG<<endl;
	} //Fin niteraciones
	if(DBGVARZ)cout<<"iiter_final="<<iiter;
	if(DBGVARZ)cout<<"   errG="<<errG<<endl;
	if (TipoCalculo==CalculoEvolucion) {
		//Estadistica (cual es mas importante)
		if(DBGVARZ) {
			cout << "<CoeffDifusion>="<<SumaDifusion /NCoeffEstacionarios;
			cout << "<SumaConveccion>="<<SumaConveccion /NCoeffEstacionarios;
			cout << "<CoeffEstacionarios>="<<SumCoeffEstacionarios /NCoeffEstacionarios;
			cout << "\t<SumCoeffTemporales>="<<SumCoeffTemporales / NCoeffEstacionarios <<endl;
		}

		Peclet=SumaConveccion/SumaDifusion;

		if (SumCoeffTemporales<2*SumCoeffEstacionarios) {
			FactorClockTime *= 0.5;

			if (glui != NULL) glui->sync_live();
		}
		if (SumCoeffTemporales>10*SumCoeffEstacionarios) {
			FactorClockTime *= 2;

			if (glui != NULL) glui->sync_live();
		}
	}

	double minF=1e20,maxF=-1e20;
	for (i=0;i<g->nVolFinito;i++) {
		if (minF>Concentracion[i]) minF=Concentracion[i];
		if (maxF<Concentracion[i]) maxF=Concentracion[i];

	}
	if(DBGVARZ) cout<<"2703: minF,maxF="<<minF<<" , "<<maxF<<endl;

	if(DBGVARZ) cout<<"CalculaConveccionDifusionDt_DominioVariableZ(...):END"<<endl;
	return(errG);
}


double CalculaFactorMallaCCRGJSM(grid3D g)
{
	int iiter,i,j,iC;
	double Factor,FactorSup,FactorInf,FactorSum=0,areamin,volmin,dijmin;
	int nFactorSum=0;
	//Resuelvo la ecuación varias veces
	FactorSup=-1e20;
	FactorInf=1e20;
	for (i=0 ; i<g.nVolFinito ; i++) {
		//Estudio el Factor en cada VolumenFinito
		double Voli=g.VolFinito[i].volumen;
		for (j=0 ; j< g.VolFinito[i].vecino.size() ;j++) {
			//recorro cada vecino==cara [i]--[j]

			Factor=g.VolFinito[i].Poligono[j].AreaPoligono * g.VolFinito[i].Poligono[j].Dab/ Voli;
			FactorSum += Factor;
			nFactorSum++;
			if (Factor>FactorSup) {
				FactorSup=Factor;
			}
			if (Factor<FactorInf) {
				FactorInf=Factor;
				areamin=g.VolFinito[i].Poligono[j].AreaPoligono;
				volmin=Voli;
				dijmin=g.VolFinito[i].Poligono[j].Dab;
				cout<< "FactorMallaCCRGJSM(Inf)="<<FactorInf<<endl;
				cout<< "Volfinito="<<i<<" de "<<g.nVolFinito<<endl;;
				cout<< "areamin="<<areamin;
				cout<< "\tdijmin="<<dijmin;
				cout<< "\tvolmin="<<volmin<<endl;
			}

		}

	}

	cout<< "FactorMallaCCRGJSM(Inf)="<<FactorInf<<endl;
	cout<< "areamin="<<areamin;
	cout<< "\tdijmin="<<dijmin;
	cout<< "\tvolmin="<<volmin<<endl;

	cout<< "FactorMallaCCRGJSM(Sup)="<<FactorSup<<endl;
	cout<< "FactorMallaCCRGJSM(Med)="<<FactorSum/nFactorSum<<endl;

	return(FactorSup);
}


double ResuelveConveccionDifusion (vector<double> &Concentracion,vector<double> &PotencialV,vector<double> &U,
		vector<double> &V,vector<double> &W,grid3D g,int niteraciones,double err0,vector<double> &TempPrevia)
{
	int iiter,i,j,iC;
	double SumaCoeffDif_F2,SumaCoeffDif,Coeff_Difusion,err,errG,nF,SumaVolumen;
	if (Concentracion.size() ==0) Concentracion.resize(g.nH3D);
	//Resuelvo la ecuación varias veces
	niteraciones=g.nH3D*50;
	if (niteraciones<30000) niteraciones=30000;
	cout <<"niteraciones="<<niteraciones<<endl;
	myfileSalida <<"niteraciones="<<niteraciones<<endl;
	double SumCoeffEstacionarios=0,SumCoeffTemporales=0;
	int    NCoeffEstacionarios=0;
	for (iiter=0 ; iiter<niteraciones ; iiter++) {
		errG=0;
		for (i=0 ; i<g.nH3D ; i++) {
			//Escribo la ecuación para cada hexahedro [i]
			SumaCoeffDif_F2=0;SumaCoeffDif=0;SumaVolumen=0;
			double PotencialV_i=PotencialV[i],PotencialV_j;
			for (j=0 ; j< g.h3D[i].vecino.size() ;j++) {
				if(1==0) cout <<"CalculaLaplacianoCero2: iiter,i,j="<<iiter<<" "<<i<<" "<<j
						<< "\tPotencialV_.size=" << PotencialV.size()
						<<"\tg.h3D[i].vecino[j]="<<g.h3D[i].vecino[j]
																   <<"\tg.h3D[i].tipo_vecino[j]="<<g.h3D[i].tipo_vecino[j]
																														<<endl;
				//recorro cada vecino==cara [i]--[j]

				if (g.h3D[i].tipo_vecino[j] == ES_BLOQUE) {
					PotencialV_j=PotencialV[ g.h3D[i].vecino[j] ];
					if (PotencialV_j>PotencialV_i)
						Coeff_Difusion = Datos_rhof * Datos_cf *(PotencialV_j-PotencialV_i) *
						g.h3D[i].Poligono[j].AreaPoligono / g.h3D[i].Poligono[j].Dab;
					else Coeff_Difusion=0;

					Coeff_Difusion += Datos_km_mu * g.h3D[i].Poligono[j].AreaPoligono / g.h3D[i].Poligono[j].Dab;

					SumaCoeffDif_F2 += Coeff_Difusion * Concentracion[ g.h3D[i].vecino[j] ] ;
					SumaCoeffDif    += Coeff_Difusion;
				} else {
					iC=g.h3D[i].vecino[j];
					double htilde_local;
					//if (g.Cara[iC].iBC==1) Se trata de Derivada normal igual a cero (nada que hacer)
					switch (g.Cara[iC].iBC) {
					case 2:
					case 3:   //Condicion Dirichlet (como un bloque con T conocida == BC2)
						PotencialV_j=g.Cara[iC].BC;
						if (PotencialV_j>PotencialV_i)
							Coeff_Difusion = Datos_rhof * Datos_cf *(PotencialV_j-PotencialV_i) *
							g.h3D[i].Poligono[j].AreaPoligono / g.h3D[i].Poligono[j].Dab;
						else Coeff_Difusion=0;

						Coeff_Difusion += Datos_km_mu * g.h3D[i].Poligono[j].AreaPoligono / g.h3D[i].Poligono[j].Dab;

						if (1==0) cout<<"g.Cara[iC].BC2="<<g.Cara[iC].BC2
								<<"\t.BC="<<g.Cara[iC].BC
								<<"\tg.Cara[iC].iBC="<<g.Cara[iC].iBC
								<<"\tiC="<<iC
								<<"\tCoeff_Difusion="<<Coeff_Difusion
								<<"\tPotencialV_i="<<PotencialV_i<<"PotencialV_j="<<PotencialV_j
								<<endl;
						SumaCoeffDif_F2 += Coeff_Difusion * g.Cara[iC].BC2;
						SumaCoeffDif    += Coeff_Difusion;
						break;

					case 11:   // Cara convectiva abajo
						htilde_local=MediaArmonica(Datos_km_mu/(g.h3D[i].Poligono[j].Dab/2),Datos_hc_inferior);

						Coeff_Difusion=htilde_local*g.h3D[i].Poligono[j].AreaPoligono;
						SumaCoeffDif_F2 += Coeff_Difusion * Datos_Tambiente;
						SumaCoeffDif    += Coeff_Difusion;
						break;
					case 12:   // Cara convectiva lateral
						htilde_local
						=MediaArmonica(Datos_km_mu/(g.h3D[i].Poligono[j].Dab/2),Datos_km2/Datos_DistanciaAlBorde);
						htilde_local=MediaArmonica(htilde_local,Datos_hc_lateral);
						Coeff_Difusion=htilde_local*g.h3D[i].Poligono[j].AreaPoligono;
						SumaCoeffDif_F2 += Coeff_Difusion * Datos_Tambiente;
						SumaCoeffDif    += Coeff_Difusion;
						break;
					}

				}
			}

			if (TipoCalculo==CalculoEvolucion) {
				///////////////////Aporte del termino temporal.....
				double CoeffT=Datos_rhom*Datos_cm*g.h3D[i].volumen/Datos_dt;

				//Estadistica (cual es mas importante)
				SumCoeffEstacionarios += SumaCoeffDif;
				SumCoeffTemporales    += CoeffT;
				NCoeffEstacionarios   ++;


				SumaCoeffDif_F2 += CoeffT*TempPrevia[i];
				SumaCoeffDif    += CoeffT;
			}
			nF = SumaCoeffDif_F2/SumaCoeffDif;
			err=fabs(nF-Concentracion[i])/(fabs(nF)+1e-20);
			if (err > errG) errG=err;
			Concentracion[i]=nF;
			//			if (F[i]==0) F[i]=0.5;
		}
		if (errG<err0) {
			if (DBG) cout<<"break: errG<"<<err0<<endl;
			break;
		}
	}
	cout<<"iiter_final="<<iiter<<endl;
	myfileSalida<<"iiter_final="<<iiter<<endl;
	if (TipoCalculo==CalculoEvolucion) {
		//Estadistica (cual es mas importante)
		cout << "<CoeffEstacionarios>="<<SumCoeffEstacionarios /NCoeffEstacionarios;
		cout << "\t<SumCoeffTemporales>="<<SumCoeffTemporales / NCoeffEstacionarios <<endl;
		myfileSalida << "<CoeffEstacionarios>="<<SumCoeffEstacionarios /NCoeffEstacionarios;
		myfileSalida << "\t<SumCoeffTemporales>="<<SumCoeffTemporales / NCoeffEstacionarios <<endl;
	}
	for (i=0 ; i<g.nH3D ; i++) {
		double A,AA,AB,AC,RA;
		double BA,BB,BC,RB;
		double CA,CB,CC,RC;
		double nx,ny,nz,DET;
		AA=0;AB=0;AC=0;RA=0;
		BA=0;BB=0;BC=0;RB=0;
		CA=0;CB=0;CC=0;RC=0;
		for (j=0 ; j< g.h3D[i].vecino.size() ;j++) {
			if (g.h3D[i].tipo_vecino[j] == ES_BLOQUE) {
				Coeff_Difusion = ( Concentracion[ g.h3D[i].vecino[j] ] - Concentracion[i] ) / g.h3D[i].Poligono[j].Dab;
			} else {
				iC=g.h3D[i].vecino[j];
				Coeff_Difusion=0;
				if (g.Cara[iC].iBC ==2 ||g.Cara[iC].iBC ==3) {
					Coeff_Difusion = ( g.Cara[iC].BC2 - Concentracion[i] ) / g.h3D[i].Poligono[j].Dab;

				}
			}
			A=g.h3D[i].Poligono[j].AreaPoligono;
			nx=g.h3D[i].Poligono[j].normal.x ;
			ny=g.h3D[i].Poligono[j].normal.y ;
			nz=g.h3D[i].Poligono[j].normal.z ;
			AA += nx*nx*A; AB+= nx*ny*A; AC+= nx*nz*A; RA += -nx*Coeff_Difusion*A;
			BA += ny*nx*A; BB+= ny*ny*A; BC+= ny*nz*A; RB += -ny*Coeff_Difusion*A;
			CA += nz*nx*A; CB+= nz*ny*A; CC+= nz*nz*A; RC += -nz*Coeff_Difusion*A;
		}
		// AA AB AC --> RA
		// BA BB BC --> RB
		// CA CB CC --> RC
		DET= AA*(BB*CC-BC*CB)+BA*(CB*AC-AB*CC)+CA*(AB*BC-BB*AC) ;
		U[i]= ( RA*(BB*CC-BC*CB)+RB*(CB*AC-AB*CC)+RC*(AB*BC-BB*AC) )/DET;
		V[i]=-( RA*(BA*CC-BC*CA)+RB*(CA*AC-AA*CC)+RC*(AA*BC-BA*AC) )/DET;
		W[i]= ( RA*(BA*CB-BB*CA)+RB*(CA*AB-AA*CB)+RC*(AA*BB-BA*AB) )/DET;
	}
	return(errG);
}

void ResetView()
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	cout <<"Antes:"<<endl;
	cout <<"Escala="<<Escala<<endl;
	int i;
	for (i=0;i<16;i++){
		MatrizRotacionGlobal[i]   =M1_MatrizRotacionGlobal[i];
		MatrizRotacionGlobalINV[i]=M1_MatrizRotacionGlobalINV[i];
	}
	cout <<"vecXEsfera=["<<vecXEsfera[0]<<","<<vecXEsfera[1]<<","<<vecXEsfera[2]<<"]"<<endl;
	cout <<"vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<"]"<<endl;
	if (gtotal != NULL) {
		double e1,e2;
		e1=8/(gtotal->xmax-gtotal->xmin);
		e2=6/(gtotal->ymax-gtotal->ymin);
		if (e1<e2)
			Escala=e1;
		else
			Escala=e2;
		vecUEsfera[0]=(vecXEsfera[0]-(gtotal->xmax+gtotal->xmin)/2)*Escala-1;
		vecUEsfera[1]=(vecXEsfera[1]-(gtotal->ymax+gtotal->ymin)/2)*Escala+1;
		vecUEsfera[2]=(vecXEsfera[2]-(gtotal->zmax+gtotal->zmin)/2)*Escala;
	}
	if (gPoly != NULL) {
		gPoly->minmax();
		double e1,e2;
		e1=8/(gPoly->xmax-gPoly->xmin);
		e2=6/(gPoly->ymax-gPoly->ymin);
		if (e1<e2)
			Escala=e1;
		else
			Escala=e2;
		vecUEsfera[0]=(vecXEsfera[0]-(gPoly->xmax+gPoly->xmin)/2)*Escala-1;
		vecUEsfera[1]=(vecXEsfera[1]-(gPoly->ymax+gPoly->ymin)/2)*Escala+1;
		vecUEsfera[2]=(vecXEsfera[2]-(gPoly->zmax+gPoly->zmin)/2)*Escala;
	}

	cout <<"Despues:"<<endl;
	cout <<"Escala="<<Escala<<endl;
	cout <<"vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<"]"<<endl;
}



void lectura_archivo_UVWZdeT_bin() {
cout<<"Leo archivo binario"<<endl;
readUVWZdeT_Binario("../../../Datos/UVWZdeT.bin")  ;
cout<<"Ok"<<endl;

pthread_mutex_unlock(&mutex_DrawRead);

sprintf(text,"%sLeido ../../../Datos/UVWZdeT.bin\n",text);if (glui != NULL) glui_edittext->set_text(text);

//Condición Inicial
Datos_dt=20*60; //Todo //Corregir
Tsimulacion=0;
Datos_km_mu=0.03;
cout<<"Tsimulacion="<<Tsimulacion<<"step_hora="<<step_hora<<endl;
}

int Iteracion_Conveccion_Difusion() {
	int i;
	DrawGraphics();
	if (U.size()==0) {

		pthread_mutex_lock(&mutex_DrawRead);

		lectura_archivo_UVWZdeT_bin();
		if (gtotal!=NULL && Datos_dt>0)
			for (i=0;i<gtotal->nVolFinito;i++) {
				F2VolumenesP[i]=F2Volumenes[i];
			}

	}


	Tsimulacion+=Datos_dt;
	Convierte_Tsimulacion_en_StepHora();
	for (int ical=0;ical<10000;ical++) {
		while (MODO_Pausa) {
			sleep(1);
		}
		sprintf(text,"Start: Una iteración %d\n",ical);if (glui != NULL) glui_edittext->set_text(text);

		DrawGraphics();
		if (DBG) cout <<"U.size()="<<U.size()<<endl;
		if (DBG) cout <<"U[step].size()="<<U[step_hora].size()<<endl;

		pthread_mutex_lock(&mutex_DrawRead);
		//CalculaConveccionDifusionDt(F2Volumenes, gtotal, 10000, 1e-6, F2VolumenesP);
		Global_error=1e-7;
		pthread_mutex_unlock(&mutex_DrawRead);
		int FlagTrataBC=0;
		double minF,maxF;
		ProyectaVolumenesAVertices(F2Volumenes,F2Nodos,&minF,&maxF,FlagTrataBC) ;

		sprintf(text,"%s minF=%f, maxF=%f,  errG=%f\n",text,minF,maxF,Global_error);if (glui != NULL) glui_edittext->set_text(text);
		//sprintf(text,"%s errG=%f\n",text,errG);//if (glui != NULL) glui_edittext->set_text(text);

		if (Global_error<1e-6) {
			Tsimulacion+=Datos_dt;
			Convierte_Tsimulacion_en_StepHora();
			for (i=0;i<gtotal->nVolFinito;i++) {
				F2VolumenesP[i]=F2Volumenes[i];
			}
		}

	}

}


void CB_keyboard(unsigned char key, int x, int y)
{
	CuentaLlamadas[__PRETTY_FUNCTION__+string(" @ ")+static_cast< std::ostringstream & >((std::ostringstream()<<__LINE__)).str()]++;

	double gLeft=-.1,gRight=.1, gBottom=-.1, gTop=.1, gNear=0, gFar=1000000;
	//cout<<"CB_keyboard()"<<endl;
	//printf("CB_keyboard: char=%d(%c), x=%d, y=%d",key,key,x,y);	cout<<endl;
	static int primerCB_keyboard=1;

	if(primerCB_keyboard) {
		sprintf(textHelp,"Esc: Exit\n");
		sprintf(textHelp,"%s   : Etapa Siguiente\n",textHelp);
		sprintf(textHelp,"%sBs : Etapa Previa Registrada\n",textHelp);
		sprintf(textHelp,"%s1  : Nada\n",textHelp);
		sprintf(textHelp,"%sA  : Add Particulas\n",textHelp);
		sprintf(textHelp,"%sB  : Bordes\n",textHelp);
		sprintf(textHelp,"%sC  : Color \n",textHelp);
		sprintf(textHelp,"%se  : Mouse== Escala \n",textHelp);
		sprintf(textHelp,"%sF  : Fullscreen \n",textHelp);
		sprintf(textHelp,"%sG  : Modo Game \n",textHelp);
		sprintf(textHelp,"%sH  : Borra tapa superior/inferior/no \n",textHelp);
	}
	primerCB_keyboard=0;
	switch(key) {
	case 27: /* ESC */

		cout <<"Calls: Cara3D::draw_caraGL()="<<Cara3D__draw_caraGL<<endl;
		cout <<"Calls: PoligonoPlano::readArchivoPol(ifstream,grid3D)="<<PoligonoPlano__readArchivoPol__ifstream_grid3D<<endl;
		cout <<"Calls: Lectura_Bahia_y_CalculosGeometricosMalla(char*)="<<Lectura_Bahia_y_CalculosGeometricosMalla_char<<endl;

/*
		itCuentaLlamadas = CuentaLlamadas.begin();
		while(itCuentaLlamadas != CuentaLlamadas.end())
		{
			cout<< setw(7) <<itCuentaLlamadas->second<<" : "<<itCuentaLlamadas->first<<endl;
			itCuentaLlamadas++;
		}
*/


		itCuentaLlamadas = CuentaLlamadas.begin();
		while(itCuentaLlamadas != CuentaLlamadas.end())
		{
//			dst[itCuentaLlamadas->second]=itCuentaLlamadas->first;

			dst.insert(pair<int, std::string> (itCuentaLlamadas -> second, itCuentaLlamadas -> first));

			itCuentaLlamadas++;
		}



		it_dst = dst.begin();
		while(it_dst != dst.end())
		{
			cout<< setw(7) <<it_dst->first <<" : "<<it_dst->second<<endl;
			it_dst++;
		}



#if 0
		std::multimap<int,std::string> dst = flip_map(CuentaLlamadas);

		// Declaring the type of Predicate that accepts 2 pairs and return a bool
		typedef std::function<bool(std::pair<std::string, int>, std::pair<std::string, int>)> Comparator;

		// Defining a lambda function to compare two pairs. It will compare two pairs using second field
		Comparator compFunctor =
				[](std::pair<std::string, int> elem1 ,std::pair<std::string, int> elem2)
				{
			return elem1.second < elem2.second;
				};

		// Declaring a set that will store the pairs using above comparision logic
		std::set<std::pair<std::string, int>, Comparator> setOfWords(
				mapOfWordCount.begin(), mapOfWordCount.end(), compFunctor);

		// Iterate over a set using range base for loop
		// It will display the items in sorted order of values
		for (std::pair<std::string, int> element : setOfWords)
			std::cout << element.first << " :: " << element.second << std::endl;

#endif

		exit(0);
		break;

	case 32: /* Espacio */
		cout<<"Lanzo Calculo_EtapaS()"<<Etapa<<endl;
		Calculo_EtapaS();
		break;
	case 8: /* BackSpace */
		cout<<"Etapa="<<Etapa<<",InicioEtapa="<<InicioEtapa<<endl;
		Etapa=InicioEtapa-1;
		Calculo_EtapaS();
		break;
	case '1': 
		//Calculo_EtapaS(true);
		break;
	case 'A': //Add particulas con el mouse
	case 'a':
		Add_Particulas=1-Add_Particulas;
		if (Add_Particulas) {
			AddMensajeRight((char *)"A=Add Particulas");
			MueveCentro=0;
			//Add_Particulas=0;
			Add_VolumenINI=0;
			Add_Voronoi=0;
			Flag_MuestraCara=0;

		} else {
			AddMensajeRight("");AddMensajeRight("");
		}
		break;
	case 'b':
	case 'B':
		ModoDibujaFrontera=1-ModoDibujaFrontera;
		if (glui != NULL) glui->sync_live();
		break;
	case 'C':
	case 'c':
		cout<<"case 'C': ColorON="<<ColorON<<"-->";
		ColorON++;
		if (ColorON==1 & F.size()!=gtotal->nV3D) ColorON++;
		//if (ColorON==1 & EtapaGlobal==ETAPA_CALCULO_T_PILA) ColorON++;

		if (ColorON==2 & F2Nodos.size()!=gtotal->nV3D) ColorON++;


		if (ColorON>3)
			ColorON=0;
		if (ColorON)    AddMensaje("C= Color on..");
		else                AddMensaje("C= Color off..");
		cout<<ColorON<<endl;

		break;

	case 'd':
		//Reset View
		ResetView();

		break;
	case 'D':
		//Reset View 2
		int i;
		Escala=M2_Escala;
		for (i=0;i<16;i++){
			MatrizRotacionGlobal[i]   =M2_MatrizRotacionGlobal[i];
			MatrizRotacionGlobalINV[i]=M2_MatrizRotacionGlobalINV[i];
		}
		cout <<"vecXEsfera=["<<vecXEsfera[0]<<","<<vecXEsfera[1]<<","<<vecXEsfera[2]<<"]"<<endl;
		cout <<"vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<"]"<<endl;
		vecDUEsfera[0]=(vecXEsfera[0]-(gtotal->xmax+gtotal->xmin)/2)*Escala;
		vecDUEsfera[1]=(vecXEsfera[1]-(gtotal->ymax+gtotal->ymin)/2)*Escala;
		vecDUEsfera[2]=(vecXEsfera[2]-(gtotal->zmax+gtotal->zmin)/2)*Escala;

		MatrizXvector4(MatrizRotacionGlobal,vecDUEsfera,vecUEsfera);

		cout <<"vecUEsfera=["<<vecUEsfera[0]<<","<<vecUEsfera[1]<<","<<vecUEsfera[2]<<"]"<<endl;
		break;
	case 'E':
	case 'e':
		glui_GrupoModoDelMouse->set_int_val(3);
		MODO_de_mover=glui_GrupoModoDelMouse->get_int_val();
		break;
	case 'F':
	case 'f':
		fullscreen = !fullscreen;
		if (fullscreen) {
			old_x = glutGet(GLUT_WINDOW_X);
			old_y = glutGet(GLUT_WINDOW_Y);
			old_width = glutGet(GLUT_WINDOW_WIDTH);
			old_height = glutGet(GLUT_WINDOW_HEIGHT);
			glutFullScreen();
		} else {
			glutReshapeWindow(old_width, old_height);
			glutPositionWindow(old_x, old_y);
		}
		break;

	case 'G':
	case 'g':
		ModoGame++;if (ModoGame>3) ModoGame=0;
		printf("ModoGame=%d\n",ModoGame);
		switch (ModoGame)
		{
		case 1:
		case 2:
		case 3:
			if (ModoGame==1) glutGameModeString("1280x1024:32@60");
			if (ModoGame==2) glutGameModeString("800x600:16@60");
			if (ModoGame==3) glutGameModeString("640x480:16@60");
			glutEnterGameMode();

			main_window = glutGetWindow();
			//		initWindow();
			glutDisplayFunc(DrawGraphics);
			glutReshapeFunc(ResizeGraphics);
			glutMouseFunc(CB_mouse);
			glutMotionFunc(CB_motion);
			glutKeyboardFunc(CB_keyboard);
			glutSpecialFunc( CB_keyboardSpecial );

			glutIdleFunc( idleevent );
			glutVisibilityFunc(visible);
			inicializacionGL();
			DrawGraphics();

			break;
		case 0:

			glutLeaveGameMode();
			main_window = glutGetWindow();
			glutDisplayFunc(DrawGraphics);
			glutReshapeFunc(ResizeGraphics);
			glutMouseFunc(CB_mouse);
			glutMotionFunc(CB_motion);
			glutKeyboardFunc(CB_keyboard);
			glutSpecialFunc( CB_keyboardSpecial );

			glutIdleFunc( idleevent );
			glutVisibilityFunc(visible);
			inicializacionGL();
			DrawGraphics();
			break;
		}
		break;

		case 'h':
		case 'H':
			DibujaSupInf++; if (DibujaSupInf>2) DibujaSupInf=0;  //0: Dibujatodo, 1:borra tapa superior, 2:borra tapa inferior
			break;
		case 'i':
		case 'I':
			ModoDibujaInterior=1-ModoDibujaInterior;
			if (glui != NULL) glui->sync_live();
			break;
		case 'j':
		case 'J':
			FlagPrintEscala=1-FlagPrintEscala;
			break;

		case 'k':
		case 'K':
			FlagPrintMatrizRotacionGlobal=1-FlagPrintMatrizRotacionGlobal;
			break;

		case 'l':
		case 'L':
		{
			R3 CentroMancha;
			CentroMancha.x=(gtotal->xmax+gtotal->xmin)/2;
			CentroMancha.y=(gtotal->ymax+gtotal->ymin)/2;
			CentroMancha.z=(gtotal->zmax+gtotal->zmin)/2;
			double RadioMancha=(gtotal->xmax-gtotal->xmin)/5;
			for (i=0;i<gtotal->nV3D;i++) {
				if( sqr(gtotal->v3D[i].x-CentroMancha.x)+sqr(gtotal->v3D[i].y-CentroMancha.y)<sqr(RadioMancha)) {
					F2Nodos[i]=1;
				} else {
					F2Nodos[i]=0;
				}
			}
			for (i=0;i<gtotal->nVolFinito;i++) {
				if( sqr(gtotal->VolFinito[i].centro.x-CentroMancha.x)+sqr(gtotal->VolFinito[i].centro.y-CentroMancha.y)<sqr(RadioMancha)) {
					F2Volumenes[i]=1;
					F2VolumenesP[i]=1;
				}else {
					F2Volumenes[i]=0;
					F2VolumenesP[i]=0;
				}
			}
			FlagEvoluciona=0;
		}
		break;

		case 'M':
		case 'm':
			MueveCentro=1-MueveCentro;
			//		gluiMueve->set_int_val(!MueveCentro);MueveCentro=gluiMueve->get_int_val();
			if (MueveCentro)    {
				AddMensajeRight((char *)"M=Mueve Centro");
				//			cout<<"E vecUEsfera="<<vecUEsfera[0]<<","<<vecUEsfera[1]<<"."<<vecUEsfera[2]<<","<<vecUEsfera[3]<<endl;

				//MueveCentro=0;
				Add_Particulas=0;
				Add_VolumenINI=0;
				Add_Voronoi=0;
				Flag_MuestraCara=0;
			}
			else  {
				AddMensajeRight("");AddMensajeRight("");
			}
			break;

		case 'n':
		case 'N':
			gluiNormales->set_int_val(!gluiNormales->int_val);
			break;
		case 'O':
		case 'o':

			if (gluiOpciones)	{
				gluiOpciones->show();

			}
			break;
		case 'P':
		case 'p':
			if (FromCB==0)
				gluiClipping->set_int_val(!ClippingON);ClippingON=gluiClipping->get_int_val();
			FromCB=0;
			/*
		ClippingON=1-ClippingON;
			 */
			if (ClippingON)    AddMensaje("P= Cli(P)ping on..");
			else                AddMensaje("P= Cli(P)ping off..");
			break;

		case 'Q':
		case 'q':
			writeTofileBMP=1-writeTofileBMP;
			iframe=-1;
			break;

		case 'R':
		case 'r':

			//		glLoadIdentity();
			//		glOrtho(gLeft, gRight, gBottom, gTop, gFar, -gFar);
			//		glGetFloatv( GL_MODELVIEW_MATRIX, (GLfloat*)MatrizRotacionGlobal );
			//		TrasponeMatriz();

			glui_GrupoModoDelMouse->set_int_val(1);MODO_de_mover=glui_GrupoModoDelMouse->get_int_val();
			break;
		case 'S':
		case 's':
			glui_GrupoModoDelMouse->set_int_val(2);MODO_de_mover=glui_GrupoModoDelMouse->get_int_val();
			break;
		case 'T':
		case 't':
			glui_GrupoModoDelMouse->set_int_val(0);MODO_de_mover=glui_GrupoModoDelMouse->get_int_val();
			break;

		case 'u':
		case 'U':

			FlagEvoluciona=1;
			break;
		case 'v':
		case 'V':
			gluiHelp->set_int_val(!MODO_MenuMENSAJES);
			MODO_MenuMENSAJES=gluiHelp->get_int_val();
			if (MODO_MenuMENSAJES)
				DrawMensajesDatos();
			break;
		case 'W':
		case 'w':
			DrawMensajesDatosGui3();
			glui3->show();
			break;
			///////////////////////////////////////////////////

		case 'X':
		case 'x':
			FlagDifusion=1-FlagDifusion;
			break;
			///////////////////////////////////////////////////

		case 'Y':
		case 'y':
			FlagConveccion=1-FlagConveccion;
			break;
			///////////////////////////////////////////////////

		case 'Z':
		case 'z':
#if defined(GLUI_GLUI_H)
			if (glui_hide)
				glui->show();
			else
				glui->hide();
			glui_hide=!glui_hide;
			//			ResizeGraphics(glutGet(GLUT_WINDOW_WIDTH),glutGet(GLUT_WINDOW_HEIGHT));
#endif
			break;

		case ',':
			Add_Voronoi=1-Add_Voronoi;
			if (Add_Voronoi) {
				AddMensajeRight((char *)",=Add Voronoi");
				MueveCentro=0;
				Add_Particulas=0;
				Add_VolumenINI=0;
				//Add_Voronoi=0;
				Flag_MuestraCara=0;

			} else {
				AddMensajeRight("");AddMensajeRight("");
			}
			break;

		case '.':
			Add_VolumenINI=1-Add_VolumenINI;
			if (Add_Voronoi) {
				AddMensajeRight((char *)",=Add VolumenINI");
				MueveCentro=0;
				Add_Particulas=0;
				//Add_VolumenINI=0;
				Add_Voronoi=0;
				Flag_MuestraCara=0;

			} else {
				AddMensajeRight("");AddMensajeRight("");
			}
			break;

		case '-':
			Flag_MuestraCara=(Flag_MuestraCara+1 % 2);
			if (Flag_MuestraCara>0) {
				if (Flag_MuestraCara==1 )AddMensajeRight((char *)",=MuestraCara,");
				if (Flag_MuestraCara==2 )AddMensajeRight((char *)",=MuestraPoligonoVoronoi,");
				MueveCentro=0;
				Add_Particulas=0;
				Add_VolumenINI=0;
				Add_Voronoi=0;

			} else {
				AddMensajeRight("");AddMensajeRight("");
			}
			break;


		case '?':
			if (gluiHelp2)	{

				TextBoxHelp->set_text(textHelp);glui_edittext->redraw_window();


				gluiHelp2->show();

			}
			break;
			///////////////////////////////////////////////////

		case '+':
			if (NumON==1 ) {
				NumEscala /=1.2;
			} else {
				char s[100];
				npasadas*=2;
				if (npasadas>maxpasadas) 	npasadas=maxpasadas;
				sprintf(s,"(+)=npasadas++=%d",npasadas);
				AddMensaje(s);
			}
			break;
		case '|':
			if (NumON==1 ) {
				NumEscala *=1.2;
			} else {
				char s[100];
				npasadas/=2;
				if (npasadas<=0) 				npasadas=1;
				sprintf(s,"(-)=npasadas--=%d",npasadas);
				AddMensaje(s);
			}
			break;
	}
	//cout<<"CB_keyboard():END"<<endl;
}

