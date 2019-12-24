#pragma once


//#include "../../EjemploOpenGL/src/grid3D.h"
//#include "../../EjemploOpenGL/src/FuncionesOpenGL.h"
#include "grid3D.h"
#include <pthread.h>




#include "FuncionesOpenGL.h"
//#include "../../Libreria-JSM/Class_Vector.h"
#include "Class_Vector.h"

#include "glui-jsm.h"

#include "Macros.h"

///Malla
int NDivZ=2,nRLC=-1,nThLC=-1;
int refinoR=0;
char FDatos[100],mensajes0[10000];
float FactorCercania=0.90,FactorAchica=0.0,FactorAchicaV=0.15,FactorZ=1;

float FactorAmbient=1,FactorDifusse=1,FactorSpecular=1,FactorEmission=1,FactorClaridad=1.7;
float FactorAmbientL1=0.35,FactorDifusseL1=0.99,FactorSpecularL1=0.99;
float FactorAmbientL2=0.35,FactorDifusseL2=0.99,FactorSpecularL2=0.99;
float GlobalFovy=30,GlobalOldFovy=30,GlobalCentros=0.05;
int BGColorR=130,BGColorG=100,BGColorB=82;

int x0Menu,y0Menu,x1Menu,y1Menu,EstadoMenu=0;

//Cuantas veces se llamó a cada funcion?
int 	Cara3D__draw_caraGL=0,
		PoligonoPlano__readArchivoPol__ifstream_grid3D=0,
		Lectura_Bahia_y_CalculosGeometricosMalla_char=0;


//Plano del clipping
//float Ax=0,By=-100,Cz=0,DD=-0.095;
int ClippingON=0,FromCB=0;
float Ax=1,By=0.17,Cz=0,DD=11000;

float TLimite=0.5;
int TLimite_if=0;

float FactorNormales=0,FactorSuavidad=1;

clock_t start_t, end_t;
double total_t;

//int CualMalla=3,nR=3,nTh=1,nZ=1;;
//int CualMalla=1,nR=5,nR2=10,nTh=5,nZ=1;double errMalla=1e-1;double Dominio_Rint=0.002; //Test 2
int CualMalla=1,nR=3,nR2=3,nTh=3,nZ=1;double errMalla=1e-1;double Dominio_Rint=0.002,Dominio_Rmed=0.075; //Test mallador
int CalculoContinuo=0,E_CalculoContinuo;

int TipoCalculo=1;
#define CalculoEstacionario 1
#define CalculoEvolucion 2

char nombre0[100];
ofstream myfileSalida,myfileVol,myfileDBG;

#define DBGMouse 0  //Mas en Macros.h


vector<int> CualesRehacer;

//CualMalla=2; // Dominio Completo 2*pi
//CualMalla=3; // Dominio Sector angulo dtheta
//CualMalla=4; // Dominio Sector angulo dtheta*10



char file_name[100];
/////
int Etapa=0,EtapaGlobal=0,EtapaGlobal2Local=0,InicioEtapaP=0,InicioEtapa=0,FinEtapas=0,Guardar=0;
int Reiterando;
static double err0;
int SaveReadMode;
///Preuba de botones
int xPos0,yPos0,xPos,yPos,iPush;
#define sqr(x) ((x)*(x))

int llamar_etapa_Siguiente_PAUSA2=0;
int NPOLY=10;
///GLUI....
#if defined(GLUI_GLUI_H)
GLUI_FileBrowser2 *FileBrowserGlui2;
GLUI_List *fb2;
GLUI *glui,*glui2,*glui3,*gluiHelp2,*gluiOpciones, *gluiEligeBC, *glui_OpcionesProgramaLevel0;
int gluiEligeBC_Visible,gluiEligeBC_w;
GLUI_Panel *gluiEligeBCpanel1;
GLUI_Panel *PanelFlimite,*PanelParticulas;

int PrintMouse=0,PrintLinea=1,PrintTiempos=1,PrintCaras=0,PrintInvert=0,PrintInicio=1,
		PrintFunciones=0,PrintParametros=0,PrintLecturas=1,PrintDZmin=0,PrintVecinos1=0;
GLUI_Checkbox *GC_PrintMouse;


char StringParametros[2000];
GLUI_TextBox *Glui3_EditParametros;

GLUI_StaticText *glui_FPS,*glui_MSG;

GLUI_Checkbox *Checkbox_particulas,*Checkbox_Veloc;
GLUI_Checkbox *gluiTranslation,*gluiRotation,*gluiMueve,*gluiNumera,*gluiHelp,*gluiNormales,*gluiBordes,*gluiInterior;
GLUI_TextBox *glui_edittext,*TextBoxHelp;
GLUI_Spinner *Spinner_particulas,*Spinner_Algunos;
int Algunos_Porcentaje=10;
vector<int> AlgunosI;
vector<int> AlgunosRND;

GLUI_StaticText *glui_porcentaje;
GLUI_Checkbox *gluiClipping;
int glui_hide=1;
void   formulario_glui();
GLUI_RadioGroup *glui_GrupoModoDelMouse;
#endif
#if defined(_GLUI_H)
GLUI *glui;
GLUI_StaticText *glui_FPS;
GLUI_Checkbox *gluiTranslation,*gluiRotation,*gluiMueve,*gluiNumera,*gluiHelp;
int glui_hide=1;
void   formulario_glui();
GLUI_RadioGroup *glui_GrupoModoDelMouse;
#endif

//Test1: pthreads

/* this variable is our reference to the second thread */
pthread_t inc_x_thread;
int x_thread=0,x_thread_p,thread_lanzado=0;
pthread_cond_t cond_calculo = PTHREAD_COND_INITIALIZER;
pthread_mutex_t mutex_calculo = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_DrawRead = PTHREAD_MUTEX_INITIALIZER;
int Inicia_calculo, Termine_calculo;



char  text0[200000] = {};
char  text[200000] = {"Hello World!"};
char  textHelp[200000] = {"Hello World!"};

double Global_error;

char  texto[200] = "*.msh3D";
///
clock_t clock0,clockF,clock2,clock2F;
clock_t clock00,clock_Time,clock_TimeP;
int Calculando=0;
float FactorClockTime=0.015,Tsimulacion=0,TsimulacionP=0,TColaP=0;
float FPS,FPS2;
int nframes=0, nframes2=0;



///MODOS
int fullscreen=0,old_x,old_y,old_width,old_height;
int NumON=0,ColorON=0,DefClippingON=0,MODO_Rotacion2=0,MODO_MenuMENSAJES=false;
int MODO_de_mover=1;
int ModoGame=0,MODO_Traslacion=0,MODO_Rotacion=1;
int MODO_Ejes=false;
int MODO_CampoVelocidades=false,MODO_Origen=false,MODO_Pausa=false;
int MODO_CampoVelocidades2=false,VelocidadEnHexahedros=false,VelocidadEnVertex=false;
int Modo_DibujaEdges=false;
int ModoDibujaNormales=false;
int ModoDibujaInterior=false,DibujaCentroVel=0;
int ModoDibujaFrontera=true;
int ModoDibujaAlgunos=false;
vector<int> CualesCarasDibuja;//Eliminar!!!
vector<int> CualesBloquesDibuja;
int Modo_DibujaCentroCaras=false;
int Modo_DibujaCentroBloques=false;
int MODO_NumeraH=false;
int MODO_NumeraFF=false;
int caso=3,CasoLectura=1,step_hora=0,stepp=0,deltat_hora=1800;
double FactorVel=0;
float stepf;
int FlagMuestraCaraSuperior=1,FlagMuestraCaraInferior=1,FlagDibujaLineas=0;
int FlagCaraSuperiorTransparente,FlagCaraSuperiorTextura,FlagCaraInferiorTextura=0;
int FlagPrintEscala=0,FlagPrintMatrizRotacionGlobal=0,Flag_MuestraCara=0;
int FlagCalculaEvolucion=0,FlagEvoluciona=0,FlagDifusion=1,FlagConveccion=0;
int ListaCaso=0;
double Peclet;

int binario=0;

int writeTofileBMP=0;
int iframe=-1,MaxFrames=300,DeboDibujarSiguiente=0;

float TimeToSave,TimeSaved,DT_toSave=0.1;


char* writedata;
char* imagendata;
char headerinfo[54];
GLuint textureID; //the array for our texture
int imagenwidth;
int imagenheight;

char* imagen2data;
GLuint textureID2; //the array for our texture
int imagen2width;
int imagen2height;

int Global_nia=1;
float factorV=10000,factorVmin=0,factorDespegaLineas=0.001;
float NumEscala=.05;
GLfloat wAngleX=0,wAngleY=0,wAngleZ=0;
float Escala=5;
GLfloat   MatrizRotacionGlobal0[16]   ={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
GLfloat   M1_MatrizRotacionGlobal[16]   ={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
GLfloat   M1_MatrizRotacionGlobalINV[16]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
GLfloat   M1_MatrizTraslacionGlobal[16] ={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};

float 	  M2_Escala=0.00358172;
GLfloat   M2_MatrizRotacionGlobal[16]= {0.169356 , -0.125896 , 0.977485 , 0 , 0.985559 , 0.021254 , -0.168017 , 0 , 0.000375712 , 0.99182 , 0.127677 , 0 , 0 , 0 , 0 , 1};
GLfloat   M2_MatrizRotacionGlobalINV[16]= {0.169356 , 0.985559 , 0.000376593 , 0 , -0.125896 , 0.0212531 , 0.99182 , 0 , 0.977485 , -0.168017 , 0.127677 , 0 , 0 , 0 , 0 , 1};

float 	  M3_Escala=0.00325611;
GLfloat   M3_MatrizRotacionGlobal[16]= {0.992515 , -0.0232533 , -0.119881 , 0 , -0.117091 , 0.0974664 , -0.988327 , 0 , 0.0346665 , 0.994966 , 0.094014 , 0 , 0 , 0 , 0 , 1};
GLfloat   M3_MatrizRotacionGlobalINV[16]= {0.992517 , -0.117092 , 0.0346663 , 0 , -0.0232536 , 0.0974663 , 0.994968 , 0 , -0.119881 , -0.988326 , 0.0940142 , 0 , 0 , 0 , 0 , 1};



GLfloat   MatrizRotacionGlobal[16]   ={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
GLfloat   MatrizRotacionGlobalINV[16]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
GLdouble   MatrizINV[16]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
GLfloat   MatrizTraslacionGlobal[16] ={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
GLfloat   vecXEsfera[4]={0,0,0,1},vecDXEsfera[4];
GLfloat   vecUEsfera[4]={0,0,0,1},vecDUEsfera[4];
GLfloat   vec1Clipping[3]={2,0,0},vec2Clipping[3]={0,0,2},vec3Clipping[3]={0,1,0},vec4Clipping[3]={0,0,0};
int width=1024,height=968;
GLfloat aspect;
int NMensajes=35;
char *Mensaje[35];
int NMensajeRight=0,MaxMensajeRight=2;
int NMensajeCenter=0,MaxMensajeCenter=5;
char *MensajeRight[2],*MensajeCenter[2];

GLuint baseF,baseBIT; 
int mode=0;
const int modeT=1,modeR=(1<<1);

int MueveCentro=0,Add_Particulas=0,Add_Voronoi=0,Add_VolumenINI=0,DibujaSupInf=0,Add_Poligono_Selected=0;
float CentroX=0,CentroY=0,CentroZ=0;
char s[100];
float VersionDatos;
//// GRID
grid3D g1,g2,*gtotal,*gPoly;

vector<vector<double> > ZdeT;

// DominioVariableZ
vector<vector<double> > FuncionF_Nodos ;
vector<vector<double> > FuncionF_Volumenes;
vector<vector<double> > FuncionDFDt_Volumenes;
vector<vector<R3> > Grad_F;
vector<vector<double> > Achica;
vector<double> AchicaPrevia,Achica_t;
int Achica_t_PrimeraVez=1;
int nprintDBGVARZ=0;
int MainMenu;


vector<vector<R3> > Coeff_H ; //Coeff_H[volumen_i][Cara_j]
vector<vector<double> > Coeff_Q ; //Coeff_Q[volumen_i][Cara_j]
vector<vector<double> > Coeff_Kn ; //Coeff_Kn[volumen_i][Cara_j]
vector<vector<double> > Coeff_Vtilde ; //Coeff_Vtilde[volumen_i][Cara_j]


vector<double > Coeff_1omega ; //Coeff_1omega[volumen_i]
vector<R3 > Coeff_2omega ; //Coeff_1omega[volumen_i]


// DominioVariableZ fin

vector<vector<double> > U,V,W;
vector<double> F;
vector<double> PotencialPsi,UU,VV,WW;
vector<double> Temperatura,UU2,VV2,WW2;
vector<double> F2Nodos;
vector<double> F2Volumenes,F2VolumenesP;

vector<double> PotencialVInfiltracion;
vector<double> TempPilaBloques,TempPilaVertices;
vector<double> TempPilaBloquesPrevia;

int DrawYES=1;
extern int npasadas;
extern  int const maxpasadas;
extern int largoCola;


void DrawMensajes();
void AddMensaje(char* newMensaje);


int   main_window;

void print_text(int x, int y, char* s) ;
void CB_keyboard(unsigned char key, int x=0, int y=0);
void Calculo_EtapaS_REAL(int inicializa=false);


#define Calculo_EtapaS printf("Calling FunctionName from %s\n",__FUNCTION__);Calculo_EtapaS_REAL

double CalculaLaplacianoCero(vector<double> &F,grid3D g,int niteraciones);

int KeyControlAltShift;

void fbreaddir2(GLUI_FileBrowser *fb, const char *d);

/// PArametros geometricos
double Dominio_Xmax=0.125;
double Dominio_Hmax=0.30;

double Dominio_Hsup=0.01;



double Datos_rhof=1000  ;//#densidad fluido agua
double Datos_rhos=1700  ;//#densidad solido
double Datos_rhom       ;//#densidad media de la pila

double Datos_phi=0.46   ;//#porosidad
double Datos_cf=4186    ; //calor especifico del fluido (agua)
double Datos_cs=1172    ; //calor especifico del solido
double Datos_cm         ; //calor especifico medio de la pila

double Datos_Tambiente=20   ;//% Temperatura ambiente
double Datos_Tinyeccion=65  ;//% Temperatura de inyeccion en borde superior
double Datos_hc_superior=5.7;
double Datos_hc_inferior=5.7;
double Datos_hc_lateral =5.7;

double Datos_KTermofilm=0.334	;//% W/(m*K)
double Datos_eTermofilm=140e-6    	;//% 140 mi

double Datos_KtEt=Datos_KTermofilm/Datos_eTermofilm;
double Datos_htilde=(Datos_hc_superior*Datos_KtEt)/(Datos_hc_superior+Datos_KtEt);

double VinyeccionLtsHr=20   ;//% 20 litros/m^2/hora
double Datos_Vinyeccion=VinyeccionLtsHr/1000/3600; //% m^3/m^2/s


double Datos_kf=0.58  ;//#conductividad fluido agua
double Datos_ks=2     ;//#conductividad solido
double Datos_km_mu=Datos_phi*Datos_kf+(1-Datos_phi)*Datos_ks;

double Datos_km2=0,Datos_DistanciaAlBorde=1;

double Datos_Kevaporacion=1e-3;
double Datos_HumedadAmbiental=0.4;
double Datos_CalorLatenteFluido=2.257e6; // Calor Latente agua: %J/kg , https://es.wikipedia.org/wiki/Calor_latente
double ErrorMax_Cada_t=1e-15;
double TiempoCalculo=0;   // Tiempo actual de calculo 0---> Tmax
double Datos_dt=1000 	      ; //dt para el calculo de evolucion [seg]
double Datos_Tmax=2*Datos_dt      ; //Tmax para el calculo de evolucion [seg]


void formulario_BC_addDatos();

void formulario_BC();
