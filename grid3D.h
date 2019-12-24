//#pragma once


#include "Class_Vector.h"
#include<vector>
#include "xmlParser.h"

#define sqr(x) ((x)*(x))

int const maxpasadas=20;

class grid3D;
class R3;

class MatrizSym3x3 {
public:
	double a,b,c,d,e,f;


	//-----
	MatrizSym3x3();
	MatrizSym3x3(const MatrizSym3x3 &M);
	MatrizSym3x3(double alfa);
	MatrizSym3x3 inversa ( );
    const MatrizSym3x3& operator= ( const MatrizSym3x3 M2);
    const R3 operator*  ( const R3& L);
    friend ostream & operator << (ostream & os, const MatrizSym3x3 & M2)
    {
    	os <<"[ "<<M2.a<< " "<<M2.b<< " "<<M2.c<< ""<<endl;
    	os <<"  "<<M2.b<< " "<<M2.d<< " "<<M2.e<< ""<<endl;
    	os <<"  "<<M2.c<< " "<<M2.e<< " "<<M2.f<< "];"<<endl;
    	return os;
    }

};



class R3{ //Puntos,vectores y Doblepuntos de \R^3
public:
	double x,y,z,L;//L se usa para trazos que unen dos puntos en los poligonos

//-----
	 R3();
	R3(double a);
	void save(ofstream &myfile);
	void read(ifstream &myfile);
	void Traslada(double dx,double dy,double dz);
	void EscalaZ(double lambda);
	void NormaUnitario();
	double Norma();
	double NormaSQR();
    const R3& operator= ( const R3& v);
    R3 operator +=(const R3& v1);

    friend R3 operator*(const R3& a,const double k)
    { R3 c; c.x=a.x*k;c.y=a.y*k;c.z=a.z*k;c.L=a.L*k;; return c; }
    friend R3 operator/(const R3& a,const double k)
    { R3 c; c.x=a.x/k;c.y=a.y/k;c.z=a.z/k;c.L=a.L/k;; return c; }
    friend R3 operator-(const R3& a,const R3& b)
        { R3 c; c.x=a.x-b.x;c.y=a.y-b.y;c.z=a.z-b.z;c.L=sqrt(sqr(c.x)+sqr(c.y)+sqr(c.z)); return c; }
    friend R3 operator+(const R3& a,const R3& b)
        { R3 c; c.x=a.x+b.x;c.y=a.y+b.y;c.z=a.z+b.z;c.L=sqrt(sqr(c.x)+sqr(c.y)+sqr(c.z)); return c; }





    friend ostream & operator << (ostream & os, const R3 & v)
    {
    	os << "[" <<v.x<<" , "<<v.y<<" , "<<v.z<<"], L="<<v.L<<endl;
    	return os;
    }

};

double ppunto(R3 a,R3 b);

double ppuntodiff(R3 a,R3 b,R3 c);

class PoligonoPlano {
public:
	double AreaPoligono,Dab;
	vector<R3> punto;
	int nPuntos;
	R3 normal;
	R3 centro;
	int CentroCalculado=0;
	grid3D *papa;
	int no,iBC;
	double nRandom;
	int signo2;
	double winX1,winY1,winX2,winY2,winX3,winY3,winZ1,winZ2,winZ3;

//------
	void save(ofstream &myfile);
	void read(ifstream &myfile);

	void drawGL();

	void draw_caraGL();

	void readArchivoPoly(ifstream &myfile,grid3D *papaL);
	void Traslada(double dx,double dy,double dz);
	void EscalaZ(double lambda);

};

class Vertex3D {
public:
	double x,y,z,z0;
	int no;
	grid3D *papa;
	R3 normalVetex;
	vector<int> soporte;

	//no se guarda
	double vz;
	int nvz;

//----
	void save(ofstream &myfile);
	void read(ifstream &myfile,grid3D *papaL);
	void Traslada(double dx,double dy,double dz);
	void EscalaZ(double lambda);
	friend ostream & operator << (ostream & os, const Vertex3D & v)
    {
    	os << "[" <<v.x<<" , "<<v.y<<" , "<<v.z<<"], z0="<<v.z0<<endl;
    	return os;
    }


};

class Cara3D {
public:
	grid3D *papa;
	int no,nvCara,nVolumenes,iBC;
	int iv[4],ih[2];
	double BC,BC2,volumen;
	R3 centro;
	R3 normalCara;
	//No se guardan....
	int signo2;
	double winX1,winY1,winX2,winY2,winX3,winY3,winZ1,winZ2,winZ3;
	double nRandom;

//----
	void save(ofstream &myfile);
	void inicializa(grid3D *pp);
	void read(ifstream &myfile,grid3D *papaL);
	void Traslada(double dx,double dy,double dz);
	void EscalaZ(double lambda);

	void drawGL();      //Impime numeror, centro, etc

	void draw_caraGL00(vector<double> &F,double minF,double maxF,int* ii);
	void draw_caraGL();  //Dibuja el poligono ==cara

};


class Tetrahedros {
	// Tetrahedros:
	//
	//  3
	//	|\ \
	//	| \  0
	//	|  / \
	//	| / \ \
	//	|/   \ \
	//  1------2

public:
//	vector<PoligonoPlano> Poligono;		// El hexahedro [i] con el g.h3D[i].vecino[j] se unen por g.h3D[i].Poligono[j]
//	vector<int> vecino;					// Lista de Hexahedros vecinos (O Cara en el borde)
//	vector<int> tipo_vecino;			// Lista del mismo largo, == ES_BOQUE o ES_CARA o ES_CARA_L2
	// Cuando g.h3D[i].tipo_vecino[j] == ES_CARA  ==> iC=g.h3D[i].vecino[j] es una Cara verdadera: g.Cara[iC]
	vector<int> dibujado;

	int iv[4];
	int icara[4];
	int no;
	double volumen;
	R3 centro;
	grid3D *papa;

//-----
//	void save(ofstream &myfile);
//	void read(ifstream &myfile,grid3D *papaL);
//	void Traslada(double dx,double dy,double dz);
//	void EscalaZ(double lambda);

//	void draw_caraGL(int ,int,int,int);
//	void draw_caraGL(vector<double> F,double,double,int ,int,int,int);
//	void draw_edgeGL(int i0,int i1);
//	void CalculaVolumen();
};



class Hexa3D {
	// Hexaedros:
	//	  2------3
	//   /|     /|
	//  6-|----7 |
	//	| |    | |
	//	| 0----|-1
	//	|/     |/
	//  4------5

public:
	vector<PoligonoPlano> Poligono;		// El hexahedro [i] con el g.h3D[i].vecino[j] se unen por g.h3D[i].Poligono[j]
	vector<int> vecino;					// Lista de Hexahedros vecinos (O Cara en el borde)
	vector<int> tipo_vecino;			// Lista del mismo largo, == ES_BOQUE o ES_CARA o ES_CARA_L2
	// Cuando g.h3D[i].tipo_vecino[j] == ES_CARA  ==> iC=g.h3D[i].vecino[j] es una Cara verdadera: g.Cara[iC]
	vector<int> dibujado;

	int iv[8];
	int icara[6];
	int no;
	double volumen;
	R3 centro;
	grid3D *papa;

//-----
	void save(ofstream &myfile);
	void read(ifstream &myfile,grid3D *papaL);
	void Traslada(double dx,double dy,double dz);
	void EscalaZ(double lambda);

	void draw_caraGL(int ,int,int,int);
	void draw_caraGL(vector<double> F,double,double,int ,int,int,int);
	void draw_edgeGL(int i0,int i1);
	void CalculaVolumen();
};


class VolumenFinito {
	// Hexaedros:
	//	  2------3
	//   /|     /|
	//  6-|----7 |
	//	| |    | |
	//	| 0----|-1
	//	|/     |/
	//  4------5

public:
	vector<PoligonoPlano> Poligono;		// El hexahedro [i] con el g.h3D[i].vecino[j] se unen por g.h3D[i].Poligono[j]
	vector<int> vecino;					// Lista de Hexahedros vecinos (O Cara en el borde)
	vector<int> tipo_vecino;			// Lista del mismo largo, == ES_BOQUE o ES_CARA o ES_CARA_L2
										// Cuando g.h3D[i].tipo_vecino[j] == ES_CARA  ==> iC=g.h3D[i].vecino[j] es una Cara verdadera: g.Cara[iC]
	vector<R3> vecino_centro;
	vector<R3> vecino_centro_caraoriginal;
	vector<R3> vecino_normal;
	vector<int> dibujado;

	int no;
	double volumen,nRandom;
	R3 centro;
	grid3D *papa;

//-----
	void save(ofstream &myfile);
	void read(ifstream &myfile,grid3D *papaL);
	void Traslada(double dx,double dy,double dz);
	void EscalaZ(double lambda);

	void draw_caraGL(int ,int,int,int);
	void draw_caraGL(vector<double> F,double,double,int ,int,int,int);
	void draw_edgeGL(int i0,int i1);
	void CalculaVolumen();
};


class TriPrisma {
	// Prisma de base Triangular:
	//    5
	//   /|\
	//  / | \
	// 3-----4
	// |  |  |
	// |  2  |
	// | / \ |
	// |/   \|
	// 0-----1
public:

	grid3D *papa;
	int no;
	int iv[6];
	double volumen;
	int icara[5];
	R3 centro;
	int CentroDibujado;

	vector<PoligonoPlano> Poligono;		// El TriPrisma [i] con el g.h3D[i].vecino[j] se unen por g.h3D[i].Poligono[j]
	vector<int> vecino;					// Lista de TriPrismas vecinos (O Cara en el borde)
	vector<int> tipo_vecino;			// Lista del mismo largo, == ES_BOQUE o ES_CARA o ES_CARA_L2
	// Cuando g.h3D[i].tipo_vecino[j] == ES_CARA  ==> iC=g.h3D[i].vecino[j] es una Cara verdadera: g.Cara[iC]
	vector<int> dibujado;

	
	
//-----
	void save(ofstream &myfile);
	void readMSH(ifstream &myfile,grid3D *papaL);
	void read(ifstream &myfile,grid3D *papaL);
	void Traslada(double dx,double dy,double dz);
	void EscalaZ(double lambda);

	void draw_caraGL(int ,int,int,int);
	void draw_caraGL(vector<double> F,double,double,int ,int,int,int);
	void draw_edgeGL(int i0,int i1);
	void CalculaVolumen();
	void DrawCentro();
};


class grid3D
{
public:
	int nH3D=0,nV3D=0,nCaras=0,nPoligonos=0,QuienGeneraPoligonos;
	double xmin,xmax,ymin,ymax,zmin,zmax;

	vector<Vertex3D> v3D;
	int nBC=0;vector<int> BC,DrawBC,lBC;
	vector<Hexa3D>   h3D;
	vector<Cara3D> Cara;
	vector<PoligonoPlano> Poligono;		// El hexahedro [i] con el g.h3D[i].vecino[j] se unen por g.h3D[i].Poligono[j]

	
	//Version >= 3
	int nTriPrisma3D=0,TriPri3DAnalizados=0,TetrahedrosAnalizados=0;
	vector<TriPrisma>   TriPrisma3D;

	//Version >= 3
	int nTetrahedros=0,nTetrahedrosAnalizados=0;
	vector<Tetrahedros>   Tetrahedro;

	int nVolFinito=0;
	vector<VolumenFinito>   VolFinito;
	vector<int>   VolFinitoSelected;
	vector<int>   Selected_Triprisma;
	int imprimir=0;
	

//----------------------

	void print_poligonos();
	void CalculaNormalVertex();
	void write(ofstream &myfile);
	void read(ifstream &myfile);
	void write(char *);
	void read(char *);



	void readPOLY(char *ifile_name);


	void readFACE(char *ifile_name);
	void readMallaXML(char *ifile_name);



	void readMSH3D(char *ifile_name);

	void minmax();

	void drawGL();
	void drawGL(vector<double> & F);

	void drawGL(vector<double> & F,double minF,double maxF );
	void drawVoronoi();
	void Voronoi_Draw_i(int iv);
	void TriPrisma_Draw_i(int iv);




	void DibujaParticula (int i,R3 & origen);
	void drawParticulas_TriPrisma();

	int EnCualTriprisma(float px,float py,float pz,int tri0);
	int EnCualTriprismaXY(float px,float py,int tri0);
	void drawVelGL(vector<double>,vector<double>,vector<double>);
	void drawVelGL2();

	void GeneraCaras(int inicia=false);	

	void DZmin(double DzMin);
	void GeneraCarasTriPri();
	void GeneraCarasTetrahedros();
	int BuscaCara3(int i0,int i1,int i2) ;

	int addBC(int iBC);
	int AddCara3(int ib,int i0,int i1,int i2,int iBC);
	int BuscaCara4(int i0,int i1,int i2,int i3) ;
	int AddCara4(int ib,int i0,int i1,int i2,int i3);
	void draw_caraGL(int ii[4]);

	void draw_caraGL00(vector<double>&F,double minF,double maxF,int ii[4]);
	void Poligonos_Generar_Version3();
	void generaPoligonos2Algunos(vector<int> &CualesRehacer);
	void Poligono_Inicial(R3 a, R3 b, PoligonoPlano &P);
	void Poligono_InicialCara(R3 a, R3 b, PoligonoPlano &P,int iC);

	void CentroCarasBloques();

	void Soporte_Vertices();
	void CalculaVolumen();
	void Particulas_Init(int i);
	void Particulas_Init_3(double posX, double posY, double posZ);


	void SeleccionaVolFinito(double px,double py,double pz);
	void SeleccionaYMuestraCara(double px,double py,double pz);

	void SeleccionaYMuestraPolygono(double px,double py,double pz);



	void cubo(int,int,int,float=1,float=1,float=1);
	void Junta(grid3D g1,grid3D g2);
	void Junta(grid3D g2);
	void Rota90Z();
	void Traslada(double dx,double dy,double dz);
	void EscalaZ(double lambda);
	grid3D(void);
	~grid3D(void);

};
#define ES_BLOQUE 1
#define ES_CARA   2
#define ES_CARA_L2   3
