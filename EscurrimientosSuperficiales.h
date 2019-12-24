
void Etapa_Calcula_Escurrimiento_Superficial()
{
	int i;
	EtapaGlobal=ETAPA_CALCULO_T_ARRIBA;EtapaGlobal2Local=Etapa;


	cout<<"------------------------------------------------------"<<endl;
	cout<<"Calculo Escurrimiento y Temperaturas  (con err0="<<err0<<")"<<endl;
	myfileSalida<<"------------------------------------------------------"<<endl;
	myfileSalida<<"Calculo Escurrimiento y Temperaturas  (con err0="<<err0<<")"<<endl;


	sprintf(text,"%sEtapa:%d:Calculo V,T Arriba (err0=%.2e)\n",text,EtapaGlobal,err0);if (glui != NULL) glui_edittext->set_text(text);



	Guardar=1;
	double minPsi,maxPsi,err,err2,errComb;
	double minTemp,maxTemp;
	cout<<"Interaciones hasta que: err<"<<err0<<endl;
	myfileSalida<<"Interaciones hasta que: err<"<<err0<<endl;

	err=CalculaMpunto(PotencialPsi,Temperatura,UU,VV,WW,*gtotal,gtotal->nH3D,err0);
	err2=CalculaTempSuperficie(Temperatura,PotencialPsi,UU2,VV2,WW2,*gtotal,gtotal->nH3D,err0);
	errComb=max(err,err2);
	U.resize(gtotal->nV3D);V.resize(gtotal->nV3D);W.resize(gtotal->nV3D);	F.resize(gtotal->nV3D);
	U.assign(gtotal->nV3D,0);V.assign(gtotal->nV3D,0);W.assign(gtotal->nV3D,0);F.assign(gtotal->nV3D,0);
	vector<int> cuantos(gtotal->nV3D);
	cuantos.assign(gtotal->nV3D,0);

	CopiaBloquesAVertices(PotencialPsi,F,&minPsi,&maxPsi,1);
	printf("minPotencialVel=%f, maxPotencialVel=%f, errPV=%.2e",minPsi,maxPsi,err);cout<<endl;
	myfileSalida<<"minPotencialVel="<<minPsi<<", maxPotencialVel="<<maxPsi<<", errPV="<<err<<endl;
	sprintf(s,"e=%.0e,FF[%.0e,%.0e]\n",errComb,minPsi,maxPsi);

	sprintf(text,"%serr=%.2e err2=%.2e\n",text,err,err2);if (glui != NULL) glui_edittext->set_text(text);


	CopiaBloquesAVertices(UU,U,&minPsi,&maxPsi,0);
	//		printf("minUU=%e, maxUU=%e, errABS=%g",minF,maxF,err);cout<<endl;
	CopiaBloquesAVertices(VV,V,&minPsi,&maxPsi,0);
	//		printf("minVV=%e, maxVV=%e, errABS=%g",minF,maxF,err);cout<<endl;
	CopiaBloquesAVertices(WW,W,&minPsi,&maxPsi,0);
	//		printf("minWW=%e, maxWW=%e, errABS=%g",minF,maxF,err);cout<<endl;




	//if (glui != NULL) glui_MSG->set_name(s);


	minTemp=maxTemp=Temperatura[0];
	for (i=0;i<gtotal->nH3D;i++) {
		if (minTemp>Temperatura[i]) minTemp=Temperatura[i];
		if (maxTemp<Temperatura[i]) maxTemp=Temperatura[i];
	}

	printf("minTemp=%f\tmaxTemp=%f (int)\terrT=%.2e",minTemp,maxTemp,err2);cout<<endl;
	myfileSalida<<"minTemp="<<minTemp<<"\tmaxTemp="<<maxTemp<<" (int)\terrT="<<err2<<endl;
	sprintf(text,"%sminT=%.2f maxT=%.2f\n",text,minTemp,maxTemp);if (glui != NULL) glui_edittext->set_text(text);

	F2Nodos.resize(gtotal->nV3D);
	CopiaBloquesAVertices(Temperatura,F2Nodos,&minTemp,&maxTemp,2);

	printf("minTemp=%f\tmaxTemp=%f (wBC)\terrT=%.2e",minTemp,maxTemp,err2);cout<<endl;
	myfileSalida<<"minTemp="<<minTemp<<"\tmaxTemp="<<maxTemp<<" (wBC)\terrT="<<err2<<endl;
	sprintf(text,"%sminT=%.2f maxT=%.2f (wBC)\n",text,minTemp,maxTemp);if (glui != NULL) glui_edittext->set_text(text);


	if (gtotal->nH3D<40){
		int j;
		cout <<"xC=[";j=0;
		for (i=0;i<gtotal->nH3D;i++) {
			cout <<gtotal->h3D[i].centro.x<<" ";j++;
			if (j>12) {cout<<"..."<<endl;j=0;}
		} cout <<"];"<<endl;

		cout <<"TC=["; j=0;
		for (i=0;i<gtotal->nH3D;i++) {
			cout <<Temperatura[i]<<" ";j++;
			if (j>12) {cout<<"..."<<endl;j=0;}
		} cout <<"];"<<endl;
	}


	err0=errComb*1e-4;
	if (err0<ErrorMax_Cada_t) err0=ErrorMax_Cada_t;
	if (errComb>ErrorMax_Cada_t)	Etapa=InicioEtapa-1;

	if (DBG) cout<<"Voy a calculaFracionVolumen(Temperatura). Etapa="<<Etapa<<endl;
	calculaFracionVolumen(Temperatura) ;

	if (DBG) cout<<"Voy Volvi calculaFracionVolumen(Temperatura). Etapa="<<Etapa<<endl;
	sprintf(text,"%sFin Etapa\n",text);	if (glui != NULL) glui_edittext->set_text(text);
}


void Etapa_Malla2D_to_Malla3D() {

	int i;
	EtapaGlobal=ETAPA_MALLADO_2D_3D;EtapaGlobal2Local=Etapa;
	//cout<<" Mallando2D-->3D "<<endl;
	myfileSalida<<"Mallando2D-->3D "<<endl;
	InicioEtapa=Etapa;


	if (glui != NULL) {
		if (PanelFlimite != NULL) PanelFlimite->disable();
		MODO_CampoVelocidades=0;
		glui->sync_live();
		Checkbox_particulas->disable();
		PanelParticulas->disable();
	}


	//sprintf(text,"%sCambio malla a dominio inferior\n",text);if (glui != NULL) glui_edittext->set_text(text);

	sprintf(text,"%sCambio CB\n",text);if (glui != NULL) glui_edittext->set_text(text);

	CualesRehacer.resize(0);
	for ( i=0 ; i<gtotal->nCaras ; i++ ) {
		if ( gtotal->Cara[i].nVolumenes == 1 && gtotal->Cara[i].iBC == 2) {
			gtotal->Cara[i].iBC = 1;    //Condicion Neumann homogeneo
			CualesRehacer.push_back(gtotal->Cara[i].ih[0]);
			//				cout<<"CualesRehacer.size()="<<CualesRehacer.size()<<endl;
			//				cout<<"CualesRehacer.end()="<<CualesRehacer.back()<<endl;
		}
	}

	gtotal->CentroCarasBloques();
	sprintf(text,"%s A la proxima: Poligonos2Algunos\n",text);if (glui != NULL) glui_edittext->set_text(text);
}


void Etapa_CompletoPoligonosNuevosNiveles()
{
	int i;
	sprintf(text,"%sGenero Poligonos2Algunos\n",text);if (glui != NULL) glui_edittext->set_text(text);

	cout<<"Genero Poligonos2Algunos: "<<CualesRehacer.size()<<endl;
	myfileSalida<<"Genero Poligonos2Algunos: "<<CualesRehacer.size()<<endl;
	gtotal->generaPoligonos2Algunos(CualesRehacer);



	sprintf(text,"%sCambio malla a dominio inferior\n",text);if (glui != NULL) glui_edittext->set_text(text);

	g2=*gtotal;
	g2.Traslada(0,0,-Dominio_Hsup);
	g2.EscalaZ(Dominio_Hmax/Dominio_Hsup/NDivZ);
	*gtotal=g2;

	sprintf(text,"%sFin Etapa\n",text);

	if (glui != NULL) glui_edittext->set_text(text);

	for (i=1;i<NDivZ;i++) {
		g2.Traslada(0,0,-Dominio_Hmax/NDivZ);
		gtotal->Junta(g2);
	}
	sprintf(text,"%sFin Etapa (listo para calcular?)\n",text);	if (glui != NULL) glui_edittext->set_text(text);
}

void Etapa_Temp3D_Arreglo_CB_3D ()
{
	int i;
	sprintf(text,"%sg.CentroCarasBloques()\n",text);	if (glui != NULL) glui_edittext->set_text(text);

	gtotal->minmax();
	gtotal->CentroCarasBloques();


	sprintf(text,"%sArreglo CB nuevas\n",text);	if (glui != NULL) glui_edittext->set_text(text);

	CualesRehacer.resize(0);
	for ( i=0 ; i<gtotal->nCaras ; i++ ) {
		if ( gtotal->Cara[i].nVolumenes == 2) {
			gtotal->Cara[i].iBC =0 ;  ///nodo interno
		} else {
			double x,y,z;
			x=gtotal->Cara[i].centro.x;
			y=gtotal->Cara[i].centro.y;
			z=gtotal->Cara[i].centro.z;
			gtotal->Cara[i].iBC = 1;    //Condicion Neumann homogeneo
			gtotal->Cara[i].BC  = 0;
			//				if (sqrt(sqr(x)+sqr(y))+1e-5>=Dominio_Rmax*cos(dTheta_med)
			//						|| z+Dominio_Hmax< 1e-5 ) {
			if ( z-1e-5 < -Dominio_Hmax ) { //Abajo
				gtotal->Cara[i].iBC = 11; //Condicion Convectiva

				int iC=i;
				if(1==0) cout<<"g.Cara[iC].BC2="<<gtotal->Cara[iC].BC2
						<<"\t.BC="<<gtotal->Cara[iC].BC
						<<"\tg.Cara[iC].iBC="<<gtotal->Cara[iC].iBC
						<<"\tiC="<<iC
						<<"\tx="<<x
						<<"\ty="<<y
						<<"\tz="<<z
						<<endl;
			}
			if ( x+1e-5 > Dominio_Xmax ) { // Adelante
				gtotal->Cara[i].iBC = 12; //Condicion Convectiva mas espesor

			}
			if ( y+1e-5 > Dominio_Xmax ) { //Derecha
				gtotal->Cara[i].iBC = 12; //Condicion Convectiva mas espesor

			}
			if (z>=-1e-5) {

				CualesRehacer.push_back(gtotal->Cara[i].ih[0]);

				gtotal->Cara[i].iBC = 2;    //Condicion Temperatura conocida
				gtotal->Cara[i].BC  = 0;
				gtotal->Cara[i].BC2  = Temperatura[ gtotal->Cara[i].ih[0] ];

				int iC=i;
				if(1==0) cout<<"g.Cara[iC].BC2="<<gtotal->Cara[iC].BC2
						<<"\t.BC="<<gtotal->Cara[iC].BC
						<<"\tg.Cara[iC].iBC="<<gtotal->Cara[iC].iBC
						<<"\tiC="<<iC
						<<"\tx="<<x
						<<"\ty="<<y
						<<"\tz="<<z
						<<endl;
			}
			int iC=i;
			if(1==0) cout<<"g.Cara[iC].BC2="<<gtotal->Cara[iC].BC2
					<<"\t.BC="<<gtotal->Cara[iC].BC
					<<"\tg.Cara[iC].iBC="<<gtotal->Cara[iC].iBC
					<<"\tiC="<<iC
					<<"\tx="<<x
					<<"\ty="<<y
					<<"\tz="<<z
					<<endl;

		}
	}

}

void Etapa_Temp3D_Calculo_temperaturas_pila_Modelo_Evolucion()
{
	int i;

	EtapaGlobal=ETAPA_CALCULO_T_PILA; EtapaGlobal2Local=Etapa;Guardar=1;

	if (Reiterando==0) {
		cout<<"============================================================"<<endl;
		myfileSalida<<"============================================================"<<endl;
		if (TipoCalculo==CalculoEvolucion) {
			double tmin=TiempoCalculo/60;
			double thoras=tmin/60;
			cout<<"Calculo para TiempoCalculo= "<<TiempoCalculo<<" seg."<<" = "<<tmin<<" min."<<" = "<<thoras<<" hrs."<<endl;
			myfileSalida<<"Calculo para TiempoCalculo= "<<TiempoCalculo<<" seg."<<" = "<<tmin<<" min."<<" = "<<thoras<<" hrs."<<endl;
		} else {
			cout<<"Calculo Estacionario:"<<endl;
			myfileSalida<<"Calculo Estacionario:"<<endl;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////
	//   Calculo de temperaturas en la pila : Modelo de Evolucion
	//
	//
	//
	//cout<<"Etapa="<<Etapa<<"IniciEtapa="<<InicioEtapa <<"EtapaGlobal2Local="<<EtapaGlobal2Local<<endl;

	cout<<"------------------------------------------------------"<<endl;
	myfileSalida<<"------------------------------------------------------"<<endl;

	printf("Etapa: Calculo T° en la Pila, con err0=%e",err0);cout<<endl;
	myfileSalida<<"Etapa: Calculo T° en la Pila, con err0="<<err0<<endl;

	if (TipoCalculo==CalculoEvolucion) {
		sprintf(text,"%sCalculo T en Pila (t=%.0f err0=%.2e)\n",text,TiempoCalculo,err0);if (glui != NULL) glui_edittext->set_text(text);
	} else {
		sprintf(text,"%sCalculo T en Pila Estacionario (err0=%.2e)\n",text,err0);if (glui != NULL) glui_edittext->set_text(text);

	}
	if (PanelFlimite != NULL)  PanelFlimite->enable();
	double err2,minTemp,maxTemp;

	printf("numero de celdas=%d",gtotal->nH3D);cout<<endl;
	myfileSalida<<"numero de celdas="<<gtotal->nH3D<<endl;
	PotencialVInfiltracion.resize(gtotal->nH3D);
	TempPilaBloques.resize(gtotal->nH3D);
	UU2.resize(gtotal->nH3D);
	VV2.resize(gtotal->nH3D);
	WW2.resize(gtotal->nH3D);
	TempPilaVertices.resize(gtotal->nV3D);
	for (i=0;i<gtotal->nH3D;i++) {
		PotencialVInfiltracion[i]=gtotal->h3D[i].centro.z*Datos_Vinyeccion;
	}

	//// Esta funcion hace el calculo de TempPilaBloques hasta que el error sea < err0
	err2=CalculaTemperaturaPilaEnTmasDt(TempPilaBloques,PotencialVInfiltracion,UU2,VV2,WW2,
			*gtotal,gtotal->nH3D,err0,TempPilaBloquesPrevia);


	F.resize(gtotal->nV3D);
	F2Nodos.resize(gtotal->nV3D);
	CopiaBloquesAVertices(TempPilaBloques,F2Nodos,&minTemp,&maxTemp,1);


	minTemp=maxTemp=TempPilaBloques[0];
	for (i=0;i<gtotal->nH3D;i++) {
		if (minTemp>TempPilaBloques[i]) minTemp=TempPilaBloques[i];
		if (maxTemp<TempPilaBloques[i]) maxTemp=TempPilaBloques[i];
	}
	//		printf("min,max TempPilaBloques=%f, %f, errABS=%g (%d)",minTemp,maxTemp,err2,gtotal->nH3D);cout<<endl;

	//		sprintf(text,"%sminT=%.2f maxT=%.2f\n errS=%.2e\n",text,minTemp,maxTemp,err2);if (glui != NULL) glui_edittext->set_text(text);

	printf("minTemp=%f\tmaxTemp=%f (int)\terrT=%.2e",minTemp,maxTemp,err2);cout<<endl;
	myfileSalida<<"minTemp="<<minTemp<<"\tmaxTemp="<<maxTemp<<" (int)\terrT="<<err2<<endl;

	sprintf(text,"%s(errT=%.2e)\n",text,err2);if (glui != NULL) glui_edittext->set_text(text);

	sprintf(text,"%sminT=%.2f maxT=%.2f (int)\n",text,minTemp,maxTemp);if (glui != NULL) glui_edittext->set_text(text);


	CopiaBloquesAVertices(TempPilaBloques,F,&minTemp,&maxTemp,2);
	CopiaBloquesAVertices(TempPilaBloques,F2Nodos,&minTemp,&maxTemp,2);



	printf("minTemp=%f\tmaxTemp=%f (wBC)\terrT=%.2e",minTemp,maxTemp,err2);cout<<endl;
	myfileSalida<<"minTemp="<<minTemp<<"\tmaxTemp="<<maxTemp<<" (wBC)\terrT="<<err2<<endl;

	sprintf(text,"%sminT=%.2f maxT=%.2f (wBC)\n",text,minTemp,maxTemp);if (glui != NULL) glui_edittext->set_text(text);


	calculaFracionVolumen(TempPilaBloques) ;

	err0=err2*1e-4;
	if (err0<ErrorMax_Cada_t) err0=ErrorMax_Cada_t;

	Reiterando=0;
	if (err2>ErrorMax_Cada_t) {
		cout<<"Reitero..."<<endl;
		myfileSalida<<"Reitero..."<<endl;
		Reiterando=1;
		Etapa=InicioEtapa-1;
	} else {
		// Termine con el calculo para este Tiempo
		if (TipoCalculo==CalculoEvolucion) {
			cout<<"Termine con el calculo para este Tiempo"<<endl;
			myfileSalida<<"Termine con el calculo para este Tiempo"<<endl;

			//cout<<"CalculoContinuo="<<CalculoContinuo<<endl;
			if (CalculoContinuo) {
				GuardarInstantanea();
				Guardar=0;
			}

			TiempoCalculo += Datos_dt ;
			if (TiempoCalculo<=Datos_Tmax) {
				//Paso al Tiempo siguiente....
				//err0=1e-10;
				for (i=0;i<gtotal->nH3D;i++) {
					TempPilaBloquesPrevia[i]=TempPilaBloques[i];
				}
				cout<<"Paso al T= "<<TiempoCalculo<<" seg"<<endl;
				myfileSalida<<"Paso al T= "<<TiempoCalculo<<" seg"<<endl;
				Etapa=InicioEtapa-1;

			} else {
				TipoCalculo=CalculoEstacionario;

				cout<<"Terminé evolucion: Comienzo caso Estacionario..."<<endl;
				myfileSalida<<"Terminé evolucion: Comienzo caso Estacionario..."<<endl;
				err0=1e-10;
				Etapa=InicioEtapa-1;
			}
		}
	}

	sprintf(text,"%sFin Etapa\n",text);

	if (glui != NULL) glui_edittext->set_text(text);

}


void calculaFracionVolumen(vector<double> &Temp) {
	int i;
	//	cout<<"calculaFracionVolumen()"<<endl;
	double FraccionVolumen=0,VolumenTotal=0;
	for (i=0;i<gtotal->h3D.size();i++) {
		VolumenTotal += gtotal->h3D[i].volumen;
		if (Temp[i]>TLimite) FraccionVolumen += gtotal->h3D[i].volumen;

	}
	//cout<<"calculaFracionVolumen()2"<<endl;
	for (i=0;i<gtotal->Cara.size();i++) {
		if (gtotal->Cara[i].iBC ==2 || gtotal->Cara[i].iBC ==3 ) {
			VolumenTotal += gtotal->Cara[i].volumen;
			if (gtotal->Cara[i].BC2>TLimite) FraccionVolumen += gtotal->Cara[i].volumen;
		}
	}

	//cout<<"Volumen con T>"<<TLimite<<" = "<<FraccionVolumen<<"\tVolumenTotal"<<VolumenTotal<<"\t(nH3D="<<gtotal->nH3D<<")"<<endl;
	FraccionVolumen /= VolumenTotal;

	if (TipoCalculo==CalculoEvolucion) {
		cout<<"Fraccion Volumen con T>"<<TLimite<<" = "<<FraccionVolumen
				<<"\tT = "<<TiempoCalculo<<endl;
		myfileSalida<<"Fraccion Volumen con T>"<<TLimite<<" = "<<FraccionVolumen
				<<"\tT = "<<TiempoCalculo<<endl;
		myfileVol<<TiempoCalculo<<" "<<FraccionVolumen<<endl;
	} else {
		cout<<"Fraccion Volumen con T>"<<TLimite<<" = "<<FraccionVolumen<<endl;
		myfileSalida<<"Fraccion Volumen con T>"<<TLimite<<" = "<<FraccionVolumen<<endl;
	}

	char s[100];
	sprintf(s,"Fraccion=%f\n",FraccionVolumen);if (glui != NULL) glui_porcentaje->set_name(s);

	//	cout<<"calculaFracionVolumen():END"<<endl;
}
