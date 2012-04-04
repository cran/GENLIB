/*! \file interface.cc
\brief Implementation toutes les fonctions d'interface de Splus

	C'est dans ce fichier que sont placer toutes les interfaces pour Splus.
	Seulement les fonctions d'ici sont caller à partir des wrappers des .ssc

\author Sébastien Leclerc
*/


#define ALLOWPRINTPROGRESS
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h> 
#include <Rcpp.h>
#include <Rcpp/as.h>
#include <RcppCommon.h>
#include <string>
#include"base.h"
#include"userInterface.h"

#include"probsimul.h"
#include"apparentement.h"
#include"congen.h"
#include"consanguinite.h"
#include"fondateur.h"
#include"outils.h"
#include"genphi.h"
#include"statanal.h"

#include "interface.h"
//#define EXPORTTYPE extern "C"  -> remplacer par RcppExport


/// Valeur textuel de l'enumeration typenoeud_t
/** \sa typenoeud_t*/
const char *stype[]=
{
	"NON-EXPLORER",
	"INUTILE",
	"NOEUD",
	"DEPART",
	"PROPOSANT",
	"GENPROPOSANTINUTILE"
};

//********************************************
//     FONCTION DE CONTROLE
//******************************************** 

RcppExport void SPLUSFlushCacheGenealogie()
{
	//Flush la cache de la genealogie et autre
	FlushGenealogie();	
	return;		
}

RcppExport void SPLUSGetTimer(SEXP sTimeInSec)
{
	//Flush la cache de la genealogie et autre
	int * TimeInSec;
	TimeInSec = INTEGER_POINTER( sTimeInSec );
	*TimeInSec = getLastTimer();
	return;		
}

RcppExport SEXP SPLUSValidateGenealogie(SEXP RGenealogie, SEXP RisValid)
{
	STARTTIMER;
	int * Genealogie, * tmp;
	Rcpp::IntegerVector gen(RGenealogie);
	Genealogie = INTEGER_POINTER(gen);
	tmp = INTEGER_POINTER(RisValid);

	*tmp = ValidateGenealogie(Genealogie);
	//for(int i=0; i<gen.size(); i++) gen[i] = Genealogie[i];
	STOPTIMER;
     return (Rcpp::List::create( Rcpp::Named("Data")    = gen, //Rcpp::wrap(tmp1),
						   Rcpp::Named("isValid") = RisValid ) );
}

/// Fonction d'interface Splus pour change le temps maximum des fonctions longue
/** \sa setCurrentMaxTime() getCurrentMaxTime() */
RcppExport void SPLUSChangeMaxProcessingTime(SEXP snewMaximum,SEXP soldMaximum)
{	
	double * newMaximum, * oldMaximum;
	newMaximum = NUMERIC_POINTER(snewMaximum);
	oldMaximum = NUMERIC_POINTER(soldMaximum);
	getCurrentMaxTime(oldMaximum);
	if (*newMaximum>=0)
		setCurrentMaxTime(*newMaximum);
	return;	
}

/// **********
//	APPARENTEMENT
// *********
#define USEALTERNATEEXCEPTION

/// Fonction d'interface Splus pour PhiMatrix
/** \sa PhiMatrix()*/
RcppExport void SPLUSPhiMatrix(SEXP sGenealogie, SEXP sProposant, SEXP sNProposant, SEXP sNiveau, SEXP sPDRetour, SEXP sPrintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * Niveau, * printit ;
	double * pdRetour;
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sProposant);

	Genealogie = INTEGER_POINTER( lGenealogie );
	proposant  = INTEGER_POINTER( lproposant );
	
	NProposant = INTEGER_POINTER(sNProposant);
	Niveau	 = INTEGER_POINTER(sNiveau);
	pdRetour	 = NUMERIC_POINTER(sPDRetour);
	printit	 = INTEGER_POINTER(sPrintit);
	
	PhiMatrix(Genealogie, proposant,*NProposant,*Niveau, pdRetour,*printit);
	STOPTIMER;
	return;
}

/// Fonction d'interface Splus pour PhiMatrix
/** \sa PhiMatrix()*/
RcppExport void SPLUSPhiMatrixMT(SEXP sGenealogie, SEXP sproposant, SEXP sNProposant, SEXP sNiveau, SEXP spdRetour, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * Niveau, * printit;
	double * pdRetour;
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sproposant);
	
	Genealogie = INTEGER_POINTER(lGenealogie);
	proposant  = INTEGER_POINTER(lproposant);
	NProposant = INTEGER_POINTER(sNProposant);
	Niveau	 = INTEGER_POINTER(sNiveau);
	pdRetour	 = NUMERIC_POINTER(spdRetour);
	printit	 = INTEGER_POINTER(sprintit);
	
	PhiMatrixMT(Genealogie, proposant, *NProposant, *Niveau, pdRetour,*printit);
	STOPTIMER;
	return;
}

/// Fonction d'interface Splus pour Phis2
/** \sa Phis2() Phis()*/
RcppExport void SPLUSPhis(SEXP sGenealogie, SEXP sproposant, SEXP sNProposant,SEXP sNiveauMin, SEXP sNiveauMax, SEXP spdRetour, 
					 SEXP sMatrixArray, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * NiveauMin, * NiveauMax, * printit;
	double * pdRetour, * MatrixArray;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sproposant);
	
	Genealogie  = INTEGER_POINTER(lGenealogie);
	proposant   = INTEGER_POINTER(lproposant);
	NProposant  = INTEGER_POINTER(sNProposant);
	NiveauMin	  = INTEGER_POINTER(sNiveauMin);
	NiveauMax	  = INTEGER_POINTER(sNiveauMax);
	pdRetour	  = NUMERIC_POINTER(spdRetour);
	MatrixArray = NUMERIC_POINTER(sMatrixArray);
	printit	 = INTEGER_POINTER(sprintit);
	
	//Peut utilise Phis en guise de comparaison  
	Phis(Genealogie, proposant,*NProposant,*NiveauMin,*NiveauMax,pdRetour,MatrixArray,*printit);
	STOPTIMER;

	return;
}

RcppExport void SPLUSPhisMT(	SEXP sGenealogie , SEXP sproposant, SEXP sNProposant, SEXP sNiveauMin, SEXP sNiveauMax, SEXP spdRetour, 
						SEXP sMatrixArray, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * NiveauMin, * NiveauMax, * printit;
	double * pdRetour, * MatrixArray;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sproposant);
	
	Genealogie  = INTEGER_POINTER(lGenealogie);
	proposant   = INTEGER_POINTER(lproposant);
	NProposant  = INTEGER_POINTER(sNProposant);
	NiveauMin	  = INTEGER_POINTER(sNiveauMin);
	NiveauMax	  = INTEGER_POINTER(sNiveauMax);
	pdRetour	  = NUMERIC_POINTER(spdRetour);
	MatrixArray = NUMERIC_POINTER(sMatrixArray);
	printit	 = INTEGER_POINTER(sprintit);
		//Peut utilise Phis en guise de comparaison  
		PhisMT(Genealogie,
		proposant,*NProposant,*NiveauMin,*NiveauMax,
		pdRetour,MatrixArray,*printit);
	STOPTIMER;

	return;
}


/// **********
//	CONSANGUINITE
// *********

/// Fonction d'interface Splus pour consan (F)
/** \sa consan()*/
RcppExport void SPLUSF(SEXP sGenealogie, SEXP sproposant, SEXP sNProposant, SEXP sNiveau, SEXP spdRetour, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * printit;
	double * pdRetour, * Niveau;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sproposant);
	Rcpp::NumericVector lpdRetour(spdRetour);
	
	Genealogie = INTEGER_POINTER(lGenealogie);
	proposant  = INTEGER_POINTER(lproposant);
	NProposant = INTEGER_POINTER(sNProposant);
	Niveau	 = NUMERIC_POINTER(sNiveau);
	pdRetour	 = NUMERIC_POINTER(lpdRetour);
	printit	 = INTEGER_POINTER(sprintit);
	
	consan(Genealogie, proposant,*NProposant,*Niveau, pdRetour,*printit);
	STOPTIMER;
	return;
}

/// Fonction d'interface Splus pour consanFs 
/** \sa consanFs()*/
RcppExport void SPLUSFS(SEXP sGenealogie, SEXP sproposant, SEXP sNProposant, SEXP sNiveauMin, SEXP sNiveauMax, SEXP spdRetour, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * printit;
	double * pdRetour, * NiveauMin, * NiveauMax;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sproposant);
	Rcpp::NumericVector lpdRetour(spdRetour);
	
	Genealogie  = INTEGER_POINTER(lGenealogie);
	proposant   = INTEGER_POINTER(lproposant);
	NProposant  = INTEGER_POINTER(sNProposant);
	NiveauMin	  = NUMERIC_POINTER(sNiveauMin);
	NiveauMax	  = NUMERIC_POINTER(sNiveauMax);
	pdRetour	  = NUMERIC_POINTER(lpdRetour);
	printit	  = INTEGER_POINTER(sprintit);
	
	consanFs(Genealogie, proposant,*NProposant,*NiveauMin,*NiveauMax, pdRetour,*printit);
	STOPTIMER;	
	return;
}

/// **********
//	OUTILS
// *********

/// Fonction d'interface Splus pour CountChild
/** \sa CountChild()*/
RcppExport void SPLUSChild(SEXP sGenealogie, SEXP splProposant,SEXP slNProposant, SEXP sretour)
{
	STARTTIMER;
	int * Genealogie, * plProposant, * lNProposant, * retour;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lplProposant(splProposant);
	
	Genealogie  = INTEGER_POINTER(lGenealogie);
	plProposant = INTEGER_POINTER(lplProposant);
	lNProposant = INTEGER_POINTER(slNProposant);
	retour	  = INTEGER_POINTER(sretour);

	CountChild(Genealogie, plProposant, *lNProposant, retour);
	STOPTIMER;
	return;
} 

/*
RcppExport  void  SPLUSTestEbranche(SEXP sGenealogie, SEXP sProposant, SEXP sNProposant, SEXP sAncetre, SEXP sNAncetre, SEXP sRetour, SEXP sTaille)
{
 Rcpp::IntegerVector lGenealogie(sGenealogie);
 Rcpp::IntegerVector lproposant(sProposant);
 Rcpp::IntegerVector lancetre(sAncetre);
 Rcpp::IntegerVector lretour(sRetour);
 
 int * Genealogie, * proposant, * ancetre, * retour, * nproposant, * nancetre, * taille ;
 Genealogie = INTEGER_POINTER( lGenealogie );
 proposant  = INTEGER_POINTER( lproposant );
 ancetre	  = INTEGER_POINTER( lancetre );
 retour	  = INTEGER_POINTER( lretour );
 
 //int lNProposant = Rcpp::as<int>(sNProposant);
 //int lNAncetre   = Rcpp::as<int>(sNAncetre);
 //int lTaille     = Rcpp::as<int>(sTaille);
 //nproposant	= &lNProposant;
 //nancetre		= &lNAncetre;
 nproposant	= INTEGER_POINTER(sNProposant);
 nancetre		= INTEGER_POINTER(sNAncetre);
 taille		= INTEGER_POINTER(sTaille); 		//&lTaille;

 int res = testEbranche(Genealogie, proposant, *nproposant, ancetre, *nancetre, retour, taille);
 return (Rcpp::List::create(	  Rcpp::Named("Data")   = lGenealogie
						, Rcpp::Named("prop")   = lproposant
						, Rcpp::Named("nProp")  = Rcpp::wrap(*nproposant)
						, Rcpp::Named("anc")    = lancetre
						, Rcpp::Named("nAnc")   = Rcpp::wrap(*nancetre)
						, Rcpp::Named("retour") = lretour
						, Rcpp::Named("taille") = Rcpp::wrap(*taille)
						));//
 return ;
}*/

/// Fonction d'interface Splus pour ebranche
/** \sa ebranche()*/
RcppExport  void  SPLUSebranche(SEXP sGenealogie, SEXP sProposant, SEXP sNProposant, SEXP sAncetre, SEXP sNAncetre, SEXP sRetour, SEXP sTaille)
{
	STARTTIMER;
	int * Genealogie, * proposant, * ancetre, * retour, * nproposant, * nancetre, * taille ;

	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sProposant);
	Rcpp::IntegerVector lancetre(sAncetre);
	Rcpp::IntegerVector lretour(sRetour);

	Genealogie = INTEGER_POINTER( lGenealogie );
	proposant  = INTEGER_POINTER( lproposant );
	ancetre	 = INTEGER_POINTER( lancetre );
	retour	 = INTEGER_POINTER( lretour );
	
	//int lNProposant = Rcpp::as<int>(sNProposant);
	//int lNAncetre   = Rcpp::as<int>(sNAncetre);
	//int lTaille     = Rcpp::as<int>(sTaille);
	//nproposant      = &lNProposant;
	//nancetre  	 = &lNAncetre;
	//taille		 = &lTaille;
	nproposant	= INTEGER_POINTER(sNProposant);
	nancetre	= INTEGER_POINTER(sNAncetre);
	taille		= INTEGER_POINTER(sTaille);
	
	ebranche(Genealogie, proposant, *nproposant, ancetre, *nancetre, retour, taille);
	STOPTIMER;
	return ;
}


/// Fonction d'interface Splus pour numeroGen
/** \sa compareGen()*/
RcppExport  void  SPLUSnumeroGen(SEXP sGenealogie, SEXP splProposant, SEXP sNProposant, SEXP sretour)		   
{

	STARTTIMER;
	int * Genealogie, * plProposant, * NProposant, * retour;
	Rcpp::IntegerVector lGenealogie (sGenealogie);
	Rcpp::IntegerVector lplProposant(splProposant);
	Rcpp::IntegerVector lretour     (sretour);

	Genealogie  = INTEGER_POINTER( lGenealogie );
	plProposant = INTEGER_POINTER( splProposant );
	retour 	  = INTEGER_POINTER( sretour );
	NProposant  = INTEGER_POINTER(sNProposant);
	
	numeroGen(Genealogie, plProposant,*NProposant, retour);
	STOPTIMER;
	return ;
}
/// Fonction d'interface Splus pour numeroGenMin
/** \sa compareGen()*/
RcppExport  void  SPLUSnumGenMin(SEXP sGenealogie, SEXP splProposant,SEXP sNProposant, SEXP sretour)		   
{

	STARTTIMER;
	int * Genealogie, * plProposant, * NProposant, * retour;
	Rcpp::IntegerVector lGenealogie (sGenealogie);
	Rcpp::IntegerVector lplProposant(splProposant);
	Rcpp::IntegerVector lretour     (sretour);

	Genealogie  = INTEGER_POINTER( lGenealogie );
	plProposant = INTEGER_POINTER( splProposant );
	retour 	  = INTEGER_POINTER( sretour );
	NProposant  = INTEGER_POINTER(sNProposant);
	
		numeroGenMin(Genealogie, plProposant,*NProposant, retour);
	STOPTIMER;

	return ;
}

/// Fonction d'interface Splus pour numeroGenMoy
/** \sa compareGen()*/
RcppExport  void  SPLUSnumGenMoy(SEXP sGenealogie, SEXP splProposant,SEXP sNProposant, SEXP sretour)		   
{

	STARTTIMER;
	int * Genealogie, * plProposant, * NProposant;
	double * retour;
	Rcpp::IntegerVector lGenealogie (sGenealogie);
	Rcpp::IntegerVector lplProposant(splProposant);
	Rcpp::IntegerVector lretour     (sretour);

	Genealogie  = INTEGER_POINTER( lGenealogie );
	plProposant = INTEGER_POINTER( splProposant );
	NProposant  = INTEGER_POINTER(sNProposant);
	retour 	  = NUMERIC_POINTER(sretour);
	
	numeroGenMoy(Genealogie, plProposant,*NProposant, retour);
	STOPTIMER;

	return ;
}
/// **********
//	CONTRIBUTION GENETIQUE
// ***********

/// Fonction d'interface Splus pour Congen
/** \sa Congen()*/
RcppExport void SPLUSConGen(SEXP sGenealogie, SEXP slProposant, SEXP sNProposant, SEXP slAncetre, SEXP sNAncetre, SEXP sdRetour, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie , * plProposant, * plAncetre, * NProposant, * NAncetre, * printit ;
	double * dRetour;
	Rcpp::IntegerVector lGenealogie(sGenealogie); // conversion automatique avec Rcpp
	Rcpp::IntegerVector lProposant(slProposant);
	Rcpp::IntegerVector lAncetre(slAncetre);
	
	Genealogie	= INTEGER_POINTER(lGenealogie);
	plProposant	= INTEGER_POINTER(lProposant);
	plAncetre		= INTEGER_POINTER(lAncetre);
	
	NProposant	= INTEGER_POINTER(sNProposant);
	NAncetre		= INTEGER_POINTER(sNAncetre);
	dRetour		= NUMERIC_POINTER(sdRetour);
	printit		= INTEGER_POINTER(sprintit);
	
	Congen(Genealogie, plProposant , *NProposant, plAncetre, *NAncetre, dRetour, *printit);
	STOPTIMER;
	return;
}

/// Fonction d'interface Splus pour Congen
/** \sa Congen()*/
RcppExport void SPLUSConGenPLUS(SEXP sGenealogie, SEXP splProposant,SEXP slNProposant, SEXP splAncetre, SEXP slNAncetre, SEXP spdSexe, 
						  SEXP spdRetour, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie , * plProposant, * plAncetre, * lNProposant, * lNAncetre, * printit ;
	double * pdRetour, * pdSexe;
	Rcpp::IntegerVector lGenealogie(sGenealogie); // conversion automatique avec Rcpp
	Rcpp::IntegerVector lProposant (splProposant);
	Rcpp::IntegerVector lAncetre   (splAncetre);
	
	Genealogie	= INTEGER_POINTER(lGenealogie);
	plProposant	= INTEGER_POINTER(lProposant);
	plAncetre		= INTEGER_POINTER(lAncetre);
	
	lNProposant	= INTEGER_POINTER(slNProposant);
	lNAncetre		= INTEGER_POINTER(slNAncetre);
	pdSexe		= NUMERIC_POINTER(spdSexe);
	pdRetour		= NUMERIC_POINTER(spdRetour);
	printit		= INTEGER_POINTER(sprintit);
	
	CongenPLUS(Genealogie, plProposant, *lNProposant, plAncetre, *lNAncetre, pdSexe, pdRetour, *printit);
	STOPTIMER;
	return;
}

/// Fonction d'interface Splus pour CongenCumul
/** \sa CongenCumul()*/
RcppExport void SPLUSCGCumul(	SEXP sGenealogie, SEXP splProposant,SEXP slNProposant, SEXP splAncetre, SEXP slNAncetre, SEXP sAncRet,
						SEXP spdRetour, SEXP spdRetourCumul, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie , * plProposant, * plAncetre, * lNProposant, * lNAncetre, * AncRet, * printit ;
	double * pdRetour, * pdRetourCumul;
	Rcpp::IntegerVector lGenealogie(sGenealogie); // conversion automatique avec Rcpp
	Rcpp::IntegerVector lProposant(splProposant);
	Rcpp::IntegerVector lAncetre(splAncetre);
	
	Genealogie	= INTEGER_POINTER(lGenealogie);
	plProposant	= INTEGER_POINTER(lProposant);
	plAncetre	= INTEGER_POINTER(lAncetre);
	AncRet		= NULL;
	lNProposant	= INTEGER_POINTER(slNProposant);
	lNAncetre	= INTEGER_POINTER(slNAncetre);
	pdRetour	= NUMERIC_POINTER(spdRetour);
	pdRetourCumul	= NUMERIC_POINTER(spdRetourCumul);
	printit		= INTEGER_POINTER(sprintit);

	CongenCumul(Genealogie, plProposant, *lNProposant, plAncetre, *lNAncetre, AncRet, pdRetour, pdRetourCumul, *printit);
	STOPTIMER
	return;
}

/// Fonction d'interface Splus pour CongenCumul
/** \sa CongenCumuldirect()*/
RcppExport void SPLUSCGCumuldirect(SEXP smatriceCG, SEXP slNProposant, SEXP splAncetre, SEXP slNAncetre, SEXP sAncRet, SEXP spdSomAnc, SEXP spdSomCumul)
{
	STARTTIMER;
	int * matriceCG, * lNProposant, * plAncetre, * lNAncetre, * AncRet;
	double * pdSomAnc, * pdSomCumul;
	Rcpp::IntegerVector lmatriceCG ( smatriceCG );
	Rcpp::IntegerVector lplAncetre ( splAncetre );
	Rcpp::IntegerVector lAncRet    ( sAncRet );
	
	matriceCG		= INTEGER_POINTER( lmatriceCG );
	plAncetre		= INTEGER_POINTER( lplAncetre );
	AncRet		= INTEGER_POINTER( lAncRet );
	lNProposant	= INTEGER_POINTER( slNProposant );
	lNAncetre		= INTEGER_POINTER( slNAncetre );
	pdSomAnc		= NUMERIC_POINTER( spdSomAnc );
	pdSomCumul	= NUMERIC_POINTER( spdSomCumul );
	 
	CongenCumuldirect(matriceCG, *lNProposant, plAncetre, *lNAncetre, AncRet,pdSomAnc,pdSomCumul);
	STOPTIMER
	return;
}

/// **********
//	DIVERS
// *********

/*FONCTION D'INTERFACE POUR SPLUS*/

/// Fonction d'interface Splus pour simul
/** \sa simul()*/
//RcppExport void SPLUSSimul(SEXP sGenealogie, SEXP sproposant, SEXP setatproposant, SEXP snproposant, SEXP sancetre, SEXP setatancetre,
RcppExport SEXP SPLUSSimul(SEXP sGenealogie, SEXP sproposant, SEXP setatproposant, SEXP snproposant, SEXP sancetre, SEXP setatancetre,
					  SEXP snancetre, SEXP snSimul, SEXP spdRetConj, SEXP spdRetSimul, SEXP spdRetProp, SEXP sprobRecomb,
					  SEXP sprobSurvieHomo, SEXP sPrintProgress)
{
	STARTTIMER;
	int * Genealogie, * proposant, * etatproposant, * nproposant, * ancetre, * etatancetre, * nancetre, * nSimul, * PrintProgress;
	double * pdRetConj, * pdRetSimul, * pdRetProp, * probRecomb;
	
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector letatproposant	( setatproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::IntegerVector letatancetre	( setatancetre );
	Rcpp::NumericVector lpdRetSimul	( spdRetSimul );
	Rcpp::NumericVector lpdRetProp	( spdRetProp );
	Rcpp::NumericVector lprobRecomb	( sprobRecomb );
	
	Genealogie	= INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER_POINTER( lproposant );
	etatproposant	= INTEGER_POINTER( letatproposant );
	ancetre		= INTEGER_POINTER( lancetre );
	etatancetre	= INTEGER_POINTER( letatancetre );
	nproposant	= INTEGER_POINTER( snproposant );
	nancetre		= INTEGER_POINTER( snancetre );
	nSimul		= INTEGER_POINTER( snSimul );
	PrintProgress	= INTEGER_POINTER( sPrintProgress );

	pdRetConj	 = NUMERIC_POINTER( spdRetConj );
	pdRetSimul = NUMERIC_POINTER( lpdRetSimul );
	pdRetProp	 = NUMERIC_POINTER( lpdRetProp );
	probRecomb = NUMERIC_POINTER( lprobRecomb );
	double probSurvieHomo  = Rcpp::as<double>(sprobSurvieHomo);

	simul(Genealogie, proposant, etatproposant, *nproposant, ancetre, etatancetre, *nancetre, *nSimul, pdRetConj, pdRetSimul,
		pdRetProp,  probRecomb, probSurvieHomo, *PrintProgress);
	STOPTIMER;
	//return;
	return Rcpp::wrap(getLastTimer());
}

/// Fonction d'interface Splus pour simulsingle
/** \sa simulsingle()*/
RcppExport void SPLUSSimulSingle(SEXP sGenealogie, SEXP sproposant, SEXP snproposant, SEXP sancetre, SEXP setatancetre, SEXP snancetre, 
						   SEXP sNSimul, SEXP spdRetour, SEXP sPrintProgress)
{
	STARTTIMER;
	int * Genealogie, * proposant, * nproposant, * ancetre, * etatancetre, * nancetre, * NSimul, * PrintProgress;
	double * pdRetour;
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector lancetre	( sancetre );
	Rcpp::IntegerVector letatancetre( setatancetre );
	
	Genealogie	= INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER_POINTER( lproposant );
	nproposant	= INTEGER_POINTER( snproposant );
	ancetre		= INTEGER_POINTER( lancetre );
	etatancetre	= INTEGER_POINTER( letatancetre );
	nancetre		= INTEGER_POINTER( snancetre );
	NSimul		= INTEGER_POINTER( sNSimul );
	PrintProgress	= INTEGER_POINTER( sPrintProgress );
	
	pdRetour		= NUMERIC_POINTER( spdRetour );
	
	simulsingle(Genealogie, proposant, *nproposant, ancetre, etatancetre, *nancetre, *NSimul, pdRetour, *PrintProgress);
	STOPTIMER;
	return;
}

/// Fonction d'interface Splus pour simulsingleFreq
/** \sa simulsingleFreq()*/
RcppExport void SPLUSSimulSingleFreq(SEXP sGenealogie, SEXP sproposant, SEXP snproposant, SEXP sancetre, SEXP setatancetre, SEXP snancetre,
							  SEXP sNSimul, SEXP spdRetour,SEXP sPrintProgress)
{
	STARTTIMER;
	int * Genealogie, * proposant, * nproposant, * ancetre, * etatancetre, * nancetre, * NSimul, * PrintProgress;
	double * pdRetour;
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::IntegerVector letatancetre	( setatancetre );
	
	Genealogie	= INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER_POINTER( lproposant );
	nproposant	= INTEGER_POINTER( snproposant );
	ancetre		= INTEGER_POINTER( lancetre );
	etatancetre	= INTEGER_POINTER( letatancetre );
	nancetre		= INTEGER_POINTER( snancetre );
	NSimul		= INTEGER_POINTER( sNSimul );
	PrintProgress	= INTEGER_POINTER( sPrintProgress );
	pdRetour		= NUMERIC_POINTER( spdRetour );
	
	simulsingleFreq(Genealogie, proposant, *nproposant, ancetre, etatancetre, *nancetre, *NSimul, pdRetour, *PrintProgress);
	STOPTIMER;
	return;
}

/// Fonction d'interface Splus pour simulsingleProb
/** \sa simulsingleProb()*/
RcppExport SEXP SPLUSSimulSingleProb( SEXP sGenealogie,SEXP sproposant, SEXP snproposant, SEXP sancetre,SEXP snancetre,
							   SEXP setatancetre,SEXP smtProb, SEXP sNSimul, SEXP sPrintProgress)
{
	STARTTIMER;
	int * Genealogie, * proposant, * nproposant, * ancetre, * etatancetre, * nancetre, * NSimul, * PrintProgress; //* mtProb, 
	
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::IntegerVector letatancetre	( setatancetre );
	Genealogie	= INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER_POINTER( lproposant );
	ancetre		= INTEGER_POINTER( lancetre );
	etatancetre	= INTEGER_POINTER( letatancetre );

	nproposant	= INTEGER_POINTER( snproposant );
	nancetre		= INTEGER_POINTER( snancetre );
	NSimul		= INTEGER_POINTER( sNSimul );
	PrintProgress	= INTEGER_POINTER( sPrintProgress );
	
	SEXP ret = simulsingleProb(Genealogie, proposant, *nproposant, ancetre, *nancetre, etatancetre, smtProb, *NSimul, *PrintProgress );
	STOPTIMER; 
	return ret;
}

/// Fonction d'interface Splus pour simulsingleFct
/** \sa simulsingleFct()*/
RcppExport SEXP SPLUSSimulSingleFct(SEXP sGenealogie, SEXP sproposant, SEXP sancetre, SEXP sancEtatAll1, SEXP sancEtatAll2, 
							 SEXP snancetre, SEXP sNSimul, SEXP sfctSousGrp, SEXP sPrintProgress)
{
	STARTTIMER;
	int * Genealogie,  * proposant, * ancetre, * ancEtatAll1, * ancEtatAll2, * nancetre, * NSimul, /* * fctSousGrp,*/ * PrintProgress;
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::IntegerVector lancEtatAll1	( sancEtatAll1 );
	Rcpp::IntegerVector lancEtatAll2	( sancEtatAll2 );
	
	Genealogie	= INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER_POINTER( lproposant );
	ancetre		= INTEGER_POINTER( lancetre );
	ancEtatAll1	= INTEGER_POINTER( sancEtatAll1 );
	ancEtatAll2	= INTEGER_POINTER( sancEtatAll2 );
	nancetre		= INTEGER_POINTER( snancetre );
	NSimul		= INTEGER_POINTER( sNSimul );
	PrintProgress	= INTEGER_POINTER( sPrintProgress );
	int nproposant = lproposant.size();
	
	SEXP ret = simulsingleFct(Genealogie, proposant, nproposant, ancetre, ancEtatAll1, ancEtatAll2, *nancetre, *NSimul, sfctSousGrp, *PrintProgress);
	STOPTIMER;
	return ret;
}

/// Fonction d'interface Splus pour prob
/** \sa prob()*/
RcppExport SEXP SPLUSProb(SEXP sGenealogie, SEXP sproposant, SEXP setatproposant,SEXP snproposant, SEXP sancetre, SEXP setatancetre, SEXP snancetre,
					 SEXP spdRetConj, SEXP spdRetSimul,SEXP sPrintProgress,SEXP sonlyConj)
{
	STARTTIMER;
	int * Genealogie, * proposant, * etatproposant, * nproposant, * ancetre, * etatancetre, * nancetre, * PrintProgress, * onlyConj;
	double * pdRetConj, * pdRetSimul;
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector letatproposant	( setatproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::IntegerVector letatancetre	( setatancetre );
	Rcpp::NumericVector lpdRetSimul	( spdRetSimul );
	
	Genealogie	= INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER_POINTER( lproposant );
	etatproposant	= INTEGER_POINTER( letatproposant );
	nproposant	= INTEGER_POINTER( snproposant );
	ancetre		= INTEGER_POINTER( lancetre );
	etatancetre	= INTEGER_POINTER( letatancetre );
	nancetre		= INTEGER_POINTER( snancetre );
	PrintProgress	= INTEGER_POINTER( sPrintProgress );
	onlyConj		= INTEGER_POINTER( sonlyConj );
	pdRetConj		= NUMERIC_POINTER( spdRetConj );
	pdRetSimul	= NUMERIC_POINTER( lpdRetSimul );
	
	SEXP ret = prob(Genealogie, proposant, etatproposant, *nproposant, ancetre, etatancetre, *nancetre, pdRetConj, pdRetSimul,
					*PrintProgress, *onlyConj);
	STOPTIMER;	
	return ret;
}


/// Fonction d'interface Splus pour CoefApparentement
/** \sa CoefApparentement()*/

RcppExport void SPLUSCoeffApparentement(SEXP sGenealogie, SEXP sproposant, SEXP snproposant, SEXP sancetre, SEXP sretour, SEXP sDuppDetection,
								SEXP sprintprogress)
{
	STARTTIMER
	int * Genealogie, * proposant, * nproposant, * ancetre, * DuppDetection, * printprogress;
	double * retour;
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::NumericVector lretour		( sretour );
	
	Genealogie	= INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER_POINTER( lproposant );
	nproposant	= INTEGER_POINTER( snproposant );
	ancetre		= INTEGER_POINTER( lancetre );
	retour		= NUMERIC_POINTER( lretour );
	DuppDetection	= INTEGER_POINTER( sDuppDetection );
	printprogress	= INTEGER_POINTER( sprintprogress );
	
	CoefApparentement(Genealogie, proposant, *nproposant, ancetre, retour, *DuppDetection, *printprogress);
	STOPTIMER;
	return; 
} 



//********************************************
//     FONCTION D'INTERFACE SPLUS CALL
//********************************************

//#ifndef MODETEST


/*! 
	\brief SPLUSCALL: Creer une genealogie binaire à l'aide d'une genealogie de forme classique No Ind, No Pere, No Mere

     Cette fonction est une interface SPLUS de la fonction CreerGenealogie()
	
	\param  SIndividu	[in] SEXP de type vecteur int contenant no individu

	\param  SPere		[in] SEXP de type vecteur int contenant no pere
								
	\param  SMere		[in] SEXP de type vecteur int contenant no mere
	
	\return Un pointeur vers un SEXP qui contient la genealogie en entier sous forme d'un long vecteur de int
	
	\remarque La genealogie composee de ind,pere,mere doit etre valide pour que la fonction s'execute correctement

	\sa CreerGenealogie()
*/
RcppExport SEXP SPLUSCALLCreerObjetGenealogie(SEXP SIndividu,SEXP SPere, SEXP SMere, SEXP SSexe)
{
	STARTTIMER
	//VARIABLE OPERATIONNEL		   conversion automatique avec Rcpp
	int * plIndividu, * plPere, * plMere, * plSexe;
//	Rcpp::IntegerVector lIndividu(SIndividu);
//	Rcpp::IntegerVector lPere(SPere);
//	Rcpp::IntegerVector lMere(SMere);
//	Rcpp::IntegerVector lSexe(0);
//	try{ lSexe = Rcpp::as<Rcpp::IntegerVector>(SSexe); } catch(std::exception &ex) { }
	
	Rcpp::IntegerVector lIndividu = Rcpp::as<Rcpp::IntegerVector>(SIndividu);
	Rcpp::IntegerVector lPere = Rcpp::as<Rcpp::IntegerVector>(SPere);
	Rcpp::IntegerVector lMere = Rcpp::as<Rcpp::IntegerVector>(SMere);
	Rcpp::IntegerVector lSexe = Rcpp::as<Rcpp::IntegerVector>(SSexe);
	
	plIndividu = INTEGER_POINTER(lIndividu);
	plPere     = INTEGER_POINTER(lPere);
	plMere     = INTEGER_POINTER(lMere);
	plSexe     = INTEGER_POINTER(lSexe);
	int lNIndividu = lIndividu.size();
	
	//Trois vecteur de meme taille
	if (lNIndividu != lPere.size() || lNIndividu != lMere.size())
	{
		ErrorHandler();
		//PROBLEM "%s","LES TROIS VECTEURS (ind,pere & mere) DOIVENT ETRE DE MEME TAILLE";
		//RECOVER(NULL_ENTRY);
	}
	
	//Verifie si le sexe est une information utile
	if (lSexe.size() != lNIndividu)	plSexe=NULL;
	
	int NombreEnfant=0; 
	int i;  

	//COMPLETE LA GENEALOGIE AVEC LES ELEMENTS MANQUANTS
	INITGESTIONMEMOIRE;	
	int* fInd	=(int*)  memalloc(lNIndividu*3,sizeof(int));	
	int* fpere	=(int*)  memalloc(lNIndividu*3,sizeof(int));	
	int* fmere	=(int*)  memalloc(lNIndividu*3,sizeof(int));
	int* fsexe 	= NULL;
	if (plSexe)
		fsexe	=(int*)  memalloc(lNIndividu*3,sizeof(int));
	CompleteGenealogie(plIndividu,plPere,plMere,plSexe,fInd,fpere,fmere,fsexe,&lNIndividu);

	//COMPTER LE NOMBRE D'ENFANT
	for(i=0;i<lNIndividu;i++) 	
	{
		if (fpere[i]!=0)	++NombreEnfant;
		if (fmere[i]!=0)	++NombreEnfant;
	}  

	//CREER UN OBJET XPLUS QUI CONTIENT TOUTE LES DONNEES
	const int TAILLESAUVEGARDE=TAILLEGENVERSION7(lNIndividu,NombreEnfant);
	int	* saveobj = new int[TAILLESAUVEGARDE];
	int * saveptr = saveobj;
//return (Rcpp::wrap( TAILLESAUVEGARDE ));

	//Creation de la genealogie
	CreerGenealogie(fInd,fpere,fmere,fsexe,lNIndividu,saveptr);
	Rcpp::IntegerVector retour(TAILLESAUVEGARDE);
	for(int i=0;i<TAILLESAUVEGARDE;i++) retour[i] = saveobj[i];
	STOPTIMER;
	for(int i=0;i<lNIndividu;i++) 
	{
		plIndividu[i] = fInd[i];
		plPere[i] = fpere[i];
		plMere[i] = fmere[i];
	}
	return (Rcpp::wrap( retour )); //saveobj;
} 

//#endif


//********************************************
//     FONCTION EXTRACTION DE GENEALOGIE
//******************************************** 

/*! 
	\brief Extraction de la genealogie sous la forme: No individu, No Pere, No mere

	La fonction Extrait la genealogie et retourne 3 vecteur soit le No Individu, le No du pere et le No de la Mere.
	
	\param genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\retval plRetIndividu [out] Vecteur de taille LengthGenealogie
								En cas de succes contient les No d'individus
	
	\retval plRetPere	[out] Vecteur de taille LengthGenealogie
								En cas de succes contient les No des peres

  	\retval plRetMere	[out] Vecteur de taille LengthGenealogie
								En cas de succes contient les No des meres

	\remark Il faut utiliser ( .C LengthGenealogie() ou SPLUS. gen.length) pour connaitre le nombre d'individu de la genealogie

	\sa LengthGenealogie()
*/
RcppExport SEXP SPLUSOutgen(SEXP Rgenealogie, SEXP RplRetIndividu, SEXP RplRetPere, SEXP RplRetMere, SEXP RplRetSexe, SEXP Rmustsort)
{
	STARTTIMER;
	//Rcpp::List par(params);
	Rcpp::IntegerVector a(Rgenealogie);
	Rcpp::IntegerVector b(RplRetIndividu);
	Rcpp::IntegerVector c(RplRetPere);
	Rcpp::IntegerVector d(RplRetMere);
	Rcpp::IntegerVector e(RplRetSexe);
	
	int * genealogie, * plRetIndividu, * plRetPere, * plRetMere, * plRetSexe, * mustsort;
	genealogie    = INTEGER_POINTER(a); //par["Data"]);
	plRetIndividu = INTEGER_POINTER(b); //par["ind"]);
	plRetPere     = INTEGER_POINTER(c); //par["pere"]);
	plRetMere     = INTEGER_POINTER(d); //par["mere"]);
	plRetSexe     = INTEGER_POINTER(e); //par["sexe"]);
	mustsort      = INTEGER_POINTER(Rmustsort);	
	
	//CREATION DU TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(genealogie, GFALSE, &lNIndividu, &Noeud); //Pas d'enfant  int test = 
	
	int havesex = (LoadNIndMasc()!=-1);
	//creation vecteur de sortie
	for(int i=0;i<lNIndividu;i++)
	{
		plRetIndividu[i]=Noeud[i].nom;
		if (Noeud[i].pere!=NULL)	plRetPere[i]=Noeud[i].pere->nom;		
		else						plRetPere[i]=0;
			
		if (Noeud[i].mere!=NULL)	plRetMere[i]=Noeud[i].mere->nom;			
		else						plRetMere[i]=0;
		//Recuperation du sexe
		if (havesex)				plRetSexe[i]=Noeud[i].sex;
		else						plRetSexe[i]=-1; //Pas utilisable
	}
	//Trie si nécessaire
	if (*mustsort)
		SortGenealogie3Vecteur(plRetIndividu,plRetPere,plRetMere,plRetSexe,lNIndividu);
	STOPTIMER;
	return (Rcpp::List::create(	Rcpp::Named("Data") = a,
							Rcpp::Named("ind")  = b,
							Rcpp::Named("father") = c,
							Rcpp::Named("mother") = d,
							Rcpp::Named("sex") = e
							));
}

/*! 
	\brief Extraction de la genealogie sous la forme: No individu, Indice Pere, Indice mere

	La fonction Extrait la genealogie et retourne 3 vecteur soit le No Individu, le Indice du pere et le Indice de la Mere.
	
	\param genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\retval plRetIndividu [out] Vecteur de taille LengthGenealogie
								En cas de succes contient les No d'individus
	
	\retval plRetPere	[out] Vecteur de taille LengthGenealogie
								En cas de succes contient les Indice des peres

  	\retval plRetMere	[out] Vecteur de taille LengthGenealogie
								En cas de succes contient les Indice des meres

	\remark Il faut utiliser ( .C LengthGenealogie() ou SPLUS. gen.length) pour connaitre le nombre d'individu de la genealogie

	\sa LengthGenealogie()
*/
RcppExport void SPLUSOutIndice(SEXP sgenealogie, SEXP splRetIndividu, SEXP splRetPere, SEXP splRetMere, SEXP splRetSexe, SEXP smustsort)
{

	STARTTIMER;			
	int * genealogie, * plRetIndividu, * plRetPere, * plRetMere, * plRetSexe, * mustsort;
	Rcpp::IntegerVector lgenealogie	( sgenealogie );
	Rcpp::IntegerVector lplRetIndividu	( splRetIndividu );
	Rcpp::IntegerVector lplRetPere	( splRetPere );
	Rcpp::IntegerVector lplRetMere	( splRetMere );
	Rcpp::IntegerVector lplRetSexe	( splRetSexe );
	
	genealogie	= INTEGER_POINTER( lgenealogie );
	plRetIndividu	= INTEGER_POINTER( lplRetIndividu );
	plRetPere 	= INTEGER_POINTER( lplRetPere );
	plRetMere 	= INTEGER_POINTER( lplRetMere );
	plRetSexe 	= INTEGER_POINTER( lplRetSexe );
	mustsort		= INTEGER_POINTER( smustsort );
	
	//CREATION DU TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(genealogie,GFALSE,&lNIndividu,&Noeud); //Pas d'enfant
	
	int havesex = (LoadNIndMasc()!=-1);
	//creation vecteur de sortie
	for(int i=0;i<lNIndividu;i++)
	{
		plRetIndividu[i]=Noeud[i].nom;
		if (Noeud[i].pere!=NULL)	plRetPere[i]=Noeud[i].pere->noind+1;		
		else						plRetPere[i]=0;
			
		if (Noeud[i].mere!=NULL)	plRetMere[i]=Noeud[i].mere->noind+1;			
		else						plRetMere[i]=0;
		//Recuperation du sexe
		if (havesex)				plRetSexe[i]=Noeud[i].sex;
		else						plRetSexe[i]=-1; //Pas utilisable
	}

	//Trie si nécessaire
	if (*mustsort)
		SortGenealogie3Vecteur(plRetIndividu,plRetPere,plRetMere,plRetSexe,lNIndividu);
	STOPTIMER;

	return;
}



//********************************************
//     FONCTIONS POUR L ANALYSE STATISTIQUE
//******************************************** 


/// Fonction d'interface Splus pour FondParGen
/** \sa initImplexe()*/
RcppExport void SPLUSFondParGen(SEXP sgenealogie, SEXP sprop, SEXP snbProp, SEXP sretour)
{
	STARTTIMER;
	int * genealogie, * prop, * nbProp, * retour;
	Rcpp::IntegerVector lgenealogie	( sgenealogie );
	Rcpp::IntegerVector lprop		( sprop );
	Rcpp::IntegerVector lretour		( sretour );
	
	genealogie = INTEGER_POINTER( lgenealogie );
	prop		 = INTEGER_POINTER( lprop );
	nbProp	 = INTEGER_POINTER( snbProp );
	retour 	 = INTEGER_POINTER( lretour );

	FondParGen(genealogie, prop, *nbProp, retour);
	STOPTIMER;

	return ;
}

