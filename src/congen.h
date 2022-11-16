/*! \file congen.h
\brief Interface des fonctions de calcul de la Contribution Genetique

Interface de toutes les fonctions en rapport avec la contribution gntique

\author Sbastien Leclerc
\contributor Jean-Franois Lefebvre
*/

#ifndef GENCONGEN
#define GENCONGEN


int Congen(int* Genealogie, int* plProposant,int lNProposant, int* plAncetre, int lNAncetre, double* pdCongen,int printprogress);

int CongenPLUS(int* Genealogie,
	int* plProposant,int lNProposant,  
	int* plAncetre, int lNAncetre, double* pdSexe, 
	double* pdCongen,int printprogress);

int CongenCumul(int* Genealogie,
	int* plProposant,int lNProposant,
	int* plAncetre, int lNAncetre,
    int* AncRet,double* pdSomAnc,double* pdSomCumul, int printprogress);

int CongenCumuldirect(int* matriceCG,
	int lNProposant,
	int* plAncetre, int lNAncetre,
    int* AncRet,double* pdSomAnc,double* pdSomCumul);

#endif






