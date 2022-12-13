/*! \file fondateur.cc
\brief Implementation des fonctions de simulation calcul de probabilite

Calcul et Analyse de diverse valeur driv du gene fondateur


\author Sbastien Leclerc 
\contributor Jean-Francois Lefebvre
*/


/// Authorise l'affichage du niveau de progression sur la sortie standard
/** 
	Si ALLOWPRINTPROGRESS est defini, les fonctions suivantes peuvent
	indiquer leur niveau de progression sur la sortie standard stdout
*/
#define ALLOWPRINTPROGRESS
#include "base.h"
#include "outils.h" 
#include "cbignum.h"
#include "hashtable.h"
#include "userInterface.h"
#include "fondateur.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <string.h> 
#include <math.h> 
#include <limits.h>
#include <vector>
#include <random>
#include <cstdlib>
#include <memory>
#include <RcppCommon.h>
#include <R.h>
//#include <Rdefines.h>

#include <Rcpp.h>
#include <Rcpp/as.h>
#include <Rcpp/Function.h>
#include <unordered_map>

#include "asa239.h"

#ifdef NEG
	#undef NEG
#endif
extern "C"
{
	#include "mpi.h"
	#include "mplogic.h"
}

// ******************************************************************** 
//
//			CONSTANTE & STRUCTURE
//
// ********************************************************************

const int MEGAOCTET = 1048576;//octet

///Liste de tableau de unsigned char.... Utilis pour representer une liste de nombre binaire de longueur variable
struct CApPath
{
	///Grand nombre binaire en format mpi reprsentant le chemin....
	mp_int num;
	///Pointeur vers l'element suivant de la liste
	CApPath *next;
};

//CONSTRUCTION ET DESTRUCTION DES PATHS
static CIndSimul** g_ExpCoeff_CheminParcouru=NULL; 
///Utilise par ExploreCoeff comme tant le dernier chemin remplis (ou le premiers)
static CApPath ** g_ExpCoeff_Path=NULL;
///Utilise par ExploreCoeff comme tant la cible de l'exploration
static CIndSimul* g_ExpCoeff_Cible=NULL;
static void FASTCALL ExploreCoeff(CIndSimul* Noeud);
static void PathDestruction(CApPath **Path,int npath);


// std::string dump_hapref(std::unordered_map<int,haplotype*> *hapRef)
// {
// 	std::stringstream out;
// 	out << "hapRef\n";
// 	for(auto h : (*hapRef)) {
// 		out << "  " << h.first << "\n";
// 		haplotype *hap = h.second;
// 		out << "    &hap:         " << hap << std::endl;
// //		out << "    hap:          " << hap->hap << "\n";
// 		out << "    pos:          " << hap->pos << "\n";
// 		out << "    fixe:         " << hap->fixe << "\n";
// 		out << "    next_segment: " << hap->next_segment << "\n";
// 	}

// 	return out.str();
// }


// ********************************************************************
//
//			PUBLIC
//
// ********************************************************************
/*! 
*/ 

class Crossovers
{
	public:		
		void Poisson_CO(const int&, double*, double*, int&, std::mt19937&, double*);
		void Poisson_ZT(const int&, double*, double*, int&, std::mt19937&, double*);
		void init_gamma(double&, double&, double&, double& );
		void Gamma_CO  (const int&, double*, double*, int&, std::mt19937&, double*);

		
 	private:
		double first_arrival[2][10000]; //Numbins = 10,000. Need to calculate first arrival times for gamma process. They are distributed differently than the interrarival times

};

void Crossovers::Poisson_CO(const int &sex, double *param, double *Morgan_len, int& NumRecomb, std::mt19937& mtgen, double *CO_array){
	double pos;			
	static std::uniform_real_distribution<> u_dist(0, 1);
	static std::poisson_distribution<int> p1_dist(param[0]);
	static std::poisson_distribution<int> p2_dist(param[1]);

	if(sex==1){
		NumRecomb = p1_dist(mtgen); //number of recombination events 
		for( int h=0; h<NumRecomb; h++ ) {
			pos = u_dist(mtgen);
			CO_array[h] = pos; 
		}
		
		std::sort(CO_array, CO_array + NumRecomb);
	}
	else{
		NumRecomb = p2_dist(mtgen); //number of recombination events 
		for( int h=0; h<NumRecomb; h++ ) {
			pos = u_dist(mtgen);
			CO_array[h] = pos; 
		}
		
		std::sort(CO_array, CO_array + NumRecomb);		
	}

}

void Crossovers::Poisson_ZT(const int &sex, double *param, double *Morgan_len, int& NumRecomb, std::mt19937& mtgen, double *CO_array){
	double chiasma_pos[40];
	double pos;
	int	   num_chiasma = 0;		
	static std::uniform_real_distribution<> u_dist(0, 1);
	static std::poisson_distribution<int> p1_dist(param[0]);
	static std::poisson_distribution<int> p2_dist(param[1]);

	if(sex==1){
		num_chiasma = p1_dist(mtgen); //number of recombination events 
		while (num_chiasma == 0){
			num_chiasma = p1_dist(mtgen);
		}
		
		for( int h=0; h<num_chiasma; h++ ) {
			pos = u_dist(mtgen);
			chiasma_pos[h] = pos; 
		}
		
	}
	else{
		num_chiasma = p2_dist(mtgen); //number of recombination events 
		while (num_chiasma == 0){
			num_chiasma = p2_dist(mtgen);
		}
		
		for( int h=0; h<num_chiasma; h++ ) {
			pos = u_dist(mtgen);
			chiasma_pos[h] = pos; 
		}
	}
	
	int counter=0;
	double u_rand;
	NumRecomb = 0;
		for (int i=0; i<num_chiasma; i++){
			u_rand = u_dist(mtgen);
			if (u_rand < 0.50){
				CO_array[counter] = chiasma_pos[i];
				counter   += 1;
				NumRecomb += 1;
			}
		}
	
	std::sort(CO_array, CO_array + NumRecomb);
}

void Crossovers::init_gamma(double& paramF,  double& paramM,  double& Morgan_LenF,  double& Morgan_LenM){
	double x;
	int wasteman = 0;

	//Number of bins is 10,000, these are the bins for reimann sum
	//Morgan_Len/10000 is the delta x of the integral
	//gamma_q is 1-CDF. We multiply by 2 because we are dividing by mean of regular distribution (which has mean 1/2)
	
	x= Morgan_LenF/10000;
	first_arrival[0][0] = 2 * (1-gammad(2*paramF*x, paramF, &wasteman)) * Morgan_LenF/10000;

	for ( int i=1; i<10000; i++){
		x = Morgan_LenF*(i+1)/10000;
		first_arrival[0][i] = 2 * (1-gammad(2*paramF*x, paramF, &wasteman)) * Morgan_LenF/10000 + first_arrival[0][i-1];
	}	

	x= Morgan_LenM/10000;
	first_arrival[1][0] = 2 * (1-gammad(2*paramM*x, paramM, &wasteman)) * Morgan_LenM/10000;

	for ( int i=1; i<10000; i++){
		x = Morgan_LenM*(i+1)/10000;
		first_arrival[1][i] = 2 *  (1-gammad(2*paramM*x, paramM, &wasteman)) * Morgan_LenM/10000 + first_arrival[1][i-1];
	}
}

void Crossovers::Gamma_CO(const int &sex, double *param, double *Morgan_len, int& NumRecomb, std::mt19937& mtgen, double *CO_array){
	double u_rand;
	double length;
	double step;
	double interrarival;
	double current_pos;
	double chiasma_pos[20];
	int Num_Chiasmata = 0;
	
	static std::uniform_real_distribution<> u_dist(0, 1);
	static std::gamma_distribution<> g1_dist(param[0],1/(2*param[0]));
	static std::gamma_distribution<> g2_dist(param[1],1/(2*param[1]));

	if (sex==1){
		length = Morgan_len[0];
	} else{
		length = Morgan_len[1];
	}

	step = length/10000;
    u_rand = u_dist(mtgen); 

	if ( u_rand > first_arrival[sex-1][9999]) NumRecomb = 0; //IF the first chiasmata is beyond the length of the chromsome then we have no recombination
	else {
		if (u_rand <= first_arrival[sex-1][0]){ 
			chiasma_pos[0] = 0.5 * step;
			Num_Chiasmata = 1;
		}
		else { //binary search to find the first position (sample between 0-1, then lookup the inverse CDF values)
			int mid, low = 0, high = 10000;
			while (high - low > 1) {
				mid = (high - low) / 2 + low;
				if (u_rand <= first_arrival[sex-1][mid])
					high = mid;
				else if (u_rand > first_arrival[sex-1][mid])
					low = mid;
			}
			Num_Chiasmata=1;
			//mid will end up being the index of the smallest value (step of reimann partial sums) we are less than or equal to
			chiasma_pos[0] = mid*step + 0.5*step;
		}

		current_pos = chiasma_pos[0];
		//After finding the first arrival time we can  use std::gamma_distribution for the rest
		if (sex==1){
			interrarival = g1_dist(mtgen);
			int counter = 1;
			while (current_pos + interrarival < Morgan_len[0]){
				Num_Chiasmata = Num_Chiasmata + 1;

				chiasma_pos[counter] = current_pos + interrarival;
				current_pos 		 = current_pos + interrarival;

				counter 	 = counter + 1;
				interrarival = g1_dist(mtgen);
			}

		}
		else{
			interrarival = g2_dist(mtgen);
			int counter = 1;
			while (current_pos + interrarival < Morgan_len[1]){
				Num_Chiasmata = Num_Chiasmata + 1;

				chiasma_pos[counter] = current_pos + interrarival;
				current_pos 		 = current_pos + interrarival;

				counter 	 = counter + 1;
				interrarival = g2_dist(mtgen);
			}
		}

		//After determining position of all chiasmata, we select each one with probability 0.5 of resolving as a crossover 
		int counter=0;
		NumRecomb = 0;
		for (int i=0; i<Num_Chiasmata; i++){
			u_rand = u_dist(mtgen);
			if (u_rand < 0.50){
				CO_array[counter] = chiasma_pos[i]/length;
				counter   += 1;
				NumRecomb += 1;
			}
		}
	}
}

void simulhaplo(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre,
						int lSimul, double* probRecomb, double* Morgan_Len, int BP_len, int model, int convert,  
						double* cm_map_FA, double* cm_map_MO, int* bp_map_FA, int* bp_map_MO, 
						std::unordered_map<int,haplotype*> *hapRef, std::string WD, int write_all_node, int seed) 
{

	try{

	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie, GTRUE, &lNIndividu, &Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);

	//CREATION OF AN ANCESTOR VECTOR *** with the reference haplotypes for the chosen ancestors. ***
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);

	//Creation des tableau
	INITGESTIONMEMOIRE;
	CIndSimul** Ordre = (CIndSimul**) memalloc(lNIndividu,sizeof(CIndSimul*));

	//Pour le sort spcial		
	int*	OrdreSaut	= (int*) memalloc(lNIndividu,sizeof(int*));				
	int NOrdre;
	
	int i;

	//Initialize all the nodes
	for(i=0;i<lNIndividu;i++)
	{
		Noeud[i].allele = 0;
		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].bFlagSort=0;
		Noeud[i].clesHaplo_1 = 0; //"0.1";
		Noeud[i].clesHaplo_2 = 0; //"0.1";
	}

	//label the nodes that are probands
	for(i=0;i<lNProposant;i++){
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
	}

	int cleFixe = 1; // haplotype keys
	
	//identifier et etiqueter les points de departs et les haplos ancetres (starting points and ancestor haplotypes)
	for(i=0;i<lNAncetre;i++)
	{		
		NoeudAnc[i]->allele = 0;
		NoeudAnc[i]->etat=GENDEPART;
	    NoeudAnc[i]->clesHaplo_1 = cleFixe++;
	    NoeudAnc[i]->clesHaplo_2 = cleFixe++;

		std::ostringstream nom;
		nom << NoeudAnc[i]->nom << ".1";
		haplotype *tmp1 = new haplotype();//[1];
		tmp1->hap  = nom.str();
		tmp1->pos  = -1;
		tmp1->fixe = 1;
		(*hapRef)[NoeudAnc[i]->clesHaplo_1] = tmp1;
		
		nom.str(std::string());
		nom << NoeudAnc[i]->nom << ".2";
		haplotype *tmp2 = new haplotype();//[1];
		tmp2->hap  = nom.str();
		tmp2->pos  = -1;
		tmp2->fixe = 1;
		(*hapRef)[NoeudAnc[i]->clesHaplo_2] = tmp2;

	}

	//identifier et marque les noeuds utile et ceux inutile a la recherche
	for(i=0;i<lNAncetre;i++)
		ExploreArbre(NoeudAnc[i]);
		
	//create the order of traversal and calculate the jumps (Jumps are unecessary, only for speeding up allele calculations, will test if it works without)
	PrepareSortPrioriteArbre(Noeud,lNIndividu);	
	NOrdre=0;

	memset(OrdreSaut,0,sizeof(int)*lNIndividu);
	for(i=0;i<lNAncetre;i++)		
		StartSortPrioriteArbre(NoeudAnc[i],Ordre,&NOrdre,OrdreSaut); // les infos de NoeudAnc sont pointes par Ordre dans "le bon ordre".

	//create mersenne twister generator to use for random  distributions 
	std::mt19937 my_rng = std::mt19937(seed);
	//initialize crossover class
	Crossovers crossovers;
  
	void (Crossovers::*SampleCO)(const int&, double*, double*, int&, std::mt19937&, double*);
	
	if 		(model==1) SampleCO=&Crossovers::Poisson_CO;
	else if (model==2) SampleCO=&Crossovers::Poisson_ZT;
	else if (model==3){
		SampleCO=&Crossovers::Gamma_CO;
		crossovers.init_gamma(probRecomb[0], probRecomb[1], Morgan_Len[0], Morgan_Len[1]);
	}

	void (*convert_dist)(int&, double*, const double&, const int&, int*, double*, int*); //(nbrecomb, CO_array, Morgan_len, bp_len, bp_map, cm_map, precision)

	if 		(convert == 0) convert_dist = &no_convert;  // no genetic/physical map is specified, then no need to scale wrt physical distance (assumed 1:1)
	else if (convert == 1) convert_dist = &convert1; 	// genetic/physical map is specified
	
	double pHap1, pHap2;

	double CO_arrayF[20]; //hold crossover positions for Father's chromosome
	double CO_arrayM[20]; //hold crossover positions for Mother's chromsome

	int BP_CO_arrayF[20]; //crossover positions in BP
	int BP_CO_arrayM[20];

	int nbRecomb1 =0;
	int nbRecomb2 =0;
	
	std::uniform_real_distribution<> u_dist(0, 1);

	//create output files
	std::string stroutHaplo = WD + "/Proband_Haplotypes.txt";
	std::ofstream outHaplo(stroutHaplo);
	outHaplo << std::dec << std::fixed<< std::setprecision(0);
	outHaplo 	<< lSimul << ";" << lNProposant << "\n";

	if(write_all_node == 1){
		std::string	stroutAllHaplo = WD + "/All_nodes_haplotypes.txt";
		std::ofstream outAllHaplo(stroutAllHaplo);
		outAllHaplo << lSimul << ";" << lNProposant << ";" << NOrdre << "\n";
		outAllHaplo << std::dec << std::fixed<< std::setprecision(0);

		for(int csimul=0;csimul<lSimul; csimul++)
		{
			int clesSim= cleFixe;
			
			for(int i=0;i<NOrdre;i++) {
				//Simulate meiosis in the parents, store the location of crossovers (in genetic distance scaled to [0,1]) in CO_array.
				(crossovers.*SampleCO)(1, probRecomb, Morgan_Len, nbRecomb1, my_rng, CO_arrayF);
				(crossovers.*SampleCO)(2, probRecomb, Morgan_Len, nbRecomb2, my_rng, CO_arrayM);

				//Now the locations of crossovers in CO_array will be converted to physical distance. 
				convert_dist(nbRecomb1, CO_arrayF, Morgan_Len[0], BP_len, bp_map_FA ,cm_map_FA, BP_CO_arrayF);
				convert_dist(nbRecomb2, CO_arrayM, Morgan_Len[1], BP_len, bp_map_MO ,cm_map_MO, BP_CO_arrayM);

				pHap1 = u_dist(my_rng); 
				makeRecombF(Ordre[i], hapRef, pHap1, nbRecomb1, BP_CO_arrayF, clesSim);
				
				pHap2 = u_dist(my_rng);
				makeRecombM(Ordre[i], hapRef, pHap2, nbRecomb2, BP_CO_arrayM, clesSim);

				if (write_all_node == 1){
					outAllHaplo <<"{"<< csimul+1 <<";"<< Ordre[i]->nom <<";" ;
					
					if(Ordre[i]->pere != NULL){
						outAllHaplo << nbRecomb1 <<",";
						if(nbRecomb1 > 0){ //Recombination event in the father
							outAllHaplo << pHap1;
							for(int j=0; j<nbRecomb1; j++){
								outAllHaplo << "," << BP_CO_arrayF[j];
							}
							outAllHaplo << ";";
						}
						else{
							outAllHaplo << pHap1 << ",0;";
						}
					}		
					else{
						outAllHaplo << "0,0,0;";
					}		

				
					if(Ordre[i]->mere != NULL){
						outAllHaplo << nbRecomb2 << ",";
						if(nbRecomb2 > 0){ //Recombination event in mother
							outAllHaplo << pHap2;
							for(int j=0;j<nbRecomb2;j++){
								outAllHaplo << "," <<  BP_CO_arrayM[j];
							}
							outAllHaplo << "}";
						}	
						else{
							outAllHaplo << pHap2 << ",0}";
						}
					}
					else{
						outAllHaplo << "0,0,0}";
					}

					haplotype* tmp = (*hapRef).find(Ordre[i]->clesHaplo_1)->second;

					int pos = tmp->pos;
					if(pos == -1) pos = BP_len;

					for( int h=0; h<2; h++ ) {
						outAllHaplo << "{" << 0 << ";" << tmp->hap << ";" << pos ; 
						while( tmp->next_segment != NULL) { 
							tmp = tmp->next_segment;
							pos = tmp->pos;
							if(pos == -1) pos = BP_len;
							outAllHaplo << ";" << tmp->hap << ";" <<  pos ;
						}
						
						outAllHaplo << "}";
						tmp = (*hapRef).find(Ordre[i]->clesHaplo_2)->second;
						pos = tmp->pos;
						if(pos == -1) pos = BP_len;
					}	
					outAllHaplo << "\n";
				}
			}
		

			for(i=0;i<lNProposant;i++){
				haplotype* tmp = (*hapRef).find(NoeudPro[i]->clesHaplo_1)->second;
				double pos = tmp->pos;
				if(pos == -1) pos = 1;

				outHaplo <<"{"<< csimul+1 <<";"<< NoeudPro[i]->nom << ";0}";

				for( int h=0; h<2; h++ ) {
					outHaplo << "{" << 0 << ";" << tmp->hap << ";" << pos ; 
			
					while( tmp->next_segment != NULL) { 
						tmp = tmp->next_segment;
						pos = tmp->pos;
						if(pos == -1) pos = BP_len;
						outHaplo << ";" << tmp->hap << ";" <<  pos;
					}
					
					outHaplo << "}";
					tmp = (*hapRef).find(NoeudPro[i]->clesHaplo_2)->second;
					pos = tmp->pos;
					if(pos == -1) pos = BP_len;
				}
				outHaplo << "\n";
			}
			//delete haplotypes from memory before next iteration of simulation
			for( int i=cleFixe; i<clesSim; i++) {
				haplotype* tmp = (*hapRef).find(i)->second; //hapKey.second;
				while(tmp->next_segment != NULL) {
					haplotype* tmp_back = tmp;
					tmp = tmp->next_segment;
					delete tmp_back;
				}
				delete tmp;
			}
			
		} // end of the for loop that goes through the # of simulations
		
		outHaplo.close();
		if(write_all_node == 1) outAllHaplo.close();

		for(int i=0; i<cleFixe; i++) {
			haplotype* tmp = (*hapRef).find(i)->second;//hapKey.second;
			while(tmp->next_segment != NULL) {
				haplotype* tmp_back = tmp;
				tmp = tmp->next_segment;
				delete tmp_back;
			}
		delete tmp;
		}
	}

	else{
	//simulation loop
		for(int csimul=0;csimul<lSimul; csimul++)
		{
			int clesSim= cleFixe;
			
			for(int i=0;i<NOrdre;i++) {
				//Simulate meiosis in the parents, store the location of crossovers (in genetic distance scaled to [0,1]) in CO_array.
				(crossovers.*SampleCO)(1, probRecomb, Morgan_Len, nbRecomb1, my_rng, CO_arrayF);
				(crossovers.*SampleCO)(2, probRecomb, Morgan_Len, nbRecomb2, my_rng, CO_arrayM);

				//Now the locations of crossovers in CO_array will be converted to physical distance. 
				convert_dist(nbRecomb1, CO_arrayF, Morgan_Len[0], BP_len, bp_map_FA ,cm_map_FA, BP_CO_arrayF);
				convert_dist(nbRecomb2, CO_arrayM, Morgan_Len[1], BP_len, bp_map_MO ,cm_map_MO, BP_CO_arrayM);

				pHap1 = u_dist(my_rng); 
				makeRecombF(Ordre[i], hapRef, pHap1, nbRecomb1, BP_CO_arrayF, clesSim);
				
				pHap2 = u_dist(my_rng);
				makeRecombM(Ordre[i], hapRef, pHap2, nbRecomb2, BP_CO_arrayM, clesSim);
			}

			for(i=0;i<lNProposant;i++){
				haplotype* tmp = (*hapRef).find(NoeudPro[i]->clesHaplo_1)->second;
				double pos = tmp->pos;
				if(pos == -1) pos = 1;

				outHaplo <<"{"<< csimul+1 <<";"<< NoeudPro[i]->nom << ";0}";

				for( int h=0; h<2; h++ ) {
					outHaplo << "{" << 0 << ";" << tmp->hap << ";" << pos ; 
			
					while( tmp->next_segment != NULL) { 
						tmp = tmp->next_segment;
						pos = tmp->pos;
						if(pos == -1) pos = BP_len;
						outHaplo << ";" << tmp->hap << ";" <<  pos;
					}
					
					outHaplo << "}";
					tmp = (*hapRef).find(NoeudPro[i]->clesHaplo_2)->second;
					pos = tmp->pos;
					if(pos == -1) pos = BP_len;
				}
				outHaplo << "\n";
			}
			//delete haplotypes from memory before next iteration of simulation
			for( int i=cleFixe; i<clesSim; i++) {
				haplotype* tmp = (*hapRef).find(i)->second; //hapKey.second;
				while(tmp->next_segment != NULL) {
					haplotype* tmp_back = tmp;
					tmp = tmp->next_segment;
					delete tmp_back;
				}
				delete tmp;
			}
		}// end of the for loop that goes through the # of simulations
		outHaplo.close();

		for(int i=0; i<cleFixe; i++) {
			haplotype* tmp = (*hapRef).find(i)->second;//hapKey.second;
			while(tmp->next_segment != NULL) {
				haplotype* tmp_back = tmp;
				tmp = tmp->next_segment;
				delete tmp_back;
			}
		delete tmp;
		}
	} 
	

	

 } catch(std::exception &ex) {
 	forward_exception_to_r(ex);
 } catch(...){
 	::Rf_error("c++ exception (unknown reason)"); 
 } 
}

void no_convert(int& nbrecomb, double* CO_array, const double& Morgan_len, const int& bp_len, int* bp_map, double* cm_map, int* BP_array){
	for(int k=0; k<nbrecomb; k++){
		BP_array[k] = (int) (CO_array[k]*bp_len);
		if (k>0){
			if(BP_array[k]==BP_array[k-1]){
				BP_array[k]=BP_array[k]+1;
			}
		}
	}

	//check if any recombinations are at same position, remove duplicate
	// if(nbrecomb>1){
	// 	double temp[20];
	// 	int j = 0;
	// 	temp[j++] = CO_array[0];
	// 	for(int k=1; k<nbrecomb; k++){
	// 		if(CO_array[k] != CO_array[k-1]){
	// 				temp[j++] = CO_array[k];
	// 		}
	// 	}
	// 	for (int k=0; k<j; k++){
	// 		CO_array[k] = temp[k];
	// 	}
	// 	nbrecomb = j;		
	// }
}

void convert1(int& nbrecomb, double* CO_array, const double& Morgan_len, const int& bp_len, int* bp_map, double* cm_map, int* BP_array){
	for (int k=0; k<nbrecomb; k++){
		double physical_distance;
		double genetic_distance = CO_array[k]*Morgan_len*100;

		int map_index = 0;
		while(genetic_distance > cm_map[map_index]) map_index++;

		physical_distance = bp_map[map_index-1]  + (bp_map[map_index] - bp_map[map_index-1])*(genetic_distance - cm_map[map_index-1])/(cm_map[map_index] - cm_map[map_index-1]);
		BP_array[k] = (int)(physical_distance);

		if (k>0){
			if(BP_array[k]==BP_array[k-1]){
				BP_array[k]=BP_array[k]+1;
			}
		}

	}
	
	// double temp[20];
	// int j = 0;
	// if(nbrecomb>1){
	// 	temp[j++] = CO_array[0];
	// 	for(int k=1; k<nbrecomb; k++){
	// 		if(CO_array[k] != CO_array[k-1]){
	// 				temp[j++] = CO_array[k];
	// 		}
	// 	}
	// 	for (int k=0; k<j; k++){
	// 		CO_array[k] = temp[k];
	// 	}
	// 	nbrecomb = j;		
	// }
}

//F for father, not female
void makeRecombF( CIndSimul *Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, int nbRecomb, int *posRecomb, int &cle )
{
    haplotype *perehap1, *perehap2;
	if(Ordre_tmp->pere != NULL){
		if (nbRecomb > 0){
			if (probHap < 0.5){
				perehap1=(*hapRef).find(Ordre_tmp->pere->clesHaplo_1)->second;
				perehap2=(*hapRef).find(Ordre_tmp->pere->clesHaplo_2)->second;
			} 
			else{
				perehap1=(*hapRef).find(Ordre_tmp->pere->clesHaplo_2)->second;
				perehap2=(*hapRef).find(Ordre_tmp->pere->clesHaplo_1)->second;
			}

			haplotype *hapChild_1 = new haplotype();
			haplotype *hapChild_deb1 = hapChild_1;
			recombine(perehap1, perehap2, hapChild_deb1, nbRecomb, posRecomb);
			Ordre_tmp->clesHaplo_1 = cle;
			(*hapRef)[cle++] = hapChild_1;
		}
		else{ //If no recombination just pass one of father's chromosomes down 
			if(probHap<0.50){
				Ordre_tmp->clesHaplo_1=Ordre_tmp->pere->clesHaplo_1;
			}
			else{
				Ordre_tmp->clesHaplo_1=Ordre_tmp->pere->clesHaplo_2;
			}
		}
	}
	else Ordre_tmp -> clesHaplo_1 = 0;
}

//M for mother, not male
void makeRecombM( CIndSimul *Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, int nbRecomb, int *posRecomb, int &cle )
{
    haplotype *merehap1, *merehap2;
 	if(Ordre_tmp->mere != NULL){
		if (nbRecomb > 0){
			if (probHap < 0.5){
				merehap1=(*hapRef).find(Ordre_tmp->mere->clesHaplo_1)->second;
				merehap2=(*hapRef).find(Ordre_tmp->mere->clesHaplo_2)->second;
			} 
			else{
				merehap1=(*hapRef).find(Ordre_tmp->mere->clesHaplo_2)->second;
				merehap2=(*hapRef).find(Ordre_tmp->mere->clesHaplo_1)->second;
			}

			haplotype *hapChild_2 = new haplotype();
			haplotype *hapChild_deb2 = hapChild_2;
			recombine(merehap1, merehap2, hapChild_deb2, nbRecomb, posRecomb);
			Ordre_tmp->clesHaplo_2 = cle;
			(*hapRef)[cle++] = hapChild_2;
		}
		else{ //If no recombination just pass one of mother's chromosomes down 
			if(probHap<0.50){
				Ordre_tmp->clesHaplo_2=Ordre_tmp->mere->clesHaplo_1;
			}
			else{
				Ordre_tmp->clesHaplo_2=Ordre_tmp->mere->clesHaplo_2;
			}
		}
	}
	else Ordre_tmp -> clesHaplo_2 = 0;
			
}


void recombine(haplotype* hapBegin, haplotype* hapEnd, haplotype* hapChild, int nbRecomb, int* posRecomb)
{
	haplotype* hap_active = hapBegin;

	for (int i=0; i < nbRecomb; i++){
		
		int position = posRecomb[i];
		// de 0 a posRecomb on prend hapBegin
		while(position > hap_active->pos && hap_active->pos != -1) {
			(*hapChild).hap          = hap_active->hap;
			(*hapChild).pos          = hap_active->pos;
			(*hapChild).fixe         = 0;
			(*hapChild).next_segment = new haplotype();//[1];
			hapChild                 = hapChild->next_segment;
			hap_active               = hap_active->next_segment;
		}

		// on ajoute la recomb pour hapChild
		(*hapChild).hap          = hap_active->hap;
		(*hapChild).pos          = position;
		(*hapChild).fixe         = 0;

		if(i%2 == 0){
			hap_active = hapEnd;
		}
		else hap_active = hapBegin;

		// on met le pointeur de hapEnd a la bonne place.
		while(position >= hap_active->pos && hap_active->pos != -1) hap_active = hap_active->next_segment;

		// on verifie que l'haplotype qui suit n'est pas le meme. Si oui, on met la nouvelle position.
		if(hap_active->hap == hapChild->hap){
			(*hapChild).pos          = hap_active->pos;
		}
		else{
			(*hapChild).next_segment = new haplotype();//[1];
			hapChild                 = hapChild->next_segment;
			(*hapChild).hap          = hap_active->hap;
			(*hapChild).pos          = hap_active->pos;
			(*hapChild).fixe         = 0;
		}
	}

	while(hap_active->pos != -1) {
		hap_active               = hap_active->next_segment;
		(*hapChild).next_segment = new haplotype();//[1]; 
		hapChild                 = hapChild->next_segment;
		(*hapChild).hap          = hap_active->hap;
		(*hapChild).pos          = hap_active->pos;
		(*hapChild).fixe         = 0;
	}

}

// int getNumberRec(double* probRecomb, int sex, int seed)
// {
// //  #if defined _WIN32 || defined _WIN64
// //    std::mt19937 gen(time(0));
// //  #else
// static  std::mt19937 mt_rand;
// //  #endif
//   if(sex==1) {
//    std::poisson_distribution<int> distribution (probRecomb[0]);
//    return distribution(mt_rand);
//   }
//   else {
//    std::poisson_distribution<int> distribution (probRecomb[1]);
//    return distribution(gen);
//   }
// }

// double getRandomNumber(int exponential)
// {
// 	static std::random_device gen;
// 	if(exponential == 0) {
//     //#if defined _WIN32 || defined _WIN64
//       //std::mt19937 gen(time(0));
//     //#else
//     //#endif
//     return double(gen())/double(gen.max());
	
//   }
//   else {
// //    std::default_random_engine position_generator;
//     std::exponential_distribution<double> alea_Exp(10.0);
//     return alea_Exp( gen );//position_generator );
//   }
// }

//no longer call this function in simulhaplo, do it directly in the main loop of simulhaplo
// int descendreHaplotypes(CIndSimul* Ordre_tmp, double probHap)
// {
//   if(Ordre_tmp->pere != NULL && Ordre_tmp->mere != NULL) {
//     if     (probHap < 0.25){ Ordre_tmp->clesHaplo_1 = Ordre_tmp->pere->clesHaplo_1; Ordre_tmp->clesHaplo_2 = Ordre_tmp->mere->clesHaplo_1; }
//     else if(probHap < 0.50){ Ordre_tmp->clesHaplo_1 = Ordre_tmp->pere->clesHaplo_2; Ordre_tmp->clesHaplo_2 = Ordre_tmp->mere->clesHaplo_1; }
//     else if(probHap < 0.75){ Ordre_tmp->clesHaplo_1 = Ordre_tmp->pere->clesHaplo_1; Ordre_tmp->clesHaplo_2 = Ordre_tmp->mere->clesHaplo_2; }
//     else				  { Ordre_tmp->clesHaplo_1 = Ordre_tmp->pere->clesHaplo_2; Ordre_tmp->clesHaplo_2 = Ordre_tmp->mere->clesHaplo_2; }
//   }
//   else if( Ordre_tmp->pere != NULL ) {  // Faire qq chose ici pour la mere et le pere qui est NULL ...
//     if(probHap < 0.5) Ordre_tmp->clesHaplo_1 = Ordre_tmp->pere->clesHaplo_1;
//     else			  Ordre_tmp->clesHaplo_1 = Ordre_tmp->pere->clesHaplo_2;
//     Ordre_tmp->clesHaplo_2 = 0;
//   }
//   else if( Ordre_tmp->mere != NULL ) {
//     if(probHap < 0.5) Ordre_tmp->clesHaplo_2 = Ordre_tmp->mere->clesHaplo_1;
//     else			  Ordre_tmp->clesHaplo_2 = Ordre_tmp->mere->clesHaplo_2;
//     Ordre_tmp->clesHaplo_1 = 0;
//   }
//   else { // les 2 sont NULL
//     Ordre_tmp->clesHaplo_1 = 0;
//     Ordre_tmp->clesHaplo_2 = 0;
//   }
//   return 0;
// }

//no longer use this function for simulhaplo now it instead uses makeRecombF for recombination of father's chromosomes and makeRecombM for mother
// void makeRecomb( CIndSimul *Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, double posRecomb, int &cle )
// {
//   int pereHap = 0, mereHap = 0;
//   if     (probHap < 0.25){ pereHap = 1; mereHap = 1; }
//   else if(probHap < 0.50){ pereHap = 1; mereHap = 0; }
//   else if(probHap < 0.75){ pereHap = 0; mereHap = 1; }

//   haplotype *hapPere, *hapMere;
//   if(Ordre_tmp->pere!= NULL && Ordre_tmp->mere!= NULL) {
//   if(pereHap == 0) hapPere = (*hapRef).find(Ordre_tmp->pere->clesHaplo_1)->second;
//   else             hapPere = (*hapRef).find(Ordre_tmp->pere->clesHaplo_2)->second;
//   if(mereHap == 0) hapMere = (*hapRef).find(Ordre_tmp->mere->clesHaplo_1)->second;
//   else             hapMere = (*hapRef).find(Ordre_tmp->mere->clesHaplo_2)->second;
//   }
//   else if( Ordre_tmp->pere != NULL ) {
//     if(pereHap == 0) hapPere = (*hapRef).find(Ordre_tmp->pere->clesHaplo_1)->second;
//     else             hapPere = (*hapRef).find(Ordre_tmp->pere->clesHaplo_2)->second;
//     hapMere = (*hapRef).find(0)->second;
//   }
//   else if( Ordre_tmp->mere != NULL ) {
//     hapPere = (*hapRef).find(0)->second;
//     if(mereHap == 0) hapMere = (*hapRef).find(Ordre_tmp->mere->clesHaplo_1)->second;
//     else             hapMere = (*hapRef).find(Ordre_tmp->mere->clesHaplo_2)->second;
//   }
//   else {
//     hapPere = (*hapRef).find(0)->second;
//     hapMere = (*hapRef).find(0)->second;
//   }

//   haplotype *hapChild_1 = new haplotype();
//   haplotype *hapChild_deb1 = hapChild_1;
//   recombine(hapPere, hapMere, hapChild_deb1, posRecomb);
//   Ordre_tmp->clesHaplo_1 = cle;
//   (*hapRef)[cle++] = hapChild_1;

//   haplotype *hapChild_2 = new haplotype();
//   haplotype *hapChild_deb2 = hapChild_2;
//   recombine(hapMere, hapPere, hapChild_deb2, posRecomb);  
//   Ordre_tmp->clesHaplo_2 = cle;
//   (*hapRef)[cle++] = hapChild_2;
  
// }

// These commented out functions: reconstruct, readSNPpos, ancestralseq, were for converting the results of simulhaplo into sequence data.
// This funcitonality is removed for now, we just provide the perl script to do the same thing. 

// bool reconstruct(std::string WD, const std::string &simufilename,const std::string &hapfilename, const std::string &SNPposfilename,const int &BPsize){

// 	try{

//     std::ifstream in (simufilename.c_str());
//     if(!in)
//     {
//         Rcpp::stop ("Cannot open the proband_haplotypes file ");
//     }

// 	WD += "/reconstructed_haplotypes.txt";
//     std::ofstream reconstructed(WD.c_str());
//     if(!reconstructed.is_open()){
//         Rcpp::stop("Can't open output file to write to. Check permissions of output directory");
//     }

//     std::vector<int> SNPpos = readSNPpos(SNPposfilename);
//     int numSNPs = SNPpos.size();

//     std::unordered_map <float, std::string> haploseqs;
//     ancestralseq(hapfilename, haploseqs);

//     std::string line;
//     std::getline(in,line); //waste first line
//     while (std::getline(in, line))
//     {
//         std::size_t tokenPos, tokenPos1; 
//         tokenPos=line.find(";");
//         tokenPos1=line.find(";", tokenPos+1);
//         reconstructed << line.substr(1,tokenPos-1) << " " << line.substr(tokenPos+1,tokenPos1-tokenPos-1) << " ";

        
//         tokenPos= line.find("}");   
//         tokenPos1= line.find("}",tokenPos+1);
//         std::string hap1(line.substr(tokenPos+2,tokenPos1-tokenPos-2));
//         tokenPos=line.find("}",tokenPos1+1);
//         std::string hap2(line.substr(tokenPos1+2,tokenPos-tokenPos1-2));
//         //reconstructing haplotype1
//         tokenPos=hap1.find(";");
//         tokenPos1=hap1.find(";",tokenPos+1);

//         int SNPpos_ind =0; 
//         float ancID;
//         int seg_end =0;

//         ancID = std::stof(hap1.substr(tokenPos+1,tokenPos1-tokenPos-1));
//         tokenPos1 = hap1.find(";", tokenPos1 +1);

//         while (tokenPos1 != std::string::npos){
//             tokenPos = hap1.find(";",tokenPos+1);
//             seg_end = std::stof(hap1.substr(tokenPos+1,tokenPos1-tokenPos-1))*BPsize;
//             tokenPos1=hap1.find(";",tokenPos1+1);

//             while(SNPpos[SNPpos_ind]<seg_end && SNPpos_ind < numSNPs){
//                 reconstructed << haploseqs[ancID].at(SNPpos_ind);
//                 SNPpos_ind++;
//             }
//             tokenPos = hap1.find(";", tokenPos +1);
//             ancID = std::stof(hap1.substr(tokenPos+1,tokenPos1-tokenPos-1));
//             tokenPos1 = hap1.find(";", tokenPos1+1);
//         }
        
//         if (seg_end == 0){
//             reconstructed << haploseqs[ancID] << std::endl;
//         }
//         else{
//             reconstructed << haploseqs[ancID].substr(SNPpos_ind, std::string::npos) <<std::endl;
//         } 
        
//         //reconstructing haplotype2
//         tokenPos=line.find(";");
//         tokenPos1=line.find(";", tokenPos+1);
//         reconstructed << line.substr(1,tokenPos-1) << " " << line.substr(tokenPos+1,tokenPos1-tokenPos-1) << " ";

//         tokenPos=hap2.find(";");
//         tokenPos1=hap2.find(";",tokenPos+1);

//         SNPpos_ind =0; 
//         seg_end =0;

//         ancID = std::stof(hap2.substr(tokenPos+1,tokenPos1-tokenPos-1));
//         tokenPos1 = hap2.find(";", tokenPos1 +1);
        
//         while (tokenPos1 != std::string::npos){
//             tokenPos = hap2.find(";",tokenPos+1);
//             seg_end = std::stof(hap2.substr(tokenPos+1,tokenPos1-tokenPos-1))*BPsize;
//             tokenPos1=hap2.find(";",tokenPos1+1);

//             while(SNPpos[ SNPpos_ind]<seg_end && SNPpos_ind < numSNPs){
//                 reconstructed << haploseqs[ancID].at(SNPpos_ind);
//                 SNPpos_ind++;
//             }

//             tokenPos = hap2.find(";", tokenPos +1);
//             ancID = std::stof(hap2.substr(tokenPos+1,tokenPos1-tokenPos-1));
//             tokenPos1 = hap2.find(";", tokenPos1+1);
//         }
        
//         if (seg_end == 0){
//             reconstructed << haploseqs[ancID] << std::endl;
//         }
//         else{
//             reconstructed << haploseqs[ancID].substr(SNPpos_ind, std::string::npos) << std::endl;
//         } 
//     }
    

//     in.close();
//     reconstructed.close();
//     return true;

// 	} catch(std::exception &ex) {
//  	forward_exception_to_r(ex);
//  	} catch(...){
//  	::Rf_error("c++ exception (unknown reason)"); 
//  	}
// 	return false; 
// }

// bool ancestralseq(const std::string &fileName, std::unordered_map<float, std::string> &haploseqs)
// {
//     std::ifstream in(fileName.c_str());

//     if(!in)
//     {
//         Rcpp::stop("Cannot open the hapfile");
//         return false;
//     }

//     float anc_id;
//     std::string anc_haplo;

//     while (in>>anc_id>>anc_haplo)
//     {
//         haploseqs[anc_id]=anc_haplo;
//     }

//     in.close();
//     return true;
// }

// std::vector<int> readSNPpos(const std::string &fileName){
    
// 	std::ifstream in(fileName.c_str());
    
//     if(!in)
//     {
//         Rcpp::stop("Cannot open the mapfile");
//     }

//     std::vector<int> vec(std::istream_iterator<int>(in), {});

//     in.close();
//     return vec;
// }

//functions to analyze the output of simulhaplo

int tb_digest_line(const std::string& chr_string, const int& myAnc, int& numHits, std::vector<int> &target_pos_L, std::vector<int> &target_pos_R){
	int lastpos = 0, counter = 1;
    bool push_back = false;

	std::size_t tokenPos, tokenPos1;
	tokenPos = chr_string.find(';');

    while(tokenPos != std::string::npos) {
		tokenPos1 = chr_string.find(';', tokenPos + 1);
        if (counter%2==1){ 
            if (std::stoi(chr_string.substr(tokenPos + 1, tokenPos1 - tokenPos - 1 - 2)) == myAnc){ //-2 because the segmentIDs for now have trailing ".1" or ".2" to encode the chromosome copy
                push_back = true;
                numHits ++;
            }
        }
        else{
            if (push_back){
                target_pos_L.push_back(lastpos);
                target_pos_R.push_back(std::stoi(chr_string.substr(tokenPos + 1, tokenPos1 - tokenPos - 1)));
                push_back = false;
            }
            lastpos = std::stoi(chr_string.substr(tokenPos + 1, tokenPos1 - tokenPos - 1));
        }
        counter++;
		tokenPos = tokenPos1;
    }
    return 0;
};

int tb_digest_line2(const std::string& recomb_string, int& pHap, int& nRecomb, int* RecPos){
    size_t c_tokenPos  = recomb_string.find(',');
    nRecomb = std::stoi(recomb_string.substr(0, c_tokenPos)); 

    size_t c_tokenPos1 = recomb_string.find(',', c_tokenPos+1);
    pHap = (std::stoi(recomb_string.substr(c_tokenPos+1, c_tokenPos1 - c_tokenPos - 1)));

    
    for (int k=0; k < nRecomb; k++) {
        c_tokenPos = recomb_string.find(',', c_tokenPos1 + 1);
        RecPos[k]  = std::stoi(recomb_string.substr(c_tokenPos1 + 1, c_tokenPos-c_tokenPos1-1));
        c_tokenPos1 = c_tokenPos;
    }
    return 0;
};

struct tb_hap{
    int nRec, pHap;
    int RecPos[20];
};

struct tb_ind;

struct tb_ind{
    int ID;
    tb_ind *parents[2];
    tb_hap chr[2];
};

int traceback_internal( tb_ind* curr_ind, int curr_chr, const int& myAnc, const int& Lpos, const int& Rpos, int* tb_path, int& pathlen){
//    tb_hap*  curr_hap  = &(curr_ind->chr[curr_chr]); //loop through the tb_ind tree, curr_hap = current haplotype, curr_ind = current individual
    tb_ind*  next_ind  = curr_ind->parents[curr_chr];
	
    bool keep_looping = true;
    int counter = 0;
    // log << "Lpos: " << Lpos << " Rpos: " << Rpos << "\n" << std::flush;
    while(keep_looping){
		// log << curr_ind->ID << " " << curr_ind->chr[curr_chr].pHap << " " << curr_ind->chr[curr_chr].nRec << "\n"  << std::flush;
        tb_path[counter] = next_ind->ID;
        counter++;
  
        int count_recomb = 0; //count the number of recombinations that occur before the ancestral segment (to see which parent its from

        if(curr_ind->chr[curr_chr].nRec ==0){
            curr_chr = curr_ind->chr[curr_chr].pHap;
            curr_ind = next_ind;
            next_ind = curr_ind->parents[curr_chr]; 
            //because pHap tells which of parent chromosomes is inherited
        }
        else{
            for(int k=0; k<(curr_ind->chr[curr_chr].nRec); k++){
				// log << curr_ind->chr[curr_chr].RecPos[k] << " " << std::flush;
                if(curr_ind->chr[curr_chr].RecPos[k] <= Lpos) count_recomb++;
                else if((curr_ind->chr[curr_chr].RecPos[k] > Lpos) & (curr_ind->chr[curr_chr].RecPos[k] < Rpos)){
					pathlen = counter;
                    return -9;
                }
            }
			// log << "\n" << count_recomb << "\n" << std::flush;
            if (count_recomb%2 == 1){
                curr_chr = curr_ind->chr[curr_chr].pHap;
                curr_chr = 1 - curr_chr; //if there is odd # of recombinations before our segment it comes from the other parent (not the start of the chromosome)
                curr_ind = next_ind;
                next_ind = curr_ind->parents[curr_chr]; 
            }
            else {
                curr_chr = curr_ind->chr[curr_chr].pHap;
                curr_ind = next_ind;
                next_ind = curr_ind->parents[curr_chr];
            }
        }

        if (next_ind->ID == myAnc) {
			keep_looping = false;
			tb_path[counter] = myAnc;
			++counter;
		}
        if (counter>100) {
            return -10;
            //throw exception
        }
	pathlen = counter;
    }
	return 0;
}


inline bool check_duplicate_path(std::vector<int> pathvec, int& new_pathlen, int* new_path){
	for (int i=0; i < new_pathlen; ++i ){
		if (pathvec.at(i) != new_path[i]) return false;
	}
	return true;
}

int simulhaplo_traceback(std::string& path_ANH, std::string& path_PH, int& myPro, int& myAnc, 
					std::vector<int>& indVec, std::vector<int>& mereVec, std::vector<int>& pereVec,
					std::vector<int>& resultvec1, std::vector<int>& resultvec2, std::vector<int>& resultvec3) 
{
	try{
	// std::ofstream log("log.txt");
    typedef std::unordered_map<int, std::unique_ptr<tb_ind>> tb_dict;
    tb_dict my_tb_dict;

    std::ifstream file_pro_haplo (path_PH);
    std::ifstream file_all_haplo (path_ANH);

    std::string line;
    std::getline(file_pro_haplo, line); //

    std::size_t tokenPos, tokenPos1;
    tokenPos =line.find(";");

    int numSim = std::stoi(line.substr(0,tokenPos));
    int numPro = std::stoi(line.substr(tokenPos+1));

    std::getline(file_all_haplo, line);
    tokenPos = line.find(';');
    int numSim2 = std::stoi(line.substr(0,tokenPos));
	tokenPos1 = line.find(';', tokenPos + 1);
    int numInd  = std::stoi(line.substr(tokenPos1+1));

    //MAKE tb_dict
    for (const int& ind : indVec){
		tb_ind* tb_ind_ptr = new tb_ind();
		tb_ind_ptr->ID = ind;
		my_tb_dict.emplace(ind, std::unique_ptr<tb_ind>(tb_ind_ptr)); //for every ind in ind vec initialize a dict entry
	}

	my_tb_dict.emplace(0, std::unique_ptr<tb_ind> (new tb_ind()));
    for (std::size_t i = 0, max = indVec.size(); i < max; i++){ //go through second time now to assign parents (who were assigned in first loop)
		my_tb_dict.at(indVec.at(i))->parents[0] = my_tb_dict.at(pereVec.at(i)).get();
        my_tb_dict.at(indVec.at(i))->parents[1] = my_tb_dict.at(mereVec.at(i)).get();
	}

    std::vector<std::vector<int>> unique_paths;
    unique_paths.reserve(30);

    std::string chr_string;
    for(int j=0; j<numSim; j++){
        //
		// log << "sim: " << j << "\n" << std::flush;
        int c1_numHits =0, c2_numHits = 0;
        std::vector<int> c1_target_pos_L, c2_target_pos_L, c1_target_pos_R, c2_target_pos_R;

        int tb_path[100];
        int pathlen =0 , Lpos = 0, Rpos = 0, curr_chr = 0, ProID;
        
        for(int i=0; i<numPro; i++){
            std::getline(file_pro_haplo,line);
            tokenPos  = line.find(';');
            tokenPos1 = line.find(';', tokenPos+1); 

            // what if they want to trace back from an internal node not a proband? can add it in later...

            ProID = std::stoi(line.substr(tokenPos+1, tokenPos1 - tokenPos - 1));
            if(ProID == myPro){
                tokenPos    = line.find('}');
                tokenPos1   = line.find('}', tokenPos + 1);
                chr_string  = line.substr(tokenPos+2, tokenPos1-tokenPos-2);
                //tb_digest_line will push back the position vectors, and modify numHits,
                tb_digest_line(chr_string, myAnc, c1_numHits, c1_target_pos_L, c1_target_pos_R);
                tokenPos    = line.find('}', tokenPos1 + 1);
                chr_string  = line.substr(tokenPos1+2, tokenPos-tokenPos1-2);
                tb_digest_line(chr_string, myAnc, c2_numHits, c2_target_pos_L, c2_target_pos_R);
            }
        }
        if (c1_numHits + c2_numHits == 0){//IF the proband doesn't have a segment from the specified ancestor then waste the next nInd lines of all_nodes
            for (int i=0; i<numInd; i++) file_all_haplo.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        else{
            int indID; //Read the All_nodes_haplo file and fill tb_dict
            for (int i=0; i<numInd; i++){
                std::getline(file_all_haplo, line);
				// log << line << "\n" << std::flush;
                tokenPos  = line.find(';');
                tokenPos1 = line.find(';', tokenPos + 1);
				// log << line.substr(tokenPos+1, tokenPos1-tokenPos-1) << "\n" << std::flush;
                indID = std::stoi(line.substr(tokenPos+1, tokenPos1-tokenPos-1));
                tb_ind* node = my_tb_dict.at(indID).get();

                tokenPos   = line.find(';', tokenPos1 + 1);
                chr_string = line.substr(tokenPos1+1, tokenPos - tokenPos1-1);
				// log << chr_string << "\n" << std::flush;
                tb_digest_line2(chr_string, node->chr[0].pHap, node->chr[0].nRec, node->chr[0].RecPos);
				// log << indID << " " << node->chr[0].nRec << " " << node->chr[0].pHap << "\n" << std::flush;

                tokenPos1  = line.find('}', tokenPos + 1);
                chr_string = line.substr(tokenPos + 1, tokenPos1 - tokenPos - 1);
				// log << chr_string << "\n" << std::flush;
                tb_digest_line2(chr_string, node->chr[1].pHap, node->chr[1].nRec, node->chr[1].RecPos);
				// log << indID << " " << node->chr[1].nRec << " " << node->chr[1].pHap << "\n" << std::flush;

            };
            //now the tb_dict is initialized, with the recomb history, need to do the traceback         
            // log << "c1 hits: " << c1_numHits << "\n" << std::flush;
            // log << "c2 hits: " << c2_numHits << "\n" << std::flush;

            for(int i=0; i<c1_numHits; i++){
				resultvec1.push_back(j+1);
                Lpos = c1_target_pos_L[i];
                Rpos = c1_target_pos_R[i];
				resultvec2.push_back(Rpos-Lpos);
                tb_ind* curr_ind  = my_tb_dict.at(myPro).get();
                curr_chr  = 0;

				// log << "traceback internal\n" << std::flush;
            	traceback_internal(curr_ind, curr_chr, myAnc, Lpos, Rpos, tb_path, pathlen);
				// log << "traceaback internal done\n" << std::flush;

				// Rcpp::Rcout  << "simulation: " << j + 1 << "\nchr1 segment: " << Lpos << "->" << Rpos << "\npathlen: " << pathlen << "\npath: ";
				// for(int h=0; h<pathlen; h++) Rcpp::Rcout << tb_path[h] << " ";
				// Rcpp::Rcout << "\n";
				
				bool duplicate_path = false;
				int path_count = 0;
				std::size_t veclen;

                for(std::vector<int>& pathvec : unique_paths ){
					++path_count;
					veclen = pathvec.size();
					// log << "check duplicate path\n" <<std::flush;
					if (veclen == (std::size_t)pathlen) duplicate_path = check_duplicate_path(pathvec, pathlen, tb_path);
					// log << "check duplicate done\n" <<std::flush;
					if (duplicate_path)    break;
					
				}

				if (!duplicate_path) {
					resultvec3.push_back(path_count + 1);
					unique_paths.emplace_back(std::vector<int> (tb_path, tb_path + pathlen));
				} else resultvec3.push_back(path_count);
                // }
                //check the num_unique paths to see if it exists yet
                //for (int h=0; h<unique_paths.size(); h++) {
                // inline function to check each element and return exists= index of path
                // if exists: results_path.push back(index+1), results_SimNum.push_back(simno), results.push_back(length])
                //
                // if doesnt exist: create new vector for unique_paths.emplace_back(std::vector<int> vec(tb_path,tb_path+pathlen) );
                //}
                //for 
            }

            for(int i=0; i<c2_numHits; i++){
				resultvec1.push_back(j+1);
                Lpos = c2_target_pos_L[i];
                Rpos = c2_target_pos_R[i];
				resultvec2.push_back(Rpos-Lpos);
                tb_ind* curr_ind  = my_tb_dict.at(myPro).get();
                curr_chr  = 1;
 
				// log << "traceback internal\n" << std::flush;
            	traceback_internal(curr_ind, curr_chr, myAnc, Lpos, Rpos, tb_path, pathlen);
				// log << "traceaback internal done\n" << std::flush;
				
				bool duplicate_path = false;
				int path_count = 0;
				std::size_t veclen;

                for(std::vector<int>& pathvec : unique_paths ){
					++path_count;
					veclen = pathvec.size();
					// log << "check duplicate\n" << std::flush;
					if (veclen == (std::size_t)pathlen) duplicate_path = check_duplicate_path(pathvec, pathlen, tb_path);
					// log << "check duplicate done\n" << std::flush;
					if (duplicate_path)    break;
				}

				if (!duplicate_path) {
					resultvec3.push_back(path_count + 1);
					unique_paths.emplace_back(std::vector<int> (tb_path, tb_path + pathlen));
				} else resultvec3.push_back(path_count);           
			}
        }
    }

	line = "";
	for (int i = 0, veclen = unique_paths.size(); i != veclen; ++i){
		line += ("\npath: "+ std::to_string(i+1)+ " ");
		for (const int& node : unique_paths.at(i)) line += (std::to_string(node)+ " ");
	}  
	Rcpp::message(Rcpp::wrap(line));

	return 0;
	} catch(std::exception &ex) {
 	forward_exception_to_r(ex);
 	} catch(...){
 	::Rf_error("c++ exception (unknown reason)"); 
 	} return 0;
};

//Following functions are all used for IBD_compare: add_overlap, check_overlap, digest_line, 

struct seg_positions{  
	int pos[300];
	int segID[300];
	int numseg = 0;
};

struct overlaps{
	int Lpos[50];
	int Rpos[50];
	int num_overlaps = 0;
};

//find IBD regions shared between hap1 and hap2
inline void seg_overlap(seg_positions& hap1, seg_positions& hap2, overlaps& overlap){
    for(int i=0; i<hap1.numseg; i++){
        int Lpos1 = hap1.pos[i];
        int Rpos1 = hap1.pos[i+1];
        for (int j=0; j<hap2.numseg; j++){
            int Lpos2 = hap2.pos[j];
            int Rpos2 = hap2.pos[j+1];

            if( Lpos2 >= Rpos1 ){ // no overlap
                break;
            }
            else if( Rpos2 > Lpos1){
                if( hap1.segID[i] == hap2.segID[j]){
					int num_overlaps = overlap.num_overlaps;
					int Lpos = (Lpos1 >= Lpos2) ? Lpos1 : Lpos2;
					int Rpos = (Rpos1 <= Rpos2) ? Rpos1 : Rpos2;
					if(num_overlaps == 0){
						overlap.num_overlaps = 1;
						overlap.Lpos[0] = Lpos;
						overlap.Rpos[0] = Rpos;
					}
					else{ //check if overlapping segment is perfectly adjacent to last segment, in which case don't count it as new segment, just elongate the last
						if(overlap.Rpos[num_overlaps-1] == Lpos){
							overlap.Rpos[num_overlaps-1] = Rpos;
						}
						else{
							overlap.Lpos[num_overlaps] = Lpos;
							overlap.Rpos[num_overlaps] = Rpos;
							overlap.num_overlaps = num_overlaps + 1;
						}
					}
                }
            }
        }
    }
};

//
inline bool check_overlaps(const int& Lpos1, const int& Rpos1, const int& Lpos2, const int& Rpos2){
	// then its non-overlapping
	if ((Lpos1 > Rpos2) || (Lpos2 > Rpos1)) return false;
	return true;
}

inline void digest_line(const std::string& chr_string, seg_positions& haplo){
    std::string          temp_string;                
    std::stringstream    ss(chr_string);

	int counter = 0;
    while(std::getline(ss, temp_string, ';')) {
        if (counter%2==1){
            if(temp_string.back()=='1'){
				haplo.segID[counter/2] = stoi(temp_string.substr(0, temp_string.size() - 2));
            }
            else{
                haplo.segID[counter/2] = - stoi(temp_string.substr(0, temp_string.size() - 2));
            }
        }
        else { 
            haplo.pos[counter/2] = stoi(temp_string);
        }
        counter++;
    }
    haplo.numseg = counter/2;
};

// haplotype of onne haploid copy of ind1 comparred against both copies (diploid) of ind2
//if ind2 is HBD in any spots that were IBD with the copy from ind1 thenn we want to consider this as one segment
//additionally if any of the IBD segments overlap (due to HBD region) then want to combine them into the maximum length single segment
//hx1 and hx2 are the comparisons from the haploid chromosome of ind1 to the two copies of ind 2
// h is the output structure that will hold the final IBD segments of this haploid copy of ind1 after adjusting for any instances of HBD in ind2
inline void check_HBD(const overlaps& hx1, const overlaps& hx2, overlaps& h){
	int memory_Lpos[50], memory_Rpos[50]; //list to store overlappinng HBD rregions. overlaps can have recursive-like effect where counting them as a longer segment then overlaps with a different segment, need to combine them all into one segment with maximal bounds
	int memory_count = 0;
	int Lpos1, Rpos1, Lpos2, Rpos2; //temporary positions
	int count = 0 ; //count for the umber of segments of the new 'h' output struct

	//std::vector default initializes all values, default initialization of bool is false
	std::vector<bool> HBD_hx2(hx2.num_overlaps); //we will double loop through hx1 and hx2, this will store for hx2 whether the segments have any HBD overlap or not

	// In future update can improve from the double loop, since they are in order don't have to do the full loop
	// though in most instances the number of IBD segments would be rather smmall and adding extra checks to save on having to do the full loop might add more overhead than you save
	for(int i=0; i<hx1.num_overlaps; ++i){
		Lpos1 = hx1.Lpos[i];
		Rpos1 = hx1.Rpos[i];
		bool any_overlap = false;
		for(int j=0; j<hx2.num_overlaps; ++j){
			Lpos2 = hx2.Lpos[j];
			Rpos2 = hx2.Rpos[j];
			if (Lpos2>= Rpos1) break;
			if (check_overlaps(Lpos1, Rpos1, Lpos2, Rpos2)){
				any_overlap = true;
				HBD_hx2[j] = true;
				if (memory_count == 0){
					memory_Lpos[memory_count] = (Lpos1 <= Lpos2) ? Lpos1 : Lpos2;
					memory_Rpos[memory_count] = (Rpos1 >= Rpos2) ? Rpos1 : Rpos2;
					++memory_count;
				}
				else{
					Lpos2 =  (Lpos1 <= Lpos2) ? Lpos1 : Lpos2; //reuse the earlier variable names, they'll get overwritten next loop
					Rpos2 =  (Rpos1 >= Rpos2) ? Rpos1 : Rpos2;
					
					if (check_overlaps(memory_Lpos[memory_count -1], memory_Rpos[memory_count -1], Lpos2, Rpos2)){
						memory_Lpos[memory_count -1] = (memory_Lpos[memory_count -1] <= Lpos2) ? memory_Lpos[memory_count -1] : Lpos2;
						memory_Rpos[memory_count -1] = (memory_Rpos[memory_count -1] >= Rpos2) ? memory_Rpos[memory_count -1] : Rpos2;
					}
					else{
						memory_Lpos[memory_count] = Lpos2;
						memory_Rpos[memory_count] = Rpos2;
						++memory_count;
					}
				}
			}
		};
		if (!any_overlap){
			h.Lpos[count] = Lpos1;
			h.Rpos[count] = Rpos1;
			++count;
		}
	}
	for (int i=0; i<hx2.num_overlaps; ++i){
		if(!HBD_hx2.at(i)){
			h.Lpos[count] = hx2.Lpos[i];
			h.Rpos[count] = hx2.Rpos[i];
			++count;
		}
	}
	
	for (int i=0; i<memory_count; ++i){
		h.Lpos[count] = memory_Lpos[i];
		h.Rpos[count] = memory_Rpos[i];
		++count;
	}

	h.num_overlaps = count;
}

void simulhaplo_compare_IBD(const int& pro1_ID, const int& pro2_ID, const int& BP_len, std::string& file_path, std::vector<int>& rvec1, std::vector<int>& rvec2, std::vector<double>& rvec3, std::vector<int>& rvec4) {
    try{
	// std::ofstream log("log.txt");
	std::ifstream  in (file_path);
    std::string line;
    std::getline(in, line); //

    std::size_t tokenPos, tokenPos1;
    tokenPos =line.find(";");

    int numPro, numSim;
    numSim = std::stoi(line.substr(0,tokenPos));
    numPro = std::stoi(line.substr(tokenPos+1));

    int ProID, SimNo;
    std::string chr_string;
	seg_positions pro1_1, pro1_2, pro2_1, pro2_2;

    for(int j=0; j<numSim; ++j){
		// log << "simulation: " << j << "\n";
        //loop through all probands
        for(int i=0; i<numPro; ++i){
            std::getline(in,line);
            tokenPos  = line.find(';');
            tokenPos1 = line.find(';', tokenPos+1); 
            ProID = std::stoi(line.substr(tokenPos+1, tokenPos1));

            //for probands of interest populate their haplo_dict
            // if(std::find(myPro.begin(), myPro.end(), ProID) != myPro.end()){
			if( ProID == pro1_ID){
                tokenPos    = line.find('}');
                tokenPos1   = line.find('}', tokenPos + 1);
                chr_string  = line.substr(tokenPos+2, tokenPos1-tokenPos-2);
				// log << ProID << "\n" << chr_string << "\n";
                digest_line(chr_string, pro1_1);

                tokenPos    = line.find('}', tokenPos1 + 1);
                chr_string = line.substr(tokenPos1+2, tokenPos-tokenPos1-2);
				// log << chr_string << "\n";
                digest_line(chr_string, pro1_2);
            }
			else if( ProID == pro2_ID){
				tokenPos    = line.find('}');
                tokenPos1   = line.find('}', tokenPos + 1);
                chr_string  = line.substr(tokenPos+2, tokenPos1-tokenPos-2);
				// log << ProID << "\n" << chr_string << "\n";
                digest_line(chr_string, pro2_1);

                tokenPos    = line.find('}', tokenPos1 + 1);
                chr_string = line.substr(tokenPos1+2, tokenPos-tokenPos1-2);
				// log << chr_string << "\n" << std::flush;
                digest_line(chr_string, pro2_2);
			}
        } 

		overlaps haploid_11x21, haploid_11x22, haploid_12x21, haploid_12x22;
		seg_overlap(pro1_1, pro2_1, haploid_11x21);
		seg_overlap(pro1_1, pro2_2, haploid_11x22);
		seg_overlap(pro1_2, pro2_1, haploid_12x21);
		seg_overlap(pro1_2, pro2_2, haploid_12x22);

		int n_seg = 0;
		int total_len = 0;

		//returns a dataframe with the following columns: simulNo, len_Shared_IBD, num_seg, mean_len, min_len, max_len
		//no printing to r console just do the dataframe
		overlaps haploid_11, haploid_12, haploid_21, haploid_22;
		check_HBD(haploid_11x21, haploid_11x22, haploid_11);
		check_HBD(haploid_12x21, haploid_12x22, haploid_12);
		check_HBD(haploid_11x21, haploid_12x21, haploid_21);
		check_HBD(haploid_12x22, haploid_11x22, haploid_22);
		// log << "simulation: " << j+1 << "\n" << std::flush;
		// Rcpp::Rcout << haploid_11.num_overlaps << " " << haploid_12.num_overlaps << " " << haploid_21.num_overlaps << " " << haploid_22.num_overlaps << "\n";
		if(haploid_11.num_overlaps + haploid_12.num_overlaps + haploid_21.num_overlaps + haploid_22.num_overlaps){
			line = "simulation: " + std::to_string(j+1);
			if(haploid_11.num_overlaps){
				line += "\npro. 1, chr. 1: ";
				// log << "pro. 1, chr. 1:\n";
				n_seg = n_seg + haploid_11.num_overlaps;
				for(int h=0; h < haploid_11.num_overlaps; ++h){
					line += std::to_string(haploid_11.Lpos[h]) + "->" + std::to_string(haploid_11.Rpos[h]) + "	";
					// log << haploid_11.Lpos[h] << "->" << haploid_11.Rpos[h] <<"\n";
					total_len = total_len + (haploid_11.Rpos[h] - haploid_11.Lpos[h]);
				}
			}
			if(haploid_12.num_overlaps){
				line += "\npro. 1, chr. 2: ";
				// log << "pro. 1, chr. 2:\n";
				n_seg = n_seg + haploid_12.num_overlaps;
				for(int h=0; h < haploid_12.num_overlaps; ++h){
					line += std::to_string(haploid_12.Lpos[h]) + "->" + std::to_string(haploid_12.Rpos[h]) +"	";
					// log << haploid_12.Lpos[h] << "->" << haploid_12.Rpos[h] <<"\n";
					total_len = total_len + (haploid_12.Rpos[h] - haploid_12.Lpos[h]);
				}
			}
			if(haploid_21.num_overlaps){
				line+= "\npro. 2, chr. 1: ";
				// log << "pro. 2, chr. 1:\n";
				n_seg = n_seg + haploid_21.num_overlaps;
				for(int h=0; h < haploid_21.num_overlaps; ++h){
					line += std::to_string(haploid_21.Lpos[h]) + "->" + std::to_string(haploid_21.Rpos[h]) +"	";					
					// log << haploid_21.Lpos[h] << "->" << haploid_21.Rpos[h] <<"\n";
					total_len = total_len + (haploid_21.Rpos[h] - haploid_21.Lpos[h]);			
				}
			}
			if(haploid_22.num_overlaps){
				line+= "\npro. 2, chr. 2: ";
				// log << "pro. 2, chr. 2:\n";
				n_seg = n_seg + haploid_22.num_overlaps;
				for(int h=0; h < haploid_22.num_overlaps; ++h){
					line += std::to_string(haploid_22.Lpos[h]) + "->" + std::to_string(haploid_22.Rpos[h]) +"	";					
					// log << haploid_22.Lpos[h] << "->" << haploid_22.Rpos[h] <<"\n";
					total_len = total_len + (haploid_22.Rpos[h] - haploid_22.Lpos[h]);
				}
			}
			Rcpp::message(Rcpp::wrap(line));	
		}			
		if (n_seg > 0){
			rvec1.push_back(j+1); //simulNo
			rvec2.push_back(n_seg);
			rvec3.push_back(100*(static_cast<float>(total_len))/(4*BP_len));
			rvec4.push_back(total_len/n_seg);
		}		
    }

	} catch(std::exception &ex) {
 	forward_exception_to_r(ex);
 	} catch(...){
 	::Rf_error("c++ exception (unknown reason)"); 
 	};
}

/*! 
	\brief Execute une simulation pour determiner les probabilites du passage d'un allele a une serie de proposant

	Calcule les probabilite qu'un proposants recoivent 0,1-2,2 allele a partir d'un groupe d'ancetre
	Calcule la probabilite conjointe que chaque proposant soit atteint
	Calcule la probatilite que de un..n proposant soit atteint
	Calcule la probabilie que chaque proposant soit atteint

	\param Genealogie	[in] Une genealogie construite  l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant  tudier
	\param plProEtat    [in] vecteur de taille lNProposant et representant l'etat a considerer pour chaque proposant
			<br>&nbsp; &nbsp;&nbsp; &nbsp;0: La condition est remplie si se proposant n'est pas malade 
			<br>&nbsp; &nbsp;&nbsp; &nbsp;1: La condition est remplie si le proposant recois 1-2 allele 
			<br>&nbsp; &nbsp;&nbsp; &nbsp;2: La condition est remplie si le proposant recois 2 allele 
	\param lNProposant	[in] Nombre d'lment du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant  tudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'lment du vecteur ancetre

	\param lSimul		[in] Nombre de simulation a effectuer

	\retval pdRetConj	[out] Un pointeur vers un double. 
					En cas de succes, le double represente la probabilite conjointe que la condition de chaque proposant soit remplis. 
	\retval pdRetSimul	[out] Un pointeur vers une vecteur de taille lNProposant.. 
					En cas de succes, Probabilite que la condition de chaque proposant soit remplis
	\retval pdRetProp	[out] Un pointeur vers une vecteur de taille lNProposant+1.
					En cas de succes represente la probabilite que 0,1,2..n condition soit remplis
									
	\param printprogress [in] imprime un message indiquant les progress accomplies

	\return 0 si la fonction est execut avec succs 
*/ 
int simul(int* Genealogie, int* plProposant, int* plProEtat,int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
		int lSimul, double* pdRetConj, double* pdRetSimul, double* pdRetProp, double* probRecomb, double probSurvieHomo,
		int printprogress)
{	
	try{
	//Validation genealogie
	if (lSimul<=0){
//		GENError("Number of simulation must be greater than zero");
		throw std::range_error("Number of simulation must be greater than zero");
		//GENError("Le nombre de simulation doit-tre suprieur  zero");
	}
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);

	//Creation des tableau
	INITGESTIONMEMOIRE;
	CIndSimul** Ordre		=(CIndSimul**) memalloc(lNIndividu,sizeof(CIndSimul*));	
	
	//Pour le sort spcial		
	int*		OrdreSaut	=(int*) memalloc(lNIndividu,sizeof(int*));				
	int NOrdre;
	
	//RECEUIL DES STATISTIQUES
	int *ProCompteur		=(int*)memalloc(lNProposant,sizeof(int*));
	int *NCompteur			=(int*)memalloc(lNProposant+1,sizeof(int*));
	int bConj;
	
	//initialisation
	//initrand();
	int i,j;
	int ap,am;
	//int *VecteurPosition=NULL;

	//INITIALISATION DE LA STRUCTURE DE NOEUD
	for(i=0;i<lNIndividu;i++)
	{
		Noeud[i].allele=0;	
		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].bFlagSort=0;	   				
	}
	
	//identifier et etiqueter les proposant
	for(i=0;i<lNProposant;i++)
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
					
	//identifier et etiqueter les points de departs
	for(i=0;i<lNAncetre;i++)
	{		
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele=interval(plAncEtat[i],0,2);		
	}

	//identifier et marque les noeuds utile et ceux inutile a la recherche
	for(i=0;i<lNAncetre;i++)
		ExploreArbre(NoeudAnc[i]);

	
	//creation d'un ordre d'execution et calcul des sauts
	PrepareSortPrioriteArbre(Noeud,lNIndividu);	
	NOrdre=0;
	memset(OrdreSaut,0,sizeof(int)*lNIndividu);
	for(i=0;i<lNAncetre;i++)		
		StartSortPrioriteArbre(NoeudAnc[i],Ordre,&NOrdre,OrdreSaut);
	
//	unsigned seed2 = time(0);
//	std::mt19937 gen(seed2);

//	boost::random_device rd;
//	std::random_device rd;		// **chgt IGES**
//	std::mt19937 gen(rd());
	
	#if defined _WIN32 || defined _WIN64
	  std::mt19937 gen(time(0));
	#else
//	  boost::random_device gen;
	  std::random_device gen;
	#endif
 	int nbannulee = 0;		// **chgt IGES**
 	int nbCasHomo = 0;		// **chgt IGES**

	//Simulation
	memset(ProCompteur,0,sizeof(int)*lNProposant);
	memset(NCompteur,0,sizeof(int)*(lNProposant+1));
	
	//CREATE_PROGRESS_BAR(lSimul,printprogress)
	for(int csimul=0;csimul<lSimul; csimul++)
	{
		//Par ordre du parent -> enfant
		//les 2 sorts n'ordonne pas dans le meme sens
		bool simAnnulee = false;
		for(int i=0;i<NOrdre;i++)
		{
			//croisement p/r au parent
			//double iRandom =urand();
			double iRandom = (double)gen()/(double)gen.max();
			//double iRandom = (double)rd()/(double)rd.max();		// pour le passage de l'allele.  **chgt IGES**
			
			if (Ordre[i]->pere!=NULL)	ap=	Ordre[i]->pere->allele;	// Ordre=vecteur de ptr vers les elts tries. allelePere (ap)
			else						ap=0;					// si pas de ptr, ap = 0.
			
			if (Ordre[i]->mere!=NULL)	am=	Ordre[i]->mere->allele;	// alleleMere (am) existe et on l'attribe
			else						am=0;					// existe pas donc am = 0.
			
			if (iRandom<TransGenCum[ap][am][0])
			{
				Ordre[i]->allele=0;
				//Saute un certain nombre de noeud
				j=i+OrdreSaut[i];
				while (i!=j)	Ordre[++i]->allele=0;
			}
			else
			{
				double alea    = (double)gen()/(double)gen.max();//(double)rd()/(double)rd.max(); // pour la recombinaison.**chgt IGES**
				int sex = Ordre[i]->sex;
				if(alea < probRecomb[1]) // tx femme (plus lev)
				{
				// si femme ou si prob < au taux male et inconnu (le plus petit)
				// si le morceau recombine (selon tx male et aussi sexe inconnu) OU si le morceau recombine (selon tx femelle)
				if( sex == GEN_FEM || alea < probRecomb[0] )
				{
					Ordre[i]->allele=0;
					//Saute un certain nombre de noeud
					j=i+OrdreSaut[i];
					while (i!=j)	Ordre[++i]->allele=0;
				}
				}
				else if (iRandom<TransGenCum[ap][am][1])
					Ordre[i]->allele=1;			
				else
				{
					Ordre[i]->allele=2;
					nbCasHomo++;
					//Descendance conditionnelle aux nombre d'alleles recus
					double alea2 = (double)gen()/(double)gen.max(); // 3e aleatoire pour determiner les cas homo
					if(alea2 > probSurvieHomo){ simAnnulee = true; break; }
				}
			}
		}
		if(!simAnnulee) {
			//Modification des compteurs
			bConj=0;
			for(int i=0;i<lNProposant;i++)
			{
//				if ( (plProEtat[i]==0 && NoeudPro[i]->allele==0) || (plProEtat[i]!=0 && NoeudPro[i]->allele>=plProEtat[i]))
				if ( (plProEtat[i]==0 && NoeudPro[i]->allele==0) || 
					(plProEtat[i]==1 && NoeudPro[i]->allele==1) ||
					(plProEtat[i]==2 && NoeudPro[i]->allele==2) ||
					(plProEtat[i]==3 && NoeudPro[i]->allele>=1) )
				{
					++ProCompteur[i];
					++bConj;
				}
			}
			++NCompteur[bConj];
			//Barre de progress
			//INCREMENT_PROGRESS_BAR()
		}
		else { // la simulation est annulee pcq la sim ne concorde pas avec la genealogie
			csimul--;
			nbannulee++;
		}
	}

	//ECRITURE DE LA VALEUR DE RETOUR
	//double* pdRetConj,double* pdRetSimul,double* pdRetProp
	for(int i=0;i<lNProposant;i++)
	{
		//printf("%f\n", double(ProCompteur[i]));
		pdRetSimul[i]=double(ProCompteur[i])/double(lSimul);
		pdRetProp[i]=double(NCompteur[i])/double(lSimul);
	}
	pdRetProp[lNProposant]=double(NCompteur[lNProposant])/double(lSimul);
	*pdRetConj=double(NCompteur[lNProposant])/double(lSimul);
	
	//FIN
	//outrand();
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}




/*! 
	\brief Execute une ou plusieurs simulation et retourne le nombre d'allele transmit a chaque proposant (pour chaque simulation)

	Calcule un etat possible pour chaque proposant en tenant compte de chaque ancetre et son etat

	\param Genealogie	[in] Une genealogie construite  l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant  tudier
	\param lNProposant	[in] Nombre d'lment du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant  tudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'lment du vecteur ancetre
	
	\param lSimul		[in] Nombre de simulation  effectuer

	\retval pdRetour	[out] Pointeur vers un vecteur de NProposant x lSimul
								En cas de succes, ce vecteur le nombre d'allele assigne a chaque proposant pour la simulation
	
	\param printprogress imprime un message indiquant les progress accomplies

	\return 0 si la fonction est execut avec succs
*/
int simulsingle(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre, int lSimul,
			 double* pdRetour,int printprogress)
{
	try{
	//VALIDATION GENEALOGIE
	if (lSimul<=0){
//		GENError("Number of simulation must be greater than zero");
		throw std::range_error("Number of simulation must be greater than zero");
		//GENError("Le nombre de simulation doit-tre suprieur  zero");
	}
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	/**D***/
	//Creation des tableau
	INITGESTIONMEMOIRE;
	CIndSimul** Ordre		=(CIndSimul**) memalloc(lNIndividu,sizeof(CIndSimul*));	
	
	//Pour le sort spcial		
	int*		OrdreSaut	=(int*) memalloc(lNIndividu,sizeof(int*));				
	int NOrdre;
	/**F***/
	//CREATION DES TABLEAU       
	int i;
	//int *VecteurPosition=NULL;

	/**D***/
	//INITIALISATION DE LA STRUCTURE DE NOEUD
	for(i=0;i<lNIndividu;i++)
	{
		Noeud[i].allele=0;	
		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].bFlagSort=0;	   				
	}
	//INITIALISATION DE LA STRUCTURE DE NOEUD
//	for(i=0;i<lNIndividu;i++)
//		Noeud[i].allele=0;			
	/**F***/
	
	//IDENTIFIER ET ETIQUETER LES PROPOSANT
	for(i=0;i<lNProposant;i++)
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;					

	//IDENTIFIER ET ETIQUETER LES POINTS DE DEPARTS
	for(i=0;i<lNAncetre;i++)
	{
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele=interval(plAncEtat[i],0,2);		
	}

	/**D***/
	//identifier et marque les noeuds utile et ceux inutile a la recherche
	for(i=0;i<lNAncetre;i++)
		ExploreArbre(NoeudAnc[i]);

	
	//creation d'un ordre d'execution et calcul des sauts
	PrepareSortPrioriteArbre(Noeud,lNIndividu);	
	NOrdre=0;
	memset(OrdreSaut,0,sizeof(int)*lNIndividu);
	for(i=0;i<lNAncetre;i++)		
		StartSortPrioriteArbre(NoeudAnc[i],Ordre,&NOrdre,OrdreSaut);
	/**F***/

	//INITIALISATION
	//initrand();
//	unsigned seed2 = time(0);
//	std::mt19937 gen(seed2);
//	boost::random_device rd;
//	std::random_device rd;		// **chgt IGES**
//	std::mt19937 gen(rd());
	#if defined _WIN32 || defined _WIN64
	  std::mt19937 gen(time(0));
	#else
	std::random_device gen;		// **chgt IGES**
//	  boost::random_device gen;
	#endif
	int ap,am;
		
	//Partie 3: SIMULATION
	//CREATE_PROGRESS_BAR(lSimul,printprogress)
	for(int csimul=0;csimul<lSimul; csimul++)
	{
/* Utilisant la facon de faire de la fonction simul */
		for(int i=0;i<NOrdre;i++)
		{
			if (Ordre[i]->pere!=NULL) ap=	Ordre[i]->pere->allele;
       		else					ap=0;
				
			if (Ordre[i]->mere!=NULL) am=	Ordre[i]->mere->allele;
       		else					am=0;
						
			if (ap==0 && am==0)		Ordre[i]->allele=0;
			else {
				double iRandom = (double)gen()/(double)gen.max();
				if (iRandom<TransGenCum[ap][am][0])		Ordre[i]->allele=0;
				else
					if (iRandom<TransGenCum[ap][am][1]) Ordre[i]->allele=1;			
					else							 Ordre[i]->allele=2;
			}
		}
	
/* Fin de la nouvelle facon de faire */ 

/* Code remplace par la facon de faire de la fonction simul
		for(i=0;i<lNIndividu;i++)
		{	  
			if (Noeud[i].etat!=GENDEPART)
			{
				if (Noeud[i].pere!=NULL) ap=	Noeud[i].pere->allele;
	       		else					 ap=0;
				
				if (Noeud[i].mere!=NULL) am=	Noeud[i].mere->allele;
	       		else					 am=0;
						
				if (ap==0 && am==0)		Noeud[i].allele=0;
				else {
					//double iRandom =urand();
					//double iRandom = (double)gen()/(double)rd.max();
					double iRandom = (double)gen()/(double)gen.max();//(double)rd()/(double)rd.max();		// pour le passage de l'allele.  **chgt IGES**
					if (iRandom<TransGenCum[ap][am][0])		Noeud[i].allele=0;
					else
						if (iRandom<TransGenCum[ap][am][1]) Noeud[i].allele=1;			
						else								Noeud[i].allele=2;
				} //fin ap==am==0
			}//fin if noeud depart
		}// fin for pour chaque noeud
** Fin du code remplace */		

        //Modification des compteurs
		const int tmp = csimul*lNProposant;
		for(i=0;i<lNProposant;i++)
			pdRetour[tmp+i]=NoeudPro[i]->allele;		

		//Barre de progress
		//INCREMENT_PROGRESS_BAR();

	}// Fin du pour chaque simulation
	///nettoyer
	//outrand();
	return 0;
 } catch(std::exception &ex) {
 	forward_exception_to_r(ex);
 } catch(...){
 	::Rf_error("c++ exception (unknown reason)"); 
 } 
 return 0;
}
/*! 
	\brief Execute une plusieur simulation et retourne le nombre d'allele transmit a chaque proposant (pour chaque simulation)

	Calcule un etat possible pour chaque proposant en tenant en comple chaque ancetre et son etat

	\param Genealogie	[in] Une genealogie construite  l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant  tudier
	\param lNProposant	[in] Nombre d'lment du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant  tudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'lment du vecteur ancetre
	
	\param lSimul		[in] Nombre de simulation  effectuer
	\param mtProb		[in] Matrice S-PLUS: tableau des probabilits

	\retval pdRetour	[out] Pointeur vers un vecteur de NProposant x lSimul
								En cas de succes, ce vecteur le nombre d'allele assigne a chaque proposant pour la simulation
	
	\param printprogress imprime un message indiquant les progress accomplies

	\return 0 si la fonction est execut avec succs
*/

SEXP simulsingleProb(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre, int* plAncEtat, SEXP mtProb,
				 int lSimul, int printprogress)
{
	try{
	//Conversion des paramtres
	Rcpp::NumericMatrix matprob(mtProb);

	//VALIDATION GENEALOGIE
	if (lSimul<=0){
//		GENError("Number of simulation must be greater than zero");
		throw std::range_error("Number of simulation must be greater than zero");
		//GENError("Le nombre de simulations doit-tre suprieur  zero");
	}
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	//CREATION DES TABLEAUX       
	int i;
	//long *VecteurPosition=NULL;

	//INITIALISATION DE LA STRUCTURE DE NOEUD
	for(i=0;i<lNIndividu;i++)
		Noeud[i].allele=0;			
	
	//IDENTIFIER ET ETIQUETER LES PROPOSANT
	for(i=0;i<lNProposant;i++)
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;					

	//IDENTIFIER ET ETIQUETER LES POINTS DE DEPARTS
	for(i=0;i<lNAncetre;i++)
	{
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele=interval(plAncEtat[i],0,2);		
	}

	//INITIALISATION
	//initrand();
//	unsigned seed2 = time(0);
//	std::mt19937 gen(seed2);
//	boost::random_device rd;
//	std::random_device rd;		// **chgt IGES**
//	std::mt19937 gen(rd());
	#if defined _WIN32 || defined _WIN64
	  std::mt19937 gen(time(0));
	#else
//	  boost::random_device gen;
	  std::random_device gen;
	#endif
	int ap,am;
		
	//Pointeur de retour
	//CSPnumeric tmp(lNProposant*lSimul);
	Rcpp::IntegerVector tmp(lNProposant*lSimul);
	
	//Partie 3: SIMULATION
	//CREATE_PROGRESS_BAR(lSimul,printprogress)
	for(int csimul=0;csimul<lSimul; csimul++)
	{
		for(i=0;i<lNIndividu;i++)
		{	  
			if (Noeud[i].etat!=GENDEPART)
			{	
				if (Noeud[i].pere!=NULL)	ap=	Noeud[i].pere->allele;
	       		else						ap=0;
				
				if (Noeud[i].mere!=NULL)	am=Noeud[i].mere->allele;
	       		else						am=0;
						
				if (ap==0 && am==0)								Noeud[i].allele=0;
				else {
//					am += 1;
//					ap += 1;
					if (Noeud[i].sex == GEN_FEM) am += 6;
//					double iRandom =urand(); 											
					//double iRandom = (double)gen()/(double)rd.max();
					double iRandom = (double)gen()/(double)gen.max();//(double)rd()/(double)rd.max();		// pour le passage de l'allele.  **chgt IGES**
					if (iRandom<(double)matprob(ap,am))			Noeud[i].allele=0;		       
					else
						if (iRandom<(double)matprob(ap,am+=3))	Noeud[i].allele=1;			
						else									Noeud[i].allele=2;
				} //fin ap==am==0
			}//fin if noeud depart
		}// fin for pour chaque noeud
	
        //Modification des compteurs
		for(i=0;i<lNProposant;i++)
			tmp(csimul*lNProposant+i)=(int)NoeudPro[i]->allele;	 //csimul*lNProposant+i+1
				
		//Barre de progress
		//INCREMENT_PROGRESS_BAR();

	}// Fin du pour chaque simulation
	///nettoyer
	//outrand();
	return tmp; //.Detach();
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}

/*! 
	\brief Execute une plusieur simulation et retourne le nombre d'allele transmit a chaque proposant en tableau
			de frquences pour toutes les simulations.

	Calcule un etat possible pour chaque proposant en tenant en compte chaque ancetre et son etat

	\param Genealogie	[in] Une genealogie construite  l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant  tudier
	\param lNProposant	[in] Nombre d'lment du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant  tudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'lment du vecteur ancetre
	
	\param lSimul		[in] Nombre de simulations  effectuer

	\retval pdRetour	[out] Pointeur vers un vecteur de NProposant x 0-1-2 pour la frquence d'allles tramises pour toutes les
							  simulations. 
	
	\param printprogress imprime un message indiquant les progress accomplies

	\return 0 si la fonction est execut avec succs
*/
int simulsingleFreq(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
				int lSimul, double* pdRetour,int printprogress)
{
	try{
	//VALIDATION GENEALOGIE
	if (lSimul<=0) {
//		GENError("Number of simulation must be greater than zero");
		throw std::range_error("Number of simulation must be greater than zero");
		//GENError("Le nombre de simulation doit-tre suprieur  zero");
	}
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	//CREATION DES TABLEAU    
	int i;
	//int *VecteurPosition=NULL;

	//INITIALISATION DE LA STRUCTURE DE NOEUD
	for(i=0;i<lNIndividu;i++)
		Noeud[i].allele=0;		
	
	//IDENTIFIER ET ETIQUETER LES PROPOSANT
	for(i=0;i<lNProposant;i++)
		NoeudPro[i]->etat = GENPROPOSANTINUTILE;					

	//IDENTIFIER ET ETIQUETER LES POINTS DE DEPARTS
	for(i=0;i<lNAncetre;i++)
	{
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele=interval(plAncEtat[i],0,2);		
	}

	//Indice de tableau pour la variable de sortie
	const int tmp0 = 0*lNProposant; //Frquence de 0 allle
	const int tmp1 = 1*lNProposant; //Frquence de 1 allle
	const int tmp2 = 2*lNProposant; //Frquence de 2 allle

	//INITIALISATION
	//initrand(); 
//	unsigned seed2 = time(0);
//	std::mt19937 gen(seed2);
//	boost::random_device rd;
//	std::random_device rd;		// **chgt IGES**
//	std::mt19937 gen(rd());
	#if defined _WIN32 || defined _WIN64
	  std::mt19937 gen(time(0));
	#else
//	  boost::random_device gen;
	  std::random_device gen;
	#endif
	int ap,am;
		
	//Partie 3: SIMULATION
	//CREATE_PROGRESS_BAR(lSimul,printprogress)
	for(int csimul=0;csimul<lSimul; csimul++)
	{
		for(i =0;i<lNIndividu;i++)
		{	  
			if (Noeud[i].etat!=GENDEPART)
			{
				if (Noeud[i].pere!=NULL) ap = Noeud[i].pere->allele;
	       		else					 ap=0;
				
				if (Noeud[i].mere!=NULL) am = Noeud[i].mere->allele;
	       		else					 am=0;
						
				if (ap==0 && am==0)							Noeud[i].allele=0;
				else {
					//double iRandom = urand();
					//double iRandom = (double)gen()/(double)rd.max();
					double iRandom = (double)gen()/(double)gen.max();//(double)rd()/(double)rd.max();		// pour le passage de l'allele.  **chgt IGES**
					if (iRandom<TransGenCum[ap][am][0])		Noeud[i].allele=0;		       
					else
						if (iRandom<TransGenCum[ap][am][1])	Noeud[i].allele=1;			
						else								Noeud[i].allele=2;
				} //fin ap==am==0
			}//fin if noeud depart
		}// fin for pour chaque noeud

		for(i=0;i<lNProposant;i++)
		{
			if (NoeudPro[i]->allele == 0)		pdRetour[tmp0+i]+= 1;
			else if (NoeudPro[i]->allele == 1)	pdRetour[tmp1+i]+= 1;
			else								pdRetour[tmp2+i]+= 1;
		}
		//Barre de progress
		//INCREMENT_PROGRESS_BAR();
	}// Fin du pour chaque simulation
	///nettoyer
	//outrand();
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}
/*! 
	\brief Execute une plusieur simulation et retourne le nombre d'allele transmit a chaque proposant (pour chaque simulation)

	Calcule un etat possible pour chaque proposant en tenant en comple chaque ancetre et son etat

	\param Genealogie	[in] Une genealogie construite  l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant  tudier
	\param lNProposant	[in] Nombre d'lment du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant  tudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'lment du vecteur ancetre
	
	\param lSimul		[in] Nombre de simulation  effectuer

	\retval pdRetour	[out] Pointeur vers un vecteur de NProposant x lSimul
								En cas de succes, ce vecteur le nombre d'allele assigne a chaque proposant pour la simulation
	
	\param printprogress imprime un message indiquant les progress accomplies

	\return 0 si la fonction est execut avec succs
*/
SEXP simulsingleFct(int * Genealogie, int * proposant, int lproposant, int * plAncetre, int * plAncEtatAll1, int * plAncEtatAll2, int lNAncetre, int lSimul, SEXP SfctSousGrp, int printprogress)
{	
	try{
	//VALIDATION GENEALOGIE
	if (lSimul<=0){
//		GENError("Number of simulation must be greater than zero");
		throw std::range_error("Number of simulation must be greater than zero");
		//GENError("Le nombre de simulation doit-tre suprieur  zero");
	}
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul ** NoeudPro=NULL; //[i]
	LoadProposant(proposant,lproposant,&NoeudPro); //[i]

	//CREATION D'UN VECTEUR D'ANCETRES
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);

	//CREATION DES TABLEAU       
	//long *VecteurPosition=NULL;

	//INITIALISATION DE LA STRUCTURE DE NOEUD
	for(int i=0;i<lNIndividu;i++)
	{
		Noeud[i].allele2Pos[0] = 0;
		Noeud[i].allele2Pos[1] = 0;
	}

	//IDENTIFIER ET ETIQUETER LES POINTS DE DEPARTS
	for(int i=0;i<lNAncetre;i++)
	{
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele2Pos[0]=interval(plAncEtatAll1[i],0,5);	
		NoeudAnc[i]->allele2Pos[1]=interval(plAncEtatAll2[i],0,5);		
	}	

	//INITIALISATION
//	initrand();
	int lptrAp[2];
	int lptrAm[2];

	//Partie 3: SIMULATION	
	//Dclaration de la liste de rsultats de la fonction de l'utilisation pour chaque simulation
	//CSPlist resultFct;
	Rcpp::List resultFct; //(lSimul);
	Rcpp::Function f(SfctSousGrp);
	//CREATE_PROGRESS_BAR(lSimul,printprogress)

	for(int csimul=0;csimul<lSimul; csimul++)
	{
		for(int i=0;i<lNIndividu;i++)
		{	  	
			if (Noeud[i].etat!=GENDEPART)
			{	
				if (Noeud[i].pere!=NULL)
				{
					lptrAp[0] = Noeud[i].pere->allele2Pos[0];
					lptrAp[1] = Noeud[i].pere->allele2Pos[1];
				}
	       		else
					lptrAp[0] = lptrAp[1] = 0;
				if (Noeud[i].mere!=NULL)
				{
					 lptrAm[0] = Noeud[i].mere->allele2Pos[0];
					 lptrAm[1] = Noeud[i].mere->allele2Pos[1];
				}
	       		else
					lptrAm[0] = lptrAm[1] = 0;
						
				int iRandom = irand(0,1);  // utilise random_device a travers irand()
				Noeud[i].allele2Pos[0]=lptrAp[iRandom];
				iRandom = irand(0,1);
				Noeud[i].allele2Pos[1]=lptrAm[iRandom];

			}//fin if noeud depart

		}// fin for pour chaque noeud	
		/*
		//Dclaration de la liste de rsultats
		//CSPlist resultGrp;
		Rcpp::List resultGrp(igrp);
		//Pour chaque groupe
		for (int i=0;i<igrp;i++)
		{	
			SEXP grp = indSelect[i];//extraction des individus du groupe

			//CSPnumeric lind(grp);
			Rcpp::IntegerVector lind(grp);
			
			//CSPnumericMatrix ans(lind.length(), 2); //output matrix.
			Rcpp::IntegerMatrix ans(lind.size(), 2);
			
			//ans.SetRowNames(CSPcharacter(grp));
			ans.attr("rownames") = Rcpp::as<Rcpp::CharacterVector>(grp);
			
			//CSPnumeric tmp(lind.length());
			Rcpp::IntegerVector tmp (lind.size());
			
			//CSPnumeric tmp2(lind.length());
			Rcpp::IntegerVector tmp2 (lind.size());
			
			for(int j=0;j<lind.length();j++)//Pour chaque individus du groupe
			{
				 // **tmp(j+1) = NoeudPro[i][j]->allele2Pos[0];	//Rsultat dans un tableau	
				 ans(j+1, 0) = NoeudPro[i][j]->allele2Pos[0];
				 // **tmp2(j+1) = NoeudPro[i][j]->allele2Pos[1];	//Rsultat dans un tableau	
				 ans(j+1, 1) = NoeudPro[i][j]->allele2Pos[1];
			}
			// **ans.SetJthColumnDirect(0, tmp.Detach());
			
			// **ans.SetJthColumnDirect(1, tmp2.Detach());
			//En liste
			// **resultGrp.Add(ans);
			resultGrp.push_back(ans);
		}
		// **CSPfunction f(SfctSousGrp); //Call the S function.
		//Rcpp::Function f(SfctSousGrp);
		//En liste de rsultats de la fct de l'utilisateur
		// **resultFct.Add(f.Eval(resultGrp));
		resultFct.push_back(f(resultGrp));
		*/
		
		Rcpp::IntegerMatrix ans(lproposant, 2);

		Rcpp::CharacterVector rowNames(lproposant);
		for(int i=0; i<lproposant; i++) { char nomLigne [10]; /*int n = */snprintf(nomLigne, 10, "%d", proposant[i]); rowNames[i] = nomLigne; }
		//for(int i=0; i<lproposant; i++) sprintf(rowNames[i], "%d", proposant[i]); // itoa(proposant[i], &rowNames[i], 10)

		Rcpp::List dimnms = Rcpp::List::create( rowNames, //Rcpp::as<Rcpp::CharacterVector>(indSelect[i]),
										Rcpp::CharacterVector::create("1", "2"));
		ans.attr("dimnames") = dimnms;
		
		for(int j=0;j<lproposant;j++)//Pour chaque individus du groupe
		{
			 ans(j, 0) = NoeudPro[j]->allele2Pos[0]; //[i] j+1
			 ans(j, 1) = NoeudPro[j]->allele2Pos[1]; //[i] j+1
		}
		resultFct.push_back(f(ans));
		//Barre de progress
		//INCREMENT_PROGRESS_BAR();

	}// Fin du pour chaque simulation

	///nettoyer
	//outrand();
	//Retourne la liste de rsultat pour chaque simulation de la fct de l'utilisateur
	return Rcpp::wrap(resultFct); //.Detach();
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}

// *************************************************************************
//	CALCUL DE LA PROBABILITE EXACTE
// ************************************************************************* 
//const int PBARINTERVAL_PROB=9;
//100 devrais tre plus que suffisant pour occup un ordinateur moderne un bon bout de temps
//const int PROB_MAXIMUM_NORDRE = 100; //646 Taille maximale d'un double (utiliser numeric_limits<double>::max() ? )

/*! 
	\brief Evalue la probabilite exacte de transfert d'un allele a une serie de proposant

	Cette fonction dtermine la probabilit conjointe qu'un ou plusieurs proposant reoivent  un certain nombre d'allle malade.
	En prenant pour acquis qu'un ou plusieurs anctre possde un ou deux d'allle malade. 
	Pour ce faire, cette fonction calcule la valeur de toute les branches de manire  obtenir la probabilit EXACT. 
	Le temps de calcul peut-tre trs prohibitif.

	\param Genealogie	[in] Une genealogie construite  l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant  tudier
	\param plProEtat    [in] vecteur de taille lNProposant et representant l'etat a considerer pour chaque proposant
			<br>&nbsp; &nbsp;&nbsp; &nbsp;0: La condition est remplie si se proposant n'est pas malade 
			<br>&nbsp; &nbsp;&nbsp; &nbsp;1: La condition est remplie si le proposant recois 1-2 allele 
			<br>&nbsp; &nbsp;&nbsp; &nbsp;2: La condition est remplie si le proposant recois 2 allele 
	\param lNProposant	[in] Nombre d'lment du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant  tudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'lment du vecteur ancetre

	\param OrdreMaximum [in] Le nombre de noeud touche au maximum, si le nombre de noeud touche est trop grand alors la procedure s'interromp automatiquement

	\retval pdRetConj	[out] Un pointeur vers un double. 
								En cas de succes, le double represente la probabilite conjointe que la condition de chaque proposant soit remplis. 
	\retval pdRetSimul	[out] Un pointeur vers une vecteur de taille lNProposant.. 
								En cas de succes, Probabilite que la condition de chaque proposant soit remplis
							
	\return 0 si la fonction est execut avec succs
*/
SEXP prob( int* Genealogie, int* plProposant, int* plProEtat,int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
		double* pdRetConj,double* pdRetSimul,int printprogress,int onlyConj)
{
	try{
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	if(LoadProposant(plProposant,lNProposant,&NoeudPro) == -1) return Rcpp::wrap(-1);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	//Creation du tableau de Noeud, d'ancetre et proposant
	INITGESTIONMEMOIRE;			
	int NOrdre;
	CIndSimul** Ordre	=(CIndSimul**) memalloc(lNIndividu,sizeof(CIndSimul*));		
	int   *OrdreSaut	= (int*)	memalloc(lNIndividu,sizeof(int));	
	
	//receuil des statistiques
	int i;
	int ap,am;

	//creation et initialisation de la structure noeud et l'ordre
	for(i=0;i<lNIndividu;i++)
	{
		Noeud[i].allele=0;
		
		Noeud[i].prob[0]=0.;
		Noeud[i].prob[1]=0.;
		Noeud[i].prob[2]=0.;			
		Noeud[i].iind=-1; //Pas un proposant

		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].bFlagSort=0;	
	}
	
	//identifier et etiqueter les proposant
	for(i=0;i<lNProposant;i++)
	{
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
		NoeudPro[i]->iind=interval(plProEtat[i],0,2);
	}

	//identifier et etiqueter les points de departs
	for(i=0;i<lNAncetre;i++)
	{
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele=interval(plAncEtat[i],0,2);
	}

	//identifier et marque les noeuds utile et ceux inutile a la recherche
	for(i=0;i<lNAncetre;i++)
		ExploreArbre(NoeudAnc[i]);
	
	//creation d'un ordre d'execution et calcul des sauts
	PrepareSortPrioriteArbre(Noeud,lNIndividu);
	NOrdre=0;
	memset(OrdreSaut,0,sizeof(int)*lNIndividu);
	for(i=0;i<lNAncetre;i++)		
		StartSortPrioriteArbre(NoeudAnc[i],Ordre,&NOrdre,OrdreSaut);

	if (NOrdre==-1){		
//		GENError("There is no link between any ancetres and any probands");
		throw std::range_error("There is no link between any ancetres and any probands");
		//GENError("Il n'y a pas de lien entre aucun des ancetres et aucun proposant");
	}
	//Liste de tableau dpendante de l'ordre
	double *Cumul		= (double*) memalloc(NOrdre+1,sizeof(double));

	// Valeur cumulative courante
	//double *PrCumul	  = (double*) memalloc(NOrdre+1,sizeof(double));	
	//Valeur actuelle pour chaque position pour chaque allele
	//double (*PrValue)[3] = (double (*)[3]) memalloc((NOrdre+1),sizeof(double[3]));	
	
	// Debut du calcul
	int n=0;
	Cumul[0]=1;
	int iteration=0;
	const int Lastindex=NOrdre-1;
	double ConjProb=0;
	
	//Initialisation de dbut de simulation
	for(i=0;i<NOrdre;i++)	
	{
		Ordre[i]->allele=-1;		
	}
	
	//On rutilise la variable bFlagSort pour indique si la condition d'un proposant
	//est acceptable ou non (1 oui, non =0)
	int nbCritereValid=0;
	for(i=0;i<lNProposant;i++)
	{
		NoeudPro[i]->bFlagSort=0; //0: non satisfait
		if (NoeudPro[i]->etat==GENPROPOSANTINUTILE)
		{
			NoeudPro[i]->prob[0]=1.0; //Si c'est un proposant inutile fait une correction
			if (NoeudPro[i]->iind==0) 
			{
				NoeudPro[i]->bFlagSort=1;
				++nbCritereValid;
			}
		}
	}
	
	//Verifie que la dure d'execution ne sera pas trop intue
//	if (NOrdre > PROB_MAXIMUM_NORDRE){
//		GENError("Execution time is too great to launch this function call");
//		throw std::exception();
		//GENError("La dure d'excution de cet appel de fonction serait surement comparable  l'ge de l'univers");
//	}
	//PROGRESS BAR	
	//Construction du vecteur de valeur
/*	const int PBar_position = MAX(NOrdre - PBARINTERVAL_PROB,0);	
	const double maximumiteration = pow(3,PBar_position+1);

	PrCumul[0]=0.0;
	PrCumul+=1; //Pour simplifier
	for(i=0;i<=PBar_position;i++)	
	{			
		PrCumul[i]=0;
		//Calcul de valeur
		const double val= pow(3,PBar_position-i);
		PrValue[i][0]=0.0;
		PrValue[i][1]=val;
		PrValue[i][2]=val+val;
	}
	CREATE_FLOAT_PROGRESS_BAR(maximumiteration,&PrCumul[PBar_position],printprogress)
*/

#ifdef USESDEBUG
	//Informatin de dbuggage
	SXLONG iter =0;
	//printf("\nThe order is = %d\n",NOrdre);
	//printf("Maximum number of iterations used:%.15G\n",pow(3,NOrdre)/2);
	//printf("\nTaskbar order  = %d\n",PBar_position);
	//printf("Progress bar iteration number:%.15G\n",maximumiteration);
#endif

	while(n>=0)
	{
			//n niveau de l'ordre qui est en traitement
			//a nombre d'allele du noeud courant

			//pour conomiser du code
			CIndSimul& nd = *Ordre[n];
			int &a = nd.allele;

			//test
			iteration++;
			++a;			
			if (a==3)
			{
				a=-1;
				n--;
				/*if (n==PBar_position){ INCREMENT_FLOAT_PROGRESS_BAR();	}*/
				#ifdef USESDEBUG
					++iter;
				#endif
			}
			else {
				//Calcul pour la progressbar
				//if (n<=PBar_position) PrCumul[n]=PrCumul[n-1]+PrValue[n][a];					

				//Trouve le nombre d'allele des parents
				if (nd.pere!=NULL)	ap = nd.pere->allele;
				else				ap = 0;
				
				if (nd.mere!=NULL)	am=	nd.mere->allele;
				else				am=0;
				
				//ebranchage (C'est ce qui empeche de calcule atteint)
				if (TransGen[ap][am][a]!=0.0) {
					//Probabilit d'tre dans l'tat courant					
					Cumul[n+1]=Cumul[n]*TransGen[ap][am][a];
					
					//Si on dsirer la probabili individuel dcoch en dessous		
					if (nd.etat==GENPROPOSANT){
						nd.prob[a]+=Cumul[n+1];							

						//Evaluation des criteres pour probabilit conjointe
						if (nd.bFlagSort==1) {
							--nbCritereValid;
							nd.bFlagSort=0;
						}
						if ( ((nd.iind==0 && a==0) || (nd.iind!=0 && a>=nd.iind)))
						{
							++nbCritereValid;
							nd.bFlagSort=1;							
						}						

						//Si tous les proposants sont considr comme valide alors...
						//Compet dans la propabilit conjointe
						if (n==Lastindex && nbCritereValid==lNProposant)												
							ConjProb+=Cumul[n+1];

						//Passe  la sous-tape suivante	
						if (n!=Lastindex && (!onlyConj || (onlyConj && nd.bFlagSort==1) ) )	++n;
					}//fin si proposant
					else
					{
						//Passe  la sous-tape suivante	
						if (n!=Lastindex)	++n;
					}//fin else proposant
				} //Fin branchage
			}//Fin de la boucle pour iteration valide
	}//Fin de la boucle principale
	//END_FLOAT_PROGRESS_BAR();

#ifdef USESDEBUG	
	//printf("\nprogressbar value: %.15G\n",PrCumul[PBar_position]);
	//printf("Number of iterations executed:%I64d\n",iter);
#endif
	
	//CALCUL DE LA PROBABILITE CONJOINTE
	*pdRetConj=ConjProb;
		
	//CALCUL DE LA PROBABILITE INDIVIDUEL (OBSOLETE) OPTIONNELLE
	for(i=0;i<lNProposant;i++)
	{
		switch(plProEtat[i])
		{
		case 0:
			pdRetSimul[i]=NoeudPro[i]->prob[0];break; 
		case 1:
			pdRetSimul[i]=NoeudPro[i]->prob[1]+NoeudPro[i]->prob[2];break;
		case 2:
			pdRetSimul[i]=NoeudPro[i]->prob[2];break;
		}
	}	
	return Rcpp::wrap(0);
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			}catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}


// *************************************************************************
//	APPARENTEMENT MODIFIER 
// ************************************************************************* 

//Facteur de charge maximale de la table de hashage utilise par CoefApparentement
/**
	si ex= 0.8 alors la table de hashage sera 20% plus grande que le nombre maximale de possibilite
	Avec ce facteur on peut se permettre d'echange memoire/performance..

  Mais dans tous les cas COAPPMOD_HASHAGECHARGEMAXIMAL ne devrais pas etre egale a 1 sinon sa vas etre extremement lent
 */
const double COAPPMOD_HASHAGECHARGEMAXIMAL=.8;

///Mmoire utilisable au grand maximum en octet;
const long double COAPPMOD_MAXMEMORY_USABLE=(long double) 4000*MEGAOCTET; //800*MEGAOCTET;
///Nombre maximal de possibilit autoris par l'algorithme...
//COAPPMOD_NBPOSSMAX_DUPP doit-tre plus petit qu'un unsigned long (d'une bonne marge)
const long double COAPPMOD_NBPOSSMAX_DUPP=(long double)ULONG_MAX*COAPPMOD_HASHAGECHARGEMAXIMAL*.7;  //.7 = Facteur scurit
//COAPPMOD_NBPOSSMAX_DUPP doit-tre plus petit qu'un XLONG (pour la progresse bar et autre)
const long double COAPPMOD_NBPOSSMAX_NODUPP=(long double)1E16;  //SERIEUSEMENT LONG SUR QUOI QUE CE SOIT  ( avant 1E12 ) **JFL**


/*! 
	\brief Calcul de l'apparentement modifier entre un ancetre et quelques proposants

	Cette fonction calcule l'apparentement modifier d'entre un anctre et n proposant. 
	Plusieurs calcul  partir d'anctre diffrent vers le mme nombre n de proposant peuvent-tre fait lors du mme appel de fonction.

	\param Genealogie	[in] Une genealogie construite  l'aide de gen.genealogie 

	\param plInput	 [in] Vecteur representant une matrice d'entier (Selon Formule: col * NLigne + ligne)
						Matrice d'entier de n colonne,  sur chaque ligne, on retrouve n proposant  tudier p/r  un anctre. 
						<br>Chaque ligne reprsente un calcul compltement distinct
						<br>Le nombre de ligne de la matrice doit tre gale au nombre d'lment du vecteur anctre pour faire la correspondance.
						<br>Ex :  
							<br>ancetre =  Vecteur 10, 20
							<br>
							<table>
								<tr>
								<td><B>10</B></td>
								<td>1</td>
								<td>2</td>
								<td>3</td>
								</tr>
								<tr>
								<td><B>20</B></td>
								<td>4</td>
								<td>5</td>
								<td>6</td>
								</tr>
							</table>
						<br>Dans ce cas, la fonction calculera.
							<br>&nbsp;&nbsp;1- 	L'apparentement modifi des proposant 1, 2, 3 avec l'anctre 10
							<br>&nbsp;&nbsp;2- 	Et l'apparentement modifi des proposant 4, 5, 6 avec l'anctre 20
						<br><br>La valeur de retour sera un vecteur avec les deux rsultats prcdents.

	\param lNColonne [in] Nombre de Colone de la matrice (Incluant la 1e colonne des ancetres)
	\param lNLigne	 [in] Nombre de ligne de la matrice
  

	\retval pdRetVecAppMod	[out] Un pointeur vers un vecteur de double de taillel lNLigne
							En cas de success, contient la valeur de l'apparentement modifier pour chaque ligne de la matrice
		 		
	\param DuppDetection	[in] Si !=0 la detection de dupplicata dans les chemins sera activer

	\param MaxMoTableHash	[in] Nombre max de Mo que la table de hash peut utiliser avant de declenche une erreur
	  
	\param Maxcombinaison	[in] Nombre maximum de combinaision a tester avant d'imprimer un message d'erreur

	\param printprogress	[in] imprime un message indiquant les progress accomplies

	\return 0 si la fonction est execut avec succs 
	
	 \sa CoPro
*/ 
int CoefApparentement(int* Genealogie,	int* plProposant, int NProposant, int* plAncetre, double* pdRetour,int DuppDetection, int printprogress)
{
	try{
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,NProposant,&NoeudPro);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **TmpNoeudAnc=NULL;	
	LoadAncetre(plAncetre,1,&TmpNoeudAnc);
	CIndSimul *NoeudAnc=*TmpNoeudAnc;	

	//Creation tableau
    INITGESTIONMEMOIRE;

	//OPTIONNEL
	CIndSimul** Ordre		=(CIndSimul**) 	memalloc(lNIndividu,sizeof(CIndSimul*));		
	CApPath**  Path			=(CApPath**) 	memalloc(NProposant,sizeof(CApPath*));
	CApPath **Current 		=(CApPath**) 	memalloc(NProposant,sizeof(CApPath*));	
	int i;	
	
	//initialisation de la structure noeud
	for(i=0;i<lNIndividu;i++)
	{		
		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].bFlagSort=0;	   		
		Noeud[i].iind=0;
		
		Ordre[i]=NULL;
	}
			
	//#2 : CREATION DU BITFIELD OU iind EST LE NO DU BITS	
	//initialisation
	for(i=0;i<NProposant;i++)
	{
 		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
		Current[i]=NULL;
		Path[i]=NULL;		
	}
	NoeudAnc->etat=GENDEPART;

	//Trouve les noeuds utile
	ExploreArbre(NoeudAnc);

	int cbit=-1;
	for(i=0;i<lNIndividu;i++)
	{
		if (Noeud[i].etat!=GENNONEXPLORER && Noeud[i].etat!=GENINUTILE)
		{
			Noeud[i].iind=++cbit;
		}
	}
	//Nombre de noeud touch, taille de la chaine de char 8 bit
	//const long NOrdre=cbit;
	const int taille=(int) ceil(double(cbit)/MP_DIGIT_BIT); //MPI
	
	//AJUSTE LA PRECISION DU MPI
	mp_set_prec(taille); //Indique a mpi de faire des nombres de la bonne taille....
	
	//RECHERCHE DE TOUS LES CHEMINS POSSIBLES ENTRE PROPOSANT ET ANCETRE
	g_ExpCoeff_CheminParcouru=Ordre; 		
	for(int cCol=0;cCol<NProposant;++cCol)
	{						
		Path[cCol]=NULL;
		g_ExpCoeff_Path=&Path[cCol];  
		g_ExpCoeff_Cible=NoeudPro[cCol];

#ifdef USESDEBUG
		//printf("\nPath finder %d -> %d\n",NoeudAnc->nom,NoeudPro[cCol]->nom);
#endif
		ExploreCoeff(NoeudAnc);				
	}

	//EVALUATION DU NOMBRE DE POSSIBILITE
	//long double dnbposs=1;
	double dnbposs=1;
	//unsigned long NbIteration = 1;
	unsigned int NbIteration = 1;
	for(int cCol=0;cCol<NProposant;++cCol)
	{						
		int tmpcountpath=0;			
		CApPath* tmp=Path[cCol];
		while (tmp!=NULL)
		{
			++tmpcountpath;
			tmp=tmp->next;
		}
		dnbposs *= tmpcountpath;
		NbIteration *= tmpcountpath; //Peut dbord
	}
	
	//VALIDE SI LE NOMBRE DE COMBINAISON EST TROP GRAND
	const long double MaxMemoryUsed = dnbposs*(taille*sizeof(mp_digit)+sizeof(mp_int))/COAPPMOD_HASHAGECHARGEMAXIMAL; //En octet
	if (DuppDetection)
	{		
		if (dnbposs>COAPPMOD_NBPOSSMAX_DUPP)
		{
			PathDestruction(Path,NProposant);
//			GENError("Number of combination to evaluate is too great for duplicata detection\nDeactivate it if you want to continue.");
			throw std::range_error("Number of combination to evaluate is too great for duplicata detection\nDeactivate it if you want to continue.");
			//GENError("Le nombre de combinaison a valuer est trop grand pour pouvoir utilis la dtection de dupplicata\n Dactiv la dtection de dupplication si vous dsir continuer"); 
		}

		if (MaxMemoryUsed>COAPPMOD_MAXMEMORY_USABLE)
		{
			PathDestruction(Path,NProposant);
//			GENError("Memory usage is too great for duplicata detection\nDeactivate it if you want to continue\n"
//				    "Maximum memory allowed: %lG Mo  memory needed: %lG Mo\n\n", 
//					COAPPMOD_MAXMEMORY_USABLE/MEGAOCTET,MaxMemoryUsed/MEGAOCTET);
			char erreur[TAILLEDESCRIPTION];
			
			double maxMemAllowed = COAPPMOD_MAXMEMORY_USABLE/MEGAOCTET;
			double memUsed = MaxMemoryUsed/MEGAOCTET;
			snprintf(erreur, TAILLEDESCRIPTION, "Memory usage is too great for duplicata detection\nDeactivate it if you want to continue\nMaximum memory allowed: %f Mo  memory needed: %f Mo\n\n",
				maxMemAllowed, memUsed);
			//sprintf(erreur, "Memory usage is too great for duplicata detection\nDeactivate it if you want to continue\nMaximum memory allowed: %Lf Mo  memory needed: %Lf Mo\n\n",
				// (COAPPMOD_MAXMEMORY_USABLE/MEGAOCTET),(MaxMemoryUsed/MEGAOCTET));
			throw std::range_error(erreur);
			//GENError("La quantite de mmoire utilis par la dtection de dupplicata est trop importante\n Dactiv la dtection de dupplication si vous dsir continuer"
					 //"\nTaille maximal permise : %lG Mo   Mmoire demand : %lG Mo\n\n",	 COAPPMOD_MAXMEMORY_USABLE/MEGAOCTET,MaxMemoryUsed/MEGAOCTET);
		}
	}
	else
	{
		if (dnbposs>COAPPMOD_NBPOSSMAX_NODUPP)
		{
			PathDestruction(Path,NProposant);
//			GENError("Execution time of this function call is too great");
			throw std::range_error("Execution time of this function call is too great");
			//GENError("La dure d'excution de cet appel de fonction serait surement comparable  l'ge de l'univers");
		}
	}

	//validation de solution triviale
	*pdRetour=0.0;//Initialise a zero la valeur de retour
	if (NbIteration==0)
	{
		PathDestruction(Path,NProposant);
		return 0;  //Operation russi mais le rsultat est simplement zero
	}

	//CREATION DE LA TABLE ANTI DUPPLICATA	
	//initrand();
	HashDouble<mp_int>* ptrHashTable=NULL;

	if (DuppDetection)
		ptrHashTable= new HashDouble<mp_int>(NbIteration,COAPPMOD_HASHAGECHARGEMAXIMAL);
	if (ptrHashTable==NULL || ptrHashTable->Error())
	{
		//PathDestruction(Path,NProposant);
//		GENError( "Insufficient memory to create an anti-duplicata table\n Deactivate duplicata detection if you want to continue\n"
//				"Memory needed: %lG Mo\n", MaxMemoryUsed/MEGAOCTET);
		char erreur[TAILLEDESCRIPTION];
		double memUsed = MaxMemoryUsed/MEGAOCTET;
		snprintf(erreur, TAILLEDESCRIPTION, "Insufficient memory to create an anti-duplicata table\n Deactivate duplicata detection if you want to continue\nMemory needed: %f Mo\n",
			memUsed);
		//sprintf(erreur, "Insufficient memory to create an anti-duplicata table\n Deactivate duplicata detection if you want to continue\nMemory needed: %Lf Mo\n",
		//	 (MaxMemoryUsed/MEGAOCTET));
		throw std::range_error(erreur);
		//GENError("Mmoire insuffisante pour utilis crer une table anti-dupplicata\n Dactiv la dtection de dupplication si vous dsir continuer"
				 //"\nMmoire ncessaire : %lG Mo\n",MaxMemoryUsed/MEGAOCTET);
	}

	//SOMMATION DE LA PROBABILITE CONJOINTE DE CHAQUE CHEMIN POSSIBLE
	//l'quivalent d'un ensemble de nbproposant boucle imbriqu.
	//Chaque boucle reprsente l'ensemble des chemins partant de l'ancetre
	//et ce rendant au proposant de cette boucle
	
	//Initialisation
	for(i=0;i<NProposant;++i)
		Current[i]=NULL;; //INITIALISE LE COMPTEUR DE BOUCLE

	int n=0; //Curseur qui indique sur que chemin on est actuellement
		
	mp_int tmpResult; //La solution courante  tudier
	mp_int tmp; 
	mp_init(&tmpResult);
	mp_init(&tmp);
	
	const int LastProposant = NProposant-1;	
	//CREATE_PROGRESS_BAR(NbIteration,printprogress);	
	while(n>=0)
	{
		//AVANCE Aincremente
		if (Current[n]==NULL) Current[n]=Path[n];
		else				  Current[n]=Current[n]->next;
						
		if (Current[n]!=NULL) {
			//EST-CE LA DERNIERE BOUCLE
			if (n==LastProposant) {
				//OUI 
				
				//INCREMENT_PROGRESS_BAR();

				//genere la solution & la distance																						
				//Pour chaque proposant
				mp_zero(&tmpResult); //Remise a zero
				for(int cCol=0;cCol<NProposant;++cCol) 
				{						
					//OR pour combiner le tout
					mpl_or(&tmpResult, &(Current[cCol]->num), &tmp);
					//Remet le rsultat dans tmpResult					
					mp_exch(&tmpResult, &tmp);									
				}
				
				//Validation Anti-dupplicata											
				if (!DuppDetection ||
					 (DuppDetection && ptrHashTable->Add(tmpResult)==GTRUE) )
				{
					//Calcul et ajout de la distance au rsultat		
					int distance=0;
					mpl_num_set(&tmpResult, &distance);  
					*pdRetour += pow2(distance-1); //Me rappelle plus pourquoi -1... 
				}
			}//Derniere boucle: NON : if (n==(NProposant-1))
			else
				++n;
		} //if (Current[n]!=NULL)
		else
			--n;		
	}//fin while(n>=0)	

	//VALEUR DE RETOUR
	//*pdRetour

	//NETTOYER LES ASSIGNEMENTS DE MMOIRES			
	//outrand();
	PathDestruction(Path,NProposant);

	mp_clear(&tmpResult);
	mp_clear(&tmp);
	if (ptrHashTable)
		delete ptrHashTable;	
	
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			}catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}


// ********************************************************************
//
//			FONCTION PRIVE
//
// ********************************************************************


/*! 
	\brief (INTERNE) Detruit un tableau de liste CApPath

		Pour chaque element du vecteur Path. detruit la liste (et tout le element contenu dans celle-ci)

	\param Path		[in] Ptr vers un vecteur de pointeur vers des liste CApPath
	\param npath	[in] Taille du vecteur Path
	
*/
static void PathDestruction(CApPath **Path,int npath)
{
	//DESTRUCTION DES PATHS (A VERIFIER)
	CApPath *Current=NULL;	
	CApPath *tmp	=NULL;
	for(int j=0;j<npath;j++)
	{
		Current=Path[j];
		while(Current!=NULL)
		{
			
			tmp=Current->next;
			mp_clear(&Current->num); 
			memfreeIN(Current); 
			Current=tmp;								
		}
	}
}


///Utilise par ExploreCoeff pour reprsent le chemin parcouru
/*
static CIndSimul** g_ExpCoeff_CheminParcouru=NULL; 
///Utilise par ExploreCoeff comme tant le dernier chemin remplis (ou le premiers)
static CApPath ** g_ExpCoeff_Path=NULL;
///Utilise par ExploreCoeff comme tant la cible de l'exploration
static CIndSimul* g_ExpCoeff_Cible=NULL;*/
/*! 
	\brief (INTERNE) Explore une genealogie et construit une liste CapPath
		
		Trouve tous les chemins entre le Noeud et la cible et sa les inscrit dans une liste CApPath
		sous forme binaire mpi_int.

		La forme binaire correspond a la somme de tout les 2<<iind des noeuds implique et peut s'etendre sur un nombre
		infini d'octet

	\attention Il faut que les noeuds soit correctement etiquette (etat)
				De plus, il faut qu'un iind unique soit assigne a chaque noeud qui sont Utile
				???? important avant mp_set_prec(taille);
				g_ExpCoeff_CheminParcouru
				g_ExpCoeff_InputPath
				g_ExpCoeff_Cible
	
	\param Noeud		[in] Noeud courant ou l'exploration est rendu (Pts de depart la majorite du temps)
	\param profondeur	[in] Profondeur actuel (au depart:0 normalement)
	\param Cible		[in] Noeud Cible 
	\param Array2		[in] Pointeur vers un vecteur de cindsimul* de taille nombre de noeud max implique 
							 (utilis pour mmoris le chemin utilis)
	\param InputPath	[in] Taille du vecteur Path
	\param taille2		[in] Taille du vecteur Path
	
	  \remark Cette fonction est recursive
*/
static void FASTCALL ExploreCoeff(CIndSimul* Noeud)
{
	//Explore l'arbre et retourne tous les path utilis		
	static int profondeur=0;

	//Trace le chemin parcouru
	g_ExpCoeff_CheminParcouru[profondeur]=Noeud;

	//on vient d'atteindre la cible?
	if (Noeud==g_ExpCoeff_Cible)
	{		
		//Creation d'un nouveau chemin
		CApPath *tmp=(CApPath*) memallocIN(1,sizeof(CApPath));
		
		//Initialise le nombre contenu....		
		mp_init( &tmp->num );		
		tmp->next=NULL;

		//Complete la srie
		*g_ExpCoeff_Path=tmp;
		g_ExpCoeff_Path=&(tmp->next);	//Avance le curseur
		
		//Ajuste le nombre en consquence....
		for(int i=0;i<=profondeur;i++)
		{
#ifdef USESDEBUG
			//printf("%d,",g_ExpCoeff_CheminParcouru[i]->nom);
#endif
			mpl_bit_set(&tmp->num, g_ExpCoeff_CheminParcouru[i]->iind);
		}
#ifdef USESDEBUG
		unsigned char buffer[100];
		mp_toradix(&tmp->num, buffer, 16);
		//printf(" == %s\n",buffer);
#endif

		return;	
	}	
	else
	{
		//Non.. on continue a cherch
		Clist *current=Noeud->fils;
		if (current!=NULL)
		{
			do
			{	
				CIndSimul *tmp=current->noeud;				
				if (tmp->etat!=GENNONEXPLORER && tmp->etat!=GENINUTILE)
				{
					++profondeur;
					ExploreCoeff(tmp);
					--profondeur;
				}
				//Fils suivant
				current=current->next;				
			}
			while(current!=NULL);
		}
	}	
	return;	
}

