//PhiMatrix called from SPLUSPhiMatrix and SPLUSPHiMatrixMT in interface.cpp
//but these interface functions are not called in the R files, 
//gen.phi now calls SPLUSPhis
// SPLUSPhis calls Phis in apparentement.cpp which does not interact w/ these fx at all
//looks like these are no longer in use then
int PhiMatrix(int* Genealogie, int* proposant, int NProposant,int Niveau, double* pdMatricePhi, int printprogress);

int PhiMatrixMT(int* Genealogie, int* proposant, int NProposant,int Niveau, double* pdMatricePhi, int printprogress);


