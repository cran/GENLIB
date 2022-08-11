/***************************
	Fichier qui contient les fonctions tr�s sp�cifique au hardware
\contributor Jean-Fran�ois Lefebvre
  ***********************************/
#ifndef GENLIBHAL
#define GENLIBHAL
#define WIN32_LEAN_AND_MEAN 

#if defined _WIN32 || defined _WIN64
 #include <Rcpp.h>
 #include <Windows.h>
#else
 #include <R_ext/RS.h>
#endif

/*FONCTION CLASSIQUE*/
int processorCount();
int thetime();

//FONCTION POUR LES SEMAPHORES
//leave CSema exposed even for OSX, Kinship4 struct requires a CSema argument
struct Opa_Cema;
typedef Opa_Cema *CSema;

#ifndef __APPLE__
int CSema_init(CSema& Semaphore,int compteinitial);
int CSema_wait(CSema& Semaphore);
int CSema_post(CSema& Semaphore);
int CSema_destroy(CSema& Semaphore);
#endif
//FONCTION POUR LE MULTITHREAD

#if defined _WIN32 || defined _WIN64
	#define THREADRETURN DWORD WINAPI  	//unsigned int
	#define THREADONEXIT return 0;
#else	
	#define THREADRETURN void*
	#define THREADONEXIT return NULL;
#endif

#ifndef __APPLE__
struct Opa_Thread;
typedef Opa_Thread *Cthread;
int  Cthread_create(Cthread& thread, THREADRETURN (*start_routine)(void*),void* arg);
void Cthread_destroy(Cthread& thread);
void Cthread_exit();
void Cthread_join(Cthread& thread);
#endif


#endif


