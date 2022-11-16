/*! \file outilanal.cc
\brief Library Genlib: outils pour les fonctions d'analyse statistique (statanal.cc)


\author Claire Gires

*/

#include "stdlib.h"
#include "outilanal.h"
#include <string.h>

Tuple::Tuple()
{
	ind = NULL;     
	tab.clear();
}


Tuple::Tuple( const Tuple& in)
{
	(*this)=in;
}


Tuple::~Tuple()
{
}

/**	\brief	Tuple::addtab  incrmente le nombre correspondant  la valeur val dans tab

	\param	val		[in] variable de type int

	\return	void
	
*/
void Tuple::addtab(int val)
{
	tab[val]++;
}

/**	\brief	retourne l'attribut tab

	\return	l'adresse de l'attribut tab

	
*/
const liste& Tuple::gettab() const
{
	return tab;
}


/**	\brief	modifie l'attribut ind

	\param	ref	[in] le nouvel individu associ au Tuple

	\return	void

	
*/
void Tuple::setNoeud(CIndSimul* ref)
{
	ind = ref;
}

/**	\brief	retourne le pointeur sur l individu associe au Tuple


	\return	CIndSimul*

	
*/
CIndSimul* Tuple::getNoeud()
{
	return ind;
}


/**	\brief	remise  zro

	Cette fonction remet  zero les autres attributs que ind, et donne  ind la valeur de ref

	\param	ref	[in] le nouvel individu associ au Tuple

	\return	void

	
*/
void Tuple::clear(CIndSimul* ref)
{
	ind = ref;
	tab.clear();
}

/**	\brief	comparaison lexicographique de deux Tuples

	teste si le Tuple courrant est plus petit que celui pass en paramtre au sens lexicographique.

	\param	lhs	[in] Premier tuple (gauche) a comparer
	\param	rhs	[in] Deuxieme tuple (droite) a comparer

	\return	int !=0 Si lhs < rhs
*/
int operator<(const Tuple& lhs, const Tuple& rhs)
{
	return lhs.tab < rhs.tab;	
}

/**	\brief	The Tuple::operator= oprateur de copie

	\param	original	[in] a parameter of type const Tuple&, Tuple  copier

	\return	Tuple&

	
*/
Tuple& Tuple::operator=(const Tuple& original)
{
	ind=original.ind;
	tab.clear();
	liste::const_iterator it = original.tab.begin();
	while(it!=original.tab.end())
	{
		tab[(*it).first]=(*it).second;
		it++;
	}
	return *this;
}


/**	\brief	fonction  qui teste l'galit entre deux Tuple

	deux Tuple sont gaux si leurs attributs tab sont gaux.
	(rappelons que tab est une map)
	Les individus auxquels ils font rfrence peuvent tre diffrents

	\param	t1	a parameter of type const Tuple&
	\param	t2	a parameter of type const Tuple&

	\return	int

	
*/
int operator==(const Tuple& t1,const Tuple& t2)
{
	return t1.tab == t2.tab;
}


/**	\brief	Fonction de comparaison de tuples pour qsort

	\param	e1	a parameter of type const void*
	\param	e2	a parameter of type const void*

	\return	int WINCDECL

*/
int WINCDECL QSORTtuple(const void* e1,const void* e2)
{
	Tuple& t1 = *((Tuple*)e1);
	Tuple& t2 = *((Tuple*)e2);

	if(t1==t2)
		return 0;
	if(t1<t2)
		return -1;
	else return 1;
}

