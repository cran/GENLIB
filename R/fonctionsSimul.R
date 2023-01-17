# 1 - gen.prob					-> garde article
# 2 - gen.simuProb				-> garde article
#     gen.simuhaplo
#     gen.simuhaplo_traceback
#     gen.simuhaplo_compare
# 3 - gen.simuSample			-> garde article
# 4 - gen.simu2All				-> garde article
# 5 - gen.simuSampleFreq			-> garde article
# 6 - gen.simuSet	 			-> garde article

#gen.prob = function(gen, pro, statePro, ancestors, stateAncestors, onlyConj = F, print.it = F, named = T)
#{
#	if(!is(gen, "GLgen"))
#		stop("Invalid parameter: gen must an instance of GLgen (see gen.genealogy)")
#		#stop("Parametre invalide: gen doit etre une instance de GLgen (voir gen.genealogy)")
#	if(!is(pro, "numeric") || !is(statePro, "numeric"))
#		stop("Invalid parameter: pro and statePro must be numeric vectors")
#		#stop("Parametre invalide: proband et statePro doivent etre des vecteurs numeriques")
#	if(!is(ancestors, "numeric") || !is(stateAncestors, "numeric"))
#		stop("Invalid parameter: ancestors and stateAncestors must be numeric vectors")
#		#stop("Parametre invalide: ancestors et stateAncestors doivent etre des vecteurs numeriques")
#	if(length(pro) != length(statePro))
#		stop("Invalid parameter: pro and statePro must have the same size")
#		#stop("Parametre invalide: Proband et statePro doivent etre de meme taille")
#	if(length(ancestors) != length(stateAncestors))
#		stop("Invalid parameter: ancestors and stateAncestors must have the same size")
#		#stop("Parametre invalide: ancestors et stateAncestors doivent etre de meme taille")
#	if(!is(print.it, "logical") || !is(named, "logical") || !is(onlyConj, "logical"))
#		stop("Invalid parameter:  print.it, named and onlyConj must be logical values")
#		#stop("Parametre invalide: print.it,named et onlyConj doivent etre des valeurs logiques")
#	#Tableau de resultat 
#	tmpconj <- double(1)
#	tmpsimul <- double(length(pro))
#	#EXPORTTYPE void SPLUSProb(long* Genealogie, 
#	#long* proband, long* stateproband,long* nproband, 
#	#long* ancestor, long* stateancestor, long* nancestor, 
#	#double* Retour,long *PrintProgress)
#	retour <- .Call("SPLUSProb",  gen@.Data, pro, statePro, length(pro), ancestors, as.integer(stateAncestors), length(ancestors),
#					tmpconj, tmpsimul, print.it, onlyConj, specialsok = T)
#	#print(paste("retour:",retour))
#	if(retour != -1) 
#   {if(!onlyConj) {
#		v <- list(conj = tmpconj, simul = tmpsimul)
#		#if(named)
#			names(v$simul) <- c(pro)
#	}
#	else v <- tmpconj
#	if(print.it) {
#		argument <- c(deparse(substitute(gen)), deparse(substitute(pro)), deparse(substitute(statePro)), deparse(substitute(
#			ancestors)), deparse(substitute(stateAncestors)))
#		header.txt <- paste("\n   ***   Appel \340 : gen.prob (", argument[1], ",", argument[2], ",", argument[3], ",", 
#			argument[4], ",", argument[5], ")  ***\n\n", sep = "")
#		cat(header.txt)
#	}
#	v
#	}
#}
#print.it = F, 
gen.simuProb = function(gen, pro, statePro, ancestors, stateAncestors, simulNo=5000, probRecomb=c(0,0), probSurvival=1.0)
{
	if(!is(gen, "GLgen"))
		stop("Invalid parameter: gen must be an instance of Glgen (see gen.genealogy)")
	if(!is(pro, "numeric") || !is(statePro, "numeric"))
		stop("Invalid parameter: pro and statePro must be numeric vectors")
	if(!is(ancestors, "numeric") || !is(stateAncestors, "numeric"))
		stop("Invalid parameter: ancestors and stateAncestors must be numeric vectors")
	if(length(pro) != length(statePro))
		stop("Invalid parameter: pro and statePro must have the same size")
	if(length(ancestors) != length(stateAncestors))
		stop("Invalid parameter: ancestors and stateAncestors must have the same size")
	if(simulNo <= 0)
		stop("Invalid parameter: simulNo must be greater than zero")
#	if(!is(named, "logical"))
#		stop("Invalid parameter: named must be logical values")
	if( !(sum(statePro>2) %in% c(0,1,2)) || !(sum(stateAncestors>2) %in% c(0,1,2)))
		warning("Some ancestors or probands states are not 0, 1 or 2. The ones below 0 will be set to 0 and the ones over 2 will be set to 2")
	
	#Tableau de resultat 
	tmpconj <- double(1)
	tmpsimul <- double(length(pro))
	tmppro <- double(length(pro) + 1)
	#extern "C" void SPLUSSimul(long* Genealogie, 
	#long* proband, long* stateproband,long* nproband, 
	#long* ancestor, long* stateancestor, long* nancestor, 
	#long* nSimul, double* pdRetConj,double* pdRetSimul,double* pdRetPro,long *PrintProgress)
	.Call("SPLUSSimul", gen@.Data, pro, statePro, length(pro), ancestors, as.integer(stateAncestors), length(ancestors),
						as.integer(simulNo), tmpconj, tmpsimul, tmppro, probRecomb, probSurvival, FALSE, specialsok = T)
	v <- list(joint = tmpconj, marginal = tmpsimul, Nprobands = tmppro)
	#if(named) {
		names(v$marginal) <- c(pro)
		names(v$Nprobands) <- c(0:length(pro))
	#}
#	if(print.it) {
#		argument <- c(deparse(substitute(gen)), deparse(substitute(pro)), deparse(substitute(statePro)), deparse(substitute(
#			ancestors)), deparse(substitute(stateAncestors)), simulNo)
#		header.txt <- paste("\n   ***   Calls : gen.simuProb (", argument[1], ",", argument[2], ",", argument[3], ",",
#			argument[4], ",", argument[5], ",", argument[6], ")  ***\n\n", sep = "")
#		cat(header.txt)
#	}
	v
}
#print.it = F, 

gen.simuHaplo = function (gen, pro=NULL, ancestors=NULL, simulNo = 1, model =1, model_params, cM_len, 
							BP_len, physical_map_Mo = NULL, physical_map_Fa = NULL, seed= 0, all_nodes =0, outDir = getwd()){
	if(!is(gen, "GLgen"))
		stop("Invalid parameter: gen must be an instance of Glgen (see gen.genealogy)")
	if(simulNo <= 0)
		stop("Invalid parameter: simulNo must be greater than zero")
	if(is.null(ancestors))
		ancestors = gen.founder(gen)
	if(is.null(pro))
		pro=gen.pro(gen)
	BP_len = as.integer(BP_len)


	if(length(model_params) != 2)
		stop("model_params must be a vector of length 2")
	if(!is(cM_len, "numeric"))
		stop("Invalid parameter: MorganLen must be a numeric vector")
	if(length(cM_len) != 2)
		stop("cM_len must be a vector of length 2")
	#comparisons to NULL don't produce boolean value	
	#if(Reconstruction==1 & (Hapfile==NULL | Mapfile==NULL))
		#stop("If reconstruction is set to 1 must specify the hap and map files")
	# if(Reconstruction==1 & BP==0)
	# 	stop("If reconstruction is set to 1, you must specify the size of the segment in BP")
	# if(Reconstruction==1 & (is.null(Hapfile) | is.null(Mapfile)))
	# 	stop("If reconstruction is set to 1, you must provide a hap file and map file")		

	if(seed==0)
		seed=as.integer(Sys.time())
	
	if(is.null(physical_map_Fa) & is.null(physical_map_Mo)){
		message("No map function specified to convert genetic distance to physical. Assumed constant along length of chromosome")
		convert = 0
		bp_map_FA = 0
		cm_map_FA = 0
		bp_map_MO = 0
		cm_map_MO = 0
	}
	else{
		#need to check that the maps are valid
		convert = 1
		if(!("BP" %in% colnames(physical_map_Fa)) || !("cM" %in% colnames(physical_map_Fa)))
			stop("column names 'BP', or 'cM' not in physical_map_Fa. Names are case-sensitive")
		if(!("BP" %in% colnames(physical_map_Mo)) || !("cM" %in% colnames(physical_map_Mo)))
			stop("column names 'BP', or 'cM' not in physical_map_Fa. Names are case-sensitive")
		bp_map_FA = as.integer(physical_map_Fa$BP)
		cm_map_FA = as.numeric(physical_map_Fa$cM)
		bp_map_MO = as.integer(physical_map_Mo$BP)
		cm_map_MO = as.numeric(physical_map_Mo$cM)

		if((bp_map_FA[1]!=0) || (cm_map_FA[1]!=0))
			stop("first element of BP and cM columns should be 0")
		if((bp_map_MO[1]!=0) || (cm_map_MO[1]!=0))
			stop("first element of BP and cM columns should be 0")
		if((bp_map_FA[length(bp_map_FA)]!=BP_len) || (cm_map_FA[length(bp_map_FA)]!=cM_len[1]))
			stop("last element of BP and cM columns should be BP_len, and cM_len, respectively")
		if((bp_map_MO[length(bp_map_MO)]!=BP_len) || (cm_map_MO[length(bp_map_MO)]!=cM_len[2]))
			stop("last element of BP and cM columns should be BP_len, and cM_len, respectively")
	}
	
	message("seed: ", seed, "\n")

	.Call("SPLUSSimulHaplo", gen@.Data, pro, length(pro), ancestors, length(ancestors), as.integer(simulNo), model_params, cM_len/100, as.integer(model), 
			as.integer(convert), as.integer(BP_len), as.integer(bp_map_FA), cm_map_FA, as.integer(bp_map_MO), cm_map_MO, 
			outDir, as.integer(all_nodes), as.integer(seed), package="GENLIB")

	if(all_nodes == 0)
		message("output file:  ", outDir, "/Proband_Haplotypes.txt\n")
	else
		message("output files:\n", outDir, "/All_nodes_haplotypes.txt\n", outDir, "/Proband_Haplotypes.txt\n")

}

gen.simuHaplo_traceback = function(gen, proID, ancestorID, all_nodes_path, proband_haplotypes_path)
{
	t <- gen.genout(gen)
	indVec <- t$ind
	fatherVec <- t$father
	motherVec <- t$mother

	message("input file paths:\n", all_nodes_path, "\n", proband_haplotypes_path, "\n")

	x = .Call("SPLUSSimulHaplo_traceback", as.integer(proID), as.integer(ancestorID), indVec, fatherVec, motherVec, all_nodes_path, proband_haplotypes_path, package="GENLIB")
	#need to pass in the vectors that will hold the results
	return(x)
}

gen.simuHaplo_IBD_compare = function (proID_1, proID_2, BP_len, proband_haplotypes_path)
{
	x = .Call("SPLUSSimulHaplo_IBD_compare", as.integer(proID_1), as.integer(proID_2), as.integer(BP_len), proband_haplotypes_path)
	return(x)
}

gen.simuHaplo_convert = function (dir = getwd())
{
	x = .Call("SPLUSSimulHaplo_convert", dir)
}
gen.simuSample = function(gen, pro, ancestors, stateAncestors, simulNo = 5000)#, named = T)
{
	if(!is(gen, "GLgen"))
		stop("Invalid parameter: gen must be an instance of Glgen (see gen.genealogy)")
	if(!is(pro, "numeric"))
		stop("Invalid parameter: pro must be numeric vectors")
	if(!is(ancestors, "numeric") || !is(stateAncestors, "numeric"))
		stop("Invalid parameter: ancestors and stateAncestors must be numeric vectors")
	if(length(ancestors) != length(stateAncestors))
		stop("Invalid parameter: ancestors and stateAncestors must have the same size")
	if(simulNo <= 0)
		stop("Invalid parameter: simulNo must be greater than zero")
#	if(!is(named, "logical"))
#		stop("Invalid parameter: named must be logical values")
	if( !(sum(stateAncestors>2) %in% c(0,1,2)) )
		warning("Some ancestors states are not 0, 1 or 2. The ones below 0 will be set to 0 and the ones over 2 will be set to 2")

	#Valeur de retour  
	tmpsimul <- double(length(pro) * simulNo)
	#extern "C" void SPLUSSimulSingle(long* Genealogie, 
	#long* proband, long* nproband, 
	#long* ancestor, long* stateancestor, long* nancestor, 
	#double* pdRetour,long* PrintProgress)
	#Call de la fonction en C
	.Call("SPLUSSimulSingle", gen@.Data, pro, length(pro), ancestors, as.integer(stateAncestors), length(ancestors), 
							  as.integer(simulNo), tmpsimul, FALSE, specialsok = TRUE)
	#Mise en matrice
	dim(tmpsimul) <- c(length(pro), simulNo)
	#if(named)
		dimnames(tmpsimul) <- list(pro, NULL)
#	if(print.it) {
#		argument <- c( deparse(substitute(gen)), deparse(substitute(pro)), deparse(substitute(ancestors)), 
#					deparse(substitute(stateAncestors)), simulNo )
#		header.txt <- paste("\n   ***   Calls : gen.simuSample (", argument[1], ",", argument[2], ",", argument[3], ",",
#						argument[4], ",", argument[5], ")  ***\n\n", sep = "")
#		cat(header.txt)
#	}
	drop(tmpsimul)
}

#gen.simu2All = function(gen, pro, ancestors, stateAncAll1, stateAncAll2, simulNo = 5000, f = f <- function(x) return(x), 
#						print.it = F, named = T)
#{
#	if(!is(gen, "GLgen"))
#		stop("Invalid parameter: gen must be an instance of Glgen (see gen.genealogy)")
#		#stop("Parametre invalide: gen doit etre une instance de GLgen (voir gen.genealogy)")
#	if(!is(pro, "numeric") & !is(pro, "GLgroup"))
#		stop("Invalid parameter: pro must be numeric vectors")
#		#stop("Parametre invalide: proband doit etre un vecteur numerique")
#	if(!is(ancestors, "numeric") || !is(stateAncAll1, "numeric") || !is(stateAncAll2, "numeric"))
#		stop("Invalid parameter: ancestors, stateAncAll1 and stateAncAll2 must be numeric vectors")
#		#stop("Parametre invalide: ancestors et stateAncestors doivent etre des vecteurs numeriques")
#	if(length(ancestors) != length(stateAncAll1))
#		stop("Invalid parameter: ancestors and stateAncAll1 must have the same size")
#		#stop("Parametre invalide: ancestors et stateAncAll1 doivent etre de meme taille")
#	if(length(ancestors) != length(stateAncAll2))
#		stop("Invalid parameter: ancestors and stateAncAll2 must have the same size")
#		#stop("Parametre invalide: ancestors et stateAncAll2 doivent etre de meme taille")
#	if(simulNo <= 0)
#		stop("Invalid parameter: simulNo must be greater than zero")
#		#stop("Parametre invalide: simulNo doit etre plus grand que zero")
#	if(!is(print.it, "logical") || !is(named, "logical"))
#		stop("Invalid parameter: print.it and named must be logical values")
#		#stop("Parametre invalide: print.it et named doivent etre des valeurs logiques")
#
#	#  J'ai aucune idee pourquoi c'est mis en liste surtout que ca cause plein de probleme cote C++ !
#	#  faique ca redevient un vecteur d'entier.
#	#if(is(pro, "numeric"))
#	#	pro = list(pro)
#	
#	#Valeur de retour  
#	#extern "C" void SPLUSSimulSingle
#	tmpsimul <- .Call("SPLUSSimulSingleFct", gen@.Data, pro, ancestors, as.integer(stateAncAll1), as.integer(stateAncAll2), 
#										length(ancestors), as.integer(simulNo),f, print.it, copy = c(F, F, F, F, F, F, F, F, F))
#	if(print.it) {
#		argument <- c(deparse(substitute(gen)), deparse(substitute(pro)), deparse(substitute(ancestors)), deparse(substitute(
#			stateAncestors)), simulNo)
#		header.txt <- paste("\n   ***   Appel a : gen.simu2All (", argument[1], ",", argument[2], ",", argument[3],
#			",", argument[4], ",", argument[5], ")  ***\n\n", sep = "")
#		cat(header.txt)
#	}
#	drop(tmpsimul)
#}
# print.it = F,
gen.simuSampleFreq = function(gen, pro, ancestors, stateAncestors, simulNo = 5000)#, named = T)
{
	if(!is(gen, "GLgen"))
		stop("Invalid parameter: gen must be an instance of Glgen (see gen.genealogy)")
	if(!is(pro, "numeric"))
		stop("Invalid parameter: pro must be numeric vectors")
	if(!is(ancestors, "numeric") || !is(stateAncestors, "numeric"))
		stop("Invalid parameter: ancestors and stateAncestors must be numeric vectors")
	if(length(ancestors) != length(stateAncestors))
		stop("Invalid parameter: ancestors and stateAncestors must have the same size")
	if(simulNo <= 0)
		stop("Invalid parameter: simulNo must be greater than zero")
#	if(!is(named, "logical"))
#		stop("Invalid parameter: named must be logical values")
	if(!(sum(stateAncestors>2) %in% c(0,1,2)))
		warning("Some ancestors states are not 0, 1 or 2. The ones below 0 will be set to 0 and the ones over 2 will be set to 2")

	#Valeur de retour  
	tmpsimul <- double(length(pro) * 3)
	#extern "C" void SPLUSSimulSingle(long* Genealogie, 
	#long* proband, long* nproband, 
	#long* ancestor, long* stateancestor, long* nancestor, 
	#double* pdRetour,long* PrintProgress)
	#Call de la fonction en C
	.Call("SPLUSSimulSingleFreq", gen@.Data, pro, length(pro), ancestors, as.integer(stateAncestors), length(ancestors), 
								as.integer(simulNo), tmpsimul,	FALSE, specialsok = TRUE)
	#Mise en matrice
	dim(tmpsimul) <- c(length(pro), 3)
	#if(named)
		dimnames(tmpsimul) <- list(pro, paste(c("Alleles transmitted 0", "Alleles transmitted 1", "Alleles transmitted 2")))
#	if(print.it) {
#		argument <- c(deparse(substitute(gen)), deparse(substitute(pro)), deparse(substitute(ancestors)), deparse(substitute(
#			stateAncestors)), simulNo)
#		header.txt <- paste("\n   ***   Calls : gen.simuSampleFreq (", argument[1], ",", argument[2], ",", argument[3],
#			",", argument[4], ",", argument[5], ")  ***\n\n", sep = "")
#		cat(header.txt)
#	}
	tmpSimul2 = data.frame(tmpsimul)
	drop(tmpSimul2)
}
#print.it = F,
gen.simuSet = function(gen, pro, ancestors, stateAncestors, probMatrix = matrix(c(c(1,0.5,0,0.5,0.25,0,0,0,0,1,1,1,1,0.75,0.5,1,0.5,0),
		c(1,0.5,0,0.5,0.25,0,0,0,0,1,1,1,1,0.75,0.5,1,0.5,0)), nrow = 3, ncol = 12), simulNo = 5000)#,  named = T)
{
	if(!is(gen, "GLgen"))
		stop("Invalid parameter: gen must be an instance of Glgen (see gen.genealogy)")
	if(!is(pro, "numeric") & !is(pro, "GLgroup"))
		stop("Invalid parameter: pro must be numeric vectors")
	if(!is(ancestors, "numeric") || !is(stateAncestors, "numeric"))
		stop("Invalid parameter: ancestors and stateAncestors must be numeric vectors")
	if(length(ancestors) != length(stateAncestors))
		stop("Invalid parameter: ancestors and stateAncestors must have the same size")
	if(simulNo <= 0)
		stop("Invalid parameter: simulNo must be greater than zero")
#	if(!is(named, "logical"))
#		stop("Invalid parameter: named must be logical values")

	# verification de la matrice
	if(sum(probMatrix>1)>0) warning("One or more value in the cumulative probability matrix (probMatrix) is larger than 1.0.")
	if(sum( (probMatrix[1:3,4:6]-probMatrix[1:3,1:3])<0 ) > 0 | sum( (probMatrix[1:3,10:12]-probMatrix[1:3,7:9])<0 ) > 0 ) 
		stop("Cumulative probability of having 0 allele transmitted must always be smaller or equal to cumulative probability of having 1 allele transmitted.")

		#Valeur de retour  
	#extern "C" void SPLUSSimulSingle
	tmpsimul <- .Call("SPLUSSimulSingleProb", gen@.Data, pro, length(pro), ancestors, length(ancestors), stateAncestors, probMatrix,
									  as.integer(simulNo), FALSE, copy = c(F, F, F, F, F, F, F, F))
	dim(tmpsimul) <- c(length(pro), simulNo)
	#if(named)
		dimnames(tmpsimul) <- list(pro, NULL)
#	if(print.it) {
#		argument <- c(deparse(substitute(gen)), deparse(substitute(pro)), deparse(substitute(ancestors)), deparse(substitute(
#			stateAncestors)), simulNo)
#		header.txt <- paste("\n   ***   Calls : gen.simuSet (", argument[1], ",", argument[2], ",", argument[3],
#			",", argument[4], ",", argument[5], ")  ***\n\n", sep = "")
#		cat(header.txt)
#	}
	drop(tmpsimul)
}

