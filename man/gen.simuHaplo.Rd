\name{gen.simuHaplo}
\alias{gen.simuHaplo}
\title{Gene dropping simulations - haplotypes}
\description{Gene dropping from founders to specified probands according to specified model of meioisis. Output .txt files will have simulated haplotypes, labelled w.r.t. the ancestor of origin}
\usage{gen.simuHaplo(gen, pro = NULL, ancestors = NULL, simulNo = 1, model
                        = 1, model_params, cM_len, BP_len, physical_map_Mo =
                        NULL, physical_map_Fa = NULL, seed = 0, all_nodes = 0,
                        outDir = getwd())}

\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Will default to all probands in the provided genealogy.}
\item{ancestors}{Vector of ancestors id numbers to be included. Will default to all ancestors of probands}
\item{simulNo}{Number of simulations to run. Defaults to 1. Results from all simulations are contained in the same output text file}
\item{model}{model=1 is Poisson model, model=2 is zero-truncated Poisson, model=3 is Gamma model}
\item{model_params}{Numeric vector of length 2. Specifies sex-specific parameters for the model of meiosis. Should be in order of male, female}
\item{cM_len}{Length of chromosome (or subregion) to be simulated in centiMorgans. Integer vector of length 2, first element is cM length for males, second for females}
\item{BP_len}{Length of chromosome (or subregion) to be simulated in Base pairs. Single integer value (male and female chromosome have same BP length)}
\item{physical_map_Mo}{Optional map for females (Mother) that specifies relationship between genetic length and physical length.
    A dataframe with a column named "BP" and another named "cM". The map is the linear piecewise function specified by the points.
    The order of the pairs should be strictly ascending, the first row should be (0,0) and the last row should be (BP_len, cM_len)
    If a map is not specified then genetic length is assumed to be constant along the physical length of the region}
\item{physical_map_Fa}{Optional map for males, should be specified in the same way}
\item{seed}{integer seed for the random number generator (C++ std library Mersenne twister) for the simulation. Default is the system time}
\item{all_nodes}{all_nodes=1 will generate Optional "All_nodes_haplotypes.txt" outpu file. Default is all_nodes=0.}
\item{outDir}{The directory to write output files to. Default is the current working directory.}
}

\value{
Will generate output files "Proband_Haplotypes.txt" and "All_nodes_haplotypes.txt" in the specified output directory.
}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.simuHaplo_traceback}}
\code{\link{gen.simuHaplo_IBD_compare}}
\code{\link{gen.simuSampleFreq}}
}
\examples{

gen_table <- data.frame(ind=c(1,2,3,4,5,6,7),
                        mother=c(0,0,2,0,0,5,3),
                        father=c(0,0,1,0,0,4,6),
                        sex=c(1,2,2,1,2,1,1)
)

#convert the genealogy table into a GENLIB object
sample_gen<-gen.genealogy(gen_table) 

gen.graph(sample_gen)
# Simulate the gene dropping from founders to probamds

gen.simuHaplo(sample_gen, simulNo=5, model = 1, model_params = c(1,1), cM_len=c(100,100),
                BP_len = 100000000, outDir=tempdir())

f<-file(paste(tempdir(),"Proband_Haplotypes.txt",sep="/"),"r")
print(readLines(f,n=10))

}
