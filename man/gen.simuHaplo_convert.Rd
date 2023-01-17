\name{gen.simuHaplo_convert}
\alias{gen.simuHaplo_convert}
\title{Convert proband simulation results into sequence data given founder haplotypes}
\description{Convert the "Proband_Haplotypes.txt" output file produced by gen.simulHaplo into sequence data according to input sequence data for founders.}
\usage{gen.simuHaplo_convert(dir)
}
\arguments{
\item{dir}{path to the directory containing all the required input files. The function expects to find the "Proband_Haplotypes.txt" file,
a file named "founders.hap" which contains the sequence data for all founder chromosomes, and another file called "founders.map" which 
contains the BP position of all included markers in the specified founder sequences.}
}
\details{
See \url{https://github.com/R-GENLIB/simuhaplo_functions} for examples on how to format the founder sequence input files.}
\value{creates a new file in the specified directory "Proband_Genotypes.txt"}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.simuHaplo}}
}
\examples{

#gen_table <- data.frame(ind=c(1,2,3,4,5,6,7,8),
#                        mother=c(0,0,2,0,0,5,3,3),
#                        father=c(0,0,1,0,0,4,6,6),
#                        sex=c(1,2,2,1,2,1,1,2))

#convert the genealogy table into a GENLIB object
#sample_gen<-gen.genealogy(gen_table) 


# Simulate the gene dropping from founders to probamds
#gen.simuHaplo(sample_gen, simulNo=5, model = 1, model_params = c(1,1), cM_len=c(100,100),
#                BP_len = 100000000, outDir=tempdir())

#after creating a "founders.map" and a "founders.hap" file in the same directory then can ran
#gen.simulHaplo_convert(dir=tempdir())
}
