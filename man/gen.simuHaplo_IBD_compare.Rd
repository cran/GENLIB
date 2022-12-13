\name{gen.simuHaplo_IBD_compare}
\alias{gen.simuHaplo_IBD_compare}
\title{Compare proband haplotypes for IBD sharing}
\description{Analyze the simulated haplotypes produced by gen.simuHaplo and identify segments of IBD sharing between the diploid haplotypes of a pair of probands}
\usage{gen.simuHaplo_IBD_compare(proID_1, proID_2, BP_len, proband_haplotypes_path)
}

\arguments{
\item{proID_1}{integer ID of a proband}
\item{proID_2}{integer ID of the second proband}
\item{BP_len}{The chromosome BP length used in the simulation}
\item{proband_haplotypes_path}{String that gives the location of the "Proband_haplotypes.txt" output file produced by gen.simuHaplo}
}

\value{
returns a dataframe with the following columns: "simulNo", "n_seg", "pIBD", "mean_seg_len". Each row is a simulation wherein the pair of probands had some IBD sharing.
The "simulNo" columns gives the simulation number, "n_seg" is the number of non-contiguous segments shared between the pair (divide this value by 2 for an average value for the pair), 
"pIBD" is the percent of the simulated chromosome shared IBD between the pair, and "mean_seg_len" is the mean length of IBD segments.
Additionally, the function will print the exact location of each IBD segment to the R console. This can be toggled by disabling messages.
}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.simuHaplo}}
\code{\link{gen.simuHaplo_traceback}}
}
\examples{

gen_table <- data.frame(ind=c(1,2,3,4,5,6,7,8),
                        mother=c(0,0,2,0,0,5,3,3),
                        father=c(0,0,1,0,0,4,6,6),
                        sex=c(1,2,2,1,2,1,1,2)
)

#convert the genealogy table into a GENLIB object
sample_gen<-gen.genealogy(gen_table) 

gen.graph(sample_gen)
# Simulate the gene dropping from founders to probamds

gen.simuHaplo(sample_gen, simulNo=5, model = 1, model_params = c(1,1), cM_len=c(100,100),
                BP_len = 100000000, outDir=tempdir())

IBD_compare <- gen.simuHaplo_IBD_compare(7, 8, 100000000, 
                paste(tempdir(), "Proband_Haplotypes.txt", sep="/"))
}
