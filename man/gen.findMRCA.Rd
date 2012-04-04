\name{gen.findMRCA}
\alias{gen.findMRCA}
\title{Finding most recent common ancestors, MRCAs}
\description{Returns MRCAs of pairs of specified individuals and the distance (number of meioses) between individuals through the MRCAs.}
\usage{gen.findMRCA(gen, individuals, NbProcess=detectCores()-1)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{individuals}{Vector of individual id numbers for which to find the MRCAs. Required}
\item{NbProcess}{Number of processes to use when running this function. Default=detectCores()-1}
}
\value{returns a matrix}
\note{This function uses the parallel programming functions of foreach, snow and doSNOW for launching several processes.}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.founder}}
\code{\link{gen.findFounders}}
\code{\link{gen.findDistance}}
\code{\link{gen.find.Min.Distance.MRCA}}
}
\examples{
#data(geneaJi) 
#genJi<-gen.genealogy(geneaJi) 
#gen.findMRCA(genJi, individuals=c(1,29), NbProcess = 1)

data(genea140) 
#gen140<-gen.genealogy(genea140) 
#gen2<-gen.branching(gen140,pro=c(409033,408728)) 
#all_founders<-gen.founder(gen2) #615 founders in genealogy of #409033 and #408728
#length(all_founders)

all_commonFounders<-gen.findFounders(gen140, individuals=c(409033,408728), NbProcess = 1)
length(all_commonFounders) #127 founders common to #409033 and #408728
#22 most recent common ancestors for #409033 and #408728
MRCA_2ind<-gen.findMRCA(gen140, individuals=c(409033,408728), NbProcess = 1)

#MRCA_2loy<-gen.findMRCA(gen140, individuals=c(409437,409422), NbProcess = 1)
#6 most recent common ancestors for #409437 and #409422
#genMRCA<-gen.branching(gen140, pro=c(409437,409422),ancestors=as.numeric(colnames(MRCA_2loy)))
#gen.graph(genMRCA)
}
\keyword{manip}
