\name{gen.founder}
\alias{gen.founder}
\title{Get founder id numbers}
\description{Returns the id numbers of the founders. Founders are defined as the individuals without parents in the genealogy (i.e., mother id=0, father id=0).}
\usage{gen.founder( gen, check = "1", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{check}{If check=1 (default), all parameters are validated.}
\item{...}{Additionnal arguments to be passed to methods.}
}
\value{A vector of integer}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.pro}}
\code{\link{gen.half.founder}}
\code{\link{gen.parent}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.founder(genJi) #There are 6 founders

data(genea140) 
gen140<-gen.genealogy(genea140)
founder140<-gen.founder(gen140)  
length(founder140) #There are 7399 founders

}
\keyword{manip}
