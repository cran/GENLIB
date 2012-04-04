\name{gen.half.founder}
\alias{gen.half.founder}
\title{Get half-founder id numbers}
\description{Returns the id numbers of the half-founders. Half-founders are defined as the individuals with only one founder parent in the genealogy.}
\usage{gen.half.founder( gen, check = "1", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{check}{If check=1 (default), all parameters are validated.}
\item{...}{Additionnal arguments to be passed to methods.}
}
\value{returns a vector of integer}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.pro}}
\code{\link{gen.founder}}
\code{\link{gen.parent}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.half.founder(genJi)  #There are 2 half-founders
}
\keyword{manip}
