\name{gen.nomen}
\alias{gen.nomen}
\title{Number of men}
\description{Returns the number of men included in the genealogy}
\usage{gen.nomen( gen, check = "1")}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{check}{If check=1 (default), all parameters are validated.}
}
\value{returns a vector of integer}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.nowomen}}
\code{\link{gen.nochildren}}
\code{\link{gen.noind}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.nomen(genJi)

}
\keyword{Utilities}
