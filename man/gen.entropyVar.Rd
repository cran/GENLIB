\name{gen.entropyVar}
\alias{gen.entropyVar}
\title{Variance of entropy index}
\description{Returns the variance of the entropy index}
\usage{gen.entropyVar( gen, pro = "0", type = "MEAN", check = "1", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 0, which will select all individuals without children.}
\item{type}{If type="MEAN" (default), the average of entropy variances (over all probands) is returned.
	If type="IND", the variance of the entropy is calculated for each specified proband.}
\item{check}{If check=1 (default), all parameters are validated.}
\item{...}{Additionnal arguments to be passed to methods.}
}
\value{A data frame with only one numeric value when type is "MEAN".
 	  When type is "IND", the number of rows equals the number of probands specified.}
\seealso{
\code{\link{gen.gc}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.entropy(genJi, type="IND")
gen.entropyVar(genJi, type="IND")
gen.entropyVar(genJi, type="MEAN")

data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.entropy(gen140)
gen.entropy(gen140, pro=c(454422, 676521, 677273, 717537, 717634, 717709, 868572))
gen.entropyVar(gen140, pro=c(454422, 676521, 677273, 717537, 717634, 717709, 868572))
gen.entropyVar(gen140, pro=c(454422, 676521, 677273, 717537, 717634, 717709, 868572), type="MEAN")

#gen.entropy(gen140, pro=c(454422, 676521, 677273, 717537, 717634, 717709, 868572), type="IND")
#gen.entropyVar(gen140, pro=c(454422, 676521, 677273, 717537, 717634, 717709, 868572), type="IND")
}
\references{
Cazes P, Cazes MH. (1996) Comment mesurer la profondeur genealogique d'une ascendance? Population (French Ed) 51:117-140.
Kouladjian K. (1986) Une mesure d'entropie genealogique. Chicoutimi, SOREP, Document III-C-43.
De Brakaeleer M, Bellis G. (1994) Genealogies et reconstitutions de familles en genetique humaine. Dossiers et Recherches, no 43, INED, Paris.
}
\keyword{manip}

