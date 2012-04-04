\name{GLgroupe-class}
\docType{class}
\alias{GLgroupe-class}
\alias{initialize,GLgroupe-method}
\alias{[,GLgroupe,ANY,ANY,ANY-method}

\title{Class \code{"GLgroupe"}}
\description{
 Object representing a set of proband in different groups.\cr
 Each element of this list is a group named with the proband number.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GLgroupe", 'list')}.
}
\section{Extends}{
Class \code{"\linkS4class{GLgen}"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Data=list)}: ... }
    \item{[}{\code{signature(object = "GLgroupe",ANY,ANY,ANY)}:... }
 }
}

\examples{
showClass("GLgroupe")
}
\keyword{classes}
