\name{Classes of GLgroupe handling}
\Rdversion{1.1}
\docType{class}
\alias{GLmultiPhiGroup-class}
\alias{GLmultiPhiGroupSingle-class}
\alias{GLmultiFGroup-class}
\alias{GLmultiFGroupSingle-class}
\alias{GLCGMatrixGroupSingle-class}
\alias{[<-,GLCGMatrixGroupSingle,ANY,ANY,ANY-method}
\alias{[,GLCGMatrixGroupSingle,ANY,ANY,ANY-method}
\alias{[<-,GLmultiFGroup,ANY,ANY,ANY-method}
\alias{[,GLmultiFGroup,ANY,ANY,ANY-method}
\alias{Dim,GLmultiFGroup-method}
\alias{[<-,GLmultiFGroupSingle,ANY,ANY,ANY-method}
\alias{[,GLmultiFGroupSingle,ANY,ANY,ANY-method}
\alias{Dim,GLmultiFGroupSingle-method}
\alias{[<-,GLmultiPhiGroupSingle,ANY,ANY,ANY-method}
\alias{[,GLmultiPhiGroupSingle,ANY,ANY,ANY-method}
\alias{Dim,GLmultiPhiGroupSingle-method}
\alias{[<-,GLmultiPhiGroup,ANY,ANY,ANY-method}
\alias{[,GLmultiPhiGroup,ANY,ANY,ANY-method}
\alias{Dim,GLmultiPhiGroup-method}

\title{Classes \code{"GLmultiPhiGroup"} \code{"GLmultiPhiGroupSingle"} \code{"GLmultiFGroup"} \code{"GLmultiFGroupSingle"} \code{"GLCGMatrixGroupSingle"} }
\description{
 Objects used to carry information mostly between invisible functions.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \cr
\code{new("GLmultiPhiGroup", MatriceCG, groupe = 'GLgroupe', grindice = 'list')} \cr
\code{new("GLmultiPhiGroupSingle", MatriceCG, groupe = 'GLgroupe', grindice = 'list')} \cr
\code{new("GLmultiFGroup", MatriceCG, groupe = 'GLgroupe', grindice = 'list')} \cr
\code{new("GLmultiFGroupSingle", MatriceCG, groupe = 'GLgroupe', grindice = 'list')} \cr
\code{new("GLCGMatrixGroupSingle", MatriceCG, groupe = 'GLgroupe', grindice = 'list')}

}
\section{Slots}{
  \describe{
    GLCGMatrixGroupSingle, GLmultiFGroupSingle, GLmultiPhiGroupSingle
    \item{\code{.Data}:}{Object of class \code{"matrix"} }
    \item{\code{groupe}:}{Object of class \code{"GLgroupe"} }
    \item{\code{grindice}:}{Object of class \code{"list"} containing the indices of the probands of 'groupe' }
    
    GLmultiPhiGroup :
    \item{\code{.Data}:}{Object of class \code{"GLmultiMatrix"} }
    \item{\code{groupe}:}{Object of class \code{"GLgroupe"} }
    \item{\code{grindice}:}{Object of class \code{"list"} containing the indices of the probands of 'groupe' }
    
    GLmultiFGroup :
    \item{\code{.Data}:}{Object of class \code{"GLmultiVector"} }
    \item{\code{groupe}:}{Object of class \code{"GLgroupe"} }
    \item{\code{grindice}:}{Object of class \code{"list"} containing the indices of the probands of 'groupe' }
  }
}
\section{Extends}{
Class \code{"\linkS4class{matrix}"}, from data part.
Class \code{"\linkS4class{array}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{structure}"}, by class "matrix", distance 3.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 4, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{[<-}{\code{signature(x = "GLCGMatrixGroupSingle", i = "ANY", j = "ANY", value = "ANY")}: ... }
    \item{[}{\code{signature(x = "GLCGMatrixGroupSingle", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{groupe}{\code{signature(x = "GLCGMatrixGroupSingle")}: ... }
	 }
}
\examples{
showClass("GLCGMatrixGroupSingle")
showClass("GLmultiFGroup")
showClass("GLmultiFGroupSingle")
showClass("GLmultiPhiGroupSingle")
showClass("GLmultiPhiGroup")
}
\keyword{classes}
