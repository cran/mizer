% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveParams.R
\name{saveParams}
\alias{saveParams}
\alias{readParams}
\title{Save a MizerParams object to file, and restore it}
\usage{
saveParams(params, file)

readParams(file)
}
\arguments{
\item{params}{A MizerParams object}

\item{file}{The name of the file or a connection where the MizerParams object
is saved to or read from.}
}
\value{
NULL invisibly
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
\code{saveParams()} saves a MizerParams object to a file. This can then be
restored with \code{readParams()}.
}
\details{
Issues a warning if the model you are saving relies on some custom functions.
Before saving a model you may want to set its metadata with
\code{\link[=setMetadata]{setMetadata()}}.
}
