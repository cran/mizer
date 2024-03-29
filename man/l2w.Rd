% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{l2w}
\alias{l2w}
\alias{w2l}
\title{Length-weight conversion}
\usage{
l2w(l, params)

w2l(w, params)
}
\arguments{
\item{l}{Lengths in cm. Either a single number or a vector with one number
for each species.}

\item{params}{A species parameter data frame or a MizerParams object.}

\item{w}{Weights in grams. Either a single number or a vector with one number
for each species.}
}
\value{
A vector with one entry for each species. \code{l2w()} returns a vector
of weights in grams and \code{w2l()} returns a vector of lengths in cm.
}
\description{
Convert between length and weight using the relationship
\deqn{w = a l^b}
where \code{a} and \code{b} are taken from the species parameter data frame.
}
\details{
This is useful for converting a length-based species parameter to a
weight-based species parameter.

If any \code{a} or \code{b} parameters are missing the default values \code{a = 0.01} and
\code{b = 3} are used for the missing values.
}
\concept{helper}
