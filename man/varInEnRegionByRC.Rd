\name{varInEnRegionByRC}
\alias{varInEnRegionByRC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Variants that have aligned read counts exceeding the cutoff
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{varInEnRegionByPeaks} is used in the \code{\link[INCA]{varInEnrichedRegion}} function. It provides indices of variants that have aligned read counts exceeding the provided cutoff.
}
\usage{
varInEnRegionByRC(vdata, counts, cutoff, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vdata}{
%%     ~~Describe \code{vdata} here~~
  A data table of variants, with the first three columns indicating the genomic positions of the variants in the order of \var{Chr}, \var{Start}, and \var{End}.
}
  \item{counts}{
%%     ~~Describe \code{counts} here~~
  A data frame or a list of Rle objects. If \code{counts} is a data frame, it must contain three columns in the order of chromosomes, lengths, and aligned read counts. For more details, refer to \code{\link[INCAData]{summarizeReadCounts}}. If \code{counts} is a list of Rle objects, each list element should be named after the corresponding chromosome, and \var{values} should represent aligned read counts. For more details, refer to \code{\link[S4Vectors]{Rle}}.
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
  A numeric value.
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
  A logical value specifying whether or not the computation should be done in parallel. Default is \code{FALSE}. DO NOT set \code{parallel} to \code{TRUE} unless the multicore parallel backend is registered for multiple cores, for example, \code{doMC::registerDoMC()}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A vector of indices of variants that have aligned read counts exceeding the provided cutoff.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Jingqi Duan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
