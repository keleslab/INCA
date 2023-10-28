\name{scoreAllelicEffect}
\alias{scoreAllelicEffect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Assign scores to variants on allelic effects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoreAllelicEffect(variants, CellLine1, CellLine2, CellNames, window = 15, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variants}{
%%     ~~Describe \code{variants} here~~
  A data frame of variants containing columns of \var{Chr}, \var{Start}, and \var{End}. If \code{variants} includes \var{INCAscore}, then the score for allelic effect will be added to \var{INCAscore}. If \code{variants} includes \var{AltCellLine}, then the score See more details for \var{AltCellLine} in  \code{\link[INCA]{alignCellGenoToVar}}.
}
  \item{CellLine1}{
%%     ~~Describe \code{CellLine1} here~~
  A list of lists
}
  \item{CellLine2}{
%%     ~~Describe \code{CellLine2} here~~
}
  \item{CellNames}{
%%     ~~Describe \code{CellNames} here~~
}
  \item{window}{
%%     ~~Describe \code{window} here~~
  A non-negative integer representing the half-window size, where \eqn{newStart = Start - window} and \eqn{newEnd = End + window}. Default is 15.
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
epg1 = list(list(peaks=peak1, threshold=0.5),
            list(counts=rc1.1, peaks=peak1[signalValue_IDR>0,], threshold=0.8),
            list(counts=rc1.2, peaks=peak1[signalValue_IDR>0,], threshold=0.8))

epg2 = list(list(peaks=peak2, threshold=0.5),
            list(counts=rc2.1, peaks=peak2[signalValue_IDR>0,], threshold=0.8),
            list(counts=rc2.2, peaks=peak2[signalValue_IDR>0,], threshold=0.8))

variants = scoreAllelicEffect(variants, epg1, epg2, c('HepG2','K562'), parallel=TRUE)
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
