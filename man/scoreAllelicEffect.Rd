\name{scoreAllelicEffect}
\alias{scoreAllelicEffect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Score of variants on allelic effects
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
}
  \item{CellLine1}{
%%     ~~Describe \code{CellLine1} here~~
}
  \item{CellLine2}{
%%     ~~Describe \code{CellLine2} here~~
}
  \item{CellNames}{
%%     ~~Describe \code{CellNames} here~~
}
  \item{window}{
%%     ~~Describe \code{window} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
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
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (variants, CellLine1, CellLine2, CellNames, window = 15,
    parallel = FALSE)
{
    variants = as.data.table(variants)
    enriched = matrix(0, ncol = 2, nrow = nrow(variants))
    enriched[varInEnrichedRegion(variants, CellLine1, window,
        parallel), 1] = 1
    enriched[varInEnrichedRegion(variants, CellLine2, window,
        parallel), 2] = 1
    variants[, "EnrichedRegion"] = as.character(factor(paste(enriched[,
        1], enriched[, 2], sep = ":"), levels = c("0:0", "0:1",
        "1:0", "1:1"), labels = c("Neither", CellNames, "Both")))
    if ("AltCellLine" \%in\% colnames(variants)) {
        ALT = variants[, AltCellLine] \%in\% CellNames
    }
    else {
        ALT = TRUE
    }
    if (!"INCAscore" \%in\% colnames(variants)) {
        variants[, "INCAscore"] = 0
    }
    variants[ALT & EnrichedRegion \%in\% CellNames, `:=`(INCAscore,
        INCAscore + 1)]
    return(variants)
  }
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
