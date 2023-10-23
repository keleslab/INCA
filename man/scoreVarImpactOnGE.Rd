\name{scoreVarImpactOnGE}
\alias{scoreVarImpactOnGE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Score variants on its impact on gene expression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoreVarImpactOnGE(variants, CellLine1, CellLine2, CellNames, window = 15, parallel = FALSE)
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
    variants[, `:=`(Gene, gsub("\\s+", "", Gene))]
    variants[Gene == "", "Gene"] = "-"
    variants = variants[, .(Gene = unlist(strsplit(Gene, "[^[:alnum:]\\-\\_\\.]",
        perl = TRUE))), by = eval(grep("^Gene$", colnames(variants),
        invert = TRUE, value = TRUE))]
    if (!"list" \%in\% class(CellLine1)) {
        CellLine1 = list(CellLine1)
    }
    if (!"list" \%in\% class(CellLine2)) {
        CellLine1 = list(CellLine2)
    }
    de = data.table(mapDEResultToGenes(variants, CellLine1[[1]],
        CellNames[1]), mapDEResultToGenes(variants, CellLine2[[1]],
        CellNames[2]))
    variants[, "DESignificant"] = as.character(factor(apply(de[,
        grep("significant", colnames(de)), with = FALSE], 1,
        paste, collapse = ":"), levels = c("NA:NA", "0:0", "0:1",
        "1:0", "1:1"), labels = c(NA, "Neither", CellNames, "Both")))
    if ("AltCellLine" \%in\% colnames(variants)) {
        ALT = variants[, AltCellLine] \%in\% CellNames
    }
    else {
        ALT = TRUE
    }
    if (!"INCAscore" \%in\% colnames(variants)) {
        variants[, "INCAscore"] = 0
    }
    variants[ALT & DESignificant \%in\% CellNames, `:=`(INCAscore,
        INCAscore + 1)]
    if (length(CellLine1) > 1 & length(CellLine2) > 1) {
        enriched = matrix(0, ncol = 2, nrow = nrow(variants))
        enriched[varInEnrichedRegion(variants, CellLine1[-1],
            window, parallel), 1] = 1
        enriched[varInEnrichedRegion(variants, CellLine2[-1],
            window, parallel), 2] = 1
        variants[, "DEBinding"] = as.character(factor(paste(enriched[,
            1], enriched[, 2], sep = ":"), levels = c("0:0",
            "0:1", "1:0", "1:1"), labels = c("Neither", CellNames,
            "Both")))
        variants[ALT & DESignificant \%in\% CellNames & DESignificant ==
            DEBinding, `:=`(INCAscore, INCAscore + 1)]
    }
    else if ("EnrichedRegion" \%in\% colnames(variants)) {
        variants[ALT & DESignificant \%in\% CellNames & DESignificant ==
            EnrichedRegion, `:=`(INCAscore, INCAscore + 1)]
    }
    else {
        message("No `peak` or `counts` files provided for assessing binding at the variant location.\n        Max score for RBP-SNV impact on gene expression is 1.")
    }
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
