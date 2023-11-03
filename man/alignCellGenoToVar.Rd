\name{alignCellGenoToVar}
\alias{alignCellGenoToVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Align genotypes of two cell lines with variants
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
alignCellGenoToVar(variants, CellLine1, CellLine2, CellNames)
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
function (variants, CellLine1, CellLine2, CellNames)
{
    variants = as.data.table(variants)
    if (!"Zyg" \%in\% colnames(variants)) {
        variants[, "Zyg"] = 1
        message("Zygosity not provided. Assume all variants are heterozygous.")
    }
    colnames(CellLine1) = c("Chr", "Pos", "Ref", paste("C1",
        c("Alt", "Alt2", "Zyg"), sep = "_"))
    CellLine1 = CellLine1[variants[, .(Chr, Pos, Ref)], on = .(Chr = Chr,
        Pos = Pos, Ref = Ref)][, -c(1:3)]
    colnames(CellLine2) = c("Chr", "Pos", "Ref", paste("C2",
        c("Alt", "Alt2", "Zyg"), sep = "_"))
    CellLine2 = CellLine2[variants[, .(Chr, Pos, Ref)], on = .(Chr = Chr,
        Pos = Pos, Ref = Ref)][, -c(1:3)]
    type = paste(apply(CellLine1, 1, function(x) any(!is.na(x)) *
        1), apply(CellLine2, 1, function(x) any(!is.na(x)) *
        1), sep = ":")
    cdata = cbind(variants[, .(Alt, Zyg)], CellLine1, CellLine2,
        AltCellLine = "Neither")
    cdata[type == "1:0" & Alt == C1_Alt & C1_Alt2 == "", "AltCellLine"] = CellNames[1]
    cdata[type == "0:1" & Alt == C2_Alt & C2_Alt2 == "", "AltCellLine"] = CellNames[2]
    cdata[type == "1:1" & Alt == C1_Alt & Alt == C2_Alt, "AltCellLine"] = "Both"
    cdata[AltCellLine == "Both" & Zyg == C1_Zyg & Zyg != C2_Zyg,
        "AltCellLine"] = CellNames[1]
    cdata[AltCellLine == "Both" & Zyg == C2_Zyg & Zyg != C1_Zyg,
        "AltCellLine"] = CellNames[2]
    cdata[C1_Alt2 != "", `:=`(C1_Alt, paste(C1_Alt, C1_Alt2,
        sep = ","))]
    cdata[C2_Alt2 != "", `:=`(C2_Alt, paste(C2_Alt, C2_Alt2,
        sep = ","))]
    cdata = cdata[, .(C1_Alt, C1_Zyg, C2_Alt, C2_Zyg, AltCellLine)]
    colnames(cdata) = gsub("C1", CellNames[1], colnames(cdata))
    colnames(cdata) = gsub("C2", CellNames[2], colnames(cdata))
    variants = cbind(variants, cdata)
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
