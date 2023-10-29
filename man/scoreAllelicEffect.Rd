\name{scoreAllelicEffect}
\alias{scoreAllelicEffect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Scores for the allelic effect of variants
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Quantify whether the variant exhibits allelic variation on RBP binding between the cell lines.
}
\usage{
scoreAllelicEffect(variants, CellLine1, CellLine2, CellNames=paste('Cell type',1:2), window=15, parallel=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variants}{
%%     ~~Describe \code{variants} here~~
  A data frame of variants containing columns of \var{Chr}, \var{Start}, and \var{End}.
}
  \item{CellLine1}{
%%     ~~Describe \code{CellLine1} here~~
  A list of lists following the \code{input} format from \code{\link[INCA]{varInEnrichedRegion}}, containing datasets related to cell line/type 1.
}
  \item{CellLine2}{
%%     ~~Describe \code{CellLine2} here~~
  A list of lists following the \code{input} format from \code{\link[INCA]{varInEnrichedRegion}}, containing datasets related to cell line/type 2.
}
  \item{CellNames}{
%%     ~~Describe \code{CellNames} here~~
  A length-2 character vector of cell line/type names.
}
  \item{window}{
%%     ~~Describe \code{window} here~~
  A non-negative integer representing the half-window size, where \eqn{newStart = Start - window} and \eqn{newEnd = End + window}. Default is 15.
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
  A logical value specifying whether or not the computation should be done in parallel. Default is \code{FALSE}. DO NOT set \code{parallel} to \code{TRUE} unless the multicore parallel backend is registered, for example, \code{doMC::registerDoMC()}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The allelic effect score (\var{AEscore}) quantifies whether the variant resides in a eCLIP-seq enriched region in one cell line but not the other and, hence, indicating allelic variation between the cell lines. This score is set to 1 if (1) the variant of interest resides within a binding region of the RBP as defined by an IDR peak in one cell line, but not in the other cell line, or a value less than 1 (typically \var{threshold} value) if (2) the variant of interest resides in a region of significant enrichment with respect to the normalized read counts, with or without repsect to the control sample, in one cell line, but not the other cell line.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A new data frame incorporating the contents of \code{variants} along with the following changes:
  \itemize{
    \item \var{EnrichedRegion} representing the cell line/type in which the variants resides within an RNA-binding protein (RBP) binding region is added.
    \item \var{AEscore} representing the allelic effect of variants on RBP binding is added.
    \item \var{AEscore} is added to \var{INCAscore}. If \code{variants} does not already contain \var{INCAscore}, \var{INCAscore} is added and set to 0 for all variants initially. If \code{variants} contains \var{AltCellLine}, then \var{AEscore} is only added to variants matched to the genotype in one of the cell lines indicated in \var{AltCellLine}. For more details on \var{AltCellLine}, refer to \code{\link[INCA]{alignCellGenoToVar}}.
  }
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
directory = 'https://raw.github.com/jduan607/INCA/master'
variants = fread(file.path(directory, 'GWAS', 'final_data.txt'))

# Peak
peak1 = fread(file.path(directory,'ENCODE_eCLIP/PeakSignals','HNRNPK_HepG2_PeakSignals.txt.gz'))
peak2 = fread(file.path(directory,'ENCODE_eCLIP/PeakSignals','HNRNPK_K562_PeakSignals.txt.gz'))

# Normalized aligned read counts
exp1 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_HepG2_NormRC_rep1.txt.gz'))
exp2 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_HepG2_NormRC_rep2.txt.gz'))

epg1 = list(list(peaks=peak1, threshold=1),
            list(counts=exp1, peaks=peak1[signalValue_IDR>0,], threshold=c(0.1,0.8)))

epg2 = list(list(peaks=peak2, threshold=1),
            list(counts=exp2, peaks=peak2[signalValue_IDR>0,], threshold=c(0.1,0.8)))

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
