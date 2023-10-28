\name{varInEnrichedRegion}
\alias{varInEnrichedRegion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Variants that reside in an enriched region.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{varInEnrichedRegion} is used in the \code{\link[INCA]{scoreAllelicEffect}} and \code{\link[INCA]{scoreVarImpactOnGE}} functions. It provides indices of variants that reside in an enriched region.
}
\usage{
varInEnrichedRegion(vdata, input, window = 15, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vdata}{
%%     ~~Describe \code{vdata} here~~
  A data frame of variants containing columns of \var{Chr}, \var{Start}, and \var{End}.
}
  \item{input}{
%%     ~~Describe \code{input} here~~
  A list of lists, where each inner list is expected to contain elements in either of the following formats:
  \enumerate{
    \item \itemize{
      \item \var{peaks}: A data frame of peaks, with the first three columns indicating the  genomic positions of the peaks in the order of chromosome, start, and end.
      \item \var{threshold}: Optional. If \var{threshold} is provided, then \code{selectPeaksBySignals} will be used to select IDR peaks and non-IDR peaks based on the provided \var{threshold} from \var{peaks}. For more details, refer to \code{\link[INCA]{selectPeaksBySignals}}.
  }
    \item \itemize{
      \item \var{counts}: A data frame or a list of Rle objects. If \code{counts} is a data frame, it must contain three columns in the order of chromosomes, lengths, and aligned read counts. For more details, refer to \code{\link[INCAData]{summarizeReadCounts}}. If \code{counts} is a list of Rle objects, each list element should be named after the corresponding chromosome, and \var{values} should represent aligned read counts. For more details, refer to \code{\link[S4Vectors]{Rle}}.
      \item \var{cutoff}: A numeric value of cutoff for aligned read counts. If both \var{cutoff} and \var{peaks} are not provided, the default cutoff of 30 will be used.
      \item \var{peaks}: Optional. If \var{cutoff} is not provided and \var{peaks} is provided, then \code{computeRCQuantile} will be used to calculate the cutoff.
      \item \var{threshold}: Optional. If \code{computeRCQuantile} is used to calculate the cutoff, the cutoff will be based on either the provided \var{threshold} or the default threshold of 0.5. For more details, refer to \code{\link[INCA]{computeRCQuantile}}.
    }
  }}
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
If \var{peaks} is provided in \code{input}, variants are considered to be located in an enriched region if they overlap with the specified peaks.
If \var{counts} is provided in \code{input}, variants are considered to be located in an enriched region if they have aligned read counts exceeding the provided \var{cutoff} or the default cutoff of 30.
\code{input} accepts list of lists, allowing enriched regions to be determined using either \var{peaks} or \var{counts} and files from multiple replicates.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of indices, with each element indicating variants that reside in an enriched region defined in the corresponding element of the \code{input}.
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
# Example variants
variants = data.table::fread('https://raw.github.com/jduan607/INCA/master/GWAS/final_data.txt')[,.(Chr,Start,End)]

# Peaks
peak = data.table::fread('https://raw.github.com/jduan607/INCA/master/ENCODE_eCLIP/PeakSignals/HNRNPK_HepG2_PeakSignals.txt.gz')
idr_peak = peak[signalValue_IDR>0, ]

# Normalized aligned read counts
exp1 = data.table::fread('https://raw.github.com/jduan607/INCA/master/ENCODE_eCLIP/NormRC/HNRNPK_HepG2_NormRC_Rep1.txt.gz')
exp2 = data.table::fread('https://raw.github.com/jduan607/INCA/master/ENCODE_eCLIP/NormRC/HNRNPK_HepG2_NormRC_Rep2.txt.gz')

# IDR peaks only
ls1 = list(peaks=idr_peak)

# IDR peaks and some non-IDR peaks selectPeaksBySignals()
ls2 = list(peaks=peak, threshold=0.5)

# Read counts and cutoff
ls3 = list(counts=exp1, cutoff=35)
ls4 = list(counts=exp2, cutoff=40)

# Read counts and cutoff from computeRCQuantile()
ls5 = list(counts=exp1, peaks=idr_peak, threshold=0.5)
ls6 = list(counts=exp2, peaks=idr_peak, threshold=0.5)

indices1 = varInEnrichedRegion(variants, list(ls1,ls3,ls4), window=15, parallel=FALSE)

indices2 = varInEnrichedRegion(variants, list(ls2,ls5,ls6), window=15, parallel=FALSE)
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
