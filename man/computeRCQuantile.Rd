\name{computeRCQuantile}
\alias{computeRCQuantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate the quantile of aligned read counts
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculate the quantile at the provided threshold of aligned read counts for variants within the specified peaks.
}
\usage{
computeRCQuantile(counts, peaks, threshold = 0.5, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
%%     ~~Describe \code{counts} here~~
  A data frame or a list of Rle objects. If \code{counts} is a data frame, it must contain three columns in the order of chromosomes, lengths, and aligned read counts. For more details, refer to \code{\link[INCAData]{summarizeReadCounts}}. If \code{counts} is a list of Rle objects, each list element should be named after the corresponding chromosome, and \var{values} should represent aligned read counts. For more details, refer to \code{\link[S4Vectors]{Rle}}.
}
  \item{peaks}{
%%     ~~Describe \code{peaks} here~~
  A data frame of peaks, with the first three columns indicating the genomic positions of the peaks in the order of chromosome, start, and end.
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
  A vector of numeric values in [0,1]. Default is 0.5.
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
  A logical value specifying whether or not the computation should be done in parallel. Default is \code{FALSE}. DO NOT set \code{parallel} to \code{TRUE} unless the multicore parallel backend is registered, for example, \code{doMC::registerDoMC()}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{counts} can accept either raw or normalized aligned read counts from the experiment, either with or without consideration of the control sample. It is recommended to use Irreproducible Discovery Rate (IDR) peaks for \code{peaks}. \code{computeRCQuantile} calculates the quantile of the aligned read counts for variants within the specified peaks based on the probability provided in \code{threshold}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A numeric value of the quantile correponding to the given threshold.
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
count = data.table::fread('https://raw.github.com/jduan607/INCA/master/ENCODE_eCLIP/NormRC/HNRNPK_HepG2_NormRC_Rep1.txt.gz')
peak = fread('https://raw.github.com/jduan607/INCA/master/ENCODE_eCLIP/PeakSignals/HNRNPK_HepG2_PeakSignals.txt.gz')[signalValue_IDR>0,]
cutoff = computeRCQuantile(count, peak, quantile=0.8, parallel=FALSE)
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
