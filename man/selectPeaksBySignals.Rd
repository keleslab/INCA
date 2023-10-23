\name{selectPeaksBySignals}
\alias{selectPeaksBySignals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Select IDR peaks and non-IDR peaks with signal values exceeding the threshold
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Select Irreproducible Discovery Rate (IDR) peaks and non-IDR peaks with signal values exceeding the quantile threshold of the signal values of the reported IDR peaks.
}
\usage{
selectPeaksBySignals(peaks, threshold = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{
%%     ~~Describe \code{peaks} here~~
A data frame with rows of peaks. \code{peaks} should contain columns indicating genomic positions of the peaks (chromosome, start, and end) and columns of signal values from peaks calls of replicate experiments and IDR. For details, refer to \link[INCAData]{summarizePeakSignals}.
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
A numeric value in [0,1]. Default is 0.5.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
A peak is considered an IDR peak if its value in the \var{IDR} column (with a name containing "IDR") is positive. The \var{IDR} column can be either binary or continuous, as long as positive values indicate IDR peaks. Next, for each \var{replicate} column (with names containing "signal"), two steps are followed: (i) the cutoff is determined as the quantile (provided in \code{threshold}) of signal values from this replicate of IDR peaks, and (ii) the non-IDR peaks with signal values from this replicate exceeding the cutoff are selected. Finally, the selected non-IDR peaks from each replicate and non-IDR peaks are extracted.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data.table consisting of IDR peaks and non-IDR peaks with signal values exceeding the quantile threshold, containing the same columns as \code{peaks}.
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
peak = fread('https://raw.github.com/jduan607/INCA/master/ENCODE_eCLIP/PeakSignals/HNRNPK_HepG2_PeakSignals.txt.gz')
new_peak = selectPeaksBySignals(peak, threshold=0.5)
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
