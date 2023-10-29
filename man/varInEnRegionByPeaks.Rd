\name{varInEnRegionByPeaks}
\alias{varInEnRegionByPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Variants that overlap with the specified peaks
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{varInEnRegionByPeaks} is used in the \code{\link[INCA]{varInEnrichedRegion}} function. It quantifies whether variants overlap with the specified peaks.
}
\usage{
varInEnRegionByPeaks(vdata, peaks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vdata}{
%%     ~~Describe \code{vdata} here~~
  A data table of variants, with the first three columns indicating the genomic positions of the variants in the order of \var{Chr}, \var{Start}, and \var{End}.
}
  \item{peaks}{
%%     ~~Describe \code{peaks} here~~
  A data frame of peaks, with the first three columns indicating the genomic positions of the peaks in the order of chromosome, start, and end.
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
A binary vector with 1 indicating variants that overlap with the specified peaks and 0 indicating variants that do not overlap.
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
