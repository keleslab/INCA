\name{compareRCToControl}
\alias{compareRCToControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate aligned read counts with respect to control
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculate aligned read counts in the experiment by comparing it the a control sample, through either subtraction or fold change.
}
\usage{
compareRCToControl(experiment, control, method = 1, pseudocount = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experiment}{
%%     ~~Describe \code{experiment} here~~
  A data frame or a list of Rle objects, containing the information of the experiment. If \code{experiment} is a data frame, it must contain three columns in the order of chromosomes, lengths, and aligned read counts. For more details, refer to \code{\link[INCAData]{summarizeReadCounts}}. If \code{experiment} is a list of Rle objects, each list element should be named after the corresponding chromosome, and \var{values} should represent aligned read counts. For more details, refer to \code{\link[S4Vectors]{Rle}}.
}
  \item{control}{
%%     ~~Describe \code{control} here~~
  A data frame or a list of Rle objects, similar to \code{experiment}, but containing information of the control sample.
}
  \item{method}{
%%     ~~Describe \code{method} here~~
  Option 1 or 2. Option 1 indicates subtraction method, \eqn{experiment - control}. Option 2 indicates fold change method, \eqn{(experiment + pseudocount) / (control + pseudocount)}.
}
  \item{pseudocount}{
%%     ~~Describe \code{pseudocount} here~~
  A positive numeric value for the pseudocount in the fold change calculation to avoid division by zero. Default is 0.1.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Despite both raw and normalized aligned read counts being accepted for \code{experiment} and \code{control}, the aligned read counts should be normalized to make them comparable.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of Rle objects, with each list element named after its corresponding chromosome. \var{values} represent aligned read counts in the experiment with respect to the control sample.
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
exp = data.table::fread('https://raw.github.com/jduan607/INCA/master/ENCODE_eCLIP/NormRC/HNRNPK_HepG2_NormRC_Rep1.txt.gz')
ctrl = data.table::fread('https://raw.github.com/jduan607/INCA/master/ENCODE_eCLIP/NormRC/HNRNPK_HepG2_NormRC_Ctrl.txt.gz')
rc = compareRCToControl(exp, ctrl, method=1)
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
