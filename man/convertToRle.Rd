\name{convertToRle}
\alias{convertToRle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Convert data frame(s) to a list of Rle objects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Convert a data frame (or a list contains data frames and Rle objects) to a list of Rle objects.
}
\usage{
convertToRle(input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
  A data frame (or a list contains data frames and Rle objects) with three columns in the order of chromosomes, lengths, and values. For details on the lengths and values, refer to \code{\link[S4Vectors]{Rle}}.
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
A list of Rle objects, with each list element named after its corresponding chromosome.
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
data = data.table::fread('https://raw.github.com/jduan607/INCA/master/ENCODE_eCLIP/NormRC/HNRNPK_HepG2_NormRC_Rep1.txt.gz')
data_rle = convertToRle(data)
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
