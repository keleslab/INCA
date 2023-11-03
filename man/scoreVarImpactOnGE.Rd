\name{scoreVarImpactOnGE}
\alias{scoreVarImpactOnGE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Scores for the RBP-variant impact on the gene expression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Score the variant harboring genes to assess whether the gene might be a target of the RBP, binding of which the variant might be impacting.
}
\usage{
scoreVarImpactOnGE(variants, CellLine1, CellLine2, CellNames=paste('Cell type',1:2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variants}{
%%     ~~Describe \code{variants} here~~
A data frame of variants containing columns of \var{Chr}, \var{Start}, and \var{End}.
}
  \item{CellLine1}{
%%     ~~Describe \code{CellLine1} here~~
A list with elements in the following order:
  \enumerate{
    \item A data frame of the differential expression analysis results. The first column is \var{gene}. The second column should contain p-values or q-values from a differential test assessing the significance of changes in expression for genes between control and knockdown conditions.
    \item A list with the following elements:
    \itemize{
      \item \var{peaks}: A data frame of peaks, with the first three columns indicating the genomic positions of the peaks in the order of chromosome, start, and end.
      \item \var{threshold}: Optional. If \var{threshold} is provided, then \code{selectPeaksBySignals} will be used to select IDR peaks and non-IDR peaks based on the provided \var{threshold} from \var{peaks}. For more details, refer to \code{\link[INCA]{selectPeaksBySignals}}.
    }

  }
}
  \item{CellLine2}{
%%     ~~Describe \code{CellLine2} here~~
 The same format as \code{CellLine1} with datasets related to cell line/type 2.
}
  \item{CellNames}{
%%     ~~Describe \code{CellNames} here~~
  A length-2 character vector of cell line/type names.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Leveraging the differential expression analysis results of RNA-seq experiments of RNA-binding protein (RBP) knockdowns by shRNA and wild type conditions, the RBP-SNV gene impact score evaluates whether (i) whether a gene’s expression changes significantly following RBP knockdown in one cell line but not the other, and (ii) whether the cell line with significant gene expression change, but not the other cell line, has an RBP binding at the variant location. When both conditions hold, it broadly implies that the encoded transcript is a binding-target of the RBP with expression consequences in that cell line, leading to a \var{DEscore} of 2. If only the former condition holds, the \var{DEscore} is 1.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A new data frame incorporating the contents of \code{variants} along with the following changes:
  \itemize{
    \item \var{Gene} is separated such that each row contains one gene.
    \item \var{DEgene} representing the cell line/type that exhibits significant changes in gene expression upon RNA-binding protein (RBP) knockdown is added.
    \item \var{DEbinding} representing the cell line/type in which genes harbor RBP binding sites is added.
    \item \var{DEscore} representing RBP-SNV impact on gene expression is added.
    \item \var{DEscore} is added to \var{INCAscore}. If \code{variants} does not already contain \var{INCAscore}, \var{INCAscore} is added and set to 0 for all variants initially. If \code{variants} contains \var{AltCellLine}, then \var{DEscore} is only added to variants matched to the genotype in one of the cell lines indicated in \var{AltCellLine}. For more details on \var{AltCellLine}, refer to \code{\link[INCA]{alignCellGenoToVar}}.
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

# Differential expression analysis
deg1 = fread(file.path(directory,'ENCODE_shRNA/DEG','HNRNPK_HepG2_DEG.txt.gz'))
deg2 = fread(file.path(directory,'ENCODE_shRNA/DEG','HNRNPK_K562_DEG.txt.gz'))

epg1 = list(deg1[,.(gene,q_value)], # first element must be DE
            list(peaks=peak1, threshold=1))
epg2 = list(deg2[,.(gene,q_value)], # first element must be DE
            list(peaks=peak2, threshold=1))

variants = scoreVarImpactOnGE(variants, epg1, egp2, c('HepG2','K562'))
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
