\name{GiG.Seq.Analysis}
\alias{GiG.Seq.Analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genomic sequence nucleotidic run analyser.
}
\description{
A function to analyse a genomic sequence in relationship with its nucleotidic composition and organisation, concretely regarding the G and C-runs within.
}
\usage{
GiG.Seq.Analysis(Name, Sequence, DNA = TRUE, Complementary = TRUE, Nucleotides = c("G", "C"), BulgeSize = 1, Density = 1e+05, byDensity = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Name}{
character, name of the DNA or RNA sequence to analyse.
}
  \item{Sequence}{
character, DNA or RNA sequence to analyse composed of the nucleotide arrangement.
}
  \item{DNA}{
logical, controls if the sequence is DNA or RNA. The factory-fresh default is \code{TRUE} assuming the sequence is DNA.
}
  \item{Complementary}{
logical, controls if the Complementary strand should be created and analyzed. The factory-fresh default is \code{TRUE}.
}
  \item{Nucleotides}{
character or vector of characters, nucleotide that composes the runs to analyse. The factory-fresh default is \code{c("G", "C")} to analyse both G and C-runs. Any other nucleotide (or letter) can be imputed but will be ignored.
}
  \item{BulgeSize}{
integer, number of acceptable non-run nucleotides to exist within runs. The factory-fresh default is \code{1}.
}
  \item{Density}{
integer, constant to calculate density results. The factory-fresh default is \code{100000}, returning results found per 100000 nucleotides. Only pertinent if \code{byDensity = TRUE}.
}
  \item{byDensity}{
logical, should the results be returned as a density. Calculated as Result.Density = (Density*Results)/(total genomic length). If set to \code{FALSE}, it will return counts. Density allows genomic size-indepedent comparisons. The factory-fresh default is \code{TRUE}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The result of GiG.Seq.Analysis is a one row data.frame with the summary of the genomic sequence.
}
\section{Column Meanings}{
\emph{Name}{: name of summary, given by the Name input.}

\emph{DNA}{: genome type, given by DNA input. If \code{DNA == TRUE}, the genome is DNA, else it is assumed to be RNA.}

\emph{Length}{: Length of the genomic sequence.}

\emph{Complementary}{: genome type, given by Complementary input. If \code{Complementary == TRUE}, the genome is double stranded, and both strands have been analysed.}

\emph{G\%seq}{: percentage of the genomic sequence which is G.}

\emph{C\%seq}{: percentage of the genomic sequence which is C.}

\emph{A\%seq}{: percentage of the genomic sequence which is A.}

\emph{UT\%seq}{: percentage of the genomic sequence which is U or T.}

\emph{N\%seq}{: percentage of the genomic sequence which is N.}

\emph{G2}{: Number of **perfect** G-runs identified with lengths betwen 2 and 5 using method 1 of G4-iM Grinder's algorithm, (GG, GGG, GGGG and GGGGG). Returned as counts or density. }

\emph{G3}{: Number of **perfect** G-runs identified with lengths betwen 3 and 5 using method 1 of G4-iM Grinder's algorithm, (GGG, GGGG and GGGGG). Returned as counts or density. }

\emph{G2X}{: Number of **perfect and imperfect** G-runs identified with lengths betwen 2 and 5 (excluding the bulges) using method 1 of G4-iM Grinder's algorithm. Returned as counts or density.}

\emph{G3X}{: Number of **perfect and imperfect** G-runs identified with lengths betwen 3 and 5 (excluding the bulges) using method 1 of G4-iM Grinder's algorithm. Returned as counts or density.}

\emph{C2}{: Number of **perfect** C-runs identified with lengths betwen 2 and 5 using method 1 of G4-iM Grinder's algorithm. (CC, CCC, CCCC and CCCCC). Returned as counts or density. }

\emph{C3}{: Number of **perfect** C-runs identified with lengths betwen 3 and 5 using method 1 of G4-iM Grinder's algorithm. (CCC, CCCC and CCCCC). Returned as counts or density. }

\emph{C2X}{: Number of **perfect and imperfect** C-runs identified with lengths betwen 2 and 5 (excluding the bulges) using method 1 of G4-iM Grinder's algorithm. Returned as counts or density.}

\emph{C3X}{: Number of **perfect and imperfect** C-runs identified with lengths betwen 3 and 5 (excluding the bulges) using method 1 of G4-iM Grinder's algorithm. Returned as counts or density.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Efres Belmonte-Reche
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
### Creating random nucleotidic sequence
Seq <- paste0(sample(c("G", "C", "T", "A", "N"), 10000, prob = c(1,1,0.6,0.6,0.01), replace = T), collapse = "")

### Analysing sequence
Rs <- GiG.Seq.Analysis(Name = "RandomSeq", Sequence = Seq, DNA = TRUE, Complementary = TRUE, Nucleotides = c("G", "C"), BulgeSize = 1, Density = 1e+05, byDensity = TRUE)

### Adding another analysis
Seq <- paste0(sample(c("G", "C", "T", "A", "N"), 10000, prob = c(1,1,0.6,0.6,0.01), replace = T), collapse = "")
Rs[2,] <- GiG.Seq.Analysis(Name = "RandomSeq2", Sequence = Seq, DNA = TRUE, Complementary = TRUE, Nucleotides = c("G", "C"), BulgeSize = 1, Density = 1e+05, byDensity = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

