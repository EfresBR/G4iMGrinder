\name{GiGList.Analysis}
\alias{GiGList.Analysis}
\title{
Summary generator for a G4-iM Grinder result list (GiGList)
}
\description{
\code{GiGList.Analysis} produces a one-row summary of the results from a \code{G4iMGrinder} output list (referred to as a “GiGList”). It aggregates counts/densities of putative quadruplex or i-Motif structures, applying optional filters on minimum score, frequency, and structure length. The resulting summary facilitates high-level comparisons across different genomes or experimental conditions.
}
\usage{
GiGList.Analysis(
  GiGList,
  iden,
  ScoreMin = c(20, 40),
  FreqMin = 10,
  LengthMin = 50,
  Density = 100000,
  byDensity = TRUE
)
}
\arguments{
  \item{GiGList}{
    \code{list}. A \code{G4iMGrinder} result list, typically named \emph{GiGList}, to be analyzed or summarized.
  }
  \item{iden}{
    \code{character}. Identification tag or label for the DNA/RNA source of the \code{GiGList}.
  }
  \item{ScoreMin}{
    \code{integer} or vector of \code{integer} values. One or more minimum score thresholds to filter structures.
    - For G-quadruplex searches, any structure with a score \code{>= ScoreMin[i]} passes the filter.
    - For i-Motif searches, the score sign is reversed (multiplied by -1) to be comparable with G-quadruplex thresholds.
    A separate column is created for each value in \code{ScoreMin}.
    Defaults to \code{c(20, 40)}, representing “medium” (≥ 20) and “high” (≥ 40) formation probability for G-quadruplexes (or ≤ -20 and ≤ -40 for i-Motifs).
  }
  \item{FreqMin}{
    \code{integer} or vector of \code{integer} values. One or more minimum frequency thresholds to filter structures based on how often they appear in the genome.
    A separate column is created for each threshold in \code{FreqMin}.
    Defaults to \code{10}, i.e., returning results with \code{frequency >= 10}.
  }
  \item{LengthMin}{
    \code{integer} or vector of \code{integer} values. One or more minimum length thresholds, applied only to results from Method 3 (M3a/M3b).
    A separate column is created for each threshold in \code{LengthMin}.
    Defaults to \code{50}, i.e., returning \code{length >= 50}.
  }
  \item{Density}{
    \code{integer}. Scaling factor for returning densities rather than raw counts. Defaults to \code{100000} (results per 100,000 nucleotides if \code{byDensity = TRUE}).
  }
  \item{byDensity}{
    \code{logical}. If \code{TRUE}, the function reports densities, computed as \code{(Density * count) / genomeLength}. If \code{FALSE}, raw counts are returned instead. Defaults to \code{TRUE}.
  }
}
\details{
Filters on score, frequency, and length can be applied simultaneously, producing separate result columns for each filter. Known quadruplex-forming or known NOT-to-form-quadruplex sequences (if included in the \code{GiGList}) are also reported. Additionally, columns indicating the percentage of “unique” structures (frequency = 1) are computed as a measure of overall structural redundancy.
}
\references{
Belmonte-Reche, E. and Morales, J. C. (2019). G4-iM Grinder: when size and frequency matter. G-Quadruplex, i-Motif and higher order structure search and analysis tool. \emph{NAR Genomics and Bioinformatics}, 2. DOI: 10.1093/nargab/lqz005

\url{https://academic.oup.com/nargab/article/2/1/lqz005/5576141}
}
\value{
A one-row \code{data.frame} summarizing the \code{GiGList}, with potential columns such as:
\item{Name}{Character, derived from the name assigned to \code{GiGList}.}
\item{iden}{Character, the \code{iden} input labeling the sequence source.}
\item{Length}{Integer, total genomic length from the \code{GiGList}.}
\item{SeqG, SeqC}{Numeric, percentage of G or C in the analyzed sequence(s).}
\item{nM2a, nM2b, nM3a, nM3b}{Count or density of structures found by Methods 2a, 2b, 3a, or 3b respectively (if those methods were used).}
\item{...S\{|x|\}, ...F\{|y|\}, ...L\{|z|\}}{Suffixes denoting columns filtered by score (\code{x}), frequency (\code{y}), or length (\code{z}).}
\item{...KTFQ, ...KNTFQ}{Suffixes denoting structures that match known-to-form or known-NOT-to-form quadruplex sequences.}
\item{...UniqPercent}{Columns indicating the percentage of unique structures within each method.}
\item{Config}{A character column summarizing the configuration of analysis or filter settings.}
}
\section{Column Meanings}{
  \describe{
    \item{\emph{Name}}{Name of the G4iMGrinder result list, inherited from \code{GiGList}.}
    \item{\emph{iden}}{Identification tag, matching \code{iden} argument.}
    \item{\emph{Length}}{Total length of the sequence(s) as computed by \code{G4iMGrinder}.}
    \item{\emph{SeqG}, \emph{SeqC}}{Percentages of G and C nucleotides in the genome.}
    \item{\emph{nM2a}, \emph{nM2b}, \emph{nM3a}, \emph{nM3b}}{Counts or densities of structures detected by each method (2a, 2b, 3a, 3b).}
    \item{\emph{.S|X|}, \emph{.F|X|}, \emph{.L|X|}}{Suffixes appended for each filter: “S” for score, “F” for frequency, “L” for length. \emph{X} indicates the numeric threshold (e.g., 20, 10, or 50).}
    \item{\emph{.KTFQ}, \emph{.KNTFQ}}{Columns indicating, for each method, how many structures are known-to-form quadruplexes or known-NOT-to-form quadruplexes (counts or densities).}
    \item{\emph{.UniqPercent}}{Percentage of unique results (\code{frequency=1}) within each method.}
    \item{\emph{Config}}{A text string describing the filters and configurations used in the summary.}
  }
}
\author{
Efres Belmonte-Reche
}
\section{Warning}{
Column names depend on the analysis conditions (methods used, filter thresholds, etc.). Combining or binding multiple \code{GiGList.Analysis} summaries requires consistent parameters. Results may not align if different G4iMGrinder settings are used.
}
\seealso{
\code{\link{G4iMGrinder}} for generating the \code{GiGList}, and other related utilities.
}
\examples{
# Example usage:

# Suppose 'Rs' is a G4iMGrinder result list (GiGList) for a DNA G-quadruplex search
Rs <- G4iMGrinder(Name = "TestGenome", Sequence = "ACGT...")

# Summarize with default filters
Summary <- GiGList.Analysis(
  GiGList = Rs,
  iden = "Parasite"
)

print(Summary)
}
