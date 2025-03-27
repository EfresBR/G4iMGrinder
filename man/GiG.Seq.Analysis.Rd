\name{GiG.Seq.Analysis}
\alias{GiG.Seq.Analysis}
\title{
Genomic sequence nucleotide run analyzer
}
\description{
\code{GiG.Seq.Analysis} examines a DNA or RNA sequence (and optionally its complementary strand) to identify and count runs of Guanines and Cytosines (G and C). It detects both perfect runs (e.g., GGGG) and imperfect runs (with bulges) of sizes ranging from 2 (e.g., GG; G-run of size 2) to 4 (e.g.,GGGG), in an non-overlapping way. The function can return raw counts or densities (per a specified length, e.g., per 100,000 nucleotides) for direct comparisons across sequences of different lengths. G and C runs are searched independently.
}
\usage{
GiG.Seq.Analysis(
  Name,
  Sequence,
  DNA = TRUE,
  Complementary = TRUE,
  Density = 1e+05,
  byDensity = TRUE
)
}
\arguments{
  \item{Name}{
    \code{character}. Name or identifier of the DNA/RNA sequence to analyze.
  }
  \item{Sequence}{
    \code{character}. The nucleotide sequence to analyze (A, C, G, T, U, or N).
  }
  \item{DNA}{
    \code{logical}. Indicates whether \code{Sequence} is DNA (\code{TRUE}) or RNA (\code{FALSE}). Defaults to \code{TRUE}.
  }
  \item{Complementary}{
    \code{logical}. If \code{TRUE}, the complementary strand is created and analyzed alongside the original strand. Defaults to \code{TRUE}.
  }
  \item{Density}{
    \code{integer}. The scaling factor for returning densities rather than raw counts. Defaults to \code{100000}, so results are reported per 100,000 nucleotides if \code{byDensity = TRUE}.
  }
  \item{byDensity}{
    \code{logical}. If \code{TRUE}, results are returned as densities, i.e., \code{(Density * runCounts) / sequenceLength}. If \code{FALSE}, raw run counts are returned. Defaults to \code{TRUE}.
  }
}
\details{
By default, this function specifically looks for runs of “G” and “C,” counting perfect runs (e.g., G2, G3, G4) and imperfect runs. Runs are analyzed in a sequential, non-overlapping manner. For example, a run of “GGGG” is counted only as G4, not as G4 plus any subset runs like G2 or G3.
}
\value{
A one-row \code{data.frame} summarizing the run analysis:
}
\section{Column Meanings}{
  \describe{
    \item{\emph{Name}}{Identifier for the analyzed sequence, matching the \code{Name} argument.}
    \item{\emph{DNA}}{Logical: \code{TRUE} if the sequence was treated as DNA; \code{FALSE} if RNA.}
    \item{\emph{Length}}{Total length of \code{Sequence}.}
    \item{\emph{Complementary}}{\code{TRUE} if the complementary strand was analyzed as well, otherwise \code{FALSE}.}
    \item{\emph{G\%seq}, \emph{C\%seq}, \emph{A\%seq}, \emph{UT\%seq}, \emph{N\%seq}}{Percentages of each nucleotide type within \code{Sequence}. (U and T are combined under UT\%seq if \code{DNA=FALSE}.)}
    \item{\emph{G2}, \emph{G3}, \emph{G2X}, \emph{G3X}}{Counts or densities of perfect (\emph{G2}, \emph{G3}) and perfect+imperfect (\emph{G2X}, \emph{G3X}) G-runs of length 2 to 4, according to the method 1 approach in G4-iM Grinder.}
    \item{\emph{C2}, \emph{C3}, \emph{C2X}, \emph{C3X}}{Counts or densities of perfect (\emph{C2}, \emph{C3}) and perfect+imperfect (\emph{C2X}, \emph{C3X}) C-runs of length 2 to 4, according to the same method 1 approach.}
  }
}
\references{
Belmonte-Reche, E. and Morales, J. C. (2019). G4-iM Grinder: when size and frequency matter. G-Quadruplex, i-Motif and higher order structure search and analysis tool. \emph{NAR Genomics and Bioinformatics}, 2. DOI: 10.1093/nargab/lqz005

\url{https://academic.oup.com/nargab/article/2/1/lqz005/5576141}
}
\author{
Efres Belmonte-Reche
}
\note{
Any additional notes or implementation details can be placed here.
}
\seealso{
\code{\link{G4iMGrinder}} for broader G4/i-Motif detection and scoring.
}
\examples{
# Creating a random nucleotide sequence of length 10,000
Seq <- paste0(
  sample(
    c("G", "C", "T", "A", "N"),
    10000,
    prob = c(1, 1, 0.6, 0.6, 0.01),
    replace = TRUE
  ),
  collapse = ""
)

# Running the analysis with default parameters
Rs <- GiG.Seq.Analysis(
  Name = "RandomSeq",
  Sequence = Seq,
  DNA = TRUE,
  Complementary = TRUE,
  byDensity = TRUE
)

# Analyzing a second sequence and storing results in the same data frame
Seq2 <- paste0(
  sample(
    c("G", "C", "T", "A", "N"),
    10000,
    prob = c(1, 1, 0.6, 0.6, 0.01),
    replace = TRUE
  ),
  collapse = ""
)
Rs[2, ] <- GiG.Seq.Analysis(
  Name = "RandomSeq2",
  Sequence = Seq2,
  DNA = TRUE,
  Complementary = TRUE,
  Density = 1e5,
  byDensity = TRUE
)

print(Rs)
}

