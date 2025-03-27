\name{G4iMGrinder}
\alias{G4iMGrinder}
\title{
Detect and analyze potential G-quadruplexes, i-Motifs, and higher-order structures in DNA or RNA sequences
}
\description{
G4iM Grinder is a flexible search engine and characterization tool designed to detect and analyze sequences capable of forming G-quadruplexes (Potential Quadruplex Sequences, PQSs), i-Motifs (Potential i-Motif Sequences, PiMS), or other higher-order quadruplex-like structures in DNA or RNA. It provides multiple “methods” to search for these motifs, allowing extensive configurability. Users can tailor the search to match specific criteria and search within their results for known quadruplex-forming or non-quadruplex-forming sequences. The results include both raw findings and frequency-weighted summaries.
}
\references{
Belmonte-Reche, E. and Morales, J. C. (2019). G4-iM Grinder: when size and frequency matter. G-Quadruplex, i-Motif and higher order structure search and analysis tool. \emph{NAR Genomics and Bioinformatics}, 2. DOI: 10.1093/nargab/lqz005

\url{https://academic.oup.com/nargab/article/2/1/lqz005/5576141}
}
\author{
Efres Belmonte-Reche
}
\note{
M1 stands for Method 1; M2 stands for Method 2; M3 stands for Method 3.
}
\usage{
G4iMGrinder(
  Name,
  Sequence,
  DNA = TRUE,
  Complementary = TRUE,
  RunComposition = "G",
  BulgeSize = 1,
  MaxIL = 3,
  MaxRunSize = 5,
  MinRunSize = 3,
  MinNRuns = 4,
  MaxNRuns = 0,
  MaxPQSSize = 33,
  MinPQSSize = 15,
  MaxLoopSize = 10,
  MinLoopSize = 0,
  LoopSeq = c("G", "T", "A", "C"),
  Method2 = TRUE,
  Method3 = TRUE,
  G4hunter = TRUE,
  cGcC = FALSE,
  PQSfinder = TRUE,
  Bt = 14,
  Pb = 17,
  Fm = 3,
  Em = 1,
  Ts = 4,
  Et = -19,
  Is = -16,
  Ei = 1,
  Ls = 1,
  ET = 1,
  WeightParameters = c(0.5, 0.5, 0),
  FreqWeight = 0,
  KnownQuadruplex = TRUE,
  KnownNOTQuadruplex = FALSE,
  RunFormula = FALSE,
  NCores = 1,
  Verborrea = TRUE
)
}
\arguments{
  \item{Name}{
    \code{character}. Name of the DNA or RNA sequence under analysis.
  }
  \item{Sequence}{
    \code{character}. The nucleotide sequence to be examined. Must be composed of valid DNA or RNA bases.
  }
  \item{DNA}{
    \code{logical}. Indicates whether \code{Sequence} is DNA (\code{TRUE}) or RNA (\code{FALSE}). Defaults to \code{TRUE}.
  }
  \item{Complementary}{
    \code{logical}. If \code{TRUE}, the complementary strand is generated and analyzed in parallel. Defaults to \code{TRUE}.
  }
  \item{RunComposition}{
    \code{character}. Nucleotide(s) used to define the “runs.” Typically \code{"G"} for G-quadruplex search or \code{"C"} for i-Motif search. Defaults to \code{"G"}.
  }
  \item{BulgeSize}{
    \code{integer}. Number of allowed non-\code{RunComposition} nucleotides within a run (used by M1). Defaults to \code{1}.
  }
  \item{MaxIL}{
    \code{integer}. Total number of additional nucleotides allowed between runs (used by M2). Defaults to \code{3}.
  }
  \item{MaxRunSize}{
    \code{integer}. Maximum length of a run. Defaults to \code{5} (used by M2).
  }
  \item{MinRunSize}{
    \code{integer}. Minimum length of a run. Defaults to \code{3} (used by M1).
  }
  \item{MinNRuns}{
    \code{integer}. Minimum number of runs required to form a structure. Defaults to \code{4} (used by M2 and M3).
  }
  \item{MaxNRuns}{
    \code{integer}. Maximum number of runs that compose a structure. Defaults to \code{0}, which disables an upper limit for run count (used by M2).
  }
  \item{MaxPQSSize}{
    \code{integer}. Maximum total length of a putative quadruplex structure. Defaults to \code{33} (used by M2).
  }
  \item{MinPQSSize}{
    \code{integer}. Minimum total length of a putative quadruplex structure. Defaults to \code{15} (used by M2 and M3).
  }
  \item{MaxLoopSize}{
    \code{integer}. Maximum number of nucleotides allowed in each loop (used by M2 and M3). Defaults to \code{10}.
  }
  \item{MinLoopSize}{
    \code{integer}. Minimum number of nucleotides allowed in each loop (used by M2 and M3). Defaults to \code{0}.
  }
  \item{LoopSeq}{
    \code{character} vector. Defines the nucleotide(s) or pattern(s) to measure or highlight within detected structures. Defaults to \code{c("G", "T", "A", "C")}.
  }
  \item{Method2}{
    \code{logical}. If \code{TRUE}, enables Method 2 (M2), which searches for size-defined structures and computes frequency (M2A and M2B). Defaults to \code{TRUE}.
  }
  \item{Method3}{
    \code{logical}. If \code{TRUE}, enables Method 3 (M3), which searches for size-unrestricted structures and computes frequency (M3A and M3B). Defaults to \code{TRUE}.
  }
  \item{G4hunter}{
    \code{logical}. If \code{TRUE}, applies the G4Hunter scoring system. Defaults to \code{TRUE}.
  }
  \item{cGcC}{
    \code{logical}. If \code{TRUE}, applies the cGcC scoring system (valid for RNA). Defaults to \code{FALSE}.
  }
  \item{PQSfinder}{
    \code{logical}. If \code{TRUE}, applies an adaptation of the PQSfinder scoring system. Defaults to \code{TRUE}.
  }
  \item{Bt}{
    \code{integer}. Tetrad stacking bonus for PQSfinder calculations. Defaults to \code{14}.
  }
  \item{Pb}{
    \code{integer}. Inter-Loop penalization constant for PQSfinder calculations. Defaults to \code{17}.
  }
  \item{Fm}{
    \code{integer}. Loop length penalization constant for PQSfinder calculations. Defaults to \code{3}.
  }
  \item{Em}{
    \code{integer}. Loop length exponential constant for PQSfinder calculations. Defaults to \code{1}.
  }
  \item{Ts}{
    \code{integer}. Tetrad supplement constant for PQSfinder calculations. Defaults to \code{4}.
  }
  \item{Et}{
    \code{integer}. Inter-Loop supplement constant for PQSfinder calculations. Defaults to \code{-19}.
  }
  \item{Is}{
    \code{integer}. Loop supplement constant for PQSfinder calculations. Defaults to \code{-16}.
  }
  \item{Ei}{
    \code{integer}. Tetrad exponential constant for PQSfinder calculations. Defaults to \code{1}.
  }
  \item{Ls}{
    \code{integer}. Inter-Loop exponential constant for PQSfinder calculations. Defaults to \code{1}.
  }
  \item{ET}{
    \code{integer}. Total formula exponential constant for PQSfinder calculations. Defaults to \code{1}.
  }
  \item{WeightParameters}{
    \code{numeric} vector of length 3. Weights for combining \code{G4hunter}, \code{PQSfinder}, and \code{cGcC} scores (in that order). Defaults to \code{c(0.5, 0.5, 0)}, producing an average of the first two.
  }
  \item{FreqWeight}{
    \code{numeric}. Weight factor for incorporating structure frequency in the final score (relevant to M2B and M3B). Defaults to \code{0}.
  }
  \item{KnownQuadruplex}{
    \code{logical}. If \code{TRUE}, matches results against known sequences that have been shown to form G-quadruplex or i-Motif in vitro. Defaults to \code{TRUE}.
  }
  \item{KnownNOTQuadruplex}{
    \code{logical}. If \code{TRUE}, matches results against known sequences shown \emph{not} to form quadruplexes. Defaults to \code{FALSE}.
  }
  \item{RunFormula}{
    \code{logical}. If \code{TRUE}, calculates and reports a symbolic formula for each detected PQS. Defaults to \code{FALSE}.
  }
  \item{NCores}{
    \code{integer}. Number of cores to use for parallel processing. Defaults to \code{1}.
  }
  \item{Verborrea}{
    \code{logical}. If \code{TRUE}, prints verbose messages about progress. Defaults to \code{TRUE}.
  }
}
\details{
}
\value{
A \code{list} containing:
\item{Configuration}{A \code{data.frame} of the parameters used in the run.}
\item{FunTime}{A \code{data.frame} with timing information for each step.}
\item{PQSM2a}{A \code{data.frame} of the M2A (size-defined) results, if \code{Method2 = TRUE}.}
\item{PQSM2b}{A \code{data.frame} of the M2B (frequency-weighted) results, if \code{Method2 = TRUE}.}
\item{PQSM3a}{A \code{data.frame} of the M3A (unrestricted-size) results, if \code{Method3 = TRUE}.}
\item{PQSM3b}{A \code{data.frame} of the M3B (frequency-weighted) results, if \code{Method3 = TRUE}.}
}
\section{Column Meanings}{
  \describe{
    \item{\emph{Start}}{Integer. Start position in \code{Sequence} (for M2A/M3A).}
    \item{\emph{Finish}}{Integer. End position in \code{Sequence} (for M2A/M3A).}
    \item{\emph{Freq}}{Integer. Frequency of occurrence (for M2B/M3B).}
    \item{\emph{Runs}}{Integer. Number of runs (e.g., G-runs in G-quadruplex).}
    \item{\emph{IL}}{Integer. Number of bulges or irregularities.}
    \item{\emph{mRun}}{Numeric. Average run size.}
    \item{\emph{Sequence}}{Character. The identified motif sequence.}
    \item{\emph{Length}}{Integer. Total length of the identified structure.}
    \item{\emph{Strand}}{Character. Indicates “+” (original) or “–” (complementary) strand, if \code{Complementary = TRUE}.}
    \item{\emph{G4Hunter}}{Numeric. Score assigned by the G4Hunter algorithm (if \code{G4hunter = TRUE}).}
    \item{\emph{pqsfinder}}{Numeric. Score from the PQSfinder adaptation (if \code{PQSfinder = TRUE}).}
    \item{\emph{cGcC}}{Numeric. Score from the cGcC algorithm (if \code{cGcC = TRUE}).}
    \item{\emph{Score}}{Numeric. Combined overall score, integrating all selected scoring methods plus frequency weighting.}
    \item{\emph{Conf.Quad.Seqs}}{Character. Known quadruplex-forming sequences detected, with counts. DNA hits have “*” after the count; RNA hits have “^”.}
    \item{\emph{Conf.NOT.Quad.Seqs}}{Character. Known non-quadruplex sequences detected, with counts. DNA hits have “*” after the count; RNA hits have “^”.}
  }
}
\examples{
  library(G4iMGrinder)

  # Example: retrieve a DNA sequence and run basic G4 search
  if (!require("seqinr")) {
    install.packages("seqinr")
    library(seqinr)
  }

  Name <- "LmajorESTs"
  Sequence <- paste0(
    read.fasta(
      file = url("http://tritrypdb.org/common/downloads/release-36/Lmajor/fasta/TriTrypDB-36_Lmajor_ESTs.fasta"),
      as.string = TRUE, legacy.mode = TRUE, seqonly = TRUE,
      strip.desc = TRUE, seqtype = "DNA"
    ),
    collapse = ""
  )

  # G-quadruplex search on DNA
  resultDNA <- G4iMGrinder(Name = Name, Sequence = Sequence)

  # G-quadruplex search on RNA (with cGcC scoring)
  resultRNA <- G4iMGrinder(Name = Name, Sequence = Sequence, DNA = FALSE, cGcC = TRUE)

  # i-Motif search in DNA
  resultIMotif <- G4iMGrinder(Name = Name, Sequence = Sequence, RunComposition = "C")

  # Customized search with bulge allowance and larger loop sizes ## More bulges and smaller G-runs (GG) increases significantly computation time
  resultCustom <- G4iMGrinder(
    Name = Name,
    Sequence = Sequence,
    BulgeSize = 2,
    MaxLoopSize = 20,
    MaxIL = 10
  )

  # Viewing results
  View(resultDNA$PQSM2a)  # M2A results
  View(resultDNA$PQSM2b)  # M2B results (with frequency weighting)
  View(resultDNA$PQSM3a)  # M3A results (unrestricted-size search)
  View(resultDNA$PQSM3b)  # M3B results (unrestricted-size with frequency weighting)
}
\keyword{datagen}
\keyword{utilities}
