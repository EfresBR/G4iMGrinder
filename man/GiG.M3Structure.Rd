\name{GiG.M3Structure}
\alias{GiG.M3Structure}
\title{
Potential Higher-Order Quadruplex Sequence (PHOQS) arrangement analyzer for Method 3A (M3A) candidates
}
\description{
\code{GiG.M3Structure} examines a specified PHOQS (Potential Higher-Order Quadruplex Sequence) candidate identified by G4-iM Grinder's Method 3A (M3A). It searches for possible sub-unit conformations (using M2A results) that fit within the PHOQS region, attempting to identify the most probable arrangement of quadruplex sub-units. Three analytical approaches are available:

\itemize{
  \item \strong{HRA} (Highest Ranking Assembly): sequentially seats the sub-units with the strongest scores, preferring known-to-form sequences.
  \item \strong{RAH} (Ranking Assembly by Highest Score): explores all possible sub-unit seatings, returning the conformation with the highest overall mean score.
  \item \strong{RAnH} (Ranking Assembly by Normalized Highest Score): explores all possible sub-unit seatings, returning the conformation with the highest “normalized” mean score (the mean score multiplied by the fraction of the PHOQS region actually occupied by PQS sub-units).
}
}
\usage{
GiG.M3Structure(
  GiGList,
  M3ACandidate,
  MAXite
)
}
\arguments{
  \item{GiGList}{
    \code{list}. A \code{G4iMGrinder} result list (GiGList) containing M2A and M3A data frames. Typically produced by \code{G4iMGrinder}.
  }
  \item{M3ACandidate}{
    \code{integer}. The row number in \code{GiGList\$PQSM3a} corresponding to the PHOQS candidate of interest.
  }
  \item{MAXite}{
    \code{integer}. Number of iterations for the random seat-allocation algorithm. A higher number ensures more exhaustive exploration of potential sub-unit conformations, especially for large PHOQS or those with many candidate sub-units. For instance, 10,000 iterations may be sufficient for a ~100-nt region with ~30 potential sub-units.
  }
}
\details{
By combining the M2A (Method 2A) sub-units with a candidate from M3A (Method 3A), \code{GiG.M3Structure} locates all conformations that could explain how multiple quadruplex sub-units (PQS) might form a single higher-order structure. The iterative random approach attempts various ways of “seating” sub-units within the PHOQS, then records the best conformations under each scoring approach (HRA, RAH, and RAnH).

Ensure \code{MAXite} is large enough to sample all or most conformations in the PHOQS region; otherwise, some arrangements may be missed.
}
\references{
Belmonte-Reche, E. and Morales, J.C. (2019). G4-iM Grinder: when size and frequency matter. G-Quadruplex, i-Motif and higher order structure search and analysis tool. \emph{NAR Genomics and Bioinformatics}, 2. \doi{10.1093/nargab/lqz005}.

\url{https://academic.oup.com/nargab/article/2/1/lqz005/5576141}
}
\author{
Efres Belmonte-Reche
}
\value{
A \code{list} describing the analysis of the selected PHOQS candidate:
\item{M2}{A \code{data.frame} of the potential sub-units (\code{M2A} results) located within the PHOQS region. An error is returned if no sub-units match the target region.}
\item{M3}{A \code{data.frame} containing one row from \code{GiGList\$PQSM3a}, representing the selected PHOQS candidate.}
\item{Potential.Arrangements}{A \code{data.frame} enumerating all discovered sub-unit conformations (arrangements) compatible with the PHOQS, found via random seat allocation repeated \code{MAXite} times.}
\item{Best.Arrangements}{A \code{list} containing up to three entries—\code{HRA}, \code{RAH}, and \code{RAnH}—each storing the best conformation(s) discovered for the given PHOQS, depending on the selected scoring method.}
}
\section{Column Meanings}{
  \describe{
    \item{\emph{RES}}{Reference for sub-units (by row name) from \code{GiGList\$PQSM2A} used in the arrangement.}
    \item{\emph{nPQS}}{Integer. Number of sub-units forming a particular conformation.}
    \item{\emph{MeanScore}}{Numeric. Mean of the sub-unit scores for that conformation.}
    \item{\emph{PQSLenghPercent}}{Numeric. Percentage of the PHOQS region occupied by the sub-units.}
    \item{\emph{idenPQS}}{Reference for sub-units in the \code{M2} data frame that form the conformation.}
    \item{\emph{nMS}}{Numeric. The “normalized mean score” = \code{MeanScore} * (\code{PQSLenghPercent}/100). Higher values indicate a high average sub-unit score and a larger fraction of the PHOQS occupied by them.}
  }
}
\examples{
# Example usage:

# 1. Running G4iMGrinder on a DNA sequence
Rs <- G4iMGrinder(Name = "ExampleSeq", Sequence = "ACGT...")

# 2. Analyzing the first PHOQS in the M3A results (row 1 of Rs$PQSM3a),
#    performing 10,000 seat-allocation iterations
firstPHOQS <- GiG.M3Structure(
  GiGList = Rs,
  M3ACandidate = 1,
  MAXite = 10000
)

# Check the best arrangements found under HRA, RAH, and RAnH
print(firstPHOQS$Best.Arrangements)
}
