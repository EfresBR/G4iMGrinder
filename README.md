# **G4-iM Grinder**

**G4-iM Grinder** is a fast, robust, and highly adaptable algorithm capable of locating, identifying, qualifying, and quantifying DNA and RNA potential quadruplex structures, such as G-quadruplexes, i-Motifs, and their higher-order variants.

Read the open-access paper on the algorithm:  
[G4-iM Grinder: when size and frequency matter. G-Quadruplex, i-Motif and higher order structure search and analysis tool](https://bit.ly/3j2UUjC)

Read about the application of the algorithm to SARS-CoV-2 and the entire Virus Realm in another open-access paper:  
[Potential G-quadruplexes and i-Motifs in the SARS-CoV-2](https://bit.ly/GiGOriginalData)

The results from both manuscripts can be found in the “Results” section below.

<div style="text-align:center">
<img src="https://www.biorxiv.org/content/biorxiv/early/2019/06/17/532382/F1.large.jpg?width=800&height=600&carousel=1" height="500" />
</div>

Please report bugs, problems, or feature requests in the Issues section.

---

## **Quadruplexes and G4-iM Grinder**

<details>

<summary>What are quadruplexes?</summary>

**[G-quadruplexes (G4s)](https://en.wikipedia.org/wiki/G-quadruplex):**  
G4s are DNA or RNA sequences rich in guanine, where four guanine bases can associate through Hoogsteen hydrogen bonding to form a square planar structure called a guanine tetrad (G-tetrad or G-quartet). Two or more G-tetrads can then stack on top of each other to form a stable G4. Unimolecular G4s occur naturally in telomeric regions and various transcriptional regulatory regions.

**[C-quadruplexes or i-Motifs (iM)](https://en.wikipedia.org/wiki/I-motif_DNA):**  
iMs are quadruplex structures formed by cytosine-rich DNA or RNA, analogous to G4s formed by guanine-rich sequences. C-rich DNA regions frequently appear in gene-regulatory portions of the genome. iMs have been experimentally observed in human cells and may play roles in cell reproduction. They also have potential applications in nanotechnology due to their pH sensitivity, serving as biosensors, nanomachines, and molecular switches.

</details>

<details>

<summary>Searching for quadruplexes?</summary>

Quadruplexes have drawn significant attention in recent years due to evidence of their functional roles across many living organisms, yet their precise formation mechanisms remain under investigation. To pinpoint potential structures, \emph{in silico} predictions rely on known \emph{in vitro} paradigms. Loops, tetrad count, run imperfections, and flanking genomic regions all appear to influence quadruplex topology and dynamics.

**G4-iM Grinder (GiG)** provides:  
1. A **search engine** that locates all possible candidates matching user-defined criteria (G-runs, C-runs, loops, etc.).  
2. A **qualification engine** that ranks or filters results by their probability of forming an actual quadruplex or i-Motif, their frequency in the genome, their overlap with known quadruplex sequences, and more.

</details>

### **Important Terms & Abbreviations**

1. **Quadruplex**: A genomic sequence that forms a G4 or iM.  
2. **G4**: A sequence confirmed to form a G4 \emph{in vitro}.  
3. **iM**: A sequence confirmed to form an iM \emph{in vitro}.  
4. **PQS**: A \emph{Potential} G4 Sequence detected \emph{in silico} but not yet confirmed \emph{in vitro}.  
5. **PiMS**: A \emph{Potential} iM Sequence detected \emph{in silico} but not yet confirmed \emph{in vitro}.

---

## **G4-iM Grinder’s Updates & News**

**Latest version: 1.6.5 (03-2025).**  
- Rewrote `GiGList.Analysis` to fix a bug counting runs with bulges; Changed it so it returns non-overlapping results.  
- Updated documentation.
- Updated DDBB to version 2.6.

<details>
<summary>Change Log</summary>

**For Version 1.6.4**  
- Adapted and further optimized `GiG.df.GenomicFeatures`.

**For Version 1.6.1**  
- Adapted and further optimized `GiG.df.GenomicFeatures`.  
- Changed `GiGList.Analysis` to accept vectors instead of single numerals in its parameters.  
- Improved results summaries in G4-iM Grinder.  
- Increased efficiency for DNA/RNA sequence handling.  
- Refined known quadruplex sequences detection to handle a growing database of confirmed G4/iM sequences.  
- Added Biostrings and biomartr dependencies.  
- G4-iM Grinder version and database version are now saved in the configuration data frame of each result.  
- Added a function to analyze genome runs (`GiG.Seq.Analysis`).  
- Added a function to analyze biological landmarks.  
- Streamlined package loading with checks for correct R version and installed dependencies.

**For Version 1.5.95**  
- Fixed a bug in the PQSfinder algorithm.  
- Changed how known G4s and iMs are detected:  
  - DNA hits get an asterisk (\*).  
  - RNA hits get a circumflex (^).  
  - Example: \`GUK1 (1*)\` or \`42.HIRA (WT) (1^)\`.
</details>

**Version 2.6 of G4-iM Grinder’s database (GiG.DB)**  
- Includes currently confirmed G4/iM sequences in SARS-CoV-2 plus over **3300** other confirmed or non-confirmed quadruplex sequences from the literature. Acknowledgments to Vasco Peixoto for introducing 400 of these sequences he found in literature.

<details>
<summary>Database Details</summary>

The **GiG.DB** includes:
1. **BioInformatic** dataframe: each entry is a nucleotide sequence from scientific studies, containing info on whether it forms a quadruplex, DNA or RNA type, etc.  
2. **Refs** dataframe: references (DOIs, PubMed IDs, etc.) for each BioInformatic entry.  
3. **BioPhysical** dataframe: T\textsubscript{m}, pH, and ion conditions for select sequences.

If you find errors or missing sequences, please open a GitHub issue at **EfresBR/G4iMGrinder**.

</details>

---

## **G4-iM Grinder’s Results**

### **SARS-CoV-2** Reference Genome Analysis

The reference genome (GCF\_009858895.2) was analyzed with a “lax” quadruplex-configuration. Results are offered in multiple formats:

- [\*.RDS](http://bit.ly/3drDhdM) (recommended for R users)  
- [\*.gff3](http://bit.ly/3ue7Edi)  
- [\*.xlsx](http://bit.ly/3s4SDc3)  

These files include all positions of PQS and PiMS found, their conservation, and any confirmed G4/iMs in SARS-CoV-2.

For **VARIANTS** found in other SARS-CoV-2 lineages and clades (not in the reference genome), see [\*.RDS](http://bit.ly/37rutk0) for GISAID-based variant data.

**VIRUS REALM** and reference genomes data were also analyzed under a “lax” quadruplex configuration.

<details>
<summary>Analytical & Raw Data</summary>

**ANALYTICAL DATA** ([Analysis.RData](http://bit.ly/3qxZH0v)):  
1. \`Analysis.Coronaviridae.fam\` – Summaries via \`GiGList.Analysis\` for the Coronaviridae family.  
2. \`Analysis.Virus.realm\` – Summaries via \`GiGList.Analysis\` for the entire virus realm.  
3. \`Baltimore.C\` – Classification tables.

**RAW DATA** ([Virus.Results.RDS](http://bit.ly/3sanTqc), ~2.4 GB):  
A large \`list\` grouping each virus family’s results (PQS, PiMS sub-lists). Methods 2A (PQSM2A) and 3 are included.

**[GISAID.refs.rar](http://bit.ly/3s5X4n4)** – References for the 17,312 SARS-CoV-2 genomes from the GISAID database.
</details>


### **G4-iM Grinder’s Results in Humans and Other Pathogens**

Results for humans and 49 pathogenic species using default G4-iM Grinder parameters (published in [the original article](https://bit.ly/3j2UUjC)) are available via [this link](https://1drv.ms/u/s!AvVGQg2rNIwDgTeth6qclA8Rz5UM?e=gmEI1a).

<details>
<summary>Database & Analysis Notes</summary>

- These analyses use **GiG.DB** V.2.5 (03-2020), which includes **2851** known-to-form / known-NOT-to-form quadruplexes.  
- ~312,072 M2A results contain at least one confirmed G4.  
- ~160,054 M2A results contain at least one confirmed iM.  

Four RData files store these results:

1. `Human.PQS.032020.RData`  
2. `Human.PiMS.032020.RData`  
3. `NonHuman.PQS.032020.RData`  
4. `NonHuman.PiMS.032020.RData`

Genomes:  
- Human genome: hg38, GRCh38.p12 (Sanger, May 2019).  
- Non-human genomes: see Section 9 of the supplementary material of the original article.

</details>

---

## **Installation and User Guide**

<details>
<summary>A. Package prerequisites</summary>

**G4-iM Grinder** is hosted at GitHub: `EfresBR/G4iMGrinder`. It requires R ≥ 4.0.0 and several CRAN/Bioconductor packages:

```r
pck <- c(
  "stringr", "stringi", "plyr", "seqinr", "stats", "parallel", 
  "doParallel", "beepr", "stats4", "devtools", "dplyr", 
  "BiocManager", "tibble"
)

foo <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}
foo(pck)
BiocManager::install(c("BiocGenerics", "S4Vectors", "Biostrings", "biomartr", "IRanges"), 
                     ask = FALSE, update = TRUE)
```

</details>

<details>
<summary>B. Package installation and loading</summary>

```r
devtools::install_github("EfresBR/G4iMGrinder")
library(G4iMGrinder)
```

</details>

<details>
<summary>C. Installation issues</summary>

Common pitfalls:

1. Missing dependencies  
2. R < 4.0.0  

Use the script below to verify:

```r
pck <- c("BiocGenerics", "S4Vectors", "stringr", "stringi", "plyr", 
         "seqinr", "stats", "parallel", "doParallel", "beepr", 
         "stats4", "devtools", "dplyr", "BiocManager", "biomartr", 
         "Biostrings")

FailFoo <- function(x){
  Info <- "Package dependencies FAILED: not installed -> "
  count <- 0
  for( i in x ){
    if( ! require( i , character.only = TRUE, quietly = TRUE ) ){
      Info <- paste0(Info, i, " ")
      count <- count +1
    }
  }
  if(count == 0){
    print("Package dependencies PASSED.")
  } else {
    print(Info)
  }
  AAA <- R.version
  if(as.numeric(AAA$major) == 4){
    if(as.numeric(AAA$minor) >= 0){
      print("R version PASSED (>= 4.0)")
    } else {
      print("R version FAILED. Update R to >= 4.0")
    }
  } else {
    print("R version FAILED. Update R to >= 4.0")
  }
}

FailFoo(pck)
```

Expected result:

```
[1] "Package dependencies PASSED."
[1] "R version PASSED (>= 4.0)"
```

If these tests pass but installation fails, please open an Issue with the full error trace.

</details>

<details>
<summary>D. (NEW) Running a genomic pre-analysis</summary>

Use `GiG.Seq.Analysis` to measure genome-wide run composition, returning a data frame of relevant features:

```r
loc <- url("http://tritrypdb.org/common/downloads/release-36/Lmajor/fasta/TriTrypDB-36_Lmajor_ESTs.fasta")
Sequence <- paste0(
  seqinr::read.fasta(file = loc, as.string = TRUE, legacy.mode = TRUE, 
                     seqonly = TRUE, strip.desc = TRUE), 
  collapse = ""
)

Pre_Rs <- GiG.Seq.Analysis(
  Name = "LmajorESTs",
  Sequence = Sequence,
  DNA = TRUE,
  Complementary = TRUE
)
```

</details>

<details>
<summary>E. Running a G4-iM Grinder analysis</summary>

```r
loc <- url("http://tritrypdb.org/common/downloads/release-36/Lmajor/fasta/TriTrypDB-36_Lmajor_ESTs.fasta")
Sequence <- paste0(
  seqinr::read.fasta(file = loc, as.string = TRUE, 
                     legacy.mode = TRUE, seqonly = TRUE, 
                     strip.desc = TRUE), 
  collapse = ""
)

Rs <- G4iMGrinder(Name = "LmajorESTs", Sequence = Sequence)
Rs2 <- G4iMGrinder(
  Name = "LmajorESTs",
  Sequence = Sequence,
  BulgeSize = 2,   
  MaxIL = 10,
  MaxLoopSize = 20
)
```

</details>

<details>
<summary>F. G4-iM Grinder’s variables</summary>

<img src="images/Variable.jpg" align="middle" height="1000" />
*(Additional PQSfinder parameters may be modified as needed.)*
</details>

<details>
<summary>G. Summarizing G4-iM Grinder results</summary>

Use `GiGList.Analysis` to consolidate results. For example:

```r
ResultTable <- GiGList.Analysis(GiGList = Rs, iden = "Predefined")
ResultTable[2, ] <- GiGList.Analysis(GiGList = Rs2, iden = "ForceLimit")
```

</details>

<details>
<summary>I. Potential Higher-Order Analysis</summary>

Method 3A (M3A) detects Potential Higher-Order Quadruplex Sequences (PHOQS). Use `GiG.M3Structure` to identify sub-unit conformations:

```r
N <- as.numeric(rownames(Rs$PQSM3a[Rs$PQSM3a$Length == max(Rs$PQSM3a$Length), ][1]))

Longest_PHOQS <- GiG.M3Structure(
  GiGList = Rs,
  M3ACandidate = N,
  MAXite = 10000
)
```

</details>

<details>
<summary>L. Searching for i-Motifs</summary>

Setting `RunComposition = "C"` targets i-Motif sequences:

```r
Rs_iM1 <- G4iMGrinder(
  Name = "LmajorESTs",
  Sequence = Sequence,
  RunComposition = "C"
)
```

</details>

<details>
<summary>M. Notes on the search engine</summary>

**G4-iM Grinder** locates overlapping or nested results that match user-defined parameters. For instance, a sequence with multiple short G-runs can generate several possible PQS overlapping one another. By default, perfect runs (e.g., `GGG`) are prioritized over slightly imperfect runs (e.g., `GCGG`) to maintain performance and to highlight the sequences most likely to form stable quadruplexes.

</details>

---

**Enjoy exploring G-quadruplexes and i-Motifs with G4-iM Grinder!**

