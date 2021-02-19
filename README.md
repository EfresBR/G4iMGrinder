# **G4-iM Grinder**
G4-iM Grinder is a fast, robust and highly adaptable algorithm. It is capable of locating, identifying, qualifying and quantifying DNA and RNA potential quadruplex structures, such as G-quadruplex, i-Motifs and their higher order versions.

Read the open-access paper on the algorithm: [G4-iM Grinder: when size and frequency matter. G-Quadruplex, i-Motif and higher order structure search and analysis tool](https://bit.ly/3j2UUjC)

Soon to be published: [Exploring G- and C-quadruplex structures as potential targets for the severe acute respiratory syndrome coronavirus 2](https://bit.ly/3j1jFMP) (Preprint version).

The Results of both manuscripts can be found in the result section.  

<div style="text-align:center"><img src="https://www.biorxiv.org/content/biorxiv/early/2019/06/17/532382/F1.large.jpg?width=800&height=600&carousel=1" height="500" /></div>

Please report in the issue section if you find bugs, problems or wish to see some features added.

###				**Quadruplexes and G4-iM Grinder**

*What are quadruplexes?*
<details>

[**G-quadruplexes (G4s)**](https://en.wikipedia.org/wiki/G-quadruplex): G4s are DNA or RNA sequences that are rich in guanine, where four guanine bases can associate through Hoogsteen hydrogen bonding to form a square planar structure called a guanine tetrad (G-tetrad or G-quartet), and two or more guanine tetrads (from G-tracts, continuous runs of guanine) can stack on top of each other to form a G4,stabilized by monovalent cation/s. The unimolecular versions occur naturally near the ends of the chromosomes, better known as the telomeric regions, and in transcriptional regulatory regions of multiple genes, both in microbes and across vertebrates including oncogenes in humans.

<div style="text-align:center"><img src ="https://upload.wikimedia.org/wikipedia/commons/thumb/5/5c/G-quadruplex.svg/600px-G-quadruplex.svg.png" height="150"/></div>

[**C-quadruplexes or i-Motifs (iM)**](https://en.wikipedia.org/wiki/I-motif_DNA): iMs are four-stranded quadruplex structures formed by cytosine-rich DNA or RNA, similar to the G-quadruplex structures that guanine-rich genomes forms. C-rich DNA regions are common in gene regulation portions of the genome. Recently, iM were discovered in human cells and were shown to play a role in cell reproduction. iM have potential applications in nanotechnology and nanomedicine, because size is more than 1 nm and less than 100 nm due to their unique pH sensitivity and have been used as biosensors, nanomachines, and molecular switches.

<div style="text-align:center"><img src ="https://external-content.duckduckgo.com/iu/?u=https%3A%2F%2Fupload.wikimedia.org%2Fwikipedia%2Fcommons%2Fthumb%2F5%2F5d%2FRNA_i-motif.gif%2F220px-RNA_i-motif.gif&f=1&nofb=1"  height="150"/></div>

[...adapted from wikipedia...]</br>

</details>
</br>

*Searching for quadruplexes*?
<details>
</br>
Quadruplexes have been studied abundantly in the last few years due to the increasing evidence of their functional role in many living organisms, yet the natural properties by which they form and work are very much unknown. To identify new structures, in silico predictions are based on in vitro verified paradigms. Loops, tetrad number, run imperfections and the flanking regions of the structures all seem to play important roles in the topology and dynamics of these secondary structures.

G4-iM Grinder (GiG) is an algorithm that incorporates a quadruplex search engine and several quadruplex qualifiers functions to locate these genomic structures.
1. **GiG's quadruplex search engine** was designed to be very tolerant as to allow the detection of all possible candidates that fulfil the user-established criteria.
2. **GiG's quadruplex qualification engine** was designed to filter the results of the search-engine. The filters are modular. They can include: scoring filters (which calculate the probability of the sequence forming an actual quadruplex or i-Motif), location in the genome filters, presence of a desired pattern filter, frequency of appearance in the genome filters, and presence of already known-to-form quadruplexes.

</br>
<div style="text-align:center"><img src ="https://www.researchgate.net/profile/Efres_Belmonte-Reche/publication/336374922/figure/fig1/AS:812376744013825@1570697117586/A-G4-iM-Grinders-workflow-when-RunComposition-G-to-find-PQS-and-PHOQS-in-a-genome.ppm" height="600"/>
</div>
</details>
</br>

**Abbreviations and other important information**
1. **Quadruplex**: 			A genomic sequence that forms an G-Quadruplex (G4) or an i-Motif (iM).
2. **G4s**: 						A sequence which has been confirmed to form a G4 *in vitro*.
3. **iM**:							A sequence which has been confirmed to form an iM *in vitro*.
4. **PQS**: 						A Potential G4 Sequence detected with an *in silico* algorithm, which has **NOT** been confirmed *in vitro*, at least initially. It requires of biophysical or biochemical methods to confirm.
5. **PiMS**: 						A Potential iM Sequence detected with an *in silico* algorithm, which has **NOT** been confirmed *in vitro*, at least initially. It requires of biophysical or biochemical methods to confirm.

###       **G4-iM Grinder's Updates & News**
G4-iM Grinder latest version: 1.6.1 (02-2021).
<details>
Changes:

For Version 1.6.1:
* Adapted and further optimized `GiG.df.GenomicFeatures`.
* Changed `GiGList.Analysis` to accept vectors instead of just single numerals in its parameters. Changed also the results returned, with better summaries of the G4-iM Grinder function.
* Changed the concept of DNA and RNA sequences in the `G4-iM Grinder` main function and other related functions to be more efficient.
* The function to find confirmed quadruplex sequences was modified to be more efficient with RAM. This is to prevent problems with an growing database of confirmed quadruplex sequences.
* Added Biostrings and biomartr dependencies. Added the packages to the package loading function.
* G4-iM Grinder version and the G4-iM Grinder database version are now saved in the configuration dataframe of each result.
* Added a function to analyze the characteristic and runs of a genome (`GiG.Seq.Analysis`).
* Added a function to analyze the biological landmarks affected by the potential quadruplex results (`GiG.df.GenomicFeatures`).
* Changed how packages are loaded so they are silent when doing so. The function will also check if all dependencies are installed and R version is at least 4.0. If any of these fail, and error will be returned asking the user to fix the problem/s before proceeding with the G4-iM Grinder analysis.

For Version 1.5.95:

1. Fixed bug in PQSfinder algorithm, which incorrectly punctuated structures.
2. Changed the way known G4s and i-Motifs structures are detected. It will now detect both DNA and RNA confirmed sequences within the results.
a. If the confirmed sequence is DNA, the results will include an asterisk (&ast;).
b. If the result sequence is RNA, it will include a circumflex (^).
c. Example: If the GUK1 DNA quadruplex was detected within the results one time, the Conf.Quad.Seqs column will state: GUK1 (1&ast;). If the 42.HIRA (WT) RNA quadruplex was detected within the results one time, the Conf.Quad.Seqs column will state: 42.HIRA (WT) (1^)
</details>
</br>


Version 2.51 of G4-iM Grinder's database (GiG.DB) is also here (09-2020). Now it includes the currently confirmed G4 and iM sequences known-to-form in the SARS-CoV-2, and over **2800 sequences** other confirmed to form (or not to form) tetraplex structures in literature.
<details>

The GiG.DB within the G4-iM Grinder package includes:

**I. BioInformatic dataframe:**

	1. Each entry is a nucleotidic sequence published in a scientific journal in relationship with its capability of forming quadruplex structures.
	2. Each entry includes 	
		A. Quadruplex		TRUE for  forming quadruplex, FALSE for NOT
		B. Genome		DNA or RNA
		C. Nucleotide		G or C, for G4s or i-Motif respectively
		D. Name			value must be unique
		E. Sequence		value must be unique
		F. Length		Length of Sequence
		G. Tm			Nº of biophysical results associated to the entry (within Biophysical dataframe)
	3. Currently there are a total of 2851 entries.
		A. 2141 form tetraplex	 and 710  dont;
		B. 283  are i-Motifs 	 and 2568 are G4s;
		C. 1858 are DNA 	 and 993  are RNA.
	4. Sequences which end in -ReV- are the reverse sequences of other entries.
		For example 	
			1. Name1 	GGTGGTGG|TTT|GG
			2. Name1-ReV- 	GG|TTT|GGTGGTGG

**II. Refs dataframe:**

	1. Each entry is the literature reference for an BioInformatic dataframe entry.
	2. Each entry includes:
		A. name 		value must be unique; Name of BioInformatic entry
		B. DOI			DOI identificator; for example: 10.1093/nargab/lqz005
		C. Pubmed	 	PubmebID identificator (PMID); for example:	29109402
		D. comments		Extra information, normally citing information
					For example: Nucleic Acids Res., 45, 7487–7493.
	3. Currently there are a total of 2851 entries.

**III. BioPhysical dataframe:**

	1. Each entry is a Biophysical result found for a particular BioInformatic entry.
	2. Data includes Tm (ºC), pH, Concentrations of sequence (uM), K+ (mM) and Na+ (mM), and the found topology.
	3. Currently there are 153 entries.

Comments: If you find an error within GiG.DB or want to include other sequences, please open an issue request in Github, **"EfresBR/G4iMGrinder"**.

</details>
</br>

###       **G4-iM Grinder's Results**

Please report in the issue section if any links are broken.

**SARS-CoV-2** reference genome (GCF_009858895.2) analysis using G4-iM Grinder and a lax quadruplex-configuration of parameters, can be found in different formats.
* **[.RDS](http://bit.ly/3drDhdM)** (R, recommended)
* **[.gff3](http://bit.ly/3ue7Edi)**
* **[.xlsx](http://bit.ly/3s4SDc3)**  (excel)</br>
They includes all the positions of the PQS and PiMS found, the conservation, and the G4s and iMs already confirmed in the manuscript and elsewhere.

**VARIANTS** of **SARS-CoV-2** found in other lineages and clades and not in the reference genome, can be found in **[.RDS](http://bit.ly/37rutk0)** format. Each entry has the (GISAID) ID, name, Lineage and Clade of the SARS-COV-2 genome where it was located. </br>
Full GiG Results of the **VIRUS REALM** and genome **REFERENCES**, analyzed using a lax quadruplex-configuration of parameters can be found next.

<details>

**ANALYTICAL DATA**:

* [Analysis.RData](http://bit.ly/3qxZH0v), is the analysis results on the raw G4-iM Grinder data. It includes 3 lists.
1.	Analysis.Coronaviridae.fam – Analysis with GiGList.Analysis function of the GiG-package of the Coronaviridae family. PQS and PiMS lists are the analysis for PQS and PiMS respectively. df.index data frame stores the identification of each genome used.
2.	Analysis.Virus.realm - Analysis with GiGList.Analysis function of the GiG-package of the entire virus realm. PQS and PiMS lists are the analysis for PQS and PiMS respectively. df.index data frame stores the identification of each genome used. Genome data frame is the analysis with the function GiG.Seq.Analysis.
3.	Baltimore.C – Baltimore Classification tables regarding each group characteristics and classification of each family into its group.


**RAW DATA**:
* [Virus.Results.RDS](http://bit.ly/3sanTqc), (2.4 Gb) includes the raw data of the G4-iM Grinder analysis on all the virus realm as a list. The list groups virus species by their families. Each species list includes a PQS and PiMS sublist. These store the composition, location, known-quadruplex sequences presence and score (amongst others) of PQS/PiMS found in each virus. The information used in this analysis was Method 2; size restricted overlapping search method (PQSM2A data.frames), although Method 3 results are also included.
* [GISAID.refs.rar](http://bit.ly/3s5X4n4), includes the references of the  17312 SARS-CoV-2 genomes downloaded from the GISAID database.



</details>
</br>



###### **SARS-CoV-2** confirmed quadruplex-related sequences

|Quadruplex	|Genome	|Nucleotide	|Name					|Sequence										|References																																																	|
|---		|---	|---		|---					|:---:											|:---:																																																		|
|TRUE		|DNA	|G			|CoVID.DNA.G4.1			|**GG**CT**GG**CAAT**GG**C**GG**				|[10.1101/2020.08.19.257493](https://bit.ly/3j1jFMP)																																						|
|TRUE		|DNA	|C			|CoVID.DNA.iM.1			|**CC**TT**CC**TTA**CCC**AGAT**CC**				|[10.1101/2020.08.19.257493](https://bit.ly/3j1jFMP)																																						|
|TRUE		|RNA	|G			|CoVID.RNA.G4.1			|**GG**CU**GG**CAAU**GG**C**GG**				|[10.1101/2020.08.19.257493](https://bit.ly/3j1jFMP)	&	[10.1002/anie.202011419](https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.202011419)	&	[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)													|
|TRUE		|RNA	|G			|CoVID.RNA.G4.2			|**GG**A**GG**A**GG**UGUUGCA**GG**				|[10.1101/2020.08.19.257493](https://bit.ly/3j1jFMP)																																						|
|TRUE		|RNA	|G			|CoVID.RNA.G4.3			|**GG**UAUGU**GG**AAA**GG**UUAU**GG**			|[10.1093/bib/bbaa114](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa114/5849532)	&		[10.1002/anie.202011419](https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.202011419)	&	[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)|
|TRUE		|RNA	|G			|CoVID.RNA.G4.4			|**GG**CUUAUA**GG**UUUAAU**GG**UAUU**GG**		|[10.1093/bib/bbaa114](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa114/5849532)	&		[10.1002/anie.202011419](https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.202011419)	&	[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)|
|TRUE		|RNA	|G			|CoVID.RNA.G4.8			|**GG**UU**GG**ACCUUU**GG**UGCA**GG**			|[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)																											|
|FALSE		|RNA	|C			|CoVID.RNA.iM.1			|**CC**UU**CC**UUA**CCC**AGAU**CC**				|[10.1101/2020.08.19.257493](https://bit.ly/3j1jFMP)																																						|
|FALSE		|DNA	|C			|SARS-DNA.iM.1			|**CC**TG**CC**TTA**CCC**AGAT**CC**				|[10.1101/2020.08.19.257493](https://bit.ly/3j1jFMP)																																						|
|FALSE		|RNA	|G			|CoVID.RNA.G4.1.MUT		|GACUGACAAUGACG									|[10.1002/anie.202011419](https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.202011419)																													|
|FALSE		|RNA	|G			|CoVID.RNA.G4.1.MUT1	|**GG**CUAGCAAU**GG**C**GG**					|[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)																											|
|FALSE		|RNA	|G			|CoVID.RNA.G4.1.MUT2	|**GG**CUAACAAU**GG**C**GG**					|[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)																											|
|FALSE		|RNA	|G			|CoVID.RNA.G4.3.MUT1	|**GG**UAUGUAGAAA**GG**UUAU**GG**				|[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)																											|
|FALSE		|RNA	|G			|CoVID.RNA.G4.3.MUT2	|**GG**UAUGUAAAAA**GG**UUAU**GG**				|[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)																											|
|FALSE		|RNA	|G			|CoVID.RNA.G4.5			|**GG**CCAU**GG**UACAUUU**GG**CUA**GG**			|[10.1002/anie.202011419](https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.202011419)																													|
|FALSE		|RNA	|G			|CoVID.RNA.G4.6			|**GG**CUUU**GG**AGACUCCGU**GG**A**GG**A**GG**	|[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)																											|
|FALSE		|RNA	|G			|CoVID.RNA.G4.7			|**GG**UAAUAAA**GG**AGCU**GG**U**GG**			|[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)																											|
|FALSE		|RNA	|G			|CoVID.RNA.G4.8.MUT1	|**GG**UUAGACCUUU**GG**UGCA**GG**				|[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)																											|
|FALSE		|RNA	|G			|CoVID.RNA.G4.8.MUT2	|**GG**UUAAACCUUU**GG**UGCA**GG**				|[10.3389/fgene.2020.587829](https://www.frontiersin.org/articles/10.3389/fmicb.2020.567317/full)																											|

|
###### **SARS-CoV-2** quadruplex-related biophysical data

|Name				|Tm(ºC)	|[K+]mM	|[Na+]mM	|Topology	|References																							|
|:---:				|:---:	|:---:	|:-:		|:---:		|:---:																								|
|CoVID.RNA.G4.1		|57.7	|100	|0			|parallel	|[10.1002/anie.202011419](https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.202011419)			|
|CoVID.RNA.G4.3		|41.2	|100	|0			|			|[10.1002/anie.202011419](https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.202011419)			|
|CoVID.RNA.G4.4		|-1		|100	|0			|unformed	|[10.1002/anie.202011419](https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.202011419)			|
|CoVID.RNA.G4.5		|-1		|100	|0			|unformed	|[10.1002/anie.202011419](https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.202011419)			|
|CoVID.RNA.G4.1.MUT	|-1		|100	|0			|unformed	|[10.1002/anie.202011419](https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.202011419)			|
|CoVID.RNA.G4.1		|54.4	|50		|0			|parallel	|[10.1101/2020.08.19.257493](https://bit.ly/3j1jFMP)												|
|CoVID.RNA.G4.2		|48.1	|50		|0			|parallel	|[10.1101/2020.08.19.257493](https://bit.ly/3j1jFMP)												|
|CoVID.RNA.G4.3		|31.5	|150	|0			|parallel	|[10.1093/bib/bbaa114](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa114/5849532)|
|CoVID.RNA.G4.4		|30.5	|150	|0			|parallel	|[10.1093/bib/bbaa114](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa114/5849532)|


</br>

G4-iM Grinder's results in **HUMANS** and **49 PATHOGENIC SPECIES** analysed using the predefined parameter configuration of the algorithm and published in the article [G4-iM Grinder: when size and frequency matter. G-Quadruplex, i-Motif and higher order structure search and analysis tool](https://bit.ly/3j2UUjC), can be found through [this link](https://1drv.ms/u/s!AvVGQg2rNIwDgTeth6qclA8Rz5UM?e=gmEI1a)

<details>

GiG.DB V.2.5 has been used to update these results (03-2020) and now include the localization of the **2851** known-to-form and known-NOT-to-form quadruplex in the database.

As of V.2.5 of GiG.DB, the total amount of results (M2A) with at least one confirmed G4 within its sequence is **312072** (236483 more than in V1.0).

As of V.2.5 of GiG.DB, the total amount of results (M2A) with at least one confirmed i-Motif within its sequence is **160054** (74171 more than in V1.0)

The 1.5 Gb .RAR compressed file hosts four RData images of the results.

1. `Human.PQS.032020.RData` for Human G-based PQS analysis
2. `Human.PiMS.032020.RData` for Human C-based PiMS analysis
3. `NonHuman.PQS.032020.RData` for non-human G-based PQS analysis
4. `NonHuman.PiMS.032020.RData` for non-human C-based PiMS analysis.


Genomes used:

	1. Human Genome - hg38, GRCh38.p12, Genome Reference Consortium Human Build 38, INSDC Assembly GCA_000001405.27 downloaded May 2019 from www.sanger.ac.uk.
	2. Non-human genomes - Please see section 9 of supplementary material of the original article for more info.

With this update, Figure 4 of G4-iM Grinder's article, which compared different tetraplex-related characteristics of each genome (including density [per 100000 nucleotides], uniqueness and Confirmed Quadruplex Sequences (CQS)) becomes:

<img src="images/Data.Analysis.V1.59, V2.5.jpg" align="middle" height="1250" />

Being the CQS columns what changes between both Figure 4s.

</details>


###       **G4-iM Grinder's Installation and User Guide**

<details>

####       **A.      Package prerequisites**

<details>

G4-iM Grinder can be downloaded from github: EfresBR/G4iMGrinder. G4-iM Grinder requires the installation of other CRAN based and Bioconductor packages.
Please, ensure all required packages are installed and R version is at least 4.0.0.
G4-iM Grinder was successfully downloaded and tested in MacOS 10.12.6, Windows 10 (x64), Ubuntu 18.04.2 (x64), Mint 19.1 (x64) and Fedora-workstation 30.
In Linux based systems, the installation of devtools may require further effort ([Check this link](https://stackoverflow.com/questions/20923209/problems-installing-the-devtools-package)).
Other OS including x86 systems have not been tested.

G4-iM Grinder has been successfully used in R 4.0.3 and R-studio 1.3.1093


```ruby

pck <- c("stringr", "stringi", "plyr", "seqinr", "stats", "parallel", "doParallel", "beepr", "stats4", "devtools", "dplyr", "BiocManager", "tibble")

#foo was written by Simon O'Hanlon Nov 8 2013.
#Thanks Simon, thanks StackOverflow and all its amazing community.

foo <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}
foo(pck)
BiocManager::install(c("BiocGenerics", "S4Vectors", "Biostrings", "biomartr", "IRanges") , ask = FALSE, update = TRUE)


```

</details>



####      **B.      Package installing and loading**

<details>

```ruby

devtools::install_github("EfresBR/G4iMGrinder")
library(G4iMGrinder)


```

</details>



####      **C.      Installation fails**

<details>

The most common reasons for failing during the installation of G4-iM Grinder are ,

1. 	Some of G4-iM Grinder's dependencies have not been installed,
2. 	R version is not at least 4.0.0

If you are having problems during installation, please, execute the following code to verify that these prerequisites are met.

```ruby


pck <- c("BiocGenerics", "S4Vectors", "stringr", "stringi", "plyr", "seqinr", "stats", "parallel", "doParallel", "beepr", "stats4", "devtools", "dplyr", "BiocManager", "biomartr", "Biostrings")

FailFoo <- function(x){
  Info <- "Package dependendies FAILED. These packages are required and are NOT installed: "
  count <- 0
  for( i in x ){
    if( ! require( i , character.only = TRUE, quietly = TRUE ) ){
      Info <- paste0(Info, i, " ")
      count <- count +1
    }
  }
  ifelse(count ==0, yes = print("Package dependencies PASSED. All required packages are installed. "),
         no = print(Info))
  AAA <- R.version
  ifelse(as.numeric(AAA$major) == 4,
         yes= ifelse(as.numeric(AAA$minor >= 0),
                     yes = print("R version requirements PASSED. R version is at least 4.0 as required."),
                     no = print("R version requirements FAILED. R needs to be updated to version >= 4.0")),
         no = print("R version requirements FAILED. R needs to be updated to version >= 4.0"))
}
FailFoo(pck)


```

The result of this code should be:

```

[1] "Package dependencies PASSED. All required packages are installed. "
[1] "R version requirements PASSED. R version is at least 4.0 as required."

```

If both the package dependencies and R version have passed the test, and still the installation fails, please, write an issue in the issue section stating the transcript of the executed commands and the full error received.

</details>



####       **D.      (NEW) Running a G4-iM Grinder pre-analysis**


<details>

Executing a genomic pre-analysis with `GiG.Seq.Analysis`.
This function can be used before a GiG analysis to determine the best search parameters to obtain quadruplex-related results.
The function’s outcome is a data frame with the most relevant genomic features, including length (in nucleotides), type of genome (DNA or ARN), strands (single or double), and G, C, T/U, A and N composition (as % of total sequence).
The function also calculates the total number of runs with different conditions (predefined parameters, bulges per run: zero and one-quantities; run lengths: two to five and three to five-length) in the genome, and returns it to the user as total counts or genomic density.
The higher the run density, the higher the probability of finding associated PQS or PiMS in the results.
```ruby

# Using a genome available online
loc <- url("http://tritrypdb.org/common/downloads/release-36/Lmajor/fasta/TriTrypDB-36_Lmajor_ESTs.fasta")
Sequence <- paste0(seqinr::read.fasta(file = loc, as.string = TRUE, legacy.mode = TRUE, seqonly = TRUE, strip.desc = TRUE), collapse = "")

# Running the pre-analysis.
require(G4iMGrinder)
Pre_Rs <- GiG.Seq.Analysis(Name = "LmajorESTs", Sequence = Sequence, DNA = TRUE, Complementary = TRUE)

```



</details>



####       **E.      Running a G4-iM Grinder analysis**


<details>

Executing a genomic G-Quadruplex analysis with G4iMGrinder function

```ruby

# Using a genome available online
loc <- url("http://tritrypdb.org/common/downloads/release-36/Lmajor/fasta/TriTrypDB-36_Lmajor_ESTs.fasta")
Sequence <- paste0(seqinr::read.fasta(file = loc, as.string = TRUE, legacy.mode = TRUE, seqonly = TRUE, strip.desc = TRUE), collapse = "")

# Executing a grind on the sequence in search of PQS
require(G4iMGrinder)
Rs  <- G4iMGrinder(Name = "LmajorESTs", Sequence = Sequence)

# Forcing the folding rule to the limit (this will take longer)
Rs2 <- G4iMGrinder(Name = "LmajorESTs", Sequence = Sequence, BulgeSize = 2,   MaxIL = 10, MaxLoopSize = 20)


```
G4-iM Grinder allows huge flexibility to adapt to any of the users requirements.

</details>



####       **F.      G4-iM Grinder's variables and their predifined values**

<details>


<img src="images/Variable.jpg" align="middle" height="1000" />
N.B. Several other parameters regarding PQSFinder are available for modification.

</details>



####       **G.      Summarizing G4-iM Grinder results**


<details>
Summarizing an analysis with GiGList.Analysis function to compare the results between genomes. This will quantify the number of results and density of each analysis. It will also give the number of results that have at least a minimum frequency, score and size. These variables can be modified. See the package documentation for more information regarding GiGList.Analysis.

```ruby

# summarizing first search
require(G4iMGrinder)
ResultTable <- GiGList.Analysis(GiGList = Rs, iden = "Predefined")

# adding the second analysis in a new row
ResultTable[2,] <- GiGList.Analysis(GiGList = Rs2, iden= "ForceLimit")


```

</details>



####       **H.      (NEW) Biological features (landmarks) affected by PQS and PiMS candidates**  

<details>

The `GiG.df.GenomicFeatures` function is suitable for determining the genomic features that share their location with (and hence may be affected by) GiG’s PQS and PiMS results.
It employs the online database connector package “biomartr” to retrieve the genomic annotations file for the sequence, with which to then match positions.
The function returns a data frame of all the matches found for the input sequences and includes different attributes (IDs, keys, relationships with other features and comments) of the matched genomic features.
Please, use the same genome associated with the annotation file.

```ruby

# Analyzing the HIV-1 virus. To do so, first lets download the genome and use it with G4-iM Grinder via the biomartr package. The virus is a ssRNA.
require(G4iMGrinder)
require(biomartr)
Sequence <- toString(read_genome(getGenome(db = "refseq", organism = "GCF_000864765.1",  reference = F)))
RsHIV <- G4iMGrinder(Name = "HIV-1", Sequence = Sequence, DNA = F, Complementary = F)

# Applying the `GiG.df.GenomicFeatures` function on the Method 2 results of the G4-iM Grinder grind.
RsHIV.GF <- GiG.df.GenomicFeatures(df = RsHIV$PQSM2a, org = "GCF_000864765.1", db = "refseq")

#Please make sure the genome analyzed with G4-iM Grinder and the annotation file are from the same organism.

```

</details>



####       **I.      Potential Higher Order Analysis**  


<details>

Executing an analysis of a higher order structure with GiG.M3Structure to analyze its potential subunit configuration. This will give all and the most interesting subunit conformations as stated in the article. See the package documentation for more information regarding GiG.M3Structure.

```ruby

# analyzing the longes PHOQS structure found in Rs$PQSM3A.

# N is the row number of the PHOQS to analyze in PQSM3a, as a numeral.
N <- as.numeric(rownames(Rs$PQSM3a[Rs$PQSM3a$Length == max(Rs$PQSM3a$Length),][1]))

require(G4iMGrinder)
Longest_PHOQS <- GiG.M3Structure(
			GiGList = Rs,
			M3ACandidate = N,
			MAXite = 10000
			)


```

</details>



####       **J.      Locating the references of Known-To-Form and Known-NOT-To-Form sequences**


<details>

Finding the reference for the Known-To-Form Quadruplex structures of an interesting Result. This procedure is the same for Known-NOT-To-Form sequences.

```ruby

# Finding the references of the known-to-form sequence 93del.
require(G4iMGrinder)
Ref93del <- GiG.DB$GiG.DB.Refs[GiG.DB$GiG.DB.Refs$Name == "93del", ]


```

</details>



####       **K.      Updating results for a pre-existing analysis**


<details>
Updating a G4-iM Grinder analysis with different variables using the GiGList.Updater function. This will avoid doing a new search analysis on the sequence and hence will be more time and resource efficient.

```ruby

# As the PHOQS structure in row 126 looks promising, we will also examine
# the Known-NOT-to-form Quadruplex of the results, quantify the % of GGG and TTA present in the sequence,
# and modify the score and frequency weight of the final score.
require(G4iMGrinder)
Rs3 <- GiGList.Updater(GiGList = Rs, KnownNOTQuadruplex = TRUE, KnownQuadruplex = TRUE,
                       LoopSeq = c("GGG", "TTA"), FreqWeight = 100, WeightParameters = c(75, 25, 0))


```

</details>



####       **L.      Grinding genomes in search of Potential i-Motif Sequences (PiMS)**


<details>

To search for potential i-Motifs in the genome we can repeat the analysis with G4iMGrinder function changing RunComposition = “C”.

```ruby

# Doing a grind in search for i-Motifs in the sequence
require(G4iMGrinder)
Rs_iM1 <- G4iMGrinder(Name = "LmajorESTs", Sequence = Sequence, RunComposition = "C")


```

</details>



####       **M.      Comments on G4-iM Grinder's Search Engine**

<details>

G4-iM Grinder locates all overlapping and nested results that fit the user-defined (or predefined if none were inserted) parameters.
For example using predefined parameters, five possible PQS (in _italics_) results will be located for the genomic sequence

> **GGGG**TTAT**GGG**TTATT**GGTGG**TTATT**GGCG**TT**GGG**

1.	_**GGGG**TTAT**GGG**TTATT**GGTGG**TTATT**GGCG**_(~~TTGGG~~)
2.	_**GGGG**TTAT**GGG**TTATT**GGTGG**TTATT**GGCG**TT**GGG**_
3.	(~~G~~)_**GGG**TTAT**GGG**TTATT**GGTGG**TTATT**GGCG**_(~~TTGGG~~)  
4.	(~~G~~)_**GGG**TTAT**GGG**TTATT**GGTGG**TTATT**GGCG**TT**GGG**_
5.	(~~GGGGTTAT~~)_**GGG**TTATT**GGTGG**TTATT**GGCG**TT**GGG**_

The only current limitation of the search engine is when a perfect (for example, **GGG**) and an imperfect (for example, **GCGG**) run coexist within the same run (for example, **GCGGG**). Although it is possible that **GCGGG** forms a run, the perfect run (**GGG**) is favored to improve computing performance and the location of more likely to form sequences. For the Genomic Sequence

> **GCGGG**TTA**GGG**TTATTT**GGG**TTA**GGG**

using predefined parameters will result in the detection of:

-	(~~GC~~)_**GGG**TTA**GGG**TTATTT**GGG**TTA**GGG**_

whilst

-	_**GCGGG**TTA**GGG**TTATTT**GGG**TTA**GGG**_

will not be detected.

Regarding frequency of the quadruplex results, Quadruplexes may actually be repeated because they form part of repetitive nucleotide sequences, including transposon families. For example, several authors have already located recurrent PQS in such repetitive elements (both human and non-human species), which depending on the location and context, may potentially grant different biological significance to the same recurrent quadruplex.

</details>



</details>
