# NTD-MANE
Is a team of young researchers primarily based at the University of Buea, Cameroon involved neglected tropical diseases (NTDs) research
for vaccines and dignostic tools development. We use immunoinformatics, proteomics, bioinformatics and immunological techniques 
to predict and investigate potential vaccine and diagnostic candidates.

## Immunoinformatics Guide (Pipeline used in our March [2019 pubplication]() )

### Selection and Retrieval of Proteins
Proteins previously reported as antigen/immunogenic were retrieved from UniProt and subjected 
to the following pipeline;

- **Signal peptide prediction** using [SignalP](http://www.cbs.dtu.dk/services/SignalP/) and 
Localization prediction using [DeepLoc](http://www.cbs.dtu.dk/services/DeepLoc/). This step 
informs us of which are potentially secretory proteins. It is important to have this 
information because secretory proteins would be better antigens/immuogens since they can be 
readily accessed by the immue system.

### Immune epitope prediction
- **Linear B-cell epitope prediction** using [BepiPred](http://www.cbs.dtu.dk/services/BepiPred/), 
[ABCpred](http://www.imtech.res.in/raghava/abcpred/), [SVMTriP](http://sysbio.unl.edu/SVMTriP/)
, and [BCPreds](http://ailab.ist.psu.edu/bcpred/). The qualities of an ideal B-cell epitope are
 important to inform our choice of epitopes. For example, they are usually 11-20 amino acids 
long (You’ll find different ranges but they are generally longer relative to T-cell epitopes 
that are usually 4-8aa longs). Again, you should understand how the scoring system of the tools
 to inform your selection of epitopes. For example, SVMTriP is score based and predicts 
tri-peptides based on similarity and propensity scores which can be used to pick best epitopes.
 You many want to settle for only epitopes that were predicted by all the tools.

- **Cytotoxic T-cell (CTL) epitope prediction** using [NetCTL 1.2](http://www.cbs.dtu.dk/services/NetCTL/).
- **Helper T-cell (HTL) epitope prediction** using [NetMHCII 2.2](http://www.cbs.dtu.dk/services/NetMHCII/).

Again, we have to know the properties of an ideal T-cell epitope. Furthermore, T-cell epitope 
prediction has to be informed by the MHC class that is predominantly involved in the immune 
response of the condition you are dealing with. For example, in our study, we limited our 
prediction of CTL epitope prediction to the MHC class I, A1 supertype. And we settled for high 
affinity HTL epitopes.

- **Construction of multi-epitope vaccine candidate:** We linked the epitopes with AAY and GPGPG
 linkers. In addition, we used an internal adjuvant in the form of a Mycobacterium tuberculosis
 ribosomal protein. Here, you just put together your sequence. 
There is no abstract technique in doing this. We essentially used MicroSoft PowerPoint for visualization.

### Accessing the properties of the chimeric peptide
- **IFN-gamma inducing epitopes** using [IFNepitope server](http://crdd.osdd.net/raghava/ifnepitope/scan.php). 
IFN-gamma is important in innate and adaptive immunity in that it stimulates macrophages an natural killer cells. 

- **Antigencity and Allergenicity predcition**. [ANTIGENpro](http://scratch.proteomics.ics.uci.edu/) and 
[Vaxijen v2.0](http://www.ddgpharmfac.net/vaxijen/VaxiJen/VaxiJen.html) predict antigenicity 
while [AllerTOP v2.0](http://www.ddg-pharmfac.net/AllerTOP) and 
[AllergenFP](http://ddg-pharmfac.net/AllergenFP/) predict allergenicity. 
We need to know if our peptide is a potential allergen as this would be an undesirable feature.

- **Physicochemical properties:** [ProtParam](http://web.expasy.org/protparam/) would give basic
 protein properties like MW, solubility, stability etc. 
[PROSO II](http://mbiljj45.bio.med.uni-muenchen.de:8888/prosoII/prosoII.seam) would also give 
information on the solubility.

- **Secondary structure Prediction** using [PSIPRED](http://bioinf.cs.ucl.ac.uk/psipred/) and 
[RaptorX](http://raptorx.uchicago.edu/StructurePropertyPred/predict/). Ranked structures are 
selected based on scoring.

- **Tertiary structure prediction** using [I-TASSER server](https://zhanglab.ccmb.med.umich.edu/I-TASSER/). 
Again, ranked structures are selected based on score. There is a huge repertoire of software 
that one can use here and these are classified by a community of protein structure prediction 
scientists called CASP.

- **Refinement of the tertiary structure** using [ModRefiner](https://zhanglab.ccmb.med.umich.edu/ModRefiner/) 
and the [GalaxyRefine server](http://galaxy.seoklab.org/cgi-bin/submit.cgi?type=REFINE).

- **Validation of tertiary structure** using [ProSA-web](https://prosa.services.came.sbg.ac.at/prosa.php). 
The [ERRAT server](http://services.mbi.ucla.edu/ERRAT/) would analyze non-bonded atom-atom 
interactions and the [RAMPAGE server](http://mordred.bioc.cam.ac.uk/~rapper/rampage.php) 
would produce a Ramachandran plot showing the proportion of the amino acids in the chimeric 
peptide that fall in energetically allowed and disallowed regions.

- **Discontinuous B-cell epitope predcition** (conformational epitopes) using [ElliPro](http://tools.iedb.org/ellipro/). 
You need to have an idea of the epitope the chimera forms when it folds. 
ElliPro produces scored epitopes which can inform you of the antigenicity/immunogenicity.

- **Molecular docking of the chimera with an appropriate immune receptor** using the [CASTp server](http://sts.bioe.uic.edu/castp/), 
[CPORT](https://milou.science.uu.nl/services/CPORT/), 
[HADDOCK 2.2](http://haddock.science.uu.nl/services/HADDOCK2.2), and 
[PRODIGY](https://nestor.science.uu.nl/prodigy/). Based on literature and previous findings, 
you should have some knowledge of the immune interaction of the condition you are dealing. 
In our case, previous findings hold that the toll-like receptor 4 (TLR4) is implicated in 
immunity to onchocerciasis. CASTp predicts binding pockets, CPORT predicts active interface 
residues which is used to inform docking. Then docking is perfomed with HADDOCK, and PRODIGY 
predicts the binding affinities.

- **In silico cloning of the chimera:** You can reverse translate the peptide to the 
corresponding nucleotide sequence using the [JCat server](http://www.prodoric.de/JCat).

- **Immune Simulation** using the [C-ImmSim server](http://150.146.2.1/C-IMMSIM/index.php). 
This software simulates the immune response based on your vaccination scheme. 
So you need to have some knowledge on the vaccination scheme. 
For example, how many injections are to be administered and the duration after each injection.

---
Here are just some of the myriad of tools available in the field on in silico vaccine design. 
You can substitute any of the tools to get a more desired outcome. 
You can as well add more tools in any one step to increase the validity of the outcome.

-------------------
Some relevant online courses and materials related to this field of research
----
- Bioinformatics (DNA and Protein sequence Analysis, Statistical Analysis in Bioinfo): https://www.edx.org/micromasters/usmx-umuc-bioinformatics
- MD Simulation related: https://www.edx.org/course/atoms-materials-predictive-theory-purduex-mse550x
- Linux Fundamentals: https://1drv.ms/b/s!AoM4yKosJgwohSfRwfxia-sgIuOe
