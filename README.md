

<h1 style = "float:letf;"> SCAN_engine </h1>

<img height = 60% src = "Psophia.jpg">

### What is SCAN about?
In Biogeography, the  field that studies the spatial distribution of species and the evolution of their environments, there is a pervasive question: - Why unrelated species often have similar geographic distributions?. Following a strict distributional approach, i.e. that considers species' distributions as the main data for analysis, until now, most methods did not use the spatial congruence between these distributions as an explicit (and controlable) parameter. And, as often remembered in the literature, congruence in geographical distributions is what lead Wallace and so many other researchers to insighst about evolution and ecology of species. This has to be the main criterion for species' distribution-based analyses in biogeography.

SCAN is totally based in species-to-species spatial relationships in an enviroment of a single-layered network analysis. Shallow network relationships allows it to recover the most obvious spatial patterns of shared distributions. When higher depths in the network are explored a plethora of other spatial relationships can be recognized, including gradients of distribution, that cannot be approapriatelly described by conventional methods.

### Why is it important?

The algorithm applies objective criteria regarding the spatial properties of entities, such as raw species' distributions, but can also be extrapolated to the analysis of environments and geographic regions. Recognized patterns may vary from highly congruent, when species have very similar distributions, to loose ones, when continuous gradients combine species overlapping or replacing each other along a transitional zone. This flexibility allows the recognition of dynamics and evolution in spatial patterns, and also allows the comparison between species and regions based on natural criteria.

 <table style = "border: 0px">
  <tr  style = "display: flex">
   <td width = 50% style = "float:left; align:center; ">
    <strong>SCAN app repository</strong>
    <br>
    <p>In this GIT version (0.21) SCAN needs a special environment until the code is updated to newer versions of R (not that easy...)</p>
    <p> Until there, you have to install <a href = "https://cran.r-project.org/bin/windows/base/old/4.2.2/"> R version 4.2.2</a> in your computer</p>
    <p> Install <a href = "https://posit.co/downloads/"> RStudio</a> and set <strong>R4.2.2</strong> as the default version (tools/global options/ R general/ R version) -> C:\Program Files\R\R-4.2.2)</p>
    <p> If needed install 'renv' package in R
    <code> install.packages("renv")</code>
    <br>
    To run SCAN_engine just type in your R console:</p>
    <code>library(shiny)</code>
    <br>
    <code>runGitHub( "cassianogatto/SCAN_engine", "cassianogatto")</code>
    <br>
  </td>
  <td width = 45% style = "float:right; align:right">
    <img width = 70%  src = "scan_maps_Icterus_Amazilia.png">
  </td>
 </tr>
 </table>

The **paper** is [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0245818) !

SCAN - *Spatial Congruence Analysis* - is a network analysis bringing an intuitive framework to search for patterns defined by spatially congruent species in biogeography: it allows a deeper understanding of species spatial relationships, and brings pure __spatial congruence__ back as a central, continuous, and explicit parameter in biogeographical analyses.

My **thesis** is also available from [here](https://repositorio.inpa.gov.br/bitstream/1/39803/3/tese_cassiano.pdf) (Intro in Portuguese, Chaps 1 & 2 in English).

Summary by AI
This doctoral thesis develops and applies a novel biogeographic method, Spatial Congruence Analysis (SCAN), to investigate the distributional patterns of Amazonian birds. SCAN identifies "chorotypes," groups of species with congruent geographic ranges, by analyzing direct and indirect spatial relationships, allowing for the detection of both highly congruent areas (similar to traditional "Areas of Endemism") and distributional gradients. The application of SCAN to Amazonian avifauna reveals a complex biogeographic structure, contradicting recent studies suggesting a lack of significant structuring, and highlights the importance of considering both historical and ecological factors in shaping species distributions. The thesis emphasizes the need to move beyond simplistic vicariance models and incorporate the nuanced interplay of various factors, such as ecological gradients and dispersal routes.

### Example
The folder <strong>example</strong> brings a small sample of New World Primate distributions to practice and build your first Chorotypes!

### Perspectives...
We are applying the method to the analysis of endemic patterns of South American Birds and Primates (with collaborators). SCAN is super intuitive, allows the gattering of tons of insights about species distributions, and now is fully converted to a standard network analysis (in R). Many network tools and concepts can now be integrated to biogeographical analysis.
