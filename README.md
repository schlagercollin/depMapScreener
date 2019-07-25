<img align="left" width="80" height="80" src="https://user-images.githubusercontent.com/23715298/61896673-4f629b80-aeca-11e9-994d-6443f1165eaa.png"><h1>DepMapScreener</h1>

<hr>

# Authorship

**Collin Schlager** 

<p>
Department of Biochemistry <br>
Department of Computer Science <br>
Stanford University <br>
<a href="mailto:schlager@stanford.edu">schlager@stanford.edu</a>
 </p>

# Features

 - [x] Virtually screen cancer cell lines using the Broad Institute's DepMap data
 - [x] Run virtual knock out screens on multiple genes
 - [x] Run virtual expression screens on any percentile of the population
 - [x] Run virtual lineage screens to compare cancer types
 - [x] Fully stand-alone R-shiny application wrapped in Electron Desktop framework
 
 # Credit
 
__Inspiration__

Edmond M. Chan et al., “WRN Helicase Is a Synthetic Lethal Target in Microsatellite Unstable Cancers,” Nature 568, no. 7753 (April 1, 2019): 551–56, https://doi.org/10.1038/s41586-019-1102-x.

__DepMap and Omics Data__
 
DepMap, Broad (2019): DepMap Achilles 19Q1 Public. figshare. Fileset. doi:10.6084/m9.figshare.7655150

Robin M. Meyers, Jordan G. Bryan, James M. McFarland, Barbara A. Weir, ... David E. Root, William C. Hahn, Aviad Tsherniak. Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. Nature Genetics 2017 October 49:1779–1784. doi:10.1038/ng.3984

Cancer Cell Line Encyclopedia Consortium, and Genomics of Drug Sensitivity in Cancer Consortium. 2015. Pharmacogenomic Agreement between Two Cancer Cell Line Data Sets. Nature 528 (7580):84–87. https://doi.org/10.1038/nature15736.

Jordi Barretina, Giordano Caponigro, Nicolas Stransky, Kavitha Venkatesan, William R. Sellers, Robert Schlegel, Levi A. Garraway, et. al. 2012. The Cancer Cell Line Encyclopedia Enables Predictive Modelling of Anticancer Drug Sensitivity. Nature 483 (7391):603–7. https://doi.org/10.1038/nature11003.

__Framework__

The integration of R-Shiny into Electron uses <a href="https://portableapps.com/node/32898">R-Portable<a/> and was inspired and adapted from Katie Sasso of the Columbus Collaboratory: <a href="https://github.com/ColumbusCollaboratory/electron-quick-start">https://github.com/ColumbusCollaboratory/electron-quick-start</a>. The template was modified to install and use custom packages locally from an in-app library directory.


