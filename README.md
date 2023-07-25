# Metagenomics-data-analysis
In this project, we intend to read the contig file and perform a series of analyses on the length of the contigs. Next, the contigs will be classified using the two methods. First, we use the BusyBee website, and then we use the method that we have developed ourselves. Using the Amber tools, we report and compare their accuracy with the ground truth file, which is used as a reference. As part of our clustering approach, we first generate features for each contig based on the percentage of each nucleotide (T,G,C,A). Nucleotide percentages for contigs are collected in a data frame. A k-means clustering algorithm is then applied using categories that are similar to truth ground categories.
| Contig | A%  | C%  | G%  | T%  |
| ------ | --- | --- | --- | --- |
| Contig1| 0.3 | 0.3 | 0.2 | 0.2 |
| Contig2| 0.4 | 0.3 | 0.2 | 0.1 |
| Contig3| 0.25| 0.25| 0.25| 0.25|
| ...    | ... | ... | ... | ... |
