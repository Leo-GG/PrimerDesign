# Quick PCR primer design

This application creates PCR primers to amplify the sequence around a selected position in the **human** genome. It returns the pair of forward and reversed primers and then uses the [UCSC In-Silico PCR tool](https://genome.ucsc.edu/cgi-bin/hgPcr) to check for alignments along the genome.
The Tm temperatures are calculated using the [MeltingTemp](https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html) module from the [Biopython](https://biopython.org/) package.
You can use the deployed app at:
- ## [Streamlit Community Cloud](primerdesign.streamlit.app)
- ## [Ploomber]()
