# MDD-DwC-mapping

The script employs mapping from MDD (Mammal Diversity Database) data model to the Darwin Core (https://dwc.tdwg.org/) standard, and updates and harmonizes various versions of MDD dataset to DwC format. Hierarchical information and synonyms were stored in custom format in MDD versions, which are now extracted and stored in DwC style using 'parentNameUsageID' and 'acceptedNameUsageID'. The synonyms are further classified into 'junior synonym' and 'senior synonym', and are highlighted in 'taxononomic status' column. The original datasets contain individual rows only for species and genus level taxonomic rank.

To convert a single version of MDD:

    src/explore_data.py --input {infile} --output {outfile}

where infile is the input file containing the raw CSV MDD version.
E.g.

    mkdir dwc
    src/explore_data.py --input data/MDD_v1.10_6615species.csv \
       --output dwc/MDD_v1.10_6615species.csv

To convert all versions in the `data` directory:

    make
