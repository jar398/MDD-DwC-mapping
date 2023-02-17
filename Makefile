
all: dwc/MDD_v1_6495species_JMamm_inDwC.csv \
     dwc/MDD_v1.1_6526species_inDwC.csv \
     dwc/MDD_v1.2_6485species_inDwC.csv \
     dwc/MDD_v1.3_6513species_inDwC.csv \
     dwc/MDD_v1.31_6513species_inDwC.csv \
     dwc/MDD_v1.4_6533species_inDwC.csv \
     dwc/MDD_v1.5_6554species_inDwC.csv \
     dwc/MDD_v1.6_6557species_inDwC.csv \
     dwc/MDD_v1.7_6567species_inDwC.csv \
     dwc/MDD_v1.8_6591species_inDwC.csv \
     dwc/MDD_v1.9_6596species_inDwC.csv \
     dwc/MDD_v1.9.1_6596species_inDwC.csv \
     dwc/MDD_v1.10_6615species_inDwC.csv

dwc/%_inDwC.csv: data/%.csv
	mkdir -p dwc
	src/explore_data.py --input $< --output $@.new
	mv -f $@.new $@
