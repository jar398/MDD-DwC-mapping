
all: dwc/MDD_v1_6495species_JMamm.csv \
     dwc/MDD_v1.1_6526species.csv \
     dwc/MDD_v1.2_6485species.csv
     dwc/MDD_v1.3_6513species.csv
     dwc/MDD_v1.31_6513species.csv
     dwc/MDD_v1.4_6533species.csv
     dwc/MDD_v1.5_6554species.csv
     dwc/MDD_v1.6_6557species.csv
     dwc/MDD_v1.7_6567species.csv
     dwc/MDD_v1.8_6591species.csv
     dwc/MDD_v1.9_6596species.csv
     dwc/MDD_v1.9.1_6596species.csv
     dwc/MDD_v1.10_6615species.csv

dwc/%.csv: data/%.csv
	mkdir -p dwc
	src/explore_data.py --input $< --output $@.new
	mv -f $@.new $@
