# You can obtain the original data/ files from Zenodo (but they may be
# cached in the data/ directory of this git repo).
# dwc/ file names are chosen to work well with listtools.

HOME=.
MDDSOURCE=data
DWC=dwc

all: $(DWC)/mdd1.0-dwc.csv \
     $(DWC)/mdd1.1-dwc.csv \
     $(DWC)/mdd1.2-dwc.csv \
     $(DWC)/mdd1.3-dwc.csv \
     $(DWC)/mdd1.31-dwc.csv \
     $(DWC)/mdd1.4-dwc.csv \
     $(DWC)/mdd1.5-dwc.csv \
     $(DWC)/mdd1.6-dwc.csv \
     $(DWC)/mdd1.7-dwc.csv \
     $(DWC)/mdd1.8-dwc.csv \
     $(DWC)/mdd1.9-dwc.csv \
     $(DWC)/mdd1.10-dwc.csv

CONVERTMDD=mkdir -p $(DWC) && time python3 $(HOME)/src/explore_data.py
$(DWC)/mdd1.0-dwc.csv: $(MDDSOURCE)/MDD_v1_6495species_JMamm.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
$(DWC)/mdd1.1-dwc.csv: $(MDDSOURCE)/MDD_v1.1_6526species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
$(DWC)/mdd1.2-dwc.csv: $(MDDSOURCE)/MDD_v1.2_6485species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
$(DWC)/mdd1.3-dwc.csv: $(MDDSOURCE)/MDD_v1.3_6513species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
$(DWC)/mdd1.31-dwc.csv: $(MDDSOURCE)/MDD_v1.31_6513species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
$(DWC)/mdd1.4-dwc.csv: $(MDDSOURCE)/MDD_v1.4_6533species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
$(DWC)/mdd1.5-dwc.csv: $(MDDSOURCE)/MDD_v1.5_6554species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
$(DWC)/mdd1.6-dwc.csv: $(MDDSOURCE)/MDD_v1.6_6557species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
$(DWC)/mdd1.7-dwc.csv: $(MDDSOURCE)/MDD_v1.7_6567species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
$(DWC)/mdd1.8-dwc.csv: $(MDDSOURCE)/MDD_v1.8_6591species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
$(DWC)/mdd1.9-dwc.csv: $(MDDSOURCE)/MDD_v1.9_6596species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
# 1.9.1 is not on zenodo (I think).  Probably no reason to worry.
$(DWC)/mdd1.9.1-dwc.csv:
$(DWC)/mdd1.10-dwc.csv: $(MDDSOURCE)/MDD_v1.10_6615species.csv
	$(CONVERTMDD) --input $< --output $@.new
	mv -f $@.new $@
