RSCRIPT?=Rscript --vanilla
BCFTOOLS?=bcftools
BGZIP?=bgzip
SOFOS?=../../sofos

default: all

all: sim01

.PHONY: default all

sim01: sim01.vcf.gz sim01-samples1.sofos.csv sim01-samples2.sofos.csv \
       sim01-sample1-ac.txt sim01-sample2-ac.txt

.PHONY: sim01

# SIM01
sim01.vcf.gz:
	$(RSCRIPT) ../vcfsim.R -a 0.1 -b 1 -m 10000 -n 10 -z 514803703 | $(BGZIP) > $@

sim01-samples1.txt: sim01.vcf.gz
	$(BCFTOOLS) query -l $< | head -n 8 > $@

sim01-samples2.txt: sim01.vcf.gz
	$(BCFTOOLS) query -l $< | tail -n 2 > $@

sim01-samples1.vcf: sim01-samples1.txt sim01.vcf.gz 
	$(BCFTOOLS) view -S $^ | $(BCFTOOLS) +fill-tags -- -t AC,AN > $@

sim01-samples2.vcf: sim01-samples2.txt sim01.vcf.gz 
	$(BCFTOOLS) view -S $^ | $(BCFTOOLS) +fill-tags -- -t AC,AN > $@

sim01-samples1.sofos.csv: sim01-samples1.vcf
	$(SOFOS) -u -r -n 4 -a 0.1 -b 1 $< > $@

sim01-samples2.sofos.csv: sim01-samples2.vcf
	$(SOFOS) -u -r -n 4 -a 0.1 -b 1 $< > $@

sim01-sample1-ac.txt: sim01-samples1.vcf
	$(BCFTOOLS) query -f '%AC\n' $< > $@

sim01-sample2-ac.txt: sim01-samples2.vcf
	$(BCFTOOLS) query -f '%AC\n' $< > $@
