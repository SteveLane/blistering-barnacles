# Time-stamp: <2017-04-13 15:44:09 (slane)>
.PHONY: all models input-data output-data clean-data clean-manuscripts clobber

all: manuscripts/censored-mle.html manuscripts/censored-mle.pdf

.INTERMEDIATES: manuscripts/censored-mle.tex

# models: stan/dynamic-governance-m0.rds

input-data: data/biofouling.rds data/imputations.rds

# output-data: data/stanfit-KIWI-dynamic-governance-m0-10.rda

################################################################################
# Rules for making stan models
# stan/dynamic-governance-m0.rds: R/compile-model.R
# 	cd $(<D); \
# 	Rscript $(<F) mname=dynamic-governance-m0

################################################################################
# Make data for feeding into models and manuscript
data/imputations.rds data/biofouling.rds: scripts/data-cleaning.R \
	data-raw/samples.csv data-raw/vessel.csv
	cd $(<D); \
	Rscript $(<F) --no-save --no-restore

################################################################################
# Rules to fit models with data
# data/stanfit-KIWI-dynamic-governance-m0-10.rda: R/fit-model.R \
# 	stan/dynamic-governance-m0.rds data/data-KIWI-10.rda
# 	cd $(<D); \
# 	Rscript $(<F) mname=dynamic-governance-m0 pathway=KIWI size=10 iter=2000

################################################################################
# Rules to make manuscripts
%.html: %.Rmd data-raw/samples.csv
	cd $(<D); \
	Rscript -e "rmarkdown::render('$(<F)')" --no-save --no-restore

%.tex: %.Rnw data/biofouling.rds
	cd $(<D); \
	Rscript -e "knitr::knit('$(<F)')" --no-save --no-restore

%.pdf: %.tex
	cd $(<D); \
	latexmk -pdf $(<F)

################################################################################
# Cleaning targets
clean-data:
	cd stan/; \
	rm -f *.rds

clean-manuscripts:
	cd manuscripts/; \
	rm -rf *.aux *.bbl *.bcf *.blg *.fdb_latexmk *.fls *.lof *.log *.lot \
		*.code *.loe *.toc *.rec *.out *.run.xml *~ *.prv \
		censored-mle.tex _region_*

clobber: clean-data clean-manuscripts
	cd manuscripts/; \
	rm -rf auto/ cache/ figure/
