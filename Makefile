# Time-stamp: <2017-04-26 09:28:52 (slane)>
.PHONY: all models input-data output-data clean-models clean-manuscripts clobber

all: manuscripts/censored-mle.html manuscripts/censored-mle.pdf \
	manuscripts/model-interrogation.html

.INTERMEDIATES: manuscripts/censored-mle.tex

models: stan/censored-mle-m0.rds \
	stan/censored-mle-m1.rds \
	stan/censored-mle-m2.rds \
	stan/censored-mle-m3.rds

input-data: data/biofouling.rds data/imputations.rds

output-data: data/censored-mle-m0-scaled.rds \
	data/censored-mle-m1-scaled.rds \
	data/censored-mle-m2-scaled.rds \
	data/censored-mle-m3-scaled.rds

################################################################################
# Make data for feeding into models and manuscript
data/imputations.rds data/biofouling.rds: scripts/data-cleaning.R \
	data-raw/samples.csv data-raw/vessel.csv
	cd $(<D); \
	Rscript $(<F) --no-save --no-restore

################################################################################
# Rules for making stan models
stan/censored-mle-m0.rds: scripts/compile-model.R stan/censored-mle-m0.stan
	cd $(<D); \
	Rscript $(<F) mname=$(basename $(@F) .rds)

stan/censored-mle-m1.rds: scripts/compile-model.R stan/censored-mle-m1.stan
	cd $(<D); \
	Rscript $(<F) mname=$(basename $(@F) .rds)

stan/censored-mle-m2.rds: scripts/compile-model.R stan/censored-mle-m2.stan
	cd $(<D); \
	Rscript $(<F) mname=$(basename $(@F) .rds)

stan/censored-mle-m3.rds: scripts/compile-model.R stan/censored-mle-m3.stan
	cd $(<D); \
	Rscript $(<F) mname=$(basename $(@F) .rds)

################################################################################
# Rules to fit models with data
data/censored-mle-m0-scaled.rds: scripts/fit-model.R \
	stan/censored-mle-m0.rds data/imputations.rds
	cd $(<D); \
	Rscript $(<F) mname=$(basename $(@F) .rds) myseed=737

data/censored-mle-m1-scaled.rds: scripts/fit-model.R \
	stan/censored-mle-m1.rds data/imputations.rds
	cd $(<D); \
	Rscript $(<F) mname=$(basename $(@F) .rds) myseed=666

data/censored-mle-m2-scaled.rds: scripts/fit-model.R \
	stan/censored-mle-m2.rds data/imputations.rds
	cd $(<D); \
	Rscript $(<F) mname=$(basename $(@F) .rds) myseed=42

data/censored-mle-m3-scaled.rds: scripts/fit-model.R \
	stan/censored-mle-m3.rds data/imputations.rds
	cd $(<D); \
	Rscript $(<F) mname=$(basename $(@F) .rds) myseed=13

################################################################################
# Rules to make manuscripts
manuscripts/censored-mle.html: manuscripts/censored-mle.Rmd \
	data-raw/samples.csv
	cd $(<D); \
	Rscript -e "rmarkdown::render('$(<F)')" --no-save --no-restore

manuscripts/model-interrogation.html: manuscripts/model-interrogation.Rmd
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
clean-models:
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
