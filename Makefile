# Time-stamp: <2017-11-05 23:03:12 (overlordR)>
.PHONY: all input-data models output-VB output-MCMC \
	paper supplement \
	clean-models clean-manuscripts clobber

all: manuscripts/censored-mle.html \
	processed-data \
	manuscripts/model-interrogation.html \
	manuscripts/censored-mle.pdf

.INTERMEDIATES: manuscripts/censored-mle.tex

input-data: data/imputations.rds

models: stan/censored-mle-m0.rds \
	stan/censored-mle-m1.rds \
	stan/censored-mle-m2.rds \
	stan/censored-mle-m3.rds \
	stan/censored-mle-m0-t.rds \
	stan/censored-mle-m1-t.rds \
	stan/censored-mle-m2-t.rds \
	stan/censored-mle-m3-t.rds

output-VB: data/censored-mle-m0-var-bayes.rds \
	data/censored-mle-m1-var-bayes.rds \
	data/censored-mle-m2-var-bayes.rds \
	data/censored-mle-m3-var-bayes.rds \
	data/censored-mle-m0-t-var-bayes.rds \
	data/censored-mle-m1-t-var-bayes.rds \
	data/censored-mle-m2-t-var-bayes.rds \
	data/censored-mle-m3-t-var-bayes.rds

output-MCMC: data/censored-mle-m0.rds \
	data/censored-mle-m1.rds \
	data/censored-mle-m2.rds \
	data/censored-mle-m3.rds \
	data/censored-mle-m0-t.rds \
	data/censored-mle-m1-t.rds \
	data/censored-mle-m2-t.rds \
	data/censored-mle-m3-t.rds

# Defaults for number of multiply imputed datasets and HMC iterations if not
# passed via cmdline.
NUMMI?=50
MCITER?=4000

################################################################################
# Make data for feeding into models and manuscript
data/imputations.rds: scripts/data-cleaning.R \
	data-raw/samples.csv data-raw/vessel.csv
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) numMI=$(NUMMI)

################################################################################
# Rule for making stan models
stan/%.rds: scripts/compile-model.R stan/%.stan
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds)

################################################################################
# Rules to fit variational bayes models
data/censored-mle-%-var-bayes.rds: stan/censored-mle-%.rds \
	scripts/fit-model-vb.R data/imputations.rds
	cd scripts; \
	Rscript --no-save --no-restore fit-model-vb.R \
		mname=$(basename $(<F) .rds) myseed=737 iter=$(MCITER)

################################################################################
# Rules to fit mcmc models
data/censored-mle-%.rds: stan/censored-mle-%.rds \
	scripts/fit-model.R data/imputations.rds
	cd scripts; \
	Rscript --no-save --no-restore fit-model.R \
		mname=$(basename $(<F) .rds) myseed=737 iter=$(MCITER)

data/censored-mle-%-t.rds: stan/censored-mle-%-t.rds \
	scripts/fit-model.R data/imputations.rds
	cd scripts; \
	Rscript --no-save --no-restore fit-model.R \
		mname=$(basename $(<F) .rds) myseed=737 iter=$(MCITER)

################################################################################
# Rules to process data (add dependencies later).
data/looic-compare.rds: scripts/post-process-compare.R output-VB
	cd $(<D); \
	Rscript --no-save --no-restore $(<F)

data/looic-compare-full.rds: scripts/post-process-compare-full.R output-MCMC
	cd $(<D); \
	Rscript --no-save --no-restore $(<F)

data/diffs.rds: scripts/post-process-t.R data/censored-mle-m3-t.rds
	cd $(<D); \
	Rscript --no-save --no-restore $(<F)

################################################################################
# Rules to make manuscripts
manuscripts/censored-mle.html: manuscripts/censored-mle.Rmd \
	data-raw/samples.csv
	cd $(<D); \
	Rscript --no-save --no-restore -e "rmarkdown::render('$(<F)')"

manuscripts/model-interrogation.html: manuscripts/model-interrogation.Rmd \
	output-VB
	cd $(<D); \
	Rscript --no-save --no-restore -e "rmarkdown::render('$(<F)')"

manuscripts/censored-mle.tex: manuscripts/censored-mle.Rnw \
	data/biofouling.rds data/imputations.rds $(ROBUST-PROC-DATA)
	cd $(<D); \
	Rscript --no-save --no-restore -e "knitr::knit('$(<F)')"

manuscripts/censored-mle-supplement.tex: manuscripts/censored-mle-supplement.Rnw \
	data/biofouling.rds data/imputations.rds $(PROC-DATA)
	cd $(<D); \
	Rscript --no-save --no-restore -e "knitr::knit('$(<F)')"

%.pdf: %.tex
	cd $(<D); \
	latexmk -pdf $(<F)

# phony rule to make paper from included figures
paper: manuscripts/censored-mle.Rnw data/biofouling.rds
	cd $(<D); \
	Rscript --no-save --no-restore -e "knitr::knit('$(<F)')"; \
	latexmk -pdf $(<F:Rnw=tex)

supplement: manuscripts/censored-mle-supplement.Rnw data/biofouling.rds
	cd $(<D); \
	Rscript --no-save --no-restore -e "knitr::knit('$(<F)')"; \
	latexmk -pdf $(<F:Rnw=tex)

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
