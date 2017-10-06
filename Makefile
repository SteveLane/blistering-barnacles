# Time-stamp: <2017-10-06 15:23:47 (slane)>
.PHONY: all input-data models output-data \
	ROBUST-PROC-DATA PROC-DATA robust-processed-data processed-data \
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
	stan/censored-mle-m4.rds \
	stan/censored-mle-m0-robust.rds \
	stan/censored-mle-m1-robust.rds \
	stan/censored-mle-m2-robust.rds \
	stan/censored-mle-m3-robust.rds \
	stan/censored-mle-m4-robust.rds

output-data: data/censored-mle-m0.rds \
	data/censored-mle-m1.rds \
	data/censored-mle-m2.rds \
	data/censored-mle-m3.rds \
	data/censored-mle-m4.rds \
	data/censored-mle-m0-robust.rds \
	data/censored-mle-m1-robust.rds \
	data/censored-mle-m2-robust.rds \
	data/censored-mle-m3-robust.rds \
	data/censored-mle-m4-robust.rds

ROBUST-PROC-DATA = graphics/obs-hist.pdf \
	graphics/imp-days1.pdf \
	graphics/imp-trips.pdf \
	graphics/imp-paint.pdf \
	graphics/plM1boat-robust.pdf \
	graphics/plM1paint-robust.pdf \
	graphics/plM3boat-robust.pdf \
	graphics/plM3paint-robust.pdf \
	graphics/plM4Type-robust.pdf \
	graphics/plSummary-robust.pdf \
	data/looic-robust.rds \
	data/diffs.rds

PROC-DATA = graphics/plM1boat.pdf \
	graphics/plM1paint.pdf \
	graphics/plM3boat.pdf \
	graphics/plM3paint.pdf \
	graphics/plM4Type.pdf \
	graphics/plSummary.pdf \
	data/looic.rds

robust-processed-data: $(ROBUST-PROC-DATA)

processed-data: $(PROC-DATA)

# Defaults for number of multiply imputed datasets and HMC iterations if not
# passed via cmdline.
NUMMI?=50
MCITER?=2000

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
# Rules to fit models with data
data/censored-mle-%.rds: scripts/fit-model.R stan/censored-mle-%.rds \
	data/imputations.rds
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds) \
		myseed=737 iter=$(MCITER)

################################################################################
# Rules to process data (add dependencies later).
$(ROBUST-PROC-DATA): scripts/post-process-robust.R robust-output-data input-data
	cd $(<D); \
	Rscript --no-save --no-restore $(<F)

$(PROC-DATA): scripts/post-process.R output-data
	cd $(<D); \
	Rscript --no-save --no-restore $(<F)

################################################################################
# Rules to make manuscripts
manuscripts/censored-mle.html: manuscripts/censored-mle.Rmd \
	data-raw/samples.csv
	cd $(<D); \
	Rscript --no-save --no-restore -e "rmarkdown::render('$(<F)')"

manuscripts/model-interrogation.html: manuscripts/model-interrogation.Rmd \
	robust-output-data
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
