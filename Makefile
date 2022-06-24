results.html: results.Rmd simulation_results.csv
	Rscript -e "rmarkdown::render('$<')"
simulation_results.csv: julia/simulation.jl data_grid.csv
	julia $<
data_grid.csv: R/r2j.R R/data_grid.R
	Rscript -e "source('$<')"