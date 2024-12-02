.PHONY: data install tables template calibrate
CALIBRATION = calibrate_params.jl calibration_utils.jl experiments/config.jl data/impvol_data.jld2
EQULIBRIUM = utils.jl equilibrium.jl experiments/config.jl
COLUMNS = actual kappa1972 nosectoral nosectoral_kappa1972
CES = CES0.5 CES1.5 
TABLES = $(CES) baseline china_1972 no_china no_io_linkages labor_adjustment trade_imbalance theta2 theta8 rho002 rho0005
.PRECIOUS: $(foreach table,$(TABLES),$(foreach column,$(COLUMNS),experiments/$(table)/$(column)/results.jld2))

# default number of Julia threads to use. otherwise `make tables PROCS=12`
PROCS = 8
JULIA = julia --project -p$(PROCS)

tables: $(foreach table,1 2 3 4left 4right 5left 5center 5right 6left 6right 7 8left 8right,output/table$(table).csv) 
ces_tables: $(foreach table,$(CES),experiments/$(table)/output_table.csv) experiments/baseline/output_table.csv 

# this takes too long to run, only run if explicitly asked `make S500`
S500: experiments/S500/output_table.csv

calibrate: $(foreach table,$(TABLES),experiments/$(table)/common_parameters.jld2) 

admissible_eos: $(wildcard experiments/CES/*/results.jld2)
experiments/CES/%/results.jld2: experiments/CES/%/common_parameters.jld2 experiments/CES/scenario.jl
	cd experiments/CES && $(JULIA) scenario.jl $(subst experiments/CES/,,$<) 

experiments/CES/2.0/common_parameters.jld2: experiments/CES/init_parameters.jl $(CALIBRATION) 
	cd experiments/CES && $(JULIA) init_parameters.jl

experiments/%/common_parameters.jld2: experiments/%/init_parameters.jl $(CALIBRATION) 
	cd $(dir $@) && $(JULIA) init_parameters.jl

define run_experiment
experiments/$(1)/%/results.jld2: $(EQULIBRIUM) experiments/$(1)/common_parameters.jld2 experiments/$(1)/%/scenario.jl experiments/$(1)/%/change_parameters.jl 
	@echo " + Compiling '$$@'"
	cd $$(dir $$@) && $(JULIA) scenario.jl > errors.log 2>&1
endef

$(foreach experiment,$(TABLES) S500,$(eval $(call run_experiment,$(experiment))))

experiments/%/output_table.csv: $(foreach column,$(COLUMNS),experiments/%/$(column)/results.jld2) output.jl table.jl
	$(JULIA) table.jl $(dir $@)

data: data/impvol_data.jld2
data/impvol_data.jld2: read_data.jl data/*.csv data/*.txt
	$(JULIA) read_data.jl

template: scenario_template.jl
	find . -name "scenario.jl" -exec cp scenario_template.jl {} \; 

# install the Julia package dependencies
install: install.jl
	$(JULIA) install.jl

# copy tables to match the order in the paper
output/table1.csv: experiments/baseline/output_table.csv
	cp $< $@
	cp $@ output/1499.csv
output/table2.csv: output/volatility_by_decade.csv
	cp $< $@
output/table3.csv: experiments/trade_imbalance/output_table.csv
	cp $< $@
output/table4left.csv: experiments/theta2/output_table.csv
	cp $< $@
output/table4right.csv: experiments/theta8/output_table.csv
	cp $< $@
output/table5left.csv: experiments/rho0005/output_table.csv
	cp $< $@
output/table5center.csv: experiments/labor_adjustment/output_table.csv
	cp $< $@
output/table5right.csv: experiments/rho002/output_table.csv
	cp $< $@
output/table6left.csv: experiments/CES0.5/output_table.csv
	cp $< $@
output/table6right.csv: experiments/CES1.5/output_table.csv
	cp $< $@
output/table7.csv: experiments/no_io_linkages/output_table.csv
	cp $< $@
output/table8left.csv: experiments/no_china/output_table.csv
	cp $< $@
output/table8right.csv: experiments/china_1972/output_table.csv
	cp $< $@
