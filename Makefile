.PHONY: data install tables template calibrate
CALIBRATION = calibrate_params.jl calibration_utils.jl experiments/config.jl data/impvol_data.jld2
EQULIBRIUM = utils.jl equilibrium.jl experiments/config.jl
COLUMNS = actual kappa1972 nosectoral nosectoral_kappa1972
CES = CES0.5 CES1.5 
TABLES = $(CES) baseline china_1972 no_china no_io_linkages labor_adjustment trade_imbalance theta2 theta8 rho002 rho0005
.PRECIOUS: $(foreach table,$(TABLES),$(foreach column,$(COLUMNS),experiments/$(table)/$(column)/results.jld2))

# default number of Julia threads to use. otherwise `make tables PROCS=12`
PROCS = 2

tables: $(foreach table,$(TABLES),experiments/$(table)/output_table.csv) 
ces_tables: $(foreach table,$(CES),experiments/$(table)/output_table.csv) experiments/baseline/output_table.csv 

# this takes too long to run, only run if explicitly asked `make S500`
S500: experiments/S500/output_table.csv

calibrate: $(foreach table,$(TABLES),experiments/$(table)/common_parameters.jld2) 

admissible_eos: $(wildcard experiments/CES/*/results.jld2)
experiments/CES/%/results.jld2: experiments/CES/%/common_parameters.jld2 experiments/CES/scenario.jl
	cd experiments/CES && julia scenario.jl $(subst experiments/CES/,,$<) 

experiments/CES/2.0/common_parameters.jld2: experiments/CES/init_parameters.jl $(CALIBRATION) 
	cd experiments/CES && julia init_parameters.jl

experiments/%/common_parameters.jld2: experiments/%/init_parameters.jl $(CALIBRATION) 
	cd $(dir $@) && julia init_parameters.jl

define run_experiment
experiments/$(1)/%/results.jld2: $(EQULIBRIUM) experiments/$(1)/common_parameters.jld2 experiments/$(1)/%/scenario.jl experiments/$(1)/%/change_parameters.jl 
	@echo " + Compiling '$$@'"
	cd $$(dir $$@) && julia -p$(PROCS) scenario.jl > errors.log 2>&1
endef

$(foreach experiment,$(TABLES) S500,$(eval $(call run_experiment,$(experiment))))

experiments/%/output_table.csv: $(foreach column,$(COLUMNS),experiments/%/$(column)/results.jld2) output.jl table.jl
	julia table.jl $(dir $@)

data: data/impvol_data.jld2
data/impvol_data.jld2: read_data.jl data/*.csv data/*.txt
	julia read_data.jl

template: scenario_template.jl
	find . -name "scenario.jl" -exec cp scenario_template.jl {} \; 

# install the Julia package dependencies
install: install.jl
	julia install.jl
