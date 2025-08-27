.PHONY: data install tables template calibrate download clean transform pipeline
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

# Data pipeline targets
pipeline: download clean transform

download: download-pwt download-unna download-wdi download-comtrade

download-pwt: input/pwt/pwt71.zip input/pwt/pwt56.zip

input/pwt/pwt71.zip input/pwt/pwt56.zip:
	julia --project code/create/download/pwt.jl

download-unna: temp/unna_download.done

temp/unna_download.done:
	julia --project code/create/download/unna.jl
	@mkdir -p temp && touch $@

download-wdi: input/wdi/WDI_csv_2015_10.zip

input/wdi/WDI_csv_2015_10.zip:
	julia --project code/create/download/wdi.jl

download-comtrade: temp/comtrade_download.done

temp/comtrade_download.done:
	julia --project code/create/download/comtrade.jl
	@mkdir -p temp && touch $@

clean: clean-pwt clean-unna clean-wdi clean-comtrade

clean-pwt: temp/pwt71_clean.csv temp/pwt56_former_clean.csv

temp/pwt71_clean.csv temp/pwt56_former_clean.csv: input/pwt/pwt71.zip input/pwt/pwt56.zip
	julia --project code/create/clean/pwt.jl

clean-unna: temp/unna_clean.csv

temp/unna_clean.csv: temp/unna_download.done
	julia --project code/create/clean/unna.jl

clean-wdi: temp/wdi_clean.csv

temp/wdi_clean.csv: input/wdi/WDI_csv_2015_10.zip
	julia --project code/create/clean/wdi.jl

clean-comtrade: temp/comtrade_clean.csv

temp/comtrade_clean.csv: temp/comtrade_download.done
	julia --project code/create/clean/comtrade.jl

transform: output/gross_output.csv output/value_added.csv output/trade_flows.csv output/price_indices.csv

output/gross_output.csv output/value_added.csv: temp/pwt71_clean.csv temp/unna_clean.csv temp/wdi_clean.csv
	julia --project code/create/transform/assemble_output_va.jl

output/trade_flows.csv: temp/comtrade_clean.csv
	julia --project code/create/transform/trade_flows.jl

output/price_indices.csv: temp/pwt71_clean.csv temp/wdi_clean.csv
	julia --project code/create/transform/price_indices.jl

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
