# parameters that govern counterfactual

ImpvolEquilibrium.remove_shock!(parameters, :global_sectoral_shock_njs)
ImpvolEquilibrium.remove_shock!(parameters, :idiosyncratic_shock_njs)

## kappa remains at 1972 level
for t=1:parameters[:T]
	parameters[:kappa_mnjt][:,:,:,t] = parameters[:kappa_mnjt][:,:,:,1]
end
