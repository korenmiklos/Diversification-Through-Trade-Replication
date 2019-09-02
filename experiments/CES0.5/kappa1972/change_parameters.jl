# parameters that govern counterfactual

## kappa remains at 1972 level
for t=1:parameters[:T]
	parameters[:kappa_mnjt][:,:,:,t] = parameters[:kappa_mnjt][:,:,:,1]
end
