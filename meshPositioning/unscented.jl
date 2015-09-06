
module Unscented
function unscented(pos_var, pos_mu, receiver_pos, C_bt, alpha)

	n = 2
	
	
	chi = zeros(n,2n)
	weights = zeros(1,2n)
	
	kappa = 3-n
	
	
	
	
	chi0 = pos_mu
	weight0 = kappa/(n + kappa)
	for i = 1:n
		chi[:,i] = pos_mu + sqrtm((n+kappa)*pos_var)[:,i]
		weights[i] = 1/(2*(n + kappa))
		chi[:,n+i] = pos_mu - sqrtm((n+kappa)*pos_var)[:,i]
		weights[n+i] = 1/(2*(n + kappa))
	end
	

	d0 = sqrt((chi0[1] - receiver_pos[1])^2 + (chi0[2] - receiver_pos[2])^2 )
	Y0 = C_bt - 20*log10(d0) - alpha*d0
	Y = zeros(1,2n)
	for i = 1:2n
		di = sqrt((chi[1,i] - receiver_pos[1])^2 + (chi[2,i] - receiver_pos[2])^2 )
		Y[i] = C_bt -20*log10(di) -alpha*di 	
	end
	rri_mu = weight0*Y0 + sum(Y.*weights)
	rri_var = weight0*(Y0-rri_mu)^2 + sum(weights.*(Y - rri_mu).^2) 


	return rri_var, rri_mu

end
export unscented

#pos_var = [1 0; 0 1]
#pos_mu = [1, 1]
#receiver_pos = [2, 2]
#C_bt = -30
#alpha = 0.5


#rri_var, rri_mu  = unscented(pos_var, pos_mu, receiver_pos, C_bt, alpha)

#println(rri_var)
#println(rri_mu)
end

