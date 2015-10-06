
close

include("unscented.jl")
using Unscented


network_on = true
make_plots = true



alpha = 0.5
dVer = 2
meas_var = 16
walk_var = 4
U = 10



# xstart = 108
# ystart = 10
xstart = rand(1,U)*10 + 20
ystart = rand(1,U) + 5

C = -[30 29 30 29 33 30 30 29 33 30 33 30 30 29 33 30 33 30 30 29 30 29 33 30 33 30 30 29 33 30 33 30 30 29 30 29 30 29 30 29]
C_bt = -30

T = 10; #time steps
N = 300; #particles

# AP:s
beaconPositions = [   58.0000   80.0000   99.0000   72.0000  106.5000 109.5000   84.0000   95.5000  108.0000  129.0000 0   19.0000   35.5000    50.5000   14.0000   22.0000 51.0000   32.0000   47.5000   69.0000 ;
                      41.5000   40.0000   40.0000   26.0000   32.0000 21.0000   10.0000   10.0000   10.0000   10.0000 30.5000   30.5000   30.5000   27.5000   21.0000 9.5000    9.5000         0         0       0]

nbr_aps = size(beaconPositions, 2)

#
w = zeros(1,N)

x = zeros(2, N, U) #Particle states
x = x + sqrt(walk_var) * randn(size(x))



for u in 1:U
    x[1,:,u] = x[1,:,u] + xstart[u]
    x[2,:,u] = x[2,:,u] + ystart[u]
end


z = zeros(2,nbr_aps)
x_next = zeros(2,N)
truePos = zeros(2, T, U)

for u = 1:U
    truePos[1,1,u] = xstart[u]
    truePos[2,1,u] = ystart[u]
end

for u = 1:U
    for t = 2:T
        truePos[:,t,u] = truePos[:,t-1,u] + sqrt(walk_var)*randn(2,1)
    end
end

x_est = zeros(2, T,U)
x_est[:,1,:] = truePos[:,1,:]
d_u = zeros(U)
z_ou = zeros(1,U)


#
tic()
for t = 2:T
    
    for u = 1:U

        d = sqrt(dVer^2 + (beaconPositions[1,:]-truePos[1,t,u]).^2 + (beaconPositions[2,:]-truePos[2,t,u]).^2)

        # Simulate data:
        for i = 1:nbr_aps
            z[1,i] = C[2*i-1] - 20*log10(d[i]) - alpha*d[i] + sqrt(meas_var)*randn()
            z[2,i] = C[2*i] - 20*log10(d[i]) - alpha*d[i] + sqrt(meas_var)*randn()
        end
        ########################
	########################
        if network_on
            
            
            for ou = 1:U
                d_u[ou] = sqrt( (truePos[1,t-1,ou]-truePos[1,t-1,u]) ^2  + (truePos[2,t-1,ou]-truePos[2,t-1,u]) ^2)
		#d_u[ou] = sqrt( (truePos[1,t,ou]-truePos[1,t,u]) ^2  + (truePos[2,t,ou]-truePos[2,t,u]) ^2)
		z_ou[ou] = C_bt - 20*log10(d_u[ou]) - alpha*d_u[ou] + sqrt(meas_var)*randn()
            end
           
            

        end
	##################################
	##################################
            

        for i = 1:N

            # Weighting:
            lnw = 0;
            d_p = sqrt(dVer^2 + (beaconPositions[1,:]-x[1,i,u]).^2 + (beaconPositions[2,:]-x[2,i,u]).^2)

            for k = 1:nbr_aps

                # 5 GHz:
                pObs = C[2*k-1] - 20*log10(d_p[k]) - alpha*d_p[k]       
                if z[1,k] < -80
                    z[1,k] = pObs
                end
                lnw = lnw -(z[1,k] - pObs)^2/(2*meas_var)


                # 2.4 GHz:
                pObs = C[2*k] - 20*log10(d_p[k]) - alpha*d_p[k]    
                if z[2,k] < -80
                    z[2,k] = pObs; 
                end
                lnw = lnw -(z[2,k] - pObs)^2/(2*meas_var)
            end

            ###################
	#########################
            if network_on
                # Between users:
                d_u_p = [];
		if u == 1
		        for ou = 1:U
		            if ou!=u
				if d_u[ou] < 15
					pos_var = cov(x[:,:,ou],vardim=2)
					pos_mu = [ mean(x[1,:,ou]), mean(x[2,:,ou]) ]
					receiver_pos = [x[1,i,u], x[2,i,u]]

					bt_rri_var, bt_rri_mu = unscented(pos_var, pos_mu, receiver_pos, C_bt, alpha)
					#bt_rri_var, bt_rri_mu = unscented(pos_var, truePos[:,t,ou], receiver_pos, C_bt, alpha)

		               	 	#d_u_p = [d_u_p ; sqrt( (truePos[1,t-1,ou]-x[1,i,u]) ^2 + (truePos[2,t-1,ou]-x[2,i,u]) ^2)]       
			
		               		#pObs = C_bt - 20*log10(d_u_p[ou]) - alpha*d_u_p[ou]
		                	lnw = lnw -(z_ou[ou] - bt_rri_mu)^2/(2*(bt_rri_var + meas_var))
		           	end          
		            end

		        end
		end
                #for ou = 1:U-1
                #    if d_u[ou] < 15
                 #       pObs = C_bt - 20*log10(d_u_p[ou]) - alpha*d_u_p[ou]
                 #       lnw = lnw -(z_ou[ou] - pObs)^2/(2*meas_var)
                #   end
               # end
            end
           #####################
	#######################
            w[i] = lnw

        end

        w = w - minimum(w)
        w = exp(w)
        w = w/sum(w)
        # Resample:
        for i = 1 : N

            temp = find(rand() .<= cumsum(w,2))
            temp = temp[1]
            x_next[:,i] = x[:,temp,u]
        end
        x_est[:,t,u] = mean(x_next,2)
        x[:,:,u] = x_next
        

        # State update:
        x[:,:,u] = x[:,:,u] + sqrt(walk_var)*randn(size(x[:,:,u]))
    end
    
    print("t = ")
    println(t)
end
toc()

meanError = mean(mean(mean(abs(truePos-x_est))))
maxError = maximum(maximum(maximum(abs(truePos-x_est))))
print("Average error = ")
println(meanError)
print("Maximum error = ")
println(maxError)
print("First user mean = ")
println(mean(mean(abs(truePos[:,:,1]-x_est[:,:,1]))))


#using Winston
#using PyPlot
if make_plots


	figure()
	plot(beaconPositions[1,:],beaconPositions[2,:],"c^")
	hold(true)
	for u = 1:U
	    plot(truePos[1,:,u], truePos[2,:,u], "b-")
	    plot(x_est[1,:,u], x_est[2,:,u], "r-")
	end
	xlim(-10, 150)
	ylim(-10, 60)
	savefig("figure1.eps")

	figure()
	plot(beaconPositions[1,:],beaconPositions[2,:],"g^")
	hold(true)
	plot(truePos[1,:,1], truePos[2,:,1],"b-")
	plot(x_est[1,:,1], x_est[2,:,1], "r-")
	xlim(-10, 150)
	ylim(-10, 60)
	savefig("figure2.eps")
end




















