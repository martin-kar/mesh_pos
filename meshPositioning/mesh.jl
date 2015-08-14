const network_on = true
const make_plots = false
const alpha = 0.5
const dVer = 2
const meas_var = 16
const walk_var = 4
const U = 100
# xstart = 108
# ystart = 10
xstart = rand(1,U)*100 + 20
ystart = rand(1,U)*40 + 5

C = -[30 29 30 29 33 30 30 29 33 30 33 30 30 29 33 30 33 30 30 29 30 29 33 30 33 30 30 29 33 30 33 30 30 29 30 29 30 29 30 29]
C_2 = -[30,30,33,30,33,33,30,33,33,30,30,33,33,30,33,33,30,30,30,30]
C_2 = reshape(C_2, 1, length(C_2))
C_5 = -[29,29,30,29,30,30,29,30,30,29,29,30,30,29,30,30,29,29,29,29]
C_5 = reshape(C_5, 1, length(C_5))

const C_bt = -30

const T = 10; #time steps
const N = 300; #particles

# AP:s
const beaconPositions = [   58.0000   80.0000   99.0000   72.0000  106.5000 109.5000   84.0000   95.5000  108.0000  129.0000 0   19.0000   35.5000    50.5000   14.0000   22.0000 51.0000   32.0000   47.5000   69.0000 ;
                      41.5000   40.0000   40.0000   26.0000   32.0000 21.0000   10.0000   10.0000   10.0000   10.0000 30.5000   30.5000   30.5000   27.5000   21.0000 9.5000    9.5000         0         0       0]

const nbr_aps = size(beaconPositions, 2)

#
w = zeros(1,N)

x = zeros(2, N, U); #Particle states
x = x + sqrt(walk_var) * randn(size(x))

x[1,:,:] = ones(N,1)*xstart
x[2,:,:] = ones(N,1)*ystart

z_2 = zeros(nbr_aps,U)
z_5 = zeros(nbr_aps,U)
x_next = zeros(2,N)
truePos = zeros(2, T, U)

randWalks = zeros(2,T,U)
randWalks[:,:,1:convert(Int64,U/2)] = cumsum(sqrt(walk_var)*randn(2,T,convert(Int64,U/2)),2)
truePos[1,:,:] = ones(T,1)*xstart
truePos[2,:,:] = ones(T,1)*ystart
truePos = truePos + randWalks

x_est = zeros(2, T,U)
x_est[:,1,:] = truePos[:,1,:]
z_ou = zeros(U,U)

const ones_U = ones(U,1)
const ones_UU = ones(U,U)
const ones_aps = ones(nbr_aps, 1)
const ones_N = ones(N,1);

const C_2_U = C_2'*ones_U'
const C_5_U = C_5'*ones_U'
const C_2_N = C_2'*ones_N'
const C_5_N = C_5'*ones_N'
const C_bt_UU = C_bt*ones_UU
const C_bt_NU = C_bt*ones(N,U);

#
tic()
for t = 2:T
    
    xDistSqr = (ones_U*beaconPositions[1,:] - (ones_aps*reshape(truePos[1,t,:],1,U))').^2
    yDistSqr = (ones_U*beaconPositions[2,:] - (ones_aps*reshape(truePos[2,t,:],1,U))').^2
    d = sqrt(dVer^2 + xDistSqr + yDistSqr)'
    
    z_2[:,:] = C_2_U - 20*log10(d) - alpha*d + sqrt(meas_var)*randn(nbr_aps,U)
    z_5[:,:] = C_5_U - 20*log10(d) - alpha*d + sqrt(meas_var)*randn(nbr_aps,U)
    if network_on
        d_x = ones_U*reshape(truePos[1,t,:],1,U)
        d_y = ones_U*reshape(truePos[2,t,:],1,U)
        d_u = sqrt((d_x - d_x').^2 + (d_y - d_y').^2)
        z_ou = C_bt_UU - 20*log10(d_u) - alpha*d_u + sqrt(meas_var)*randn(U,U)
    end
        
    for u = 1:U

        # Weighting:
            
        d_x2 = ((ones_N*beaconPositions[1,:])' - ones_aps*x[1,:,u]).^2
        d_y2 = ((ones_N*beaconPositions[2,:])' - ones_aps*x[2,:,u]).^2
        d_p = sqrt(dVer^2 + d_x2 + d_y2)
            
	# Should not have Gaussian noise term here.
        pObs_2 = C_2_N - 20*log10(d_p) - alpha*d_p
        pObs_5 = C_2_N - 20*log10(d_p) - alpha*d_p
            
        lnw = - (z_2[:,u]*ones_N' - pObs_2).^2/(2*meas_var) - (z_5[:,u]*ones_N' - pObs_5).^2/(2*meas_var)

        if network_on
            # Between users:
            d_x2 = (ones_N*squeeze(truePos[1,t-1,:],1) - x[1,:,u]'*ones(1,U)).^2
            d_y2 = (ones_N*squeeze(truePos[2,t-1,:],1) - x[2,:,u]'*ones(1,U)).^2
            d_u_p = sqrt(d_x2 + d_y2)
            
            pObs = C_bt_NU - 20*log10(d_u_p) - alpha*d_u_p
            luw = -(ones_N*z_ou[u,:] - pObs).^2/(2*meas_var)
            luw[:,u] = 0.0
            w = sum(lnw,1) + sum(luw,2)'
	else
	    w = sum(lnw,1)
        end
        w = w - minimum(w)   
        w = exp(w)
        w = w./sum(w)
        csum = squeeze(cumsum(w,2),1)
	# Resampling:
        for i = 1:N
            temp = searchsortedfirst(csum, rand())
            x_next[:,i] = x[:,temp,u]
        end
        x_est[:,t,u] = mean(x_next,2)
        x[:,:,u] = x_next
        
        # State update:
        x[:,:,u] = x[:,:,u] + sqrt(walk_var) * randn(size(x[:,:,u]))
    end
end
toc()

meanError = mean(mean(mean(abs(truePos-x_est))))
maxError = maximum(maximum(maximum(abs(truePos-x_est))))
print("Average error = ")
println(meanError)
print("Maximum error = ")
println(maxError)

#print("x_est = ")
#println(x_est)

#print("truePos = ")
#println(truePos)


if make_plots
	using PyPlot
	#using Winston

	PyPlot.figure()
	PyPlot.plt.plot(beaconPositions[1,:],beaconPositions[2,:],"c^")
	hold(true)
	for u = 1:U
	    if u < U/2
		PyPlot.plt.plot(truePos[1,:,u], truePos[2,:,u], "b-")
	    else
		PyPlot.plt.plot(truePos[1,:,u], truePos[2,:,u], "b+")
	    end
	    PyPlot.plt.plot(x_est[1,:,u], x_est[2,:,u], "r-")
	end
	xlim(-10, 150)
	ylim(-10, 60)
	savefig("figure1.eps")

	PyPlot.figure()
	PyPlot.plt.plot(beaconPositions[1,:],beaconPositions[2,:],"g^")
	hold(true)
	PyPlot.plt.plot(truePos[1,:,1], truePos[2,:,1],"b-")
	PyPlot.plt.plot(x_est[1,:,1], x_est[2,:,1], "r-")
	xlim(-10, 150)
	ylim(-10, 60)
	savefig("figure2.eps")
end
