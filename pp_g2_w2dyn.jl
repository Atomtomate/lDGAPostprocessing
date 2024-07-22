 using HDF5

function read_worm_comp(file_path)
	g2worm = Array{Array{ComplexF64, 3},1}(undef,6)
	f = h5open(file_path,"r")
	g2worm[1] = read(f["worm-last/ineq-001/g4iw-worm/00001/value"])
	g2worm[2] = read(f["worm-last/ineq-001/g4iw-worm/00004/value"])
	g2worm[3] = read(f["worm-last/ineq-001/g4iw-worm/00007/value"])
	g2worm[4] = read(f["worm-last/ineq-001/g4iw-worm/00010/value"])
	g2worm[5] = read(f["worm-last/ineq-001/g4iw-worm/00013/value"])
	g2worm[6] = read(f["worm-last/ineq-001/g4iw-worm/00016/value"])	
	return g2worm
end

function read_1_part_gf(f_dmft)
    f = h5open(f_dmft,"r")
    g_raw =  read(f["dmft-last/ineq-001/giw/value"])
    g0_raw =  read(f["dmft-last/ineq-001/g0iw/value"])
    giw = 0.5*(g_raw[:,1]+g_raw[:,2])[:]
    g0iw = 0.5*(g0_raw[:,1]+g0_raw[:,2])[:]
    return g0iw, giw
end

function calc_2_part_gf(g2worm)
	gupup = 0.5 .* (g2worm[1] .+ g2worm[6])
	gupdown = 0.5 .* (g2worm[2] .+ g2worm[5])
	gupdownbar = 0.5 .* (g2worm[3] .+ g2worm[4])
	return gupup, gupdown, gupdownbar
end

function calc_m_d_gf(gupup,gupdown)
	g_m = gupup .- gupdown
	g_d = gupup .+ gupdown
	return g_m, g_d
end

function calc_chi(beta,g1,g2)
	chiout = g2
	om_0 = floor(Int64,length(g2[:,1,1])/2+0.5) #cast to int
	n_nu = length(g2[1,:,1])
	for j in 1:n_nu
		for k in 1:n_nu
			chiout[om_0,j,k] -= beta*g1[floor(Int64,j)] *g1[floor(Int64,k)] 
		end
	end
	return chiout
end

function write_compatible_file_2part(chiupup,chiupdown,output_path)
	n_om = length(chiupup[:,1,1])
	n_nu = length(chiupup[1,:,1])
	open(output_path,"w+") do f
		for i in 0:n_om*n_nu*n_nu-1
			nu2_count = i%n_nu+1
			nu1_count = floor(Int64,i/n_nu)%n_nu+1
			om_count = floor(Int64,i/n_nu/n_nu)+1
			write(f,string(floor(Int64,om_count-1-(n_om-1)/2))*"\t"*string(floor(Int64,nu1_count-1-n_nu/2))*"\t"*string(floor(Int64,nu2_count-1-n_nu/2))*"\t"*string(real.(chiupup[om_count,nu1_count,nu2_count]))*"\t"*string(imag.(chiupup[om_count,nu1_count,nu2_count]))*"\t"*string(real.(chiupdown[om_count,nu1_count,nu2_count]))*"\t"*string(imag.(chiupdown[om_count,nu1_count,nu2_count]))*"\n")
		end
	end
end

function read_chi_phys_ph(file_path)
	chi_phys_in = read(h5open(file_path,"r")["stat-001/ineq-001/ntau-n0/value"])
	chi_phys_out = Array{Float64,2}(undef,length(chi_phys_in[:,1,1,1,1]),4)
	chi_phys_out[:,1] = chi_phys_in[:,1,1,1,1]
	chi_phys_out[:,2] = chi_phys_in[:,1,1,2,1]
	chi_phys_out[:,3] = chi_phys_in[:,2,1,1,1]
	chi_phys_out[:,4] = chi_phys_in[:,2,1,2,1]
	return chi_phys_out
end

function read_chi_phys_ph_worm(file_path)
	f = h5open(file_path,"r")
	chi_ph_uu = 0.5*(read(f["worm-last/ineq-001/p2iw-worm/00001/value"])+read(f["worm-last/ineq-001/p2iw-worm/00016/value"]))
	chi_ph_ud = 0.5*(read(f["worm-last/ineq-001/p2iw-worm/00004/value"])+read(f["worm-last/ineq-001/p2iw-worm/00013/value"]))	
	return chi_ph_uu, chi_ph_ud
end

function read_chi_phys_pp(file_path)
	chi_pp = Array{Array{ComplexF64, 1},1}(undef,4)
	f = h5open(file_path,"r")
	chi_pp[1] = read(f["worm-last/ineq-001/p2iwpp-worm/00004/value"])
	chi_pp[2] = read(f["worm-last/ineq-001/p2iwpp-worm/00007/value"])
	chi_pp[3] = read(f["worm-last/ineq-001/p2iwpp-worm/00010/value"])
	chi_pp[4] = read(f["worm-last/ineq-001/p2iwpp-worm/00013/value"])
	return chi_pp
end

function convert_chi_asympt(omrange,nfill,beta,chi_phys_tau)
	taurange = length(chi_phys_tau[:,1])-1
	chi_phys_om_upup =  zeros(ComplexF64,omrange)
	chi_phys_om_updown =  zeros(ComplexF64,omrange)
	for om in 0:omrange-1
		for t in 0:taurange-1
			chi_phys_om_upup[om+1] += exp(2*pi*1im*t*om/taurange)*(chi_phys_tau[t+1,1]+chi_phys_tau[t+1,4])*0.5/taurange
			chi_phys_om_updown[om+1] += exp(2*pi*1im*t*om/taurange)*(chi_phys_tau[t+1,2]+chi_phys_tau[t+1,3])*0.5/taurange
		end
	end
	chi_phys_om_upup[1] = chi_phys_om_upup[1]-nfill^2/4
	chi_phys_om_updown[1] = chi_phys_om_updown[1]-nfill^2/4
	return beta^3  .*chi_phys_om_upup, beta^3  .* chi_phys_om_updown
end

function calc_chi_pp(chi_pp)
	chiupdown = 0.5 .* (chi_pp[1] .+ chi_pp[4])
	chiupdownbar = 0.5 .* (chi_pp[2] .+ chi_pp[3])
	return chiupdown, chiupdownbar
end

function write_compatible_file_asympt(chiupup,chiupdown,chipp,output_path,beta)
	n_om = length(chiupup[:])
	open(output_path,"w+") do f
		for i in 0:n_om-1
			om_freq = 2*pi*i/beta
			write(f,string(om_freq)*"\t"*string(real.(chiupup[i+1]))*"\t"*string(imag.(chiupup[i+1]))*"\t"*string(real.(chiupdown[i+1]))*"\t"*string(imag.(chiupdown[i+1]))*"\t"*string(real.(chipp[i+1]))*"\t"*string(imag.(chipp[i+1]))*"\n")
		end
	end
end

function write_compatible_files_g1(g0iw,giw,beta,output_path_g0,output_path_g)
    n_nu = length(g0iw)
	open(output_path_g0,"w+") do f
		for i in 0:n_nu-1
            nu_freq = 2*pi*(i+1)/beta
            write(f,string(nu_freq)*"\t"*string(real.(1 ./(g0iw[i+1])))*"\t"*string(imag.(1 ./g0iw[i+1]))*"\n")
		end
	end
	open(output_path_g,"w+") do g 
		for i in 0:n_nu-1
            nu_freq = 2*pi*(i+1)/beta
            write(g,string(nu_freq)*"\t"*string(real.((giw[i+1])))*"\t"*string(imag.(giw[i+1]))*"\n")
		end
	end
end

function convert_w2dyn_g1(inputfile_g1,output_path_g,output_path_g0)
    f_g1 = h5open(inputfile_g1,"r")
	beta = read(attributes(f_g1[".config"])["general.beta"])
    close(f_g1)
    g0iw, giw = read_1_part_gf(inputfile_g1)
    write_compatible_files_g1(g0iw,giw,beta,output_path_g0,output_path_g)
end

function convert_w2dyn_vertex(input_file_g1,input_file_g2,output_path)
	f_g1 = h5open(input_file_g1,"r")
	beta = read(attributes(f_g1[".config"])["general.beta"])
	gupup, gupdown, gupdownbar = calc_2_part_gf(read_worm_comp(input_file_g2))
	g1_bare = read(f_g1["dmft-last/ineq-001/giw/value"])
	g1 = 0.5*(g1_bare[:,2,1]+g1_bare[:,1,1])[floor(Int64,length(g1_bare[:,1,1])/2)-floor(Int64,length(gupup[1,:,1])/2)+1:floor(Int64,length(g1_bare[:,1,1])/2)+floor(Int64,length(gupup[1,:,1])/2)]
	chiupup = beta .* calc_chi(beta,g1,gupup)
	chiupdown = beta .* calc_chi(beta,g1,gupdown)
	write_compatible_file_2part(chiupup,chiupdown,output_path)
end

function convert_w2dyn_chi_phys(input_file_chi_ph,input_file_chi_pp,output_path)
	ph_file = h5open(input_file_chi_ph)
	beta = read(attributes(ph_file[".config"])["general.beta"])
	nfill = read(ph_file["stat-last/ineq-001/rho1/value"])[1,1,1,1]+read(ph_file["stat-001/ineq-001/rho1/value"])[2,1,2,1]
	print(string(nfill)*"\n")
	chi_phys_tau_ph = read_chi_phys_ph(input_file_chi_ph)
	chi_pp_in = read_chi_phys_pp(input_file_chi_pp)
	omrange = length(chi_pp_in[1])
	om_0 = floor(Int64,(omrange+1)/2)
	chiupup, chiupdown = convert_chi_asympt(om_0,nfill,beta,chi_phys_tau_ph)
	chipp, chi_not_needed = beta^2*calc_chi_pp(chi_pp_in)
	write_compatible_file_asympt(chiupup,chiupdown,chipp[om_0:end],output_path,beta)
end

function convert_w2dyn_chi_phys_worm(input_file_chi_ph,input_file_chi_pp,output_path)
	ph_file = h5open(input_file_chi_ph)
	beta = read(attributes(ph_file[".config"])["general.beta"])
	nfill = read(ph_file["worm-last/ineq-001/rho1/value"])[1,1,1,1]+read(ph_file["worm-last/ineq-001/rho1/value"])[2,1,2,1]
	chi_pp_in = read_chi_phys_pp(input_file_chi_pp)
	om_0 = floor(Int64,(length(chi_pp_in[1])+1)/2)
	chiupup, chiupdown = read_chi_phys_ph_worm(input_file_chi_ph)
	chipp, chi_not_needed = beta^2*calc_chi_pp(chi_pp_in)
	write_compatible_file_asympt(chiupup[om_0:end],chiupdown[om_0:end],chipp[om_0:end],output_path,beta)
end
