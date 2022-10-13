using DelimitedFiles
using Statistics

nu = readdlm("nusse.out")
sz = size(nu,1)/10
av1 = mean(nu[end-2000:end-1000,2])
av2 = mean(nu[end-1000:end,2])
test1 = abs((av2-av1)/av2)*100

if test1>2
	println("Difference greater than 2%")
end

#println("Mean plate NU for 100 time units is $av1 and $av2 with total run time $sz")
#println("$av1  $av2")
writedlm("nusse.out", [av1 av2] )

rms = readdlm("rms_vel.out")
col2 = mean(rms[end-2500:end,2])
col3 = mean(rms[end-2500:end,3])

#println("$col2  $col3")
writedlm("avg_rms_vals", [col2 col3] )

# Mean temperature calculations
using HDF5
f=h5open("stafield_master.h5","r")
fd=h5open("stafield_data.h5","r")
names(f)
A=read(f)
Ad=read(fd)
# vxmean=Ad["vx_mean"]/A["averaging_time"]
vymean=Ad["Vr_mean"]./A["averaging_time"]
vzmean=Ad["Vz_mean"]./A["averaging_time"]
tempmean=Ad["Vth_mean"]./A["averaging_time"]
# vxrms=Ad["vx_rms"]./A["averaging_time"]
vyrms=Ad["Vr_rms"]./A["averaging_time"]
vzrms=Ad["Vz_rms"]./A["averaging_time"]
temprms=Ad["Vth_rms"]./A["averaging_time"]
using DelimitedFiles
f=readdlm("axicor.out")
zm=f[1:end-1,2]
using Plots
plotly()
plot(zm,vzrms)
plot!(zm,vyrms)
