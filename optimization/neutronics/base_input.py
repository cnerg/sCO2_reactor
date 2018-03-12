from string import Template
# This is for r=0.5, PD=1.28, c=0.00031
base_string = Template("""\
MCNP6 pin cell study radius:${r_core} cm
c Cell Card
1 1 -${fuel_rho} -1 4 -5 imp:n=1 
2 2 -1.7 (-2 1 -7 6):(-1 5 -7):(6 -1 -4) imp:n=1 
99 0 7:-6:(2 -7 6) imp:n=0 
 
c Surface Card
1 CZ ${r_core}
2 CZ ${r_refl}
4 PZ 0
5 PZ ${core_z}
6 PZ ${refl_min}
7 PZ ${refl_max}

c Data Card
c
c Burnup
burn time=0.25 0.5 9R 10 5R 100 3R 365.25 8R power=0.09941 pfrac=1 29R 
        bopt=1 14 -1 mat= 1
      omit= -1 111  66159 67163 67164 67166 68163 68165 68169 69166   
         69167 69171 69172 69173 70168 70169 70170 70171 70172 70173  
         70174 6014 7016 39087 39092 39093 40089 40097 41091 41092    
         41096 41097 41098 41099 42091 42093 70175 70176 71173 71174  
         71177 72175 72181 72182 73179 73183 74179 74181 8018 8019    
         9018 10021 12027 13026 13028 14027 14031 16031 16035 16037
         17034 17036 17038 18037 18039 20045 20047 21044 21046 21047
         21048 22045 22051 23047 23048 23049 23052 23053 23054 24049
         24051 24055 24056 25051 25052 25053 25054 25056 25057 25058
         26053 26055 26059 26060 26061 27057 27060 27061 27062 27063 
         27064 28057 28063 28065 29062 29064 29066 40088 41100 42101
         43097 43098 44097
${fuel_string}
C name: Carbon, Graphite (reactor grade)
C density = 1.7
m2
     5010 -1.8431e-07
     5011 -8.1569e-07
     6012 -9.8841e-01
     6013 -1.1584e-02
VOL ${fuel_vol} ${refl_vol}
kcode 10000 1 15 60
ksrc  0 0 1
mode n
F4:n 1
E4 1E-9 10ILOG 1E-5 140ILOG 10
print""")

base_submit = Template("""\
# This is a "normal" job.
universe = vanilla

# Files to place stdout and stderr
output = ${htc_dir}/output/process_${sym}(Cluster)_${sym}(Process).out
error = ${htc_dir}/error/process_${sym}(Cluster)_${sym}(Process).err
log = ${htc_dir}/log/process_${sym}(Cluster)_${sym}(Process).log
#+IsBuildJob = true
#requirements = (Target.HasGluster == true) && (OpSysAndVer =?= "SL6") && (IsBuildSlot == true)
requirements = (Target.HasGluster == true) && (OpSysAndVer =?= "SL6") 

# Copy environment variables from the submit node
getenv = True
environment = "PATH=/bin:/usr/bin:/usr/local/bin LD_LIBRARY_PATH= PYTHONPATH="

# Indicate if/when files should be transfered
should_transfer_files = YES
when_to_transfer_output = ON_EXIT

# Indicate the executable to be run and any other input files needed
executable = ${htc_dir}/send_install_run_mcnp.bash
transfer_input_files = ${htc_dir}/mcnp_run_tar_send.bash, ${htc_dir}/inputs/${sym}(comp_file)

# Command line arguments for the executable
arguments = ${sym}(comp_file)

# Don't send emails
notification = never

# Request resources
request_cpus = 1
request_memory = ${mem}MB
request_disk = ${disk}GB

queue comp_file from input_list.txt""")
