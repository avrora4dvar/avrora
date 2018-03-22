#COMPILER SETTINGS
#Modify 

#Directory containing header (.h) files
CPP = ../Include/
#Options for preprocessing with header files
#Default: -Mpreprocess -I$(CPP)
CPPFLAGS = -Mpreprocess -I$(CPP)
#Fortran compiler
#pgf90 (default)|mpif90 (MPI)
FC = mpif90
#Options for fortran compiler
#OpenMP: -mp (OpenMP)/-O2 (optimalization)
FFLAGS = -O3 -fast 
#Location of libraries (netcdf,LPACK)
LIB = ../Lib/libnetcdf.a -L/usr/lib64/atlas -llapack -lblas
#LIB = -L/usr/lib64/ -lnetcdff
#Directory in which the executables will be created
#NB: .o,.mod files will be created in the current directory 
#List with dependencies between .f files and modules
dependencyList = MakeDepend


#############################################################################
# DEPENDENCIES AND RULES
# Do not modify

clean_list := core *.o *.oo *.mod *.f lib*.a *.bak
.SUFFIXES: .o .f .F .F90 
%.o: %.F 
	$(FC) -c $(FFLAGS) $(CPPFLAGS) $^ $(LIB)
%.o: %.F90
	$(FC) -c $(FFLAGS) $(CPPFLAGS) $^ $(LIB)

modules :=	../Src/Modules \
		../Src/Nonlinear \
		../Src/Tangent \
		../Src/Adjoint \
		../Src/Utility \
		../Src/Cov \
                ../Tools  
#Set file types for which make must look
vpath %.F90 $(modules): 
vpath %.F $(modules)

###########################################################################
# VARIABLES
# For each executable create a list of associated .F/.F90 files

# VARIABLES concerning ad_sample_field
ad_sample_field_SRCS = ad_sample_field.F \
                       ad_sample_surf_pnt.F \
                       ad_interp_3d.F \
                       ad_demean_ssh.F \
                       ad_daily.F \
                       ad_avg_2d.F \
                       ad_interp_4d.F \
                       interp1.F \
                       nc_read.F \
                       nc_dim.F \
                       create_his_file.F \
                       zero_his_file.F \
                       unique.F 

test_SRCS =  	test.F90 \
		mod_BHM.F90
test_OBJS=$(test_SRCS:.F=.o)
test_BIN= test
$(test_BIN): $(test_OBJS)
	$(FC) $(FFLAGS) $(LIB) $(CPPFLAGS) -o $@ $^


ad_sample_field_OBJS=$(ad_sample_field_SRCS:.F=.o)
ad_sample_field_BIN= ad_sample_field
$(ad_sample_field_BIN): $(ad_sample_field_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

# VARIABLES tl_identity
tl_identity_SRCS = tl_identity.F90 \
		   mod_netcdf.F90 \
                   mod_roms.F90
tl_identity_OBJS = $(tl_identity_SRCS:.F90=.o)
tl_identity_BIN=tl_identity
$(tl_identity_BIN): $(tl_identity_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

# VARIABLES ad_identity
ad_identity_SRCS = ad_identity.F90 \
		   mod_netcdf.F90 \
                   mod_roms.F90
ad_identity_OBJS = $(ad_identity_SRCS:.F90=.o)
ad_identity_BIN=ad_identity
$(ad_identity_BIN): $(ad_identity_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)


# VARIABLES tl_sample
tl_sample_SRCS =	tl_sample.F90 \
			mod_netcdf.F90 \
			mod_interp.F90 \
                        mod_sample.F90
tl_sample_OBJS = $(tl_sample_SRCS:.F90=.o)
tl_sample_BIN=tl_sample
$(tl_sample_BIN): $(tl_sample_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

# VARIABLES ad_sample
ad_sample_SRCS =	ad_sample.F90 \
			mod_netcdf.F90 \
			mod_interp.F90 \
                        mod_sample.F90 \
                        mod_roms.F90
ad_sample_OBJS = $(ad_sample_SRCS:.F90=.o)
ad_sample_BIN=ad_sample
$(ad_sample_BIN): $(ad_sample_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)


# VARIABLES concerning tl_sample_field
tl_sample_field_SRCS = tl_sample_field.F \
                       tl_sample_surf_pnt.F \
		       tl_interp_3d.F \
                       tl_interp_4d.F \
                       tl_demean_ssh.F \
		       tl_daily.F \
                       interp1.F \
                       nc_read.F \
                       nc_dim.F \
                       unique.F \
                       tl_avg_2d.F 
tl_sample_field_OBJS = $(tl_sample_field_SRCS:.F=.o)
tl_sample_field_BIN = tl_sample_field
$(tl_sample_field_BIN): $(tl_sample_field_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)


#------------------------------------------------------------------------

# VARIABLES concerning AD:
ad_SRCS = ad_avrora.F       \
          ad_get_data.F     \
          ad_exchange.F     \
          ad_ini_out.F      \
          ad_init_ubar.F    \
          ad_init_zeta.F    \
          ad_coriolis.F     \
          ad_comb_frc.F     \
          ad_set_depth.F    \
          ad_pre_step3d.F   \
          ad_prsgrd.F       \
          ad_omega.F        \
          ad_rhs3d.F        \
          ad_rho_eos.F      \
          ad_set_data.F     \
          ad_set_massflux.F \
          ad_set_vbc.F      \
          ad_step2d.F       \
          ad_step3d_t.F     \
          ad_step3d_uv.F    \
          ad_t3dbc.F        \
          ad_t3dmix2.F      \
          ad_u2dbc.F        \
          ad_u3dbc.F        \
          ad_uv_adv.F       \
          ad_uv3dmix2.F     \
          ad_zetabc.F       \
          allocate_arrays.F \
          exchange.F        \
          fwd_fast_vars.F   \
          get_2dfld.F       \
          get_3dfld.F       \
          get_state_r.F     \
          inp_par.F         \
          omega.F           \
          read_frc_time.F   \
          read_fwd_time.F   \
          read_grid.F       \
          rho_eos.F         \
          set_bounds.F      \
          set_depth.F       \
          set_massflux.F    \
          set_state.F       \
          set_visc2.F       \
          set_weights.F     \
          mean_field.F      \
          mpi_exchange.F    \
	  mod_ocean.F       \
	  mod_parallel.F    \
	  mod_forces.F      \
	  mod_output.F      
ad_OBJS = $(ad_SRCS:.F=.o)
ad_BIN = ad_avrora
$(ad_BIN): $(ad_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)


# VARIABLES concerning TL:
tl_SRCS = tl_avrora.F       \
          tl_coriolis.F     \
          tl_get_data.F     \
          tl_initial.F      \
          tl_init_ubar.F    \
          tl_init_zeta.F    \
          tl_omega.F        \
          tl_output.F       \
          tl_set_depth.F    \
          tl_pre_step3d.F   \
          tl_prsgrd.F       \
          tl_rho_eos.F      \
          tl_rhs3d.F        \
          tl_set_data.F     \
          tl_set_massflux.F \
          tl_set_vbc.F      \
          tl_step2d.F       \
          tl_step3d_t.F     \
          tl_step3d_uv.F    \
          tl_t3dmix2.F      \
          tl_t3dbc.F        \
          tl_u2dbc.F        \
          tl_u3dbc.F        \
          tl_uv_adv.F       \
          tl_uv3dmix2.F     \
          tl_zetabc.F       \
          allocate_arrays.F \
          exchange.F        \
          fwd_fast_vars.F   \
          get_2dfld.F       \
          get_3dfld.F       \
          get_state.F       \
          inp_par.F         \
          omega.F           \
          read_frc_time.F   \
          read_fwd_time.F   \
          read_grid.F       \
          set_bounds.F      \
          set_depth.F       \
          set_massflux.F    \
          set_visc2.F       \
          set_state.F       \
          set_weights.F     \
          mean_field.F      \
	  mod_ocean.F       \
	  mod_parallel.F    \
	  mod_forces.F      \
	  mod_output.F      \
          mpi_exchange.F
tl_OBJS = $(tl_SRCS:.F=.o)
tl_BIN = tl_avrora
$(tl_BIN): $(tl_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

#---------------------------------------------------------------------------

# VARIABLES concerning ad_balance_3D
ad_balance_3d_SRCS = ad_balance_3D.F \
                     get_grddim.F \
                     ad_output_psi_T.F \
                     ad_prsgrd_SDJ.F \
                     index1d.F \
                     gauss_band.F \
                     ad_input_ini.F \
                     get_3dfld.F \
                     get_2dfld.F \
		     mod_output.F 
ad_balance_3d_OBJS = $(ad_balance_3d_SRCS:.F=.o)
ad_balance_3d_BIN = ad_balance_3D
$(ad_balance_3d_BIN): $(ad_balance_3d_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

# VARIABLES concerning COV_INI_TEMPU_UNI
cov_ini_temp_uni_SRCS = cov_ini_temp_uni.F
cov_ini_temp_uni_OBJS=$(cov_ini_temp_uni_SRCS:.F=.o)
cov_ini_temp_uni_BIN= cov_ini_temp_uni
$(cov_ini_temp_uni_BIN): $(cov_ini_temp_uni_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)


# VARIABLES concerning tl_balance_3D
tl_balance_3d_SRCS = tl_balance_3D.F \
                     tl_output_ini.F \
                     index1d.F \
                     tl_prsgrd_SDJ.F \
                     gauss_band.F \
                     get_3dfld.F \
                     get_2dfld.F \
                     get_grddim.F \
	             mod_output.F   
tl_balance_3d_OBJS = $(tl_balance_3d_SRCS:.F=.o)
tl_balance_3d_BIN=tl_balance_3D
$(tl_balance_3d_BIN): $(tl_balance_3d_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

#-------------------------------------------------------------------------

# VARIABLES concerning cov_ensO
cov_ensO_SRCS = cov_ens.F90 \
                mod_netcdf.F90 \
                mod_constants.F90 \
                mod_roms.F90 \
                mod_interp.F90
cov_ensO_OBJS=$(cov_ensO_SRCS:.F90=.o)
cov_ensO_BIN= cov_ens
$(cov_ensO_BIN): $(cov_ensO_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

#------------------------------------------------------------------------

create_pre_SRCS = create_pre.F90 \
            mod_netcdf.F90 \
            mod_BHM.F90 \
            mod_interp.F90 
create_pre_OBJS= $(create_pre_SRCS:.F90=.o)
create_pre_BIN=create_pre
$(create_pre_BIN): $(create_pre_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)


create_ini_SRCS = create_ini.F90 \
            mod_netcdf.F90 \
            mod_roms.F90\
            mod_interp.F90 
create_ini_OBJS= $(create_ini_SRCS:.F90=.o)
create_ini_BIN=create_ini
$(create_ini_BIN): $(create_ini_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

tl_interp_SRCS = tl_interp.F90 \
            mod_netcdf.F90 \
            mod_roms.F90 \
            mod_interp.F90 
tl_interp_OBJS= $(tl_interp_SRCS:.F90=.o)
tl_interp_BIN=tl_interp
$(tl_interp_BIN): $(tl_interp_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

ad_interp_SRCS = ad_interp.F90 \
            mod_netcdf.F90 \
            mod_roms.F90 \
            mod_interp.F90 
ad_interp_OBJS= $(ad_interp_SRCS:.F90=.o)
ad_interp_BIN=ad_interp
$(ad_interp_BIN): $(ad_interp_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

create_background_SRCS = create_background.F90 \
            mod_filter.F90 \
            mod_netcdf.F90 \
            mod_time.F90 \
            mod_roms.F90 \
            mod_interp.F90
create_background_OBJS= $(create_background_SRCS:.F90=.o)
create_background_BIN=create_background
$(create_background_BIN): $(create_background_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)
ini_gauss_SRCS = ini_fourier.F90 \
                 mod_netcdf.F90 \
                 mod_constants.F90 \
                 normal.F90
ini_gauss_OBJS = $(ini_gauss_SRCS:.F90=.o)
ini_gauss_BIN=ini_fourier
$(ini_gauss_BIN): $(ini_gauss_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

create_wind_SRCS = create_wind_eof.F90 \
                   mod_netcdf.F90 \
                   mod_BHM.F90
create_wind_OBJS = $(create_wind_SRCS:.F90=.o)
create_wind_BIN =create_wind
$(create_wind_BIN): $(create_wind_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

# VARIABLES concerning cov_ensO
rbcg_SRCS = RBCG.F90 \
            mod_netcdf.F90
rbcg_OBJS=$(rbcg_SRCS:.F90=.o)
rbcg_BIN= RBCG
$(rbcg_BIN): $(rbcg_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

# VARIABLES create_mean
create_mean_SRCS = create_mean.F90 \
		   mod_netcdf.F90 \
                   mod_roms.F90
create_mean_OBJS = $(create_mean_SRCS:.F90=.o)
create_mean_BIN=create_mean
$(create_mean_BIN): $(create_mean_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

# VARIABLES cov_scale
cov_scale_SRCS = cov_scale.F90 \
		  mod_netcdf.F90 \
                  mod_BHM.F90 
cov_scale_OBJS = $(cov_scale_SRCS:.F90=.o)
cov_scale_BIN=cov_scale
$(cov_scale_BIN): $(cov_scale_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)


# VARIABLES ini_scale
ini_scale_SRCS = ini_scale.F90 \
	         mod_netcdf.F90 \
                 mod_BHM.F90 
ini_scale_OBJS = $(ini_scale_SRCS:.F90=.o)
ini_scale_BIN=ini_scale
$(ini_scale_BIN): $(ini_scale_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

pcg_SRCS = PCG.F90 \
           mod_netcdf.F90 \
           mod_BHM.F90 
pcg_OBJS = $(pcg_SRCS:.F90=.o)
pcg_BIN=PCG
$(pcg_BIN): $(pcg_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)

mso_SRCS = MSO.F90 \
           mod_netcdf.F90 \
           mod_BHM.F90 
mso_OBJS = $(mso_SRCS:.F90=.o)
mso_BIN=mso
$(mso_BIN): $(mso_OBJS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LIB)



#-----------------------------------------------------------------------

#########################################################################
#  Cleaning targets.
#  Do not change
.PHONY: clean
clean:
	rm -f $(clean_list)
include $(dependencyList)