EXE=${GRASP}/bin/rmcdhf
LIBS=-L ${GRASP}/lib/ -ldvd90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/libdvd90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

# PS
# References to the co_* files, related to the continuum orbitals calculations, have been added bellow
# PS END

OBJS= \
	mpi_s.o \
	dsubrs.o \
	dsubrs_I.o \
	fco.o \
	fco_I.o \
	gco.o \
	gco_I.o \
	setcof.o \
	setcof_I.o \
	lodcsh2GG.o \
	lodcsh2GG_I.o \
	setcsl.o \
	setcsl_I.o \
	setdbg.o \
	setdbg_I.o \
	setham.o \
	setham_I.o \
	xpot.o \
	xpot_I.o \
	ypot.o \
	ypot_I.o \
	dacon.o \
	dacon_I.o \
	dampck.o \
	dampck_I.o \
	setlag.o \
	setlag_I.o \
	setmcp.o \
	setmcp_I.o \
	setmix.o \
	setmix_I.o \
	setsum.o \
	setsum_I.o \
	setxuv.o \
	setxuv_I.o \
	setxv.o \
	setxv_I.o \
	setxz.o \
	setxz_I.o \
	lagcon.o \
	lagcon_I.o \
	cofpot.o \
	cofpot_I.o \
	consis.o \
	consis_I.o \
	csfwgt.o \
	csfwgt_I.o \
	dampor.o \
	dampor_I.o \
	defcor.o \
	defcor_I.o \
	eigen.o \
	eigen_I.o \
	engoutgg.o \
	engoutgg_I.o \
	endsum.o \
	endsum_I.o \
	estim.o \
	estim_I.o \
	getaldwt.o \
	getaldwt_I.o \
	getald.o \
	getald_I.o \
	prtrsl.o \
	prtrsl_I.o \
	getoldwt.o \
	getoldwt_I.o \
	getold.o \
	getold_I.o \
	getscd.o \
	getscd_I.o \
	hmout.o \
	hmout_I.o \
	out.o \
	out_I.o \
	in.o \
	in_I.o \
	prwf.o \
	prwf_I.o \
	ispar.o \
	ispar_I.o \
	itjpo.o \
	itjpo_I.o \
	maneig.o \
	maneig_I.o \
	matrix.o \
	matrix_I.o \
	maxarr.o \
	maxarr_I.o \
	newco.o \
	newco_I.o \
	outbnd.o \
	outbnd_I.o \
	newe.o \
	newe_I.o \
	orbout.o \
	orbout_I.o \
	orthor.o \
	orthor_I.o \
	orthy.o \
	orthy_I.o \
	solve.o \
	solve_I.o \
	improv.o \
	improv_I.o \
	scf.o \
	scf_I.o \
	strsum.o \
	strsum_I.o \
	co_spherical_bessel.o \
	co_getmax.o \
	co_save_to_file.o \
	co_normalization.o \
	co_phase_shift.o \
	co_scattering_length.o \
	co_add_polarization_potential.o \
	rscfvu.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
