FYPPS=$(wildcard *.fpp)
F90S=$(patsubst %.fpp,%.f90,$(wildcard *.f90 *.fpp))
OBJS=$(patsubst %.f90,%.o,$(F90S))

MOD_OBJS=$(patsubst %.fpp,%.o,$(patsubst %.f90,%.o,$(wildcard m*.f90 m*.fpp)))

FC = gfortran

EXEC = run_code
PROGRAM = p_main
PRG_OBJ = $(PROGRAM).o

all : $(PROGRAM)

define CREATE_FYPP_RULE
$(2) : $(1)
	fypp $(FYPPFLAGS) $(1) $(2)
endef

define CREATE_FC_RULE
$(2) : $(1)
	$(FC) $(FFLAGS) -c -o $(2) $< $(1) $(LNKFLAGS)
endef

$(foreach l,$(FYPPS),$(eval $(call CREATE_FYPP_RULE,$(l),$(patsubst %.fpp,%.f90,$(l)))))
$(foreach l,$(F90S), $(eval $(call CREATE_FC_RULE,  $(l),$(patsubst %.f90,%.o,  $(l)))))

$(PROGRAM) : $(OBJS)
	$(FC) $(OBJS) $(FFLAGS) -o $(EXEC) $(LNKFLAGS)

$(PRG_OBJ) : $(MOD_OBJS)

.PHONY: all clean

clean:
	rm -rf $(OBJS) $(PROGRAM) $(EXEC) $(patsubst %.o,%.mod,$(MOD_OBJS))

# Modules depend on each other
# m_startup.o : m_parameters.o
# m_bcs.o : m_parameters.o
# m_weno.o : m_parameters.o m_bcs.o
# m_qpqs.o : m_parameters.o m_weno.o
# m_dt.o : m_parameters.o
# m_leftright.o : m_parameters.o
# m_flux.o : m_parameters.o m_qpqs.o
# m_riemann_solvers.o : m_parameters.o m_leftright.o  m_flux.o 
# m_rhs.o : m_parameters.o m_riemann_solvers.o
# m_rk.o : m_parameters.o m_qpqs.o m_rhs.o
#
#
