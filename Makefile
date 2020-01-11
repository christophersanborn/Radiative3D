CC=gcc
CPP=g++
FLAGS=-O3 -g -Wall
DEFINES=

build_dir = build
objdir = $(build_dir)

OUT_EXEC = main
objects  = geom_s2.o geom_r3.o geom_r4.o probability.o sources.o \
           scatterers.o events.o phonons.o grid.o media.o model.o \
	   rtcoef.o dataout.o scatparams.o cmdline.o user.o global.o \
           ecs.o elastic.o main.o

objects := $(addprefix $(objdir)/,$(objects))

.PHONY : default cleanall clean neat anyway directories .FORCE
default : directories $(OUT_EXEC)

directories : $(objdir)

$(objdir) :
	mkdir -p $@

##//////
## Versioning:
##
ifneq ($(wildcard .git),)
  gitstatus = $(shell git status --untracked-files=no --porcelain)
  ifneq ($(gitstatus),)
    gitstatus = , with local modifications
  endif
  revision = "\"$(shell git rev-parse --short HEAD)$(gitstatus) (git)\""
else ifneq ($(wildcard .svn),)
  branch = $(shell basename "`svn info | grep '^URL:'`")
  revision = "\"$(branch)@$(shell svnversion) (svn)\""
else
  revision =
endif
ifneq ($(revision),)
  DEFINES += -DREVISION_NUM=$(revision)
endif

##//////
## Platform-Specific Stuff:
##
machine = $(shell uname -m)
ifeq ($(machine),armv7l)
  FLAGS += -Wno-psabi   # Silence ABI warnings on Raspberry Pi
endif


##//////
## Configutation Options:
##
##   These definitions and recipes encapsulate compile-time build
##   options.
##

.PHONY : config-float config-double
.FORCE :

configopt_files = config/opt-fptype.hpp

config/opt-fptype.hpp :
	@$(MAKE) config-double
        # Default to config-double if opt-fptype does not exist

config-float :
	@echo "Configure - Floating-point type: float"
	@echo '// ***Auto-Generated***' > config/opt-fptype.hpp
	@echo '#define TYPEDEFS_H_FP_TYPE float' >> config/opt-fptype.hpp

config-double :
	@echo "Configure - Floating-point type: double"
	@echo '// ***Auto-Generated***' > config/opt-fptype.hpp
	@echo '#define TYPEDEFS_H_FP_TYPE double' >> config/opt-fptype.hpp

config-long-double :
	@echo "Configure - Floating-point type: long double"
	@echo '// ***Auto-Generated***' > config/opt-fptype.hpp
	@echo '#define TYPEDEFS_H_FP_TYPE long double' >> config/opt-fptype.hpp


##//////
## Headers, and their dependencies:
##

#      /* Fundamentals: */
#

typedefs_hpp = typedefs.hpp config/opt-fptype.hpp
raytype_hpp  = raytype.hpp
array_hpp    = array.hpp
complex_hpp  = complex.hpp  $(typedefs_hpp)
probability_hpp = probability.hpp $(typedefs_hpp)

#
#      /* Libraries: */
#

geom_base_hpp = geom_base.hpp $(typedefs_hpp)
geom_s2_hpp   = geom_s2.hpp $(geom_base_hpp)
geom_r3_hpp   = geom_r3.hpp $(geom_base_hpp)
geom_r4_hpp   = geom_r4.hpp $(geom_r3_hpp)
geom_hpp      = geom.hpp    $(geom_r3_hpp) $(geom_s2_hpp)
tensors_hpp   = tensors.hpp  $(geom_r3_hpp)
elastic_hpp   = elastic.hpp $(typedefs_hpp)
ecs_hpp       = ecs.hpp     $(geom_hpp) $(elastic_hpp)


#
#      /* Radiative3D Computation Core: */
#

rtcoef_hpp  = rtcoef.hpp  $(geom_hpp) $(raytype_hpp) $(complex_hpp)
phonons_hpp = phonons.hpp $(geom_hpp) $(raytype_hpp)
sources_hpp = sources.hpp $(geom_hpp) $(raytype_hpp) $(probability_hpp)
events_hpp  = events.hpp  $(sources_hpp) $(tensors_hpp)
scatparams_hpp = scatparams.hpp $(geom_hpp) $(elastic_hpp)
scatterers_hpp = scatterers.hpp $(sources_hpp) $(scatparams_hpp)
grid_hpp    = grid.hpp    $(ecs_hpp)
media_hpp   = media.hpp   $(raytype_hpp) $(elastic_hpp) $(geom_hpp) $(array_hpp)
model_hpp   = model.hpp   $(events_hpp) $(scatterers_hpp) $(grid_hpp) $(media_hpp)

#
#      /* Radiative3D User Interaction and Input/Output: */
#

params_hpp  = params.hpp
cmdline_hpp = cmdline.hpp $(geom_hpp)
dataout_hpp = dataout.hpp $(geom_hpp) $(raytype_hpp) $(tensors_hpp)

#
#


##//////
## Executable:
##
$(OUT_EXEC) : $(objects)
	$(CPP) -c main.cpp $(FLAGS) -DREVISION_NUM=$(revision) -o $(objdir)/main.o
	$(CPP) $^ -o $@ $(FLAGS)

anyway : $(objects)
	$(CPP) -c main.cpp $(FLAGS) $(DEFINES)
	$(CPP) $^ -o $(OUT_EXEC) $(FLAGS)


##//////
## Objects Files:
##

comd =  # (Common Build Dependencies)
        # (If these files change, rebuild everything)

$(objdir)/geom_s2.o : geom_s2.cpp $(geom_s2_hpp) $(geom_r3_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/geom_r3.o : geom_r3.cpp $(geom_r3_hpp) $(geom_s2_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/geom_r4.o : geom_r4.cpp $(geom_r4_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/probability.o : probability.cpp $(probability_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/elastic.o : elastic.cpp $(elastic_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/sources.o : sources.cpp $(sources_hpp) $(phonons_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/events.o : events.cpp $(events_hpp) $(phonons_hpp) $(dataout_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/scatparams.o : scatparams.cpp $(scatparams_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/scatterers.o : scatterers.cpp $(scatterers_hpp) $(phonons_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/phonons.o : phonons.cpp $(phonons_hpp) $(media_hpp) $(rtcoef_hpp) \
                        $(scatterers_hpp) $(dataout_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/rtcoef.o : rtcoef.cpp $(rtcoef_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/grid.o : grid.cpp $(grid_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/media.o : media.cpp $(media_hpp) $(grid_hpp) $(rtcoef_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/model.o : model.cpp $(model_hpp) $(phonons_hpp) $(dataout_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/dataout.o: dataout.cpp $(dataout_hpp) $(phonons_hpp) $(media_hpp) \
	               $(ecs_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/cmdline.o: cmdline.cpp $(cmdline_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/user.o : user.cpp  $(grid_hpp) $(comd) user_*_inc.cpp
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/global.o : global.cpp  $(ecs_hpp) $(dataout_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/ecs.o : ecs.cpp $(ecs_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@

$(objdir)/main.o : main.cpp  $(params_hpp) $(model_hpp) $(dataout_hpp) \
                   $(rtcoef_hpp) $(cmdline_hpp) $(comd)
	$(CPP) -c $< $(FLAGS) -o $@


## Cleanup:
##
cleanall :
	rm -f $(configopt_files)
	rm -f $(objects)
	rm -f $(OUT_EXEC)

clean :
	@echo "*NOTE: 'Make clean' removes object files and target exe, but"
	@echo "        does NOT remove auto-generated config files.  To remove"
	@echo "        these too, do 'make cleanall'."
	rm -f $(objects)
	rm -f $(OUT_EXEC)

neat :
	rm -f $(objects)
