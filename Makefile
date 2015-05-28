SHELL := /bin/bash
CPP := g++

DIRS := lib bin

CFLAGS := -std=c++11 -Wall -O3 -Itools -Iparts

FJ_DIR    := $(shell fastjet-config --prefix)
FJ_CFLAGS := -I$(FJ_DIR)/include
FJ_LIBS   := -L$(FJ_DIR)/lib -lfastjet

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

LHAPDF_CFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LIBS   := $(shell lhapdf-config --ldflags)

.PHONY: all misc clean

HIST_SRC := $(filter-out src/hist_weights.cc,$(wildcard src/hist_*.cc))
HIST_OBJ := $(patsubst src/%.cc,lib/%.o,$(HIST_SRC))
HIST_EXE := $(patsubst src/%.cc,bin/%,$(HIST_SRC))

all: $(DIRS) bin/inspect_bh bin/reweigh bin/plot bin/merge_parts bin/overlay bin/scale bin/cross_section $(HIST_EXE)

misc: bin/hist_weights bin/cross_section_hist bin/cross_section_bh

# directories #######################################################
$(DIRS):
	@mkdir -p $@

# tools #############################################################
lib/timed_counter.o: lib/%.o: tools/%.cc tools/%.hh
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) -c $(filter %.cc,$^) -o $@

lib/csshists.o lib/hist_range.o: lib/%.o: tools/%.cc tools/%.hh
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

# parts #############################################################
lib/BHEvent.o lib/SJClusterAlg.o lib/weight.o lib/hist_wt.o: lib/%.o: parts/%.cc parts/%.hh
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

lib/rew_calc.o: lib/%.o: parts/%.cc parts/%.hh parts/BHEvent.hh
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) $(LHAPDF_CFLAGS) -c $(filter %.cc,$^) -o $@

# main objects ######################################################
lib/inspect_bh.o lib/reweigh.o lib/plot.o lib/merge_parts.o lib/overlay.o lib/hist_weights.o lib/cross_section_hist.o lib/cross_section_bh.o lib/scale.o lib/cross_section.o: lib/%.o: src/%.cc
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

$(HIST_OBJ): lib/%.o: src/%.cc
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) $(FJ_CFLAGS) \
		-DCONFDIR="\"`pwd -P`/config\"" \
		-c $(filter %.cc,$^) -o $@

# executables #######################################################
bin/cross_section_bh bin/cross_section_hist bin/inspect_bh bin/merge_parts bin/scale: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS)

bin/reweigh: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(LHAPDF_LIBS) -lboost_program_options

bin/hist_weights: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) -Wl,--no-as-needed $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_regex

bin/plot: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options -lboost_regex

bin/overlay: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS)

bin/cross_section: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) -Wl,--no-as-needed $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(FJ_LIBS) -lboost_program_options

$(HIST_EXE): bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) -Wl,--no-as-needed $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(FJ_LIBS) -lboost_program_options -lboost_regex

# Objects' dependencies #############################################
lib/hist_wt.o: parts/weight.hh

lib/inspect_bh.o: parts/BHEvent.hh

lib/reweigh.o: tools/timed_counter.hh parts/rew_calc.hh parts/BHEvent.hh

lib/hist_weights.o: tools/csshists.hh

lib/overlay.o: tools/propmap.hh tools/hist_range.hh

lib/cross_section_bh.o: parts/BHEvent.hh

lib/cross_section.o: tools/int_range.hh tools/real_range.hh tools/timed_counter.hh parts/BHEvent.hh parts/weight.hh parts/fj_jetdef.hh

$(HIST_OBJ): tools/csshists.hh tools/int_range.hh tools/real_range.hh tools/timed_counter.hh tools/catstr.hh tools/senum.hh parts/BHEvent.hh parts/SJClusterAlg.hh parts/weight.hh parts/hist_wt.hh parts/fj_jetdef.hh

# EXE dependencies ##################################################
bin/inspect_bh: lib/BHEvent.o

bin/reweigh: lib/timed_counter.o lib/rew_calc.o lib/BHEvent.o

bin/hist_weights: lib/csshists.o

bin/overlay: lib/hist_range.o

bin/cross_section_bh: lib/BHEvent.o

bin/cross_section: lib/timed_counter.o lib/BHEvent.o lib/weight.o

$(HIST_EXE): lib/csshists.o lib/timed_counter.o lib/BHEvent.o lib/SJClusterAlg.o lib/weight.o lib/hist_wt.o

clean:
	rm -rf bin/* lib/*
