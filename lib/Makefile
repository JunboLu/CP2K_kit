f2py := /home/lujunbo/bin/deepmd-kit/bin/f2py3.9

all : dynamic_mod geometry_mod rmsd_mod statistic_mod

dynamic_mod : dynamic_mod.f90
	$(f2py) -c dynamic_mod.f90 -m dynamic_mod

geometry_mod : geometry_mod.f90
	$(f2py) -c geometry_mod.f90 -m geometry_mod

rmsd_mod : rmsd_mod.f90
	$(f2py) -c rmsd_mod.f90 -m rmsd_mod

statistic_mod : statistic_mod.f90
	$(f2py) -c statistic_mod.f90 -m statistic_mod
