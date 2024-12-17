ALL_FIGURE_NAMES=$(shell cat Differential_Geometry.figlist)
ALL_FIGURES=$(ALL_FIGURE_NAMES:%=%.pdf)

allimages: $(ALL_FIGURES)
^^I@echo All images exist now. Use make -B to re-generate them.

FORCEREMAKE:

-include $(ALL_FIGURE_NAMES:%=%.dep)

%.dep:
^^Imkdir -p "$(dir $@)"
^^Itouch "$@" # will be filled later.

figs/Differential_Geometry-figure0.pdf: 
^^Ixelatex -halt-on-error -interaction=batchmode -jobname "figs/Differential_Geometry-figure0" "\def\tikzexternalrealjob{Differential_Geometry}\input{Differential_Geometry}"

figs/Differential_Geometry-figure1.pdf: 
^^Ixelatex -halt-on-error -interaction=batchmode -jobname "figs/Differential_Geometry-figure1" "\def\tikzexternalrealjob{Differential_Geometry}\input{Differential_Geometry}"

figs/Differential_Geometry-figure2.pdf: 
^^Ixelatex -halt-on-error -interaction=batchmode -jobname "figs/Differential_Geometry-figure2" "\def\tikzexternalrealjob{Differential_Geometry}\input{Differential_Geometry}"

figs/Differential_Geometry-figure2.pdf: figs/Differential_Geometry-figure2.md5
figs/Differential_Geometry-figure3.pdf: 
^^Ixelatex -halt-on-error -interaction=batchmode -jobname "figs/Differential_Geometry-figure3" "\def\tikzexternalrealjob{Differential_Geometry}\input{Differential_Geometry}"

figs/Differential_Geometry-figure3.pdf: figs/Differential_Geometry-figure3.md5
