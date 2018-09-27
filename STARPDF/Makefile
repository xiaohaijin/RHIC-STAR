EXE=STARBranch.pdf
SRC=STARBranch.tex
BIB=STARBranch.aux

$(EXE):$(SRC)
	xelatex $(SRC)

.PHONY.:look
look:
	evince STARBranch.pdf

.PHONY.:clean
clean:
	$(RM) STARBranch.aux STARBranch.log STARBranch.pdf STARBranch.dvi STARBranch.bbl STARBranch.blg STARBranch.toc STARBranch.out STARBranch.nav STARBranch.snm STARBranch.thm STARBranch.lot STARBranch.lof
