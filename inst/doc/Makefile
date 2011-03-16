all: RNAinteractMAPK clean

RNAinteractMAPK: RNAinteractMAPK.tex
	texi2dvi --pdf RNAinteractMAPK.tex

clean:
	-$(RM) -f Figure-*
	-$(RM) -f Tab*
	-$(RM) -rf Dmel2PPMAPK
	-$(RM) -rf networks
	-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	-$(RM) -f Rplots.pdf

