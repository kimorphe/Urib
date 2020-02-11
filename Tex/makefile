saft.pdf: saft.dvi
	dvipdfmx -p a4 -r 600 saft.dvi
saft.dvi: saft.tex paths.eps length.eps
	platex saft.tex

paths.eps: paths.svgz
	inkscape -z -f paths.svgz -E paths.eps
length.eps: length.svgz
	inkscape -z -f length.svgz -E length.eps



