assignment-1.pdf: assignment-1.tex ../plots/* ../../header.tex ../../../bibliography.bib

	pdflatex assignment-1.tex
	bibtex assignment-1.aux
	pdflatex assignment-1.tex
	pdflatex assignment-1.tex
	gnome-open assignment-1.pdf &
