assignment-2.pdf: assignment-2.tex ../plots/* ../../header.tex ../../../bibliography.bib

	pdflatex assignment-2.tex
	bibtex assignment-2.aux
	pdflatex assignment-2.tex
	pdflatex assignment-2.tex
	gnome-open assignment-2.pdf &
