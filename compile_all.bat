:: @echo off

for %%f in (*.tex) do (
	pdflatex "%%f"
	bibtex "%%f"
	pdflatex "%%f"
)
:: to pause
:: pause

:: to delete some auxiliary files
del *.log *.aux *.out *.aaa *.bak