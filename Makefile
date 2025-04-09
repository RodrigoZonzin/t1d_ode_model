push: 
	git add .
	git commit -m "upado $(shell date +%d-%m-%Y) by $(USER)"
	git push

run_solve: 
	python3 codigo_solveivp.py -o resultados.pdf
