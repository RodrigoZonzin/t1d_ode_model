push: 
	git add .
	git commit -m "upado $(shell date +%d-%m-%Y) by $(USER)"
	git push
