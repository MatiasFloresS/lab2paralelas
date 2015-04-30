main:
	@clear
	@g++ -msse3 -msse4.1 -o sort sort.cpp -fopenmp -lm
	@echo "Compilación exitosa..."
	@./sort -i ../../../../opt/tpp/floats_134217728.raw -o ordenada.raw -N 134217728 -d 0 -L 2

clean:
	@echo "Limpiando..."
	@rm 	sort
	@rm ordenada.raw