SHELL=/bin/bash

environment.yml: FORCE
	@if [[ $$CONDA_DEFAULT_ENV -eq "base" ]]; then \
		conda env export -n immunopipe | grep -v "^prefix" > $@; \
	else \
		echo "Must run in base conda environment"; \
		exit 1; \
	fi

clean:
	echo "Cleaning"
	rm -f environment.yml

FORCE: ;
