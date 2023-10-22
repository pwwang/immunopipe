api:
	pipen ref --pipeline "immunopipe:Immunopipe" \
		--destdir docs/processes \
		--replace-titles "Envs=Environment Variables" \
		--hide-sections Requires \
		--hide-sections Input \
		--hide-sections Output
