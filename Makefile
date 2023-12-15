local:
	pipen ref --pipeline "immunopipe:Immunopipe" \
		--destdir docs/processes \
		--replace-titles "Envs=Environment Variables" \
		--hide-sections Requires \
		--hide-sections Input \
		--hide-sections Output \
		--replace "https://pwwang.github.io/immunopipe=.." \
		--replace "${HOME}=~" \
		--replace "${USER}=user"

api:
	pipen ref --pipeline "immunopipe:Immunopipe" \
		--destdir docs/processes \
		--replace-titles "Envs=Environment Variables" \
		--hide-sections Requires \
		--hide-sections Input \
		--hide-sections Output

.PHONY: local api
