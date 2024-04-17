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

test:
	@poetry run pytest tests

test-data:
	@curl -s https://raw.githubusercontent.com/pwwang/immunopipe-AdrienneML-2020/master/prepare-data.sh | bash /dev/stdin tests/data/prepared false

.PHONY: local api test
