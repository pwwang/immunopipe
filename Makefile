local:
	pipen ref --pipeline "immunopipe:Immunopipe" \
		--destdir docs/processes \
		--replace-titles "Envs=Environment Variables" \
		--hide-sections Requires \
		--replace "https://pwwang.github.io/immunopipe=.." \
		--replace "${HOME}=~" \
		--replace "${USER}=user"

api:
	pipen ref --pipeline "immunopipe:Immunopipe" \
		--destdir docs/processes \
		--replace-titles "Envs=Environment Variables" \
		--hide-sections Requires

test:
	@poetry run pytest tests -v

test-data:
	curl -s \
		https://raw.githubusercontent.com/pwwang/immunopipe-example/master/prepare-data.sh | \
		bash /dev/stdin tests/running/data/prepared false

.PHONY: local api test
