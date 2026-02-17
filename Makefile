local:
	pipen ref --pipeline "immunopipe:Immunopipe" \
		--destdir docs/processes \
		--replace-titles "Envs=Environment Variables" \
		--hide-sections Requires \
		--replace "https://pwwang.github.io/immunopipe=.." \
		--replace "${HOME}=~" \
		--replace "${USER}=pwwang"

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

version:
	@bash -c 'if [ -z "$(word 2,$(MAKECMDGOALS))" ]; then \
		CURRENT_VERSION=$$(grep "^__version__" immunopipe/version.py | sed "s/__version__ = \"\(.*\)\"/\1/"); \
		MAJOR=$$(echo $$CURRENT_VERSION | cut -d. -f1); \
		MINOR=$$(echo $$CURRENT_VERSION | cut -d. -f2); \
		PATCH=$$(echo $$CURRENT_VERSION | cut -d. -f3); \
		NEW_PATCH=$$((PATCH + 1)); \
		NEW_VERSION="$$MAJOR.$$MINOR.$$NEW_PATCH"; \
	else \
		NEW_VERSION="$(word 2,$(MAKECMDGOALS))"; \
	fi; \
	echo "Updating version to $$NEW_VERSION"; \
	sed -i "s/^version = .*/version = \"$$NEW_VERSION\"/" pyproject.toml; \
	sed -i "s/^__version__ = .*/__version__ = \"$$NEW_VERSION\"/" immunopipe/version.py; \
	LAST_MERGE=$$(git log --grep="Merge remote-tracking branch '"'"'origin/master'"'"' into dev" --format="%H" -n 1 2>/dev/null || echo ""); \
	if [ -z "$$LAST_MERGE" ]; then \
		COMMITS=$$(git log --pretty=format:"- %s" HEAD); \
	else \
		COMMITS=$$(git log --pretty=format:"- %s" $$LAST_MERGE..HEAD); \
	fi; \
	if [ -n "$$COMMITS" ]; then \
		tail -n +3 docs/CHANGELOG.md > docs/CHANGELOG.md.tmp; \
		head -n 2 docs/CHANGELOG.md > docs/CHANGELOG.md.new; \
		printf "\n## %s\n\n%s\n\n" "$$NEW_VERSION" "$$COMMITS" >> docs/CHANGELOG.md.new; \
		cat docs/CHANGELOG.md.tmp >> docs/CHANGELOG.md.new; \
		mv docs/CHANGELOG.md.new docs/CHANGELOG.md; \
		rm -f docs/CHANGELOG.md.tmp; \
	fi; \
	echo "Version updated to $$NEW_VERSION"'

.PHONY: local api test version
