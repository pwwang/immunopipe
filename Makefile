local:
	uv run pipen ref --pipeline "immunopipe:Immunopipe" \
		--destdir docs/processes \
		--replace-titles "Envs=Environment Variables" \
		--hide-sections Requires \
		--replace "https://pwwang.github.io/immunopipe=.." \
		--replace "${HOME}=~" \
		--replace "${USER}=pwwang"

api:
	uv run pipen ref --pipeline "immunopipe:Immunopipe" \
		--destdir docs/processes \
		--replace-titles "Envs=Environment Variables" \
		--hide-sections Requires

test:
	@uv run pytest tests -v

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

# Delete documentation versions from gh-pages using mike. Usage: make ddversion 1.0.0 2.0.0
# Supports wildcards with 'x': make ddversion 2.0.x (deletes all 2.0.*)
ddversion:
	@bash -c 'VERSIONS_INPUT="$(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))"; \
	if [ -z "$$VERSIONS_INPUT" ]; then \
		echo "Usage: make ddversion VERSION [VERSION...]"; \
		echo ""; \
		echo "Delete documentation versions from gh-pages using mike."; \
		echo ""; \
		echo "Arguments:"; \
		echo "  VERSION    Version number to delete (e.g., 1.0.0)"; \
		echo "             Supports wildcards with '\''x'\'' (e.g., 2.0.x deletes all 2.0.*)"; \
		echo ""; \
		echo "Examples:"; \
		echo "  make ddversion 1.0.0           # Delete version 1.0.0"; \
		echo "  make ddversion 1.0.0 2.0.0     # Delete versions 1.0.0 and 2.0.0"; \
		echo "  make ddversion 2.0.x           # Delete all versions starting with 2.0"; \
		echo "  make ddversion 2.x.x           # Delete all versions starting with 2"; \
		exit 0; \
	fi; \
	AVAILABLE_VERSIONS=$$(uv run mike list); \
	VERSIONS_TO_DELETE=""; \
	for VERSION_PATTERN in $$VERSIONS_INPUT; do \
		if echo "$$VERSION_PATTERN" | grep -q "x"; then \
			PATTERN=$$(echo "$$VERSION_PATTERN" | sed "s/x/.*/g"); \
			MATCHED=$$(echo "$$AVAILABLE_VERSIONS" | grep "^$$PATTERN" || true); \
			if [ -z "$$MATCHED" ]; then \
				echo "Error: No versions match pattern $$VERSION_PATTERN"; \
				echo "Available versions:"; \
				echo "$$AVAILABLE_VERSIONS"; \
				exit 1; \
			fi; \
			VERSIONS_TO_DELETE="$$VERSIONS_TO_DELETE $$MATCHED"; \
			MATCHED_DISPLAY=$$(echo "$$MATCHED" | tr "\n" "," | sed "s/,$$//"); \
			echo "Pattern $$VERSION_PATTERN matched: $$MATCHED_DISPLAY"; \
		else \
			if ! echo "$$AVAILABLE_VERSIONS" | grep -q "^$$VERSION_PATTERN$$"; then \
				echo "Error: Version $$VERSION_PATTERN not found in mike list"; \
				echo "Available versions:"; \
				echo "$$AVAILABLE_VERSIONS"; \
				exit 1; \
			fi; \
			VERSIONS_TO_DELETE="$$VERSIONS_TO_DELETE $$VERSION_PATTERN"; \
		fi; \
	done; \
	for VERSION in $$VERSIONS_TO_DELETE; do \
		echo "Deleting version $$VERSION..."; \
		uv run mike delete $$VERSION; \
	done; \
	echo "Pushing to gh-pages..."; \
	git push origin gh-pages; \
	echo "Successfully deleted versions"'

ddver: ddversion

# Prevent make from trying to interpret version arguments as targets
%:
	@:

.PHONY: local api test version ddversion ddver
