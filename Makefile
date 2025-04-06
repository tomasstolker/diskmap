.PHONY: help pypi pypi-test docs coverage test clean

help:
	@echo "pypi - submit to PyPI server"
	@echo "pypi-check - check the distribution for PyPI"
	@echo "pypi-test - submit to TestPyPI server"
	@echo "docs - generate Sphinx documentation"
	@echo "coverage - check code coverage"
	@echo "test - run unit tests"
	@echo "clean - remove artifacts"

pypi:
	python -m build
	twine upload dist/*

pypi-check:
	python -m build
	twine check dist/*

pypi-test:
	python -m build
	twine upload -r testpypi dist/*

docs:
	rm -rf docs/api
	sphinx-apidoc -o docs diskmap
	cd docs/
	$(MAKE) -C docs clean
	$(MAKE) -C docs html

coverage:
	coverage run --source=diskmap -m pytest
	coverage report -m

test:
	pytest --cov=diskmap/ --cov-report=xml

clean:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -rf {} +
	rm -rf docs/.ipynb_checkpoints
	rm -rf docs/*.fits
	rm -rf docs/*.dat
	rm -rf docs/_build/
	rm -rf docs/api
	rm -rf build/
	rm -rf dist/
	rm -rf diskmap.egg-info/
	rm -rf .pytest_cache/
	rm -f .coverage
	rm -f coverage.xml
