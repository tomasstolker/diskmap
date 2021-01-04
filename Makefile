.PHONY: help pypi pypi-test docs coverage test clean

help:
	@echo "pypi - submit package to the PyPI server"
	@echo "pypi-test - submit package to the TestPyPI server"
	@echo "docs - generate Sphinx documentation"
	@echo "coverage - check code coverage"
	@echo "test - run test cases"
	@echo "clean - remove all artifacts"

pypi:
	python setup.py sdist bdist_wheel
	twine upload dist/*

pypi-test:
	python setup.py sdist bdist_wheel
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

docs:
	sphinx-apidoc -o docs diskmap
	cd docs/
	$(MAKE) -C docs clean
	$(MAKE) -C docs html

coverage:
	coverage run --source=diskmap -m pytest
	coverage report -m

test:
	pytest --cov=diskmap/

clean:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -rf {} +
	rm -rf docs/.ipynb_checkpoints
	rm -rf docs/*.fits
	rm -rf docs/*.dat
	rm -rf docs/_build/
	rm -rf build/
	rm -rf dist/
	rm -rf diskmap.egg-info/
	rm -rf .pytest_cache/
	rm -f .coverage
