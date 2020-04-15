.PHONY: help pypi pypi-test clean

help:
	@echo "pypi - submit package to the PyPI server"
	@echo "pypi-test - submit package to the TestPyPI server"
	@echo "clean - remove all artifacts"

pypi:
	python setup.py sdist bdist_wheel
	twine upload dist/*

pypi-test:
	python setup.py sdist bdist_wheel
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

clean:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -rf {} +
	rm -rf build/
	rm -rf dist/
	rm -rf diskmap.egg-info/
	rm -f .coverage
