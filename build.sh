rm -rf build *.egg* dist cdftpy-0.0.4
python setup.py sdist bdist_wheel
twine upload --skip-existing dist/*