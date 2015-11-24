# pip uninstall molpher
python setup.py clean --all
# python setup.py build_ext
python setup.py bdist_wheel
python setup.py test
#pip install .
