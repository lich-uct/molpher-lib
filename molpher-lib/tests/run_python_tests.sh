. ./lib/set_environ.sh
ulimit -c unlimited
python3 -m unittest -v tests/tests.py
