import sys


def shorten_repr(cls, self):
    vals = super(cls, self).__repr__().split(';')
    first = vals[0]
    second = vals[1].split(' at ')[1].strip('>').strip()
    return first + ' at ' + second


