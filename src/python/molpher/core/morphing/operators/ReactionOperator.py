import molpher.swig_wrappers.core
from molpher.core import MolpherMol
from molpher.core._utils import shorten_repr
from molpher.core.morphing.operators import MorphingOperator


class ReactionOperator(molpher.swig_wrappers.core.ReactionOperator):

    def __repr__(self):
        return shorten_repr(self.__class__, self)

    def morph(self):
        ret = super(ReactionOperator, self).morph()
        if ret:
            ret.__class__ = MolpherMol
        return ret

    def getOriginal(self):
        ret = super(ReactionOperator, self).getOriginal()
        if ret:
            ret.__class__ = MolpherMol
        return ret

    @staticmethod
    def getDefaultOperators():
        ret = molpher.swig_wrappers.core.ReactionOperator.getDefaultOperators()
        if ret:
            for x in ret:
                x.__class__ = MorphingOperator
        return ret