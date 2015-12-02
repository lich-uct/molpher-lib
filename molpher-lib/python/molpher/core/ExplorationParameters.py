import molpher

class ExplorationParameters(molpher.swig_wrappers.core.ExplorationParameters):
    """
    Houses morphing parameters and can be used to initialize a `molpher.core.ExplorationTree` instance.

    :param kwargs: morphing parameters

    The names, default values and descriptions of the morphing parameters are given in the following table:

    ..  csv-table:: a title
        :header: "Name", "Default", "Description"
        :widths: 10, 10, 50

        "source", "''", "SMILES string of the source molecule."
        "target", "''", "SMILES string of the target molecule."

    ..  todo::
        Finish the table.


    :type kwargs: :py:class:`dict`

    """

    class UnknownParameterException(Exception):
        """
        Indicates that a parameter in `ExplorationParameters.params` is unknown.
        """

        pass

    def __init__(self, **kwargs):
        super(ExplorationParameters, self).__init__()
        self._SETTERS_MAP = {
            'source' : self.setSourceMol
            , 'target' : self.setTargetMol
            , 'operators' : self.setChemOperators
            , 'fingerprint' : self.setFingerprint
            , 'similarity' : self.setSimilarityCoef
            , 'weight_min' : self.setMinAcceptableMolecularWeight
            , 'weight_max' : self.setMaxAcceptableMolecularWeight
            , 'accept_min' : self.setCntCandidatesToKeep
            , 'accept_max' : self.setCntCandidatesToKeepMax
            , 'far_produce' : self.setCntMorphs
            , 'close_produce' : self.setCntMorphsInDepth
            , 'far_close_threshold' : self.setDistToTargetDepthSwitch
            , 'max_morphs_total' : self.setCntMaxMorphs
            , 'non_producing_survive' : self.setItThreshold
        }
        if kwargs:
            self._update_instance(kwargs)

    def _update_instance(self, options):
        for key in options:
            if key in self._SETTERS_MAP:
                self._SETTERS_MAP[key](options[key])
            else:
                raise self.UnknownParameterException('Unknown option: {0}'.format(key))


    @property
    def params(self):
        """
        Holds a dictionary of current parameter values for this instance.
        Also, a new dictionary of parameters can be assigned to change them.

        :return: a dictionary of parameters
        :rtype: :py:class:`dict`
        """

        return {
            'source' : self.getSourceMol().getSMILES()
            , 'target' : self.getTargetMol().getSMILES()
            , 'operators' : self.getChemOperators()
            , 'fingerprint' : self.getFingerprint()
            , 'similarity' : self.getSimilarityCoef()
            , 'weight_min' : self.getMinAcceptableMolecularWeight()
            , 'weight_max' : self.getMaxAcceptableMolecularWeight()
            , 'accept_min' : self.getCntCandidatesToKeep()
            , 'accept_max' : self.getCntCandidatesToKeepMax()
            , 'far_produce' : self.getCntMorphs()
            , 'close_produce' : self.getCntMorphsInDepth()
            , 'far_close_threshold' : self.getDistToTargetDepthSwitch()
            , 'max_morphs_total' : self.getCntMaxMorphs()
            , 'non_producing_survive' : self.getItThreshold()
        }

    @params.setter
    def params(self, options):
        self._update_instance(options)

    @property
    def is_valid(self):
        """
        Checks if this instance represents valid parameters.

        ..  todo::
            Document what exactly is checked.

        :return: validity as :py:class:`True` or :py:class:`False`
        :rtype: :py:class:`bool`
        """

        return self.valid()


