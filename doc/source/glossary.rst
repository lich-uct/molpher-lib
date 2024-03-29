.. |nbsp| unicode:: 0xA0
   :trim:

.. _glossary:

Glossary
========

..  glossary::
    :sorted:

    exploration tree
        This is the data structure :term:`Molpher` generates as it progresses through :term:`chemical space`.
        It is represented in the API as an instance of
        the :class:`~molpher.core.ExplorationTree.ExplorationTree` class.

    SAScore
        Available from :class:`molpher.core.MolpherMol` instances via the :attr:`sascore` parameter.
        It is the synthetic availability score according to Ertl et al. [3]_. The higher the score,
        the more difficult the compound should be to synthesize. Compounds with SAScore of more
        than 6 are considered synthetically inaccessible and are discarded by the :class:`~.operations.FilterMorphsOper`
        operation and :meth:`~.core.ExplorationTree.filterMorphs` method by default.

    SAScore.dat
        The data file used for computation of the synthetic feasibility scores.

    similarity measure
        A similarity measure measures similarity between two molecules using a :term:`molecular fingerprint`.

    molecular fingerprint
        Molecular fingerprints are binary strings that encode certain structural features of a molecule.

    candidate morphs
        Morphs generated from the leaves of an :term:`exploration tree`, but not yet attached to it.

    source molecule
        The root of the exploration tree. All morphs in the tree are descendants of this molecule.

    target molecule
        Target molecule of the exploration process. The algorithm attempts to generate a path
        in chemical space leading from the :term:`source molecule` to this one.

    exploration parameters
        Those define how the :term:`exploration tree` behaves (i.e. what :term:`similarity measure` to
        use and what :term:`fingerprint <molecular fingerprint>`, how the :term:`candidate morphs` are filtered or
        how often the :term:`tree <exploration tree>` is pruned). See
        the description of the :class:`~molpher.core.ExplorationData.ExplorationData` class for details.

    morph
        Any molecule generated by :term:`Molpher` or *Molpher-lib*.

    chemical operators
        Also called *morphing operators*, are a set of predefined structural
        modifications that a molecule can undergo during :term:`molecular morphing`.

    tree snapshot
        A tree snapshot is a file on disk that contains all the information needed to create
        an :term:`exploration tree` instance.
        It can be saved and loaded using the :meth:`~molpher.swig_wrappers.core.ExplorationTree.save`
        and :meth:`~molpher.core.ExplorationTree.ExplorationTree.create` methods.

    XML template
        A file in XML format that can be used as a configuration file and is loaded the same way as a tree snapshot
        (using the :meth:`~molpher.core.ExplorationTree.ExplorationTree.create` method).

    tree operation
        Any operation that manipulates the :term:`exploration tree`. Can be defined
        by implementing the :class:`molpher.core.operations.TreeOperation` interface.

    morphing iteration
        Any sequence of operations that finishes with attaching new generation of morphs to a tree.
        For example, an iteration is committed when `ExplorationTree.extend()` is called.

    morph generation
        The morphs attached to the tree upon committing a :term:`morphing iteration`.

    chemical space path
        A consecutive sequence of :term:`morphs <morph>` created by iteratively applying
        a :term:`chemical operator <chemical operators>` on leaves in the :term:`exploration tree`.

    chemical space
        Chemical space has many definitions, but in this documentation we define it as
        *space spanned by all possible (i.e. energetically stable) molecules and chemical compounds* [1]_.

    morphing parameters
        A set of restrictions and rules the :term:`exploration tree` follows when
        certain :term:`operations <tree operation>` are performed.

    non-producing molecule
        A molecule that has not produced any morphs that would
        have a value of the objective function lower than itself.

    Molpher
        Software developed as a collaboration between the `SIRET Research Group <http://siret.ms.mff.cuni.cz>`_
        (`Faculty of Mathematics and Physics, Charles University in Prague, <http://www.mff.cuni.cz>`_),
        `Laboratory of Informatics and Chemistry <http://ich.vscht.cz>`_ (`UCT Prague <http://www.vscht.cz>`_)
        and `CZ-OPENSCREEN <http://www.openscreen.cz>`_ [2]_. The main goal of the project
        is to implement an effective tool for :term:`chemical space` exploration with :term:`molecular morphing`.
        Molpher itself is a standalone program and is available from `GitHub <https://github.com/siret/Molpher>`_.

    molecular morphing
        Molecular morphing is an atom-based *de novo* method of computer-aided drug design
        and it was first implemented in the
        :term:`Molpher` program |nbsp| [2]_. It uses a set of :term:`chemical operators` that modify
        structures of compounds to 'travel'
        through :term:`chemical space` and sample certain biologically interesting areas.

    selectors
        Located in the :mod:`~molpher.core.selectors` module, these serve the purpose of selecting
        various options for the morphing algorithm. See :class:`~molpher.core.ExplorationData.ExplorationData`
        for more information on how to affect the way :term:`morphs <morph>` are generated.

    high-throughput screening
        High-throughput screening (HTS) is a specialized method of chemical biology which relies on robotics
        to conduct many experiments in a very short time (usually in the order of tens of thousands a day) with
        the goal to find compounds that have the potential to affect a biological target and, thus, could
        become leads in a drug discovery project.

    virtual screening
        Virtual screening (VS) can be regarded as a computational alternative to :term:`HTS <high-throughput screening>`.
        The goal of virtual screening is the same: to find bioactive compounds, but in this case the compounds
        are represented *in silico* and saved in a database. Various computational techniques can then be used
        to automatically probe the database for structures that could have biological activity in the real world.
        Therefore, the point of VS is to create plausible hypotheses about the usefulness of the compounds in the
        database and discard all that very likely do not satisfy the requirements for a bioactive compounds in a
        given project.

.. [1] https://en.wikipedia.org/wiki/Chemical_space
.. [2] D. Hoksza, P. Škoda, M. Voršilák, and D. Svozil, “Molpher: a software framework for systematic chemical space exploration,” Journal of Cheminformatics, vol. 6, no. 1, p. 7, Mar. 2014. DOI: `10.1186/1758-2946-6-7 <http://dx.doi.org/10.1186/1758-2946-6-7>`_
.. [3] Peter Ertl, P & Schuffenhauer, Ansgar. (2009). Estimation of Synthetic Accessibility Score of Drug-Like Molecules Based on Molecular Complexity and Fragment Contributions. Journal of cheminformatics. 1. 8. DOI:`10.1186/1758-2946-1-8 <http://dx.doi.org/10.1186/1758-2946-1-8>`_.
