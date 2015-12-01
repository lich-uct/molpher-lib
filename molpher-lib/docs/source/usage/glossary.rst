.. _glossary:

Glossary
========

The glossary of important terms used in the `documentation <index>`.

..  glossary::

    exploration tree
        Definition of the exploration tree term.

    Molpher-backend project
        The original Molpher project (located in the :file:`backend/` folder in the repository root).
        It contains the implementation of Molpher backend server.
        It needs to be compiled as a static library so that the `Molpher-lib project` can
        link against it during compilation.

    Molpher-lib project
        The implementation of Molpher-lib. It contains the :file:`include/` directory with the header files
        that define the C++ API.

    common headers
        Located in :file:`commons/` folder in the repository root, these header files are shared
        between `Molpher-backend project`, `Molpher-lib project` and the frontend.
