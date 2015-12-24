Introduction
============

Before we delve more into the technical realm and elaborate on the details of the library itself,
let us first introduce the ideas behind `Molpher` and introduce some important terminology.

What is Molpher?
----------------

In early stages of a drug discovery endeavor a large library of compunds often needs to be screened
for molecules that have potential to become a new lead in the project. This
is usually done in a laboratory setting using a technique called
`high-throughput screening (HTS) <https://en.wikipedia.org/wiki/High-throughput_screening>`_.
This can be a very effective method, but is still very limiting in the number
of compounds that can be screened when we consider the incredible vastness of `chemical space`.
Therefore, a lot of computational tools have been developed for the purpose
of prioritizing certain compounds for the laboratory experiments (such as HTS).

`Molpher` is one of those tools. In short, it generates a large computational
library of putatively active compounds by connecting known active molecules
by a `chemical space path`. Like an actual path, a `chemical space path`
has a start and a destination. In `chemical space` this translates to the concept
of a `source molecule` and a `target molecule`. `Molpher` works on the premise
that other potentially active compounds might lie on a path between
those two, if they are both active on the same target. We call the process
of searching for such a path *molecular morphing*.

What is Molpher-lib?
--------------------

`Molpher` as a software was originally directed at users (i. e. chemists) that might use it
in their day to day practice. However, further research showed that there are many
possibile implementations of the algorithm and that the exploration can be driven
by vastly different means. Thus came the idea to develop a more flexible
and extensible solution so that anyone can easily use `Molpher` to implement
their own ideas and easily make `Molpher` fit their own needs.

To make it as easy as possible for new users to get up to speed with the library,
we have prepared a :doc:`tutorial <usage/tutorial>` in the :doc:`next section <usage/index>`.