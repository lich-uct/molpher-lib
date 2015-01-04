/*
 Copyright (c) 2012 Petr Koupy

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

/*
 XXX Qt MOC hack (QTBUG-22829)

 There is a bug in Meta Object Compiler of Qt 4.8.0 (and possibly below) that
 prevents it to process certain macros in Boost 1.48 (and possibly above).
 Bug is being tracked by Qt developers, but it is unsure when this will be
 fixed. In the meantime, include this header to any header that defines
 Qt object and at the same time includes either directly or indirectly Boost
 headers. Note that this header have to be included by the very first #include
 directive in the header file.
 */

#ifdef QTMOC_HACK

#ifdef Q_MOC_RUN
#define BOOST_TT_HAS_OPERATOR_HPP_INCLUDED
#endif

#endif