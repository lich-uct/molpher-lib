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
 XXX NetBeans parser hack

 NetBeans parser tends to go into infinite loop when analyzing IDL macros
 from RCF (RCF_METHOD_*). Following macro wrapper hides problematic macros
 from NetBeans. Include this header into every IDL interface definition and
 wrap C(x) macro around RCF_METHOD_* macros.
 */

#ifdef C
#error Macro definition conflict.
#endif

#ifdef NETBEANS_HACK

#ifndef NOT_NETBEANS
#define C(x)
#else
#define C(x) x
#endif

#endif

#ifndef C
#define C(x) x
#endif
