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

#include "ReducerFactory.h"

#include "KamadaKawaiReducer.h"
#include "PcaReducer.h"

DimensionReducer *ReducerFactory::Create(DimRedSelector selector)
{
    switch (selector) {
    case DR_KAMADAKAWAI:
        return new KamadaKawaiReducer();
    case DR_PCA:
        return new PcaReducer();
    default: // use PCA as default
        return new PcaReducer();
    }
}

void ReducerFactory::Recycle(DimensionReducer *reducer)
{
    delete reducer;
}
