/*
 Copyright (c) 2012 Marek Mikes

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

#include <cassert>

#include "GlobalObjectsHolder.h"
#include "components/Bookmarks.h"

GlobalObjectsHolder gGlobalObjectsHolder;

void GlobalObjectsHolder::Add(Bookmarks *bookmarks)
{
    // bookmarks should be created only once
    assert(NULL == mBookmarks);

    mBookmarks = bookmarks;
}

void GlobalObjectsHolder::Remove(Bookmarks *bookmarks)
{
    assert(mBookmarks == bookmarks);

    mBookmarks = NULL;
}

Bookmarks *GlobalObjectsHolder::GetBookmarks() const
{
    // bookmarks should be already created
    assert(NULL != mBookmarks);

    return mBookmarks;
}
