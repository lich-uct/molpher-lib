/*
 Copyright (c) 2012 Vladimir Fiklik
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

#include <string>
#include <map>

#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include "global_types.h"

class PasswordCache
{
public:
    static void CachePassword(JobId jobId, std::string &password);
    static bool ResolvePassword(JobId jobId, std::string &password);

private:
    PasswordCache();
    PasswordCache(const PasswordCache &other); // singleton - do not implement
    PasswordCache &operator=(const PasswordCache &other); // singleton - do not implement

    static PasswordCache instance; // eager singleton instance

    typedef boost::mutex Guard;
    typedef boost::unique_lock<Guard> Lock;
    Guard mMapGuard;

    typedef std::map<JobId, std::string> PasswordMap;
    PasswordMap mPasswordMap;
};
