// ---------------------------------------------------------------------
// Copyright (C) 2015 Chris Garry
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>
// ---------------------------------------------------------------------


#ifndef PIPP_UTF8_H
#define PIPP_UTF8_H

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cstdio>
#include <string>

#if defined(__unix__) || defined(OS_OSX)
#include <sys/param.h>		// define or not BSD macro
#endif

// 64-bit fseek for various platforms
#if defined(__GLIBC__) || defined(__gnu_hurd__)
#define fseek64 fseeko64  // GNU
#define ftell64 ftello64  // GNU
#elif defined(_WIN32)
#define fseek64 _fseeki64  // Windows
#define ftell64 _ftelli64  // Windows
#else // POSIX
#define fseek64 fseeko  // OS X, DragonFly BSD, FreeBSD, OpenBSD, NetBSD, musl
#define ftell64 ftello  // OS X, DragonFly BSD, FreeBSD, OpenBSD, NetBSD, musl
#endif

#endif  // PIPP_UTF8_H
