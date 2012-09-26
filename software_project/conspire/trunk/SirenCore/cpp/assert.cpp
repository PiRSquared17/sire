/********************************************\
  *
  *  Siren - C++ metaobject library
  *
  *  Copyright (C) 2011  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#include "Siren/siren_assert.h"
#include "Siren/siren.hpp"

using namespace Siren;
using namespace Siren::Qt;

REGISTER_SIREN_CLASS( Siren::assertation_error )

assertation_error::assertation_error() : Exception()
{}

assertation_error::assertation_error(const char *condition, CODELOC_ARGS)
                : Exception( Siren::tr(
            "Assertation failed: %1").arg(condition), CODELOC_PASS_ARGS)
{}

assertation_error::assertation_error(const assertation_error &other) : Exception(other)
{}

assertation_error::~assertation_error() throw()
{}

void assertation_error::throwSelf() const
{
    throw *this;
}
