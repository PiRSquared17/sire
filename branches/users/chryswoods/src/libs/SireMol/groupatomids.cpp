/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "groupatomids.h"
#include "atomidentifier.h"
#include "residentifier.h"
#include "chainidentifier.h"
#include "segidentifier.h"
#include "cgidentifier.h"

using namespace SireMol;

//////
////// Explicitly instantiate the GroupAtomID classes
//////

template class GroupAtomID<ResID,AtomID>;
template class GroupAtomID<ChainID,AtomID>;
template class GroupAtomID<SegID,AtomID>;
template class GroupAtomID<CGID,AtomID>;
