/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIRESYSTEM_SYSIDENTIFIER_H
#define SIRESYSTEM_SYSIDENTIFIER_H

#include "sysid.h"

#include <boost/shared_ptr.hpp>

namespace SireSystem
{

/** This is a generic holder for any SysID class! 

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT SysIdentifier : public SysID
{
public:
    SysIdentifier();
    SysIdentifier(const SysID &sysid);
    SysIdentifier(const SysIdentifier &other);
    
    ~SysIdentifier();
    
    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId<SysIdentifier>() );
    }
    
    const char* what() const
    {
        return SysIdentifier::typeName();
    }
    
    SysIdentifier* clone() const
    {
        return new SysIdentifier(*this);
    }
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const SysID& base() const;
    
    SysIdentifier& operator=(const SysIdentifier &other);
    SysIdentifier& operator=(const SysID &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const SysIdentifier &other) const;
    bool operator!=(const SysIdentifier &other) const;
    
    bool operator==(const SysID &other) const;
    bool operator!=(const SysID &other) const;

private:
    /** Pointer to the SysID */
    boost::shared_ptr<SysID> d;
};

}

Q_DECLARE_METATYPE(SireSystem::SysIdentifier);

#endif
