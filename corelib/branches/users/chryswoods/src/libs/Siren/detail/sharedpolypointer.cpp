/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "sharedpolypointer.hpp"

#include "Siren/errors.h"

#include <QDebug>

using namespace Siren;
using namespace Siren::detail;

/** Destructor */
SharedPolyPointerBase::~SharedPolyPointerBase()
{}

/** Save the contents of the pointer to the stream

    \throw SireError::unknown_type
    \throw SireError::program_bug
*/
void SharedPolyPointerBase::save(QDataStream &ds, const QString &objname,
                                 const void *data)
{
    if (objname.isEmpty() or data == 0)
    {
        ds << QString::null;
    }
    else
    {
        //get the ID number of this type
        int id = QMetaType::type( objname.toLatin1().constData() );

        if ( id == 0 or not QMetaType::isRegistered(id) )
            throw Siren::unknown_type(QObject::tr(
                "The object with type \"%1\" does not appear to have been "
                "registered with QMetaType. It cannot be streamed! (%2, %3)")
                    .arg(objname).arg(id).arg(QMetaType::isRegistered(id)), CODELOC);

        //save the object type name
        ds << QString(objname);

        //use the QMetaType streaming function to save this table
        if (not QMetaType::save(ds, id, data))
            throw Siren::program_bug(QObject::tr(
                "There was an error saving the object of type \"%1\". "
                "Has the programmer remembered to add a RegisterMetaType for this class?")
                    .arg(objname), CODELOC);
    }
}

/** Read a SharedPolyPointer from the passed QDataStream

    \throw SireError::unknown_type
    \throw SireError::program_bug
*/
void* SharedPolyPointerBase::read(QDataStream &ds)
{
    //read the type name
    QString type_name;
    ds >> type_name;

    if (type_name.isNull())
    {
        //this is a null pointer
        return 0;
    }
    else
    {
        //get the type that represents this name
        int id = QMetaType::type( type_name.toLatin1().constData() );

        if ( id == 0 or not QMetaType::isRegistered(id) )
            throw Siren::unknown_type( QObject::tr(
                  "Cannot deserialise an object of type \"%1\". "
                  "Ensure that the library or module containing "
                  "this type has been loaded and that it has been registered "
                  "with QMetaType.").arg(type_name), CODELOC );

        //create a default-constructed object of this type
        void *ptr = QMetaType::construct(id,0);

        if (not ptr)
            throw Siren::program_bug( QObject::tr(
                    "Could not create an object of type \"%1\" despite "
                    "this type having been registered with QMetaType. This is "
                    "a program bug!!!").arg(type_name), CODELOC );

        //load the object from the datastream
        if ( not QMetaType::load(ds, id, ptr) )
            throw Siren::program_bug(QObject::tr(
                "There was an error loading the object of type \"%1\"")
                     .arg(type_name), CODELOC);

        return ptr;
    }
}

void SharedPolyPointerBase::throwInvalidCast(const QString &oldtype,
                                             const QString &newtype)
{
    throw Siren::invalid_cast( QObject::tr(
                "Cannot cast from a SharedPolyPointer of type \"%1\""
                "to a SharedPolyPointer of type \"%2\"")
                    .arg( oldtype, newtype ), CODELOC );
}