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

#include "SireError/exception.h"
#include "SireStream/datastream.h"

#include "SireError/errors.h"

#include "getbacktrace.h"

#include <QDataStream>

#include <QDebug>

using namespace SireError;
using namespace SireStream;

/** Construct a null exception */
exception::exception()
{}

/** Constructor.
    \param error The error message associated with this exception.
    \param place From where in the code this exception was thrown. This
                 is supplied automatically via the 'CODELOC' macro.
*/
exception::exception(QString error, QString place) : err(error), plce(place)
{
    pidstr = QObject::tr("unidentified thread");
    bt = getBackTrace();
}

/** Copy constructor */
exception::exception(const exception &other)
          : std::exception(other), err(other.err), plce(other.plce), bt(other.bt)
{}

/** Destructor */
exception::~exception() throw()
{}

/** Return a clone of this exception */
exception* exception::clone() const
{
    //get the ID number of this type
    int id = QMetaType::type( this->what() );

    if ( id == 0 or not QMetaType::isRegistered(id) )
        throw SireError::unknown_type(QObject::tr(
            "The exception with type \"%1\" does not appear to have been "
            "registered with QMetaType. It cannot be cloned! (%2, %3)")
                .arg(this->what()).arg(id).arg(QMetaType::isRegistered(id)), CODELOC);

    return static_cast<exception*>( QMetaType::construct(id,this) );
}

/** Pack this exception into a binary QByteArray - this packs the exception
    with its type information */
QByteArray exception::pack() const
{
    QByteArray data;
    
    //reserve 128K of space for the exception (should be way
    //more than enough!)
    data.reserve( 128 * 1024 );
    
    QDataStream ds(&data, QIODevice::WriteOnly);
    
    qDebug() << "Saving an exception of type" << this->what();
    
    //get the ID number of this type
    int id = QMetaType::type( this->what() );

    if ( id == 0 or not QMetaType::isRegistered(id) )
        throw SireError::unknown_type(QObject::tr(
            "The exception with type \"%1\" does not appear to have been "
            "registered with QMetaType. It cannot be streamed! (%2, %3)")
                .arg(this->what()).arg(id).arg(QMetaType::isRegistered(id)), CODELOC);

    //save the object type name
    ds << QString(this->what());

    //use the QMetaType streaming function to save this table
    if (not QMetaType::save(ds, id, this))
        throw SireError::program_bug(QObject::tr(
            "There was an error saving the exception of type \"%1\". "
            "Has the programmer added a RegisterMetaType for this exception?")
                .arg(this->what()), CODELOC);

    return data;
}

/** Unpack an exception from the raw data in the passed bytearray */
boost::shared_ptr<SireError::exception> exception::unpack(const QByteArray &data)
{
    QDataStream ds(data);
    
    //read the type of exception first
    QString type_name;
    
    ds >> type_name;
    
    qDebug() << "Loading an exception of type" << type_name;
    
    //get the type that represents this name
    int id = QMetaType::type( type_name.toLatin1().constData() );

    if ( id == 0 or not QMetaType::isRegistered(id) )
        throw SireError::unknown_type( QObject::tr(
              "Cannot deserialise an exception of type \"%1\". "
              "Ensure that the library or module containing "
              "this exception has been loaded and that it has been registered "
              "with QMetaType.").arg(type_name), CODELOC );
    
    //construct an exception of this type
    boost::shared_ptr<exception> ptr(  
                        static_cast<exception*>(QMetaType::construct(id,0)) );
                        
    //load the object from the datastream
    if ( not QMetaType::load(ds, id, ptr.get()) )
        throw SireError::program_bug(QObject::tr(
            "There was an error loading the exception of type \"%1\"")
                 .arg(type_name), CODELOC);

    return ptr;
}

/** Unpack the exception contained in the raw binary data 'data', and then
    throw that exception */
void exception::unpackAndThrow(const QByteArray &data)
{
    boost::shared_ptr<exception> e = exception::unpack(data);
    e->throwSelf();
}

/** Return a pretty string representation of all of the information stored
    in the exception, suitable for printing to the screen or to a log file. */
QString exception::toString() const throw()
{
    return QObject::tr("Exception '%1' thrown by process '%2'.\n"
                       "%3\n"
                       "Thrown from %4\n"
                       "__Backtrace__\n"
                       "%5\n"
                       "__EndTrace__")
             .arg(what()).arg(pid()).arg(why()).arg(where()).arg(bt.join("\n"));
}

/** Return the error message associated with this exception.
    \return The error message associated with this exception.
*/
QString exception::error() const throw()
{
    return err;
}

/** Return the location (function, line, file) from which this exception
    was originally thrown. The amount of information available will depend
    on the compiler used to compile the program, and come from the 'CODELOC'
    macro.
    \return Description of from where the exception was thrown.
*/
QString exception::from() const throw()
{
    return plce;
}

/** Return the function backtrace when the exception was constructed */
QStringList exception::trace() const throw()
{
    return bt;
}

/** Overloaded functions to give logical names... */
QString exception::where() const throw()
{
    return this->from();
}

/** Overloaded functions to give logical names... */
QString exception::why() const throw()
{
    return this->error();
}

/** Return the identifying string for the thread on the process that threw
    this exception. You should call the static SireError::exception::setPID
    function in each thread for the result of this call to be intelligable */
QString exception::pid() const throw()
{
    return pidstr;
}
