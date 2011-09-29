#ifndef SIREN_QT4SUPPORT_H
#define SIREN_QT4SUPPORT_H

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

// Definition of Qt classes that are required for Siren.
// These will (eventually!) be replaced with Siren's own classes,
// thereby removing the dependency on Qt

#include <QString>
#include <QStringList>
#include <QMutex>
#include <QMutexLocker>
#include <QByteArray>
#include <QDataStream>
#include <QSharedData>
#include <QSharedDataPointer>
#include <QHash>
#include <QVector>

#include <QDebug>

namespace Siren
{
    typedef QString String;
    typedef QStringList StringList;
    typedef QMutex Mutex;
    typedef QMutexLocker MutexLocker;
    typedef QByteArray ByteArray;
    typedef QDataStream DataStream;
    typedef QSharedData SharedData;
    
    template<class T>
    struct imp_shared_ptr
    {
        typedef QSharedDataPointer<T> Type;
    };
    
    template<class T>
    struct Vector
    {
        typedef QVector<T> Type;
    };

    template<class Key, class Value>
    struct Hash
    {
        typedef QHash<Key,Value> Type;
    };

    #if !defined(QT_NO_DEBUG_STREAM)
    QDebug sirenDebug() { return qDebug(); }
    #endif

} // end of namespace Siren

#endif // SIREN_QT4SUPPORT_H

