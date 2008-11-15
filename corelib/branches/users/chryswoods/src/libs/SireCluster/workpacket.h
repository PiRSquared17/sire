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

#ifndef SIRECLUSTER_WORKPACKET_H
#define SIRECLUSTER_WORKPACKET_H

#include "sireglobal.h"

#include <QSharedData>

#include "SireBase/sharedpolypointer.hpp"

SIRE_BEGIN_HEADER

namespace SireCluster
{
class WorkPacket;
class WorkPacketBase;
class ErrorPacket;
class WorkTest;
}

QDataStream& operator<<(QDataStream&, const SireCluster::WorkPacketBase&);
QDataStream& operator>>(QDataStream&, SireCluster::WorkPacketBase&);

QDataStream& operator<<(QDataStream&, const SireCluster::ErrorPacket&);
QDataStream& operator>>(QDataStream&, SireCluster::ErrorPacket&);

QDataStream& operator<<(QDataStream&, const SireCluster::WorkPacket&);
QDataStream& operator>>(QDataStream&, SireCluster::WorkPacket&);

QDataStream& operator<<(QDataStream&, const SireCluster::WorkTest&);
QDataStream& operator>>(QDataStream&, SireCluster::WorkTest&);

namespace SireCluster
{

/** This is the base class of all WorkPackets. A WorkPacket
    contains all of the code and input data for a piece of work,
    and also contains space to return the output and current
    progress 
    
    @author Christopher Woods
*/
class SIRECLUSTER_EXPORT WorkPacketBase : public QSharedData
{

friend QDataStream& ::operator<<(QDataStream&, const WorkPacketBase&);
friend QDataStream& ::operator>>(QDataStream&, WorkPacketBase&);

public:
    WorkPacketBase();
    WorkPacketBase(const WorkPacketBase &other);
    
    virtual ~WorkPacketBase();
    
    virtual bool shouldPack() const;
    virtual int approximatePacketSize() const;
    
    void abort();
    
    void runChunk();
    
    float progress() const;
    
    virtual bool hasFinished() const=0;
    
    virtual bool isError() const;
    virtual void throwError() const;
    
    bool wasAborted() const;
    
    virtual WorkPacketBase* clone() const=0;
    
    static const char* typeName()
    {
        return "SireCluster::WorkPacketBase";
    }
    
    virtual const char* what() const=0;

    template<class T>
    bool isA() const
    {
        return dynamic_cast<const T*>(this) != 0;
    }
    
    template<class T>
    const T& asA() const
    {
        return dynamic_cast<const T&>(*this);
    }

protected:
    virtual float chunk()=0;

    WorkPacketBase& operator=(const WorkPacketBase&);
    
private:
    /** The current progress of the work */
    float current_progress;
    
    /** Whether or not the work was aborted */
    bool was_aborted;
};

/** This class is the generic holder for all work packets.
    
    @author Christopher Woods
*/
class SIRECLUSTER_EXPORT WorkPacket
{

friend QDataStream& ::operator<<(QDataStream&, const WorkPacket&);
friend QDataStream& ::operator>>(QDataStream&, WorkPacket&);

public:
    WorkPacket();
    WorkPacket(const WorkPacketBase &work);
    WorkPacket(const WorkPacket &other);
    
    ~WorkPacket();
    
    WorkPacket& operator=(const WorkPacket &other);

    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId<WorkPacket>() );
    }
    
    const char* what() const
    {
        return WorkPacket::typeName();
    }

    QByteArray pack() const;
    
    static WorkPacket unpack(const QByteArray &data);

    bool shouldPack() const;

    bool isNull() const;

    void abort();

    void runChunk() throw();
    
    float progress() const;
    
    bool wasAborted() const;
    
    bool hasFinished() const;

    bool isError() const;
    void throwError() const;

    const WorkPacketBase& base() const;

    template<class T>
    bool isA() const;
    
    template<class T>
    const T& asA() const;

private:
    void setError(const SireError::exception &e) throw();

    /** Implicitly shared pointer to the work */
    SireBase::SharedPolyPointer<WorkPacketBase> d;
};

/** This is a packet that contains an error. This is returned
    if something went wrong while running a WorkPacket
    
    @author Christopher Woods
*/
class SIRECLUSTER_EXPORT ErrorPacket : public WorkPacketBase
{

friend QDataStream& ::operator<<(QDataStream&, const ErrorPacket&);
friend QDataStream& ::operator>>(QDataStream&, ErrorPacket&);

public:
    ErrorPacket();
    ErrorPacket(const SireError::exception &e);

    ErrorPacket(const ErrorPacket &other);
    
    ~ErrorPacket();
    
    ErrorPacket& operator=(const ErrorPacket &other);
    
    ErrorPacket* clone() const
    {
        return new ErrorPacket(*this);
    }
    
    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId<ErrorPacket>() );
    }
    
    const char* what() const
    {
        return ErrorPacket::typeName();
    }
    
    int approximatePacketSize() const;
    
    bool isError() const;
    void throwError() const;
    
    bool hasFinished() const;

protected:
    float chunk();
    
private:
    /** A binary representation of the error */
    QByteArray error_data;
};

/** This is a small test packet that can be used to test
    that everything is working */
class SIRECLUSTER_EXPORT WorkTest : public WorkPacketBase
{

friend QDataStream& ::operator<<(QDataStream&, const WorkTest&);
friend QDataStream& ::operator>>(QDataStream&, WorkTest&);

public:
    WorkTest();
    WorkTest(int start, int end, int step=1);
    WorkTest(const WorkTest &other);
    
    ~WorkTest();
    
    WorkTest& operator=(const WorkTest &other);
    
    WorkTest* clone() const
    {
        return new WorkTest(*this);
    }
    
    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId<WorkTest>() );
    }
    
    const char* what() const
    {
        return WorkTest::typeName();
    }
    
    int approximatePacketSize() const;
    
    bool hasFinished() const;

protected:
    float chunk();
    
private:
    /** The range to loop over */
    qint32 current, start, end, step;
};

}

Q_DECLARE_METATYPE( SireCluster::WorkPacket )
Q_DECLARE_METATYPE( SireCluster::ErrorPacket )
Q_DECLARE_METATYPE( SireCluster::WorkTest )

SIRE_EXPOSE_CLASS( SireCluster::WorkPacketBase )
SIRE_EXPOSE_CLASS( SireCluster::WorkPacket )
SIRE_EXPORT_CLASS( SireCluster::ErrorPacket )
SIRE_EXPOSE_CLASS( SireCluster::WorkTest )

SIRE_END_HEADER

#endif
