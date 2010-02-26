/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIRECLUSTER_NETWORK_COMMUNICATOR_H
#define SIRECLUSTER_NETWORK_COMMUNICATOR_H

#include <QUuid>
#include <QList>

#include <boost/function.hpp>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireCluster
{
namespace network
{

class HostInfo;
class Message;
class Envelope;

/** A Communicator provides the ability to 
    send "Message" objects to another process.
    
    @author Christopher Woods
*/
class SIRECLUSTER_EXPORT Communicator
{
public:
    static const HostInfo& getLocalInfo();
    static HostInfo getHostInfo(const QUuid &host_uid);

    static void send(const Message &message, const QUuid &recipient);
    static void send(const Message &message, const QList<QUuid> &recipients);
    
    static void broadcast(const Message &message);
  
    static void received(const Envelope &message);  
    static void received(const Message &message);
    static void received(const QByteArray &message);

    static void addNeighbour(const HostInfo &hostinfo,
                             boost::function<void(const QByteArray&)> send_function);

    static void addNeighbour(const HostInfo &hostinfo,
                             const QUuid &router,
                             boost::function<void(const QByteArray&)> send_function,
                             int njumps=1);
};

} // end of namespace network
} // end of namespace SireCluster

SIRE_END_HEADER

#endif
