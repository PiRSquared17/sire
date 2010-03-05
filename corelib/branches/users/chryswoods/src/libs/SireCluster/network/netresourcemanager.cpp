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

#include "netresourcemanager.h"
#include "communicator.h"
#include "nodemessages.h"

#include "SireCluster/resources/backend.h"

using namespace SireCluster;
using namespace SireCluster::network;
using namespace SireCluster::resources;

void NetResourceManager::init()
{}

void NetResourceManager::end()
{}

QUuid NetResourceManager::reserveResource(int expires)
{
    QHash<QUuid,quint64> msgs = Communicator::broadcast( ReserveRequest() );
    Communicator::awaitSent(msgs);

    return QUuid();
}

QUuid NetResourceManager::reserveResource(const QString &description,
                                          int expires)
{
    return QUuid();
}

QList<QUuid> NetResourceManager::reserveResources(int n, int expires)
{
    return QList<QUuid>();
}

QList<QUuid> NetResourceManager::reserveResources(const QString &description,
                                                  int n, int expires)
{
    return QList<QUuid>();
}

ActiveBackend NetResourceManager::collectReservation(const QUuid &uid)
{
    return ActiveBackend();
}

QHash<QUuid,ActiveBackend>
NetResourceManager::collectReservation(const QList<QUuid> &uids)
{
    return QHash<QUuid,ActiveBackend>();
}

void NetResourceManager::releaseReservation(const QUuid &uid)
{}

void NetResourceManager::releaseReservation(const QList<QUuid> &uids)
{}
