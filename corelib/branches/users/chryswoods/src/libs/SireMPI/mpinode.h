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

#ifndef SIREMPI_MPINODE_H
#define SIREMPI_MPINODE_H

#include <QUuid>

#include <boost/shared_ptr.hpp>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMPI
{

class MPINodes;
class MPIWorker;
class MPIPromise;

namespace detail
{
class MPINodeData;
class MPINodesData;
}

bool exec_running();

int exec(int argc, const char **argv);
void bg_exec(int argc, const char **argv);

void shutdown();

/** This class provides a handle to represent a node that
    is made available via MPI. It provides information about
    the node (number of available processors and memory etc.)
    and provides the ID necessary to allow jobs to be assigned
    to it.
    
    @author Christopher Woods
*/
class SIREMPI_EXPORT MPINode
{

friend class MPINodes;

friend class detail::MPINodeData;
friend class detail::MPINodesData;

public:
    MPINode();
    
    MPINode(const MPINode &other);

    ~MPINode();
    
    MPINode& operator=(const MPINode &other);
    
    bool operator==(const MPINode &other) const;
    bool operator!=(const MPINode &other) const;
    
    const QUuid& UID() const;

    bool isMaster() const;

    bool isNull() const;

    bool isBusy() const;
        
    bool wasAborted() const;
        
    void wait();    
    bool wait(int time);

    void returnNode();

    MPINodes communicator() const;
    int rank() const;

    MPIPromise start(const MPIWorker &worker);
    void stop();
    void abort();
    
    double getProgress();
    
    MPIPromise getInterimResult();
    MPIPromise getResult();

protected:
    MPINode(const MPINodes &communicator, int rank, bool is_master);
    MPINode(const boost::shared_ptr<detail::MPINodeData> &data);

private:
    /** Shared pointer to the data for this node - this
        is used to hide the implementation to ensure binary
        compatibility between versions of Sire that have MPI
        and those that don't */
    boost::shared_ptr<detail::MPINodeData> d;
};

}

SIRE_EXPOSE_CLASS( SireMPI::MPINode )

SIRE_END_HEADER

#endif
