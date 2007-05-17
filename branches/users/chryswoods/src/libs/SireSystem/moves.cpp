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

#include "moves.h"
#include "simsystem.h"

#include "SireMol/moleculegroupid.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireBase;
using namespace SireMol;
using namespace SireStream;

///////////
/////////// Implementation of MovesBase
///////////

static const RegisterMetaType<MovesBase> r_movesbase(MAGIC_ONLY,
                                                     "SireSystem::MovesBase");

/** Serialise to a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds, const MovesBase&)
{
    writeHeader(ds, r_movesbase, 0);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, MovesBase&)
{
    VersionID v = readHeader(ds, r_movesbase);

    if (v != 0)
        throw version_error(v, "0", r_movesbase, CODELOC);

    return ds;
}

/** Constructor */
MovesBase::MovesBase() : QSharedData()
{}

/** Destructor */
MovesBase::~MovesBase()
{}

///////////
/////////// Implementation of SameMoves
///////////

static const RegisterMetaType<SameMoves> r_samemoves;

/** Serialise to a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds,
                                          const SameMoves &samemoves)
{
    writeHeader(ds, r_samemoves, 1);

    SharedDataStream sds(ds);

    sds << samemoves.single_move << static_cast<const MovesBase&>(samemoves);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, SameMoves &samemoves)
{
    VersionID v = readHeader(ds, r_samemoves);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> samemoves.single_move >> static_cast<MovesBase&>(samemoves);
    }
    else
        throw version_error(v, "1", r_samemoves, CODELOC);

    return ds;
}

/** Null constructor - constructs a set of null moves */
SameMoves::SameMoves() : MovesBase()
{}

/** Construct a set of moves where they are all 'move' */
SameMoves::SameMoves(const Move &move)
          : MovesBase(), single_move(move)
{}

/** Copy constructor */
SameMoves::SameMoves(const SameMoves &other)
          : MovesBase(), single_move(other.single_move)
{}

/** Destructor */
SameMoves::~SameMoves()
{}

///////////
/////////// Implementation of Moves
///////////

static const RegisterMetaType<Moves> r_moves;

/** Serialise to a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds, const Moves &moves)
{
    writeHeader(ds, r_moves, 1);

    //we can only serialise the moves when they are not busy!
    Moves &non_const_moves = const_cast<Moves&>(moves);

    non_const_moves.pause();

    //serialise the moves
    SharedDataStream sds(ds);

    sds << moves.d
        << moves.sysid
        << moves.ntotal
        << moves.ncomplete;

    //we can now unpause the moves
    non_const_moves.play();

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, Moves &moves)
{
    VersionID v = readHeader(ds, r_moves);

    if (v == 1)
    {
        //destream into a new Moves object
        Moves new_moves;

        SharedDataStream sds(ds);

        sds >> new_moves.d
            >> new_moves.sysid
            >> new_moves.ntotal
            >> new_moves.ncomplete;

        //copy the new Moves object
        moves = new_moves;
    }
    else
        throw version_error(v, "1", r_moves, CODELOC);

    return ds;
}

static SharedPolyPointer<MovesBase> shared_null( new SameMoves() );

/** Null constructor */
Moves::Moves() : d(shared_null), sysid(0), ntotal(0), ncomplete(0)
{}

/** Construct from the passed moves */
Moves::Moves(const MovesBase &moves)
      : d(moves), sysid(0), ntotal(0), ncomplete(0)
{}

/** Construct from a passed Move - all moves in the set will
    be this move */
Moves::Moves(const Move &move)
      : d(new SameMoves(move)), sysid(0), ntotal(0), ncomplete(0)
{}

/** Copy constructor */
Moves::Moves(const Moves &other)
{
    *this = other;
}

/** Destructor */
Moves::~Moves()
{
    //we cannot delete this object until it has finished!
    QMutexLocker lkr(&movemutex);
}

/** Copy assignment */
Moves& Moves::operator=(const Moves &other)
{
    if (this == &other)
        return *this;

    //we need to wait for this copy to stop running (if it is)
    QMutexLocker lkr(&movemutex);

    Moves &non_const_other = const_cast<Moves&>(other);

    //ok, we are definitely not running - is 'other' running?
    if (non_const_other.movemutex.tryLock())
    {
        //neither we, nor 'other' are running - we can take
        //a straight copy
        d = other.d;

        sysid = other.sysid;
        ntotal = other.ntotal;
        ncomplete = other.ncomplete;

        non_const_other.movemutex.unlock();
    }
    else
    {
        //'other' is running - it needs to be paused to be
        //deep-copied
        non_const_other.pause();

        //take a deep copy
        d = other.d->clone();

        sysid = other.sysid;
        ntotal = other.ntotal;
        ncomplete = other.ncomplete;

        //we can restart 'other'
        non_const_other.play();
    }

    return *this;
}

/** Initialise all of the moves with the passed System
     - this is used to catch a lot of errors before the
       simulation starts, rather than while it is running
       (or even months into it running!) */
void Moves::assertCompatibleWith(QuerySystem &system) const
{
    //can only do this while the moves aren't active
    QMutexLocker lkr( const_cast<QMutex*>(&movemutex) );
    d->assertCompatibleWith(system);
}

/** Set the energy component that all of the moves in this
    set will follow */
void Moves::setEnergyComponent(const Symbol &component)
{
    //can only do this while the moves aren't active
    QMutexLocker lkr(&movemutex);

    d->setEnergyComponent(component);
}

/** Return a copy of the moves that are part of this run.
    This allows you to query how they performed during
    the simulation */
QList<Move> Moves::moves()
{
    //can only do this while the moves aren't active
    QMutexLocker lkr(&movemutex);

    return d->moves();
}

/** This internal function is used to run the
    moves from ncompleted to ntotal */
void Moves::_pvt_run(SimSystem &system, bool record_stats)
{
    if (ncomplete >= ntotal)
        //there's nothing to do
        return;

    //get a pointer to the MovesBase implementation -
    // this detaches the object from shared storage
    MovesBase *moves = d.data();

    //run the moves, one at a time
    do
    {
        //get a write lock...
        QMutexLocker lkr(&runmutex);

        //check that we have not exceeded the total
        //number of moves to run (e.g. we've finished
        //the job, or the total number of moves has
        //been set to zero because the job has been aborted)
        if (ncomplete >= ntotal)
            break;

        //get the next move to run...
        Move &nextmove = moves->nextMove();

        //apply this move to the system
        nextmove.move(system);

        if (record_stats)
            //commit the move - this increments the statistics
            system.commit();

        //increment the number of complete moves
        ++ncomplete;

    } while(true);
}

/** Run 'nmoves' moves of the system. This will run the
    moves in a way that allows the simulation to be paused
    between moves. If record_stats is true, then the
    system will be monitored at the end of each move. */
void Moves::run(SimSystem &system, quint32 nmoves, bool record_stats)
{
    //we need to ensure that only a single system is
    //being run at a time
    QMutexLocker lkr(&movemutex);

    ntotal = nmoves;
    ncomplete = 0;

    sysid = system.ID();

    //run the moves
    this->_pvt_run(system, record_stats);
}

/** Continue the moves in this set on the passed system.
    This function is used to continue a set of moves that
    were part of a snapshot */
void Moves::resume(SimSystem &system, bool record_stats)
{
    //we need to ensure that only a single system
    //is being run at a time
    QMutexLocker lkr(&movemutex);

    sysid = system.ID();

    if (ncomplete >= ntotal)
        //the simulation has finished
        return;

    //run the moves
    this->_pvt_run(system, record_stats);
}

/** Rerun the moves in this set on the passed system.
    While the same number of moves will be attempted,
    there is no guarantee that the moves will be the same,
    or in the same order as the last run. */
void Moves::rerun(SimSystem &system, bool record_stats)
{
    QMutexLocker lkr(&movemutex);

    ncomplete = 0;
    sysid = system.ID();

    if (ncomplete >= ntotal)
        return;

    this->_pvt_run(system, record_stats);
}

/** Pause the simulation - the thread running the simulation is
    blocked until the simulation is unpaused (played) */
void Moves::pause()
{
    if (not movemutex.tryLock())
    {
        //get a read-lock on the simulation lock
        runmutex.lock();

        if (ncomplete >= ntotal)
            //the simulation is complete or has been aborted
            //No need to pause it any more
            runmutex.unlock();
    }
    else
        //the simulation is not running
        movemutex.unlock();
}

/** Abort the simulation - this will stop the moves and allow
    the thread running the moves to move on. */
void Moves::stop()
{
    if (not movemutex.tryLock())
    {
        QMutexLocker lkr(&runmutex);

        //set ntotal == ncomplete - this tells the simulation to stop
        ntotal = ncomplete;
    }
    else
        //the simulation is not running
        movemutex.unlock();
}

/** Play the simulation - playing a simulation that is not paused
    is an error, and will result in program termination. */
void Moves::play()
{
    if (not movemutex.tryLock())
    {
        if (ncomplete < ntotal)
            runmutex.unlock();
    }
    else
        //the simulation is not running
        movemutex.unlock();
}

/** Return the number of moves necessary to complete one cycle
    of this set */
int Moves::count() const
{
    return d->count();
}

/** Return the percentage completion of the run
    (-1 implies not started, 100 == finished) */
int Moves::percentProgress() const
{
    //this doesn't need to be protected with
    //a mutex as it doesn't matter if the data
    //is a little out of date

    if (ntotal == 0)
        return -1;
    else
        return ncomplete * 100 / ntotal;
}
