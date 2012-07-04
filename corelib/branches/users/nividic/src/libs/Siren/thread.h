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

#ifndef SIREN_THREAD_H
#define SIREN_THREAD_H

#include <QThread>

#include "sirenglobal.h"

#include "mutex.h"
#include "waitcondition.h"

SIREN_BEGIN_HEADER

namespace Siren
{

void init_rand(int ID=0);

/** This is a thin wrapper around QThread that provides
    automatic registration of Siren functions, and also
    provides a catch for any uncaught exceptions from
    the thread. This prevents a background thread from 
    crashing the entire program if it loses an exception.
    The exception is passed to the global error handler 
    
    Note that derived classes must reimplement
    Thread::threadMain rather than QThread::run
    
    @author Christopher Woods
*/
class SIREN_EXPORT Thread : public QThread
{
public:
    Thread(QObject *parent=0);
    Thread(QString name, QObject *parent=0);
    virtual ~Thread();

    void start();

protected:
    void run();

    void signalStarted();
    
    virtual void threadMain()=0;

private:
    /** The name of this thread */
    QString thread_name;
    
    /** Mutex used to wait for the thread to start */
    Mutex start_mutex;
    
    /** Waitcondition used to wait for the thread to start */
    WaitCondition start_waiter;
    
    bool thread_is_starting;
};

}

SIREN_END_HEADER

#endif