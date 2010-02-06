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

#include "waitcondition.h"
#include "mutex.h"
#include "forages.h"

using namespace Siren;

WaitCondition::WaitCondition()
{}

WaitCondition::~WaitCondition()
{}

const unsigned long chunk = 5000;

bool WaitCondition::wait( Mutex *m, unsigned long time )
{
    if (not m)
        return false;

    if (time == ULONG_MAX )
    {
        while (for_ages())
        {
            if (w.wait(&(m->m), chunk))
                return true;
        }
        
        return false;
    }
    else
    {
        if (time <= chunk)
            return w.wait( &(m->m), time );
    
        while (for_ages())
        {
            int wait = qMin(time, chunk);
            
            if (w.wait(&(m->m), wait))
                return true;
                
            time -= wait;
            
            if (time <= 0)
                return false;
        }
        
        return false;
    }
}

bool WaitCondition::wait( QMutex *m, unsigned long time )
{
    if (not m)
        return false;

    if (time == ULONG_MAX )
    {
        while (for_ages())
        {
            if (w.wait(m, chunk))
                return true;
        }
        
        return false;
    }
    else
    {
        if (time <= chunk)
            return w.wait(m, time);
    
        while (for_ages())
        {
            int wait = qMin(time, chunk);
            
            if (w.wait(m, wait))
                return true;
                
            time -= wait;
            
            if (time <= 0)
                return false;
        }
        
        return false;
    }
}

void WaitCondition::wakeOne()
{
    w.wakeOne();
}

void WaitCondition::wakeAll()
{
    w.wakeAll();
}