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

#ifndef SIREMOVE_MPIPROMISE_H
#define SIREMOVE_MPIPROMISE_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMove
{

namespace detail
{
class MPIWorker;
}

/** This class provides the template-independent parts of the 
    the MPIPromise class
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT MPIPromiseBase
{
public:
    ~MPIPromiseBase();
    
    bool isNull() const;

    void abort();
    void stop();
    
    double progress();
    
    bool hasFinished();
    bool isError();
    bool isAborted();

    void wait();
    
    bool wait(int ms);
    
protected:
    MPIPromiseBase();
    MPIPromiseBase(const boost::shared_ptr<detail::MPIWorker> &worker);

    QByteArray resultData();
    QByteArray interimResultData();

private:
    boost::shared_ptr<detail::MPIWorker> worker_ptr;
};

template<class T>
class SIREMOVE_EXPORT MPIPromise : public MPIPromiseBase
{

friend class MPINode;  // so can call MPIWorker constructor

public:
    MPIPromise();
    MPIPromise(const MPIPromise<T> &other);
    
    ~MPIPromise();

    const T& result();
    const T& interimResult();

private:
    MPIPromise(const boost::shared_ptr<detail::MPIWorker> &worker);

    /** The current value of the result */
    boost::shared_ptr<T> current_result;
    
    /** The interim result */
    boost::shared_ptr<T> interim_result;
};

template<class T>
MPIPromise<T>::MPIPromise();

template<class T>
MPIPromise<T>::MPIPromise(const MPIPromise<T> &other);

template<class T>
MPIPromise<T>::MPIPromise(const boost::shared_ptr<detail::MPIWorker> &worker);

template<class T>
MPIPromise<T>::~MPIPromise();

template<class T>
const T& MPIPromise<T>::result()
{
    if (current_result.get() == 0)
    {
        QByteArray result_data = MPIPromiseBase::resultData();
        
        if (result_data.isEmpty())
        {
            //the job was aborted!
            current_result.reset();
            interim_result.reset();
        }
        else
        {
            QDataStream ds(result_data);
        
            current_result.reset( new T() );
        
            ds >> *current_result;
        
            interim_result = current_result;
        }
    }
    
    return *current_result;
}

template<class T>
const T& MPIPromise<T>::interimResult()
{
    if (current_result.get() == 0)
    {
        QByteArray result_data = MPIPromiseBase::interimResultData();
        
        if (result_data.isEmpty())
        {
            //the job was aborted!
            interim_result.reset();
        }
        else
        {
            QDataStream ds(result_data);
            
            interim_result.reset( new T() );
            
            ds >> *interim_result;
        }
    }
    
    return *interim_result;
}

}

SIRE_END_HEADER

#endif
